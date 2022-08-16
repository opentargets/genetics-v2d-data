#!/usr/bin/env python3
"""Generates a look-up table that includes all disease mappings present in the Genetics Portal studies."""

import argparse
import logging

import pandas as pd
from pathlib import Path

from get_therapeutic_areas import *


def main(
    studies: str, finngen_version: int, finngen_mappings: str, ukb_original_mappings: str, ukb_updated_mappings: str, output_path: str
) -> None:

    # 1. Extract mappings per data source GWAS catalog traits from study table (these do not require OT mapping)
    gwas_catalog_mappings = get_gwas_catalog_mappings(studies)
    valid_ukb = get_ukb_mappings(ukb_original_mappings, ukb_updated_mappings)
    valid_finngen = get_finngen_mappings(finngen_version, finngen_mappings)
    # Assert there are no studies with a null mapping
    for source in [gwas_catalog_mappings, valid_ukb, valid_finngen]:
        if 'proposed_efos' not in source.columns:
            null_studies = source[source['trait_reported'].isna()]
        else:
            source = source.explode('proposed_efos')
            null_studies = source[source['proposed_efos'].isna()]
        assert len(null_studies) == 0, f"Studies in {source} contain invalid mappings."
    logging.info('Disease mappings loaded and split per data source.')

    # 2. Merge updated mappings across data sources
    genetics_mappings = (
        pd.concat([valid_finngen, valid_ukb, gwas_catalog_mappings], ignore_index=True)
        # Coalesce all mappings in trait_efos
        .assign(trait_efos=lambda x: x.proposed_efos.combine_first(x.trait_efos))
        .drop('proposed_efos', axis=1)
        .explode('trait_efos')
    )
    assert len(genetics_mappings) == (
        len(gwas_catalog_mappings.explode('trait_efos'))
        + len(valid_ukb.explode('proposed_efos'))
        + len(valid_finngen.explode('proposed_efos'))
    ), "WARNING! Some mappings went missing during the merge."

    # Check everything is an ontology ID and that there are no mappings without a TA
    assert (
        genetics_mappings['trait_efos'].str.contains('\w+_\d+', regex=True).all() == False
    ), 'WARNING! There are invalid EFO IDs'
    genetics_mappings = genetics_mappings.loc[genetics_mappings['trait_efos'].str.contains('\w+_\d+', regex=True), :]

    # 3. Bring therapeutic areas
    genetics_mappings_w_ta = (
        build_therapeutic_areas(genetics_mappings)
        # A study/trait can be mapped to multiple EFOs, each with a different set of therapeutic areas.
        # All the therapeutic areas and EFOs are collected into the same column. The most significant TA
        # per study is extracted. The result of collecting these is a multidimensional array that must be flattened.
        .groupby(['study_id', 'trait_reported'], dropna=False)
        .agg({'therapeutic_areas': list, 'trait_efos': list})
        .reset_index()
    )
    genetics_mappings_w_ta['therapeutic_areas'] = genetics_mappings_w_ta['therapeutic_areas'].apply(
        lambda X: flatten_array(X)
    )

    # Extract the most relevant TA from the array
    genetics_mappings_w_ta['trait_category'] = genetics_mappings_w_ta['therapeutic_areas'].apply(
        get_prioritised_therapeutic_area
    )
    genetics_mappings_w_ta.drop('therapeutic_areas', axis=1, inplace=True)
    logging.info('EFO loaded. Therapeutic areas built.')

    assert len(
        genetics_mappings_w_ta[genetics_mappings_w_ta['trait_category'].isna()]
    ), 'WARNING! There are EFO IDs without a therapeutic area.'
    assert len(genetics_mappings_w_ta) == len(
        genetics_mappings_w_ta['study_id'].unique()
    ), 'WARNING! There are duplicated studies.'

    # Assert no studies are lost in the process of adding the TA
    assert len(genetics_mappings_w_ta.study_id.unique()) == len(
        genetics_mappings.study_id.unique()
    ), 'WARNING! Some studies were lost in the process of adding the TA.'

    # 4. Format and write output
    genetics_mappings_w_ta.to_parquet(output_path)
    logging.info(f'{output_path} successfully generated. Exiting.')


def read_input_file(path: str):
    """
    Automatically detect the format of the input data and read it into the Spark dataframe. The supported formats
    are: a single CSV file; a directory with Parquet files.
    """
    if 'parquet' in path:
        path = Path(path)
        if Path.is_dir(path):
            data_dir = Path(path)
            full_df = pd.concat(pd.read_parquet(parquet_file) for parquet_file in data_dir.glob('*.parquet'))
            return pd.read_parquet(full_df)
        if Path.is_file(path):
            return pd.read_parquet(path)
    elif 'json' in path:
        return pd.read_json(path, lines=True)
    elif 'spreadsheet' in path:
        return pd.read_csv(path.replace('/edit?usp=sharing', '/export?format=csv'))
    elif 'csv' in path:
        return pd.read_csv(path)
    else:
        raise ValueError(f'Unsupported input file format: {path}')


def get_gwas_catalog_mappings(studies: str) -> pd.DataFrame:
    """Extracts GWAS catalog trait mappings from study table (these do not require OT mapping)."""

    studies_df = read_input_file(studies)
    return studies_df[studies_df['study_id'].str.startswith('GCST')].filter(
        items=['study_id', 'trait_reported', 'trait_efos']
    )


def get_ukb_mappings(ukb_original_mappings: str, ukb_updated_mappings: str) -> pd.DataFrame:
    """
    UK Biobank trait mappings are merged from two different sources:
    1. ukbb_old_df: Initial curation done using EFO2.
    2. ukbb_new_df: Review of some studies from ukbb_old_df where the mapping is incorrect or is updated.
    A study being in both datasets means that it is either incorrect or updated:
    If updated (candidate = True) --> The new mapping is kept.
    If incorrect (candidate = False) --> The mapping is dropped.
    """
    ukbb_old_df = get_ukb_original_mappings(ukb_original_mappings)
    ukbb_new_df = get_ukb_updated_mappings(ukb_updated_mappings)

    return (
        ukbb_old_df.merge(ukbb_new_df, on=['study_id'], how='outer', indicator=True)
        .query("_merge == 'left_only' or _merge == 'both' and candidate == True")
        # Coalesce EFOs into a single column
        .assign(proposed_efos=lambda X: X.proposed_efos_y.combine_first(X.proposed_efos_x))
        .filter(items=['study_id', 'trait_reported', 'proposed_efos'])
    )


def get_ukb_updated_mappings(ukb_updated_mappings: str) -> pd.DataFrame:
    """Extracts valid UK Biobank trait mappings from the curation spreadsheet."""

    return (
        read_input_file(ukb_updated_mappings)
        .rename(columns={'traitName': 'trait_reported', 'candidateId': 'proposed_efos'})
        .query('candidate == True | current == True')
        .dropna(how='all')
        .filter(items=['study_id', 'trait_reported', 'current', 'currentEfo', 'candidate', 'proposed_efos'])
    )


def get_ukb_original_mappings(ukb_original_mappings: str) -> pd.DataFrame:
    """Extracts original UK Biobank trait mappings from the curation JSON ."""

    return (
        read_input_file(ukb_original_mappings)
        .drop('curation_confidence', axis=1)
        .rename(columns={'efos': 'proposed_efos'})
        .explode('proposed_efos')
    )


def get_finngen_mappings(finngen_version: int, finngen_mappings: str) -> pd.DataFrame:
    """
    Extracts Finngen trait mappings from the curation spreadsheet
    
    Args:
      finngen_version (int): The version of the Finngen data that you are using.
      finngen_mappings (str): The path to the Finngen trait mappings spreadsheet.
    
    Returns:
      A dataframe with the following columns:
        - study_id
        - trait_reported
        - proposed_efos
    """

    version = f'FINNGEN_R{finngen_version}_'

    return (
        read_input_file(finngen_mappings)
        .query('valid == True')
        .filter(items=['NAME', 'LONGNAME', 'efo_cls'])
        # Trim all strings to have a clean mapped id
        .apply(lambda x: x.str.strip())
        # Group data
        .groupby('NAME')
        .agg(lambda x: list(set(x)))
        .reset_index()
        .rename(columns={'NAME': 'study_name', 'LONGNAME': 'trait_reported', 'efo_cls': 'proposed_efos'})
        .explode('trait_reported')
        .assign(study_id=lambda x: version + x.study_name)
        .drop('study_name', axis=1)
    )


def build_therapeutic_areas(genetics_mappings: pd.DataFrame) -> pd.DataFrame:
    """Therapeutic areas per trait are built into the mappings table."""
    efo_tas_df = extract_therapeutic_areas_from_owl()
    return genetics_mappings.merge(efo_tas_df, left_on='trait_efos', right_on='efo_id', how='left').drop('efo_id', axis=1)


def flatten_array(arr: List) -> List:  # sourcery skip: use-contextlib-suppress
    """Flattens a bidimensional array."""
    try:
        return [i for sublist in arr for i in sublist]
    except Exception:
        pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--studies', help='Location of file with all studies.', required=True)
    parser.add_argument(
        '--finngen_mappings',
        help='URL of the spreadsheet that contains all Finngen disease mappings',
        nargs='?',
        default='https://docs.google.com/spreadsheets/d/1RRWfUTLy4TO9XmBzcbJ2wPRdda3qISjRS4PJmEdxE3k/edit?usp=sharing',
    )
    parser.add_argument(
        '--finngen_version',
        help='The version of the Finngen manifest the study table is based on.',
        required=True, type=int,
    )
    parser.add_argument(
        '--ukb_original_mappings',
        help='JSON file that contains the initial UK Biobank disease mappings. File can be found in gs://genetics-portal-input/ukb_phenotypes/ukb_efo_annotation-2021-10-25.json',
        required=True,
    )
    parser.add_argument(
        '--ukb_updated_mappings',
        help='URL of the spreadsheet that contains the updated UK Biobank disease mappings resulting from upgrading to EFO3',
        default='https://docs.google.com/spreadsheets/d/1PotmUEirkV36dh-vpZ3GgxQg_LcOefZKbyTq0PNQ6NY/edit?usp=sharing',
    )
    parser.add_argument(
        '--output', help='LUT that includes all disease mappings present in the Genetics Portal studies.', required=True
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format='%(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )

    # Report input data:
    logging.info(f'Collection of all studies path: {args.studies}')
    logging.info(f'Finngen curation spreadsheet URL: {args.finngen_mappings}')
    logging.info(f'UK Biobank initial curation file path: {args.ukb_original_mappings}')
    logging.info(f'UK Biobank curation updates spreadsheet URL: {args.ukb_updated_mappings}')
    logging.info(f'Output disease/EFO LUT: {args.output}')

    main(
        studies=args.studies,
        finngen_mappings=args.finngen_mappings,
        finngen_version=args.finngen_version,
        ukb_original_mappings=args.ukb_original_mappings,
        ukb_updated_mappings=args.ukb_updated_mappings,
        output_path=args.output,
    )
