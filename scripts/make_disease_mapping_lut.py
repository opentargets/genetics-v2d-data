#!/usr/bin/env python3
"""Generates a look-up table that includes all disease mappings present in the Genetics Portal studies."""

import argparse
import logging
from pathlib import Path, PurePath
import tempfile

import pandas as pd

from extract_therapeutic_areas import *

def main(
    study_index: str,
    finngen_spreadsheet: str,
    ukbb_old_mappings: str,
    ukbb_new_mappings: str,
    output_disease_lut: str  
) -> None:

    # 1. Extract mappings per data source GWAS catalog traits from study table (these do not require OT mapping)
    gwas_catalog_mappings = get_gwas_catalog_mappings(study_index)
    valid_ukb = get_ukbb_mappings(ukbb_old_mappings, ukbb_new_mappings)
    valid_finngen = get_finngen_mappings(finngen_spreadsheet)
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
        .assign(trait_efos=lambda x: x.proposed_efos.combine_first(x.trait_efos)).drop('proposed_efos', axis=1)
        .explode('trait_efos')
    )

    assert len(genetics_mappings) == (
        len(gwas_catalog_mappings.explode('trait_efos')) +
        len(valid_ukb.explode('proposed_efos')) +
        len(valid_finngen.explode('proposed_efos'))
    ), "WARNING! Some mappings went missing during the merge."

    # 3. Bring therapeutic areas
    genetics_mappings_w_ta = build_therapeutic_areas(genetics_mappings)
    genetics_mappings_w_ta['therapeutic_area'] = genetics_mappings_w_ta['therapeutic_areas'].apply(get_prioritised_therapeutic_area)
    logging.info('EFO loaded. Therapeutic areas built.')

    # Check everything is an ontology ID and that there are no mappings without a TA
    assert genetics_mappings_w_ta['trait_efos'].str.contains('\w+_\d+', regex=True).all() == False, 'WARNING! There are invalid EFO IDs'
    assert len(genetics_mappings_w_ta[genetics_mappings_w_ta['therapeutic_area'].isna()]), 'WARNING! There are EFO IDs without a therapeutic area.'
    genetics_mappings_w_ta = genetics_mappings_w_ta.loc[genetics_mappings_w_ta['trait_efos'].str.contains('\w+_\d+', regex=True), :]

    # 4. Format and write output
    (
        genetics_mappings_w_ta
        .groupby(['study_id', 'trait_reported', 'therapeutic_area'], dropna=False).agg(lambda x: list(set(x))).reset_index()
        .to_parquet(output_disease_lut)
    )
    logging.info(f'{output_disease_lut} successfully generated. Exiting.')

def read_input_file(path: str):
    """
    Automatically detect the format of the input data and read it into the Spark dataframe. The supported formats
    are: a single CSV file; a directory with Parquet files.
    """
    if 'parquet' in path:
        path = Path(path)
        if Path.is_dir(path):
            data_dir = Path(path)
            full_df = pd.concat(
                pd.read_parquet(parquet_file)
                for parquet_file in data_dir.glob('*.parquet')
            )
            return pd.read_parquet(full_df)
        if Path.is_file(path):
            return pd.read_parquet(path)
    elif 'json' in path:
        return pd.read_json(path, lines=True)
    else:
        return pd.read_csv(path.replace('/edit?usp=sharing', '/export?format=csv'))

def get_gwas_catalog_mappings(
    study_index: str
) -> pd.DataFrame:
    """Extracts GWAS catalog trait mappings from study table (these do not require OT mapping)."""
    
    studies_df = read_input_file(study_index)
    return (
        studies_df[studies_df['study_id'].str.startswith('GCST')]
        .filter(items=['study_id', 'trait_reported', 'trait_efos'])
    )

def get_ukbb_mappings(
    ukbb_old_mappings: str,
    ukbb_new_mappings: str
) -> pd.DataFrame:
    """
    UK Biobank trait mappings are merged from two different sources:
    1. ukbb_old_df: Initial curation done using EFO2.
    2. ukbb_new_df: Review of some studies from ukbb_old_df where the mapping is incorrect or is updated.
    A study being in both datasets means that it is either incorrect or updated:
    If updated (candidate = True) --> The new mapping is kept.
    If incorrect (candidate = False) --> The mapping is dropped.
    """
    ukbb_old_df = get_ukbb_old_mappings(ukbb_old_mappings)
    ukbb_new_df = get_ukbb_new_mappings(ukbb_new_mappings)

    return (
        ukbb_old_df.merge(ukbb_new_df, on=['study_id'], how='outer', indicator=True)
        .query("_merge == 'left_only' or _merge == 'both' and candidate == True")
        .assign(proposed_efos=lambda X: X.proposed_efos_y.combine_first(X.proposed_efos_x))
        .filter(items=['study_id', 'proposed_efos'])
    )

    
def get_ukbb_new_mappings(
    ukbb_new_mappings: str
) -> pd.DataFrame:
    """Extracts valid UK Biobank trait mappings from the curation spreadsheet."""

    return (
        read_input_file(ukbb_new_mappings)
        .rename(columns={'traitName': 'trait_reported', 'candidateId':'proposed_efos'})
        .query('candidate == True | current == True')
        .dropna(how='all')
        .filter(items=['study_id', 'trait_reported', 'current', 'currentEfo', 'candidate', 'proposed_efos'])
    )

def get_ukbb_old_mappings(
    ukbb_old_mappings: str
) -> pd.DataFrame:
    """Extracts initial UK Biobank trait mappings from the curation JSON ."""

    return (
        read_input_file(ukbb_old_mappings)
        .drop('curation_confidence', axis=1)
        .rename(columns={'efos': 'proposed_efos'})
        .explode('proposed_efos')
    )

def get_finngen_mappings(
    finngen_spreadsheet: str
) -> pd.DataFrame:
    """Extracts Finngen trait mappings from the curation spreadsheet."""
    
    return (
        read_input_file(finngen_spreadsheet)
        .query('valid == True')
        .filter(items=['NAME', 'LONGNAME', 'efo_cls'])
        # Trim all strings to have a clean mapped id
        .apply(lambda x: x.str.strip())
        # Group data
        .groupby('NAME').agg(lambda x: list(set(x))).reset_index()
        .rename(columns={
            'NAME':'study_name',
            'LONGNAME':'trait_reported',
            'efo_cls':'proposed_efos'
        })
        .explode('trait_reported')
        .assign(study_id=lambda x: 'FINNGEN_R5_' + x.study_name)
        .drop('study_name', axis=1)
    )

def build_therapeutic_areas(
    genetics_mappings: pd.DataFrame
) -> pd.DataFrame:
    """Therapeutic areas per trait are built into the mappings table."""
    
    efo_tas_df = extract_therapeutic_areas_from_owl()
    genetics_mappings_w_trait = (
        genetics_mappings
        .merge(
            efo_tas_df,
            left_on='trait_efos',
            right_on='efo_id',
            how='left')
        .drop('efo_id', axis=1)
    )

    return genetics_mappings_w_trait

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--in_studies', help='Directory of parquet files that stores the study index.', required=True)
    parser.add_argument('--in_finngen_mappings', help='URL of the spreadsheet that contains all Finngen disease mappings', nargs='?', default='https://docs.google.com/spreadsheets/d/1yrQPpsRi-mijs_BliKFZjeoxP6kGIs9Bz-02_0WDvAA/edit?usp=sharing')
    parser.add_argument('--in_ukbb_old_mappings', help='JSON file that contains the initial UK Biobank disease mappings. File can be found in gs://genetics-portal-input/ukb_phenotypes/ukb_efo_annotation.190828.json', required=True)
    parser.add_argument('--in_ukbb_new_mappings', help='URL of the spreadsheet that contains the updated UK Biobank disease mappings resulting from upgrading to EFO3', default='https://docs.google.com/spreadsheets/d/1PotmUEirkV36dh-vpZ3GgxQg_LcOefZKbyTq0PNQ6NY/edit?usp=sharing')
    parser.add_argument('--out_disease_lut', help='Parquet file that stores all disease mappings present in the Genetics Portal studies', required=True)
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format='%(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )

    # Report input data:
    logging.info(f'Study index file path: {args.in_studies}')
    logging.info(f'Finngen curation spreadsheet URL: {args.in_finngen_mappings}')
    logging.info(f'UK Biobank curation spreadsheet URL: {args.in_ukbb_new_mappings}')
    logging.info(f'UK Biobank initial curation file path: {args.in_ukbb_old_mappings}')

    main(
        args.in_studies,
        args.in_finngen_mappings,
        args.in_ukbb_old_mappings,
        args.in_ukbb_new_mappings,
        args.out_disease_lut
    )
