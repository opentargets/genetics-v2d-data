#!/usr/bin/env python3
"""Generates a look-up table that includes all disease mappings present in the Genetics Portal studies."""

import argparse
from pathlib import Path, PurePath
import tempfile

import pandas as pd

from extract_therapeutic_areas import *

def main(
    study_index: str,
    finngen_spreadsheet: str,
    ukbb_spreadsheet: str,
    output_disease_lut: str  
) -> None:
    # Load data
    studies_df = read_input_file(study_index)
    finngen_df = read_input_file(finngen_spreadsheet)
    ukbb_df = read_input_file(ukbb_spreadsheet)

    # 1. Extract mappings per data source GWAS catalog traits from study table (these do not require OT mapping)
    gwas_catalog_mappings = get_gwas_catalog_mappings(studies_df)
    valid_ukb = get_ukbb_mappings(ukbb_df)
    valid_finngen = get_finngen_mappings(finngen_df)
    # Assert there are no studies with a null mapping
    for source in [gwas_catalog_mappings, valid_ukb, valid_finngen]:
        if 'proposed_efos' not in source.columns:
            null_studies = source[source['trait_reported'].isna()]
        else:
            source = source.explode('proposed_efos')
            null_studies = source[source['proposed_efos'].isna()]
        assert len(null_studies) == 0, f"Studies in {source} contain invalid mappings."
    

    # 4. Merge updated mappings across data sources
    genetics_mappings = (
        pd.concat([valid_finngen, valid_ukb, gwas_catalog_mappings], ignore_index=True)
        # Coalesce all mappings in trait_efos
        .assign(trait_efos=lambda x: x.proposed_efos.combine_first(x.trait_efos)).drop('proposed_efos', axis=1)
        .explode('trait_efos')
        .drop_duplicates()
    )

    # 5. Bring therapeutic areas
    genetics_mappings_w_ta = build_therapeutic_areas(genetics_mappings)
    genetics_mappings_w_ta['therapeutic_area'] = genetics_mappings_w_ta['therapeutic_areas'].apply(get_prioritised_therapeutic_area)

    # Check everything is an ontology ID and that there are no mappings without a TA
    assert genetics_mappings_w_ta['trait_efos'].str.contains('\w+_\d+', regex=True).all() == False, 'WARNING! There are invalid EFO IDs'
    assert len(genetics_mappings_w_ta[genetics_mappings_w_ta['therapeutic_area'].isna()]), 'WARNING! There are EFO IDs without a therapeutic area.'
    genetics_mappings_w_ta = genetics_mappings_w_ta.loc[genetics_mappings_w_ta['trait_efos'].str.contains('\w+_\d+', regex=True), :]

    # 6. Format and write output
    output_table = (genetics_mappings_w_ta
        .groupby(['study_id', 'trait_reported', 'therapeutic_area']).agg(lambda x: list(set(x))).reset_index()
    )
    output_table.to_parquet(output_disease_lut)

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
    else:
        return pd.read_csv(path.replace('/edit?usp=sharing', '/export?format=csv'))

def get_gwas_catalog_mappings(
    studies_df: pd.DataFrame
) -> pd.DataFrame:
    """Extracts GWAS catalog traits from study table (these do not require OT mapping)."""
    
    return (
        studies_df[studies_df['study_id'].str.startswith('GCST')]
        .filter(items=['study_id', 'trait_reported', 'trait_efos'])
    )

def get_ukbb_mappings(
    ukbb_df: pd.DataFrame
) -> pd.DataFrame:
    """Extracts UK Biobank traits from the curation spreadsheet."""

    return (ukbb_df
            .query('candidate == True')
            .dropna(how='all')
            .filter(items=['study_id', 'traitName', 'candidateId'])
            # Trim potential whitespaces
            .apply(lambda x: x.str.strip())
            # Group data
            .groupby('study_id').agg(lambda x: list(set(x))).reset_index()
            .rename(columns={
                'traitName': 'trait_reported',
                'candidateId': 'proposed_efos'
            })
            .explode('trait_reported')
    )

def get_finngen_mappings(
    finngen_df: pd.DataFrame
) -> pd.DataFrame:
    """Extracts UK Biobank traits from the curation spreadsheet."""
    
    return (finngen_df
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

def write_evidence_strings(output_table: pd.DataFrame, output_file: str) -> None:
    """Exports the table to a parquet file."""
    with tempfile.TemporaryDirectory() as tmp_dir_name:
        (
            output_table.coalesce(1).write.format('parquet').mode('overwrite')
            .option('compression', 'org.apache.hadoop.io.compress.GzipCodec').save(tmp_dir_name)
        )
        parquet_chunks = [f for f in Path.iterdir(tmp_dir_name) if f.endswith('.parquet')]
        assert len(parquet_chunks) == 1, f'Expected one Parquet file, but found {len(parquet_chunks)}.'
        Path.rename(PurePath.joinpath.join(tmp_dir_name, parquet_chunks[0]), output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--in_studies', help='Directory of parquet files that stores the study index.', required=True)
    parser.add_argument('--in_finngen_mappings', help='URL of the spreadsheet that contains all Finngen disease mappings', nargs='?', default='https://docs.google.com/spreadsheets/d/1yrQPpsRi-mijs_BliKFZjeoxP6kGIs9Bz-02_0WDvAA/edit?usp=sharing')
    parser.add_argument('--in_ukbb_mappings', help='URL of the spreadsheet that contains the updated UK Biobank disease mappings resulting from upgrading to EFO3', default='https://docs.google.com/spreadsheets/d/1iTGRVPXsHizNXdDnj0on9zfURjTLP8A73mn7CNZiVw8/edit?usp=sharing')
    parser.add_argument('--out_disease_lut', help='Parquet file that stores all disease mappings present in the Genetics Portal studies', required=True)
    args = parser.parse_args()
    main(
        args.in_studies,
        args.in_finngen_mappings,
        args.in_ukbb_mappings,
        args.out_disease_lut
    )
