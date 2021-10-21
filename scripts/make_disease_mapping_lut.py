#!/usr/bin/env python3
"""Generates a look-up table that includes all disease mappings present in the Genetics Portal studies."""

import argparse
from pathlib import Path
from typing import List, Union

import pandas as pd

from extract_therapeutic_areas import *

def main(
    study_index: str,
    finngen_spreadsheet: str,
    ukbb_spreadsheet: str,
    disease_index: str,
    output_disease_lut: str  
):
    # Load data
    studies_df = read_input_file(study_index)
    finngen_df = read_input_file(
        finngen_spreadsheet.replace('/edit?usp=sharing', '/export?format=csv')
    )
    ukbb_df = read_input_file(
        ukbb_spreadsheet.replace('/edit?usp=sharing', '/export?format=csv')
    )

    # 1. Extract mappings per data source GWAS catalog traits from study table (these do not require OT mapping)
    gwas_catalog_mappings = get_gwas_catalog_mappings(studies_df)
    valid_ukb = get_ukbb_mappings(ukbb_df)
    valid_finngen = get_finngen_mappings(finngen_df)

    # 4. Join static with dinamycally imported mappings TODO: concat and not mergesl
    genetics_mappings = (gwas_catalog_mappings
        .merge(
            pd.concat([valid_finngen, valid_ukb], ignore_index=True),
            on=["study_id", "trait_reported"],
            how="outer"
        ))
    # Coalesce trait_efos to include the updated mappings 
    genetics_mappings['trait_efos'] = genetics_mappings['proposed_efos'].combine_first(genetics_mappings.trait_efos)
    genetics_mappings.drop(columns='proposed_efos', inplace=True)

    # 5. Bring therapeutic areas
    '''
    genetics_mappings_w_ta = (genetics_mappings
        .explode('trait_efos')
        # Get list of all TAs from the disease index
        .merge(
            disease_index_df.filter(items=['id', 'therapeuticAreas']),
            left_on='trait_efos',
            right_on='id',
            how='left'
        )
        .drop(columns='id')
        .explode('therapeuticAreas')
        # Group data
        .groupby('study_id').agg(lambda x: list(set(x))).reset_index()
        .explode('trait_reported')
    )
    genetics_mappings_w_ta['trait_category'] = genetics_mappings_w_ta['therapeuticAreas'].apply(lambda X: get_more_relevant_ta(X))
    genetics_mappings_w_ta = genetics_mappings_w_ta.explode('trait_category').explode('trait_efos')
    '''
    tas = extract_therapeutic_areas_from_owl(OWL_FILENAME)
    genetics_mappings_w_ta = (

    )

    # Check everything is an ontology ID
    assert genetics_mappings_w_ta['trait_efos'].str.contains('\w+_\d+', regex=True).all() == False, 'WARNING! There is at least one mapping with an invalid ID'
    genetics_mappings_w_ta = genetics_mappings_w_ta.loc[genetics_mappings_w_ta['trait_efos'].str.contains('\w+_\d+', regex=True), :]
    genetics_mappings_final = (genetics_mappings_w_ta
        .groupby('study_id').agg(lambda x: list(set(x))).reset_index()
    )
    genetics_mappings_final.to_parquet(output_disease_lut)

def read_input_file(path: str):
    """
    Automatically detect the format of the input data and read it into the Spark dataframe. The supported formats
    are: a single CSV file; a directory with Parquet files.
    """
    if 'parquet' in path:
        data_dir = Path(path)
        full_df = pd.concat(
            pd.read_parquet(parquet_file)
            for parquet_file in data_dir.glob('*.parquet')
        )
        return pd.read_parquet(full_df)
    elif 'csv' in path:
        return pd.read_csv(path)

def get_gwas_catalog_mappings(
    studies_df: pd.DataFrame
) -> pd.DataFrame:
    """Extracts GWAS catalog traits from study table (these do not require OT mapping)."""
    
    return  (
        # Drop studies without a mapping
        studies_df[studies_df['trait_efos'].str.len() > 0]
        .filter(items=['study_id', 'trait_reported', 'trait_efos'])
        .explode('trait_efos')
        # Group data
        .groupby('study_id').agg(lambda x: list(set(x))).reset_index()
        # Drop Finngen studies to later bring them all from the Finngen dataset
        #.query('study_id.str.contains("FINNGEN")==False', engine='python')
        .explode('trait_reported')
    )

def get_ukbb_mappings(
    ukbb_df: pd.DataFrame
) -> pd.DataFrame:

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
    
    return (finngen_df
        .query('valid == True')
        .filter(items=['NAME', 'LONGNAME', 'efo_cls'])
        # Trim all strings to have a clean mapped id
        .apply(lambda x: x.str.strip())
        # Group data
        .groupby('NAME').agg(lambda x: list(set(x))).reset_index()
        .rename(columns={
            'NAME':'study_id',
            'LONGNAME':'trait_reported',
            'efo_cls':'proposed_efos'
        })
        .explode('trait_reported')
        .assign(study_id=lambda x: 'FINNGEN_R5_' + x.study_id)
    )
    #valid_finngen['study_id'] = 'FINNGEN_R5_' + valid_finngen['study_id']

def get_therapeutic_areas():
    efo_owl_parsed_df = extract_therapeutic_areas_from_owl(owl_url)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--in_studies', help='Directory of parquet files that stores the study index.', default='https://docs.google.com/spreadsheets/d/1yrQPpsRi-mijs_BliKFZjeoxP6kGIs9Bz-02_0WDvAA/edit?usp=sharing', required=True)
    parser.add_argument('--in_finngen-mappings', help='URL of the spreadsheet that contains all Finngen disease mappings', required=True)
    parser.add_argument('--in_ukbb-mappings', help='URL of the spreadsheet that contains the updated UK Biobank disease mappings resulting from upgrading to EFO3', default='https://docs.google.com/spreadsheets/d/1iTGRVPXsHizNXdDnj0on9zfURjTLP8A73mn7CNZiVw8/edit?usp=sharing', required=True)
    parser.add_argument('--in_disease-index', help='Directory of parquet files that stores the OT disease index to extract the therapeutic areas', required=True)
    parser.add_argument('--out_disease-lut', help='Parquet file that stores all disease mappings present in the Genetics Portal studies', required=True)
    args = parser.parse_args()
    main(
        args.in_studies,
        args.in_finngen-mappings,
        args.in_ukbb-mappings,
        args.in_disease-index,
        args.out_disease-lut
    )