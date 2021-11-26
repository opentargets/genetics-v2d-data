#!/usr/bin/env python3
"""Generates a look-up table that includes all disease mappings present in the Genetics Portal studies."""

import argparse
from pathlib import Path
from typing import List, Union

import pandas as pd

# A dict of all therapeutic areas and their ids ranked by order of relevance
THERAPEUTIC_AREAS = {
    'therapeutic_area': [
        'cell proliferation disorder', 'infectious disease',
        'pregnancy or perinatal disease', 'animal disease',
        'disease of visual system', 'cardiovascular disease',
        'pancreas disease', 'gastrointestinal disease',
        'reproductive system or breast disease', 'integumentary system disease',
        'endocrine system disease', 'respiratory or thoracic disease',
        'urinary system disease', 'musculoskeletal or connective tissue disease',
        'disease of ear', 'immune system disease',
        'hematologic disease', 'nervous system disease',
        'psychiatric disorder', 'nutritional or metabolic disease',
        'genetic, familial or congenital disease', 'injury, poisoning or other complication',
        'phenotype', 'measurement', 'biological process'],
    'id': [
        'MONDO_0045024', 'EFO_0005741', 'OTAR_0000014',
        'EFO_0005932', 'MONDO_0024458', 'EFO_0000319',
        'EFO_0009605', 'EFO_0010282', 'OTAR_0000017',
        'EFO_0010285', 'EFO_0001379', 'OTAR_0000010',
        'EFO_0009690', 'OTAR_0000006', 'MONDO_0021205',
        'EFO_0000540', 'EFO_0005803', 'EFO_0000618',
        'MONDO_0002025', 'MONDO_0024297', 'OTAR_0000018',
        'OTAR_0000009', 'EFO_0000651','EFO_0001444',
        'GO_0008150']
}


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
    disease_index_df = read_input_file(disease_index)

    # 1. Extract static mappings from study table
    static_mappings = (
        # Drop studies without a mapping
        studies_df[studies_df['trait_efos'].str.len() > 0]
        .filter(items=['study_id', 'trait_reported', 'trait_efos'])
        .explode('trait_efos')
        # Group data
        .groupby('study_id').agg(lambda x: list(set(x))).reset_index()
        # Drop Finngen studies to later bring them all from the Finngen dataset
        .query('study_id.str.contains("FINNGEN")==False', engine='python')
        .explode('trait_reported')
    )

    # 2. Extract new updated mappings from UK Biobank
    valid_ukb = (ukbb_df
        .query('candidate == True')
        .dropna(how='all')
        .filter(items=['study_id', 'traitName', 'candidateId'])
        # Trim all strings to have a clean mapped id
        .apply(lambda x: x.str.strip())
        # Group data
        .groupby('study_id').agg(lambda x: list(set(x))).reset_index()
        .rename(columns={
            'traitName':'trait_reported',
            'candidateId':'proposed_efos'
        })
        .explode('trait_reported')
    )

    # 3. Extract mappings from Finngen
    valid_finngen = (finngen_df
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
    )
    valid_finngen['study_id'] = 'FINNGEN_R5_' + valid_finngen['study_id']

    # 4. Join static with dinamycally imported mappings
    genetics_mappings = (static_mappings
        .merge(
            pd.concat([valid_finngen, valid_ukb], ignore_index=True),
            on=["study_id", "trait_reported"],
            how="outer"
        ))
    # Coalesce trait_efos to include the updated mappings 
    genetics_mappings['trait_efos'] = genetics_mappings['proposed_efos'].combine_first(genetics_mappings.trait_efos)
    genetics_mappings.drop(columns='proposed_efos', inplace=True)

    # 5. Bring therapeutic areas
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

    # Check everything is an ontology ID
    assert genetics_mappings_w_ta['trait_efos'].str.contains('\w+_\d+', regex=True).all() == False, 'WARNING! There is at least one mapping with an invalid ID'
    genetics_mappings_w_ta = genetics_mappings_w_ta.loc[genetics_mappings_w_ta['trait_efos'].str.contains('\w+_\d+', regex=True), :]
    genetics_mappings_final = (genetics_mappings_w_ta
        .groupby('study_id').agg(lambda x: list(set(x))).reset_index()
    )
    genetics_mappings_final.to_parquet(output_disease_lut)

def read_input_file(path: str):
    """Automatically detect the format of the input data and read it into the Spark dataframe. The supported formats
    are: a single CSV file; a directory with Parquet files."""
    if 'parquet' in path:
        data_dir = Path(path)
        full_df = pd.concat(
            pd.read_parquet(parquet_file)
            for parquet_file in data_dir.glob('*.parquet')
        )
        return pd.read_parquet(full_df)
    elif 'csv' in path:
        return pd.read_csv(path)

def get_more_relevant_ta(
    tas: List[str]
) -> Union[List[str], str]:
    """Uses the index of the therapeutic areas df to get the label of the more relevant TA by selecting which one has the minimal index."""
    sorted_tas = pd.DataFrame(data=THERAPEUTIC_AREAS)
    try:
        if len(tas) > 0:
            min_index = float('inf')
            for ta in tas:
                idx = sorted_tas.index[sorted_tas['id'] == ta]
                if idx < min_index: 
                    min_index = idx
            ta = sorted_tas.iloc[min_index]["therapeutic_area"].values[0]
            return ta
    except TypeError:
        return "Uncategorised"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--in_studies', help='Directory of parquet files that stores the study index.', required=True)
    parser.add_argument('--in_finngen-mappings', help='URL of the spreadsheet that contains all Finngen disease mappings', required=True)
    parser.add_argument('--in_ukbb-mappings', help='URL of the spreadsheet that contains the updated UK Biobank disease mappings resulting from upgrading to EFO3', required=True)
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