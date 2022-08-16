"""
Processes Finngen's manifest to extract all studies and their metadata in the OTG format.
"""
# coding: utf-8

import argparse
from collections import OrderedDict
import logging

import numpy as np
import pandas as pd


def main(input_path: str, output_path: str) -> None:

    logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.INFO)

    # Read manifest
    manifest = (
        pd.read_json(input_path, orient='records').filter(
            items=['phenocode', 'phenostring', 'category', 'num_cases', 'num_controls']
        )
        # When phenostring is not provided, phenotype extracted from the phenocode
        .assign(
            phenostring=lambda df: df.apply(
                lambda row: row['phenostring'] if row['phenostring'] and row['phenostring'] != '' else row['phenocode'],
                axis=1,
            )
        )
        # Renaming columns to accomodate OTG schema:
        .rename(
            columns={
                'phenocode': 'study_id',
                'phenostring': 'trait',
                'category': 'trait_category',
                'num_cases': 'n_cases',
                'num_controls': 'n_controls',
            }
        )
    )

    logging.info(f"{input_path} has been loaded. Formatting...")

    # Format table:
    manifest['study_id'] = 'FINNGEN_R6_' + manifest['study_id']
    manifest['n_total'] = manifest['n_cases'] + manifest['n_controls']
    manifest['pmid'] = ''
    manifest['pub_date'] = '2022-01-24'
    manifest['pub_author'] = 'FINNGEN_R6'
    manifest['ancestry_initial'] = 'European=' + manifest['n_total'].astype(str)
    manifest['n_replication'] = 0
    manifest['ancestry_replication'] = ''
    manifest['pub_journal'] = ''
    manifest['pub_title'] = ''
    manifest['trait_efos'] = np.nan

    cols = OrderedDict(
        [
            ('study_id', 'study_id'),
            ('pmid', 'pmid'),
            ('pub_date', 'pub_date'),
            ('pub_journal', 'pub_journal'),
            ('pub_title', 'pub_title'),
            ('pub_author', 'pub_author'),
            ('trait', 'trait_reported'),
            ('trait_efos', 'trait_efos'),
            ('ancestry_initial', 'ancestry_initial'),
            ('ancestry_replication', 'ancestry_replication'),
            ('n_total', 'n_initial'),
            ('n_cases', 'n_cases'),
            ('n_replication', 'n_replication'),
        ]
    )

    manifest = manifest.loc[:, list(cols.keys())].rename(columns=cols)

    manifest.to_json(output_path, orient='records', lines=True)
    logging.info(f"{len(manifest)} studies have been saved in {output_path}. Exiting.")


def parse_args():
    """
    Load command line args.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', metavar="<str>", type=str, required=True)
    parser.add_argument('--output', metavar="<str>", help=("Output"), type=str, required=True)

    return parser.parse_args()


if __name__ == '__main__':

    # Parse args
    args = parse_args()

    main(input_path=args.input, output_path=args.output)
