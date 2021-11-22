# coding: utf-8

from collections import OrderedDict

import argparse
import numpy as np
import pandas as pd


def main():

    # Parse args
    args = parse_args()

    # read manifest
    FINNGEN_manifest = (pd.read_json(args.in_manifest, lines=True).rename(
            columns={
                'phenocode': 'study_id',
                'phenostring': 'trait',
                'category': 'trait_category',
                'num_cases': 'n_cases',
                'num_controls': 'n_controls'
            }
        )
    )

    keep_columns = [
        'study_id',
        'trait',
        'trait_category',
        'n_cases',
        'n_controls'
    ]
    FINNGEN_manifest = FINNGEN_manifest[keep_columns]

    # Format table:
    FINNGEN_manifest['study_id'] = 'FINNGEN_R5_' + FINNGEN_manifest['study_id']
    FINNGEN_manifest['n_total'] = FINNGEN_manifest['n_cases'] + FINNGEN_manifest['n_controls']
    FINNGEN_manifest['pmid'] = ''
    FINNGEN_manifest['pub_date'] = '2021-5-11'
    FINNGEN_manifest['pub_author'] = 'FINNGEN_R5'
    FINNGEN_manifest['ancestry_initial'] = 'European=' + FINNGEN_manifest['n_total'].astype(str)
    FINNGEN_manifest['n_replication'] = 0
    FINNGEN_manifest['ancestry_replication'] = ''
    FINNGEN_manifest['pub_journal'] = ''
    FINNGEN_manifest['pub_title'] = ''
    FINNGEN_manifest['trait_efos'] = np.nan

    cols = OrderedDict([
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
        ('n_replication', 'n_replication')
        ])

    FINNGEN_manifest = FINNGEN_manifest.loc[:, list(cols.keys())].rename(columns=cols)

    # Write
    FINNGEN_manifest.to_json(args.outf, orient='records', lines=True)

def parse_args():
    """
    Load command line args.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_manifest', metavar="<str>", type=str, required=True)
    parser.add_argument('--outf', metavar="<str>", help=("Output"), type=str, required=True)
    
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()

#    FINNGEN_manifest_path="/home/xg1/r5_finngen.json"
#    FINNGEN_EFO_path="/home/xg1/finngen_df5_efo_mapping.lastest.nonull.tsv"
# "/home/xg1/genetics-v2d-data/tmp/210526/merged_study_table.old.json"
# "/home/xg1/genetics-v2d-data/tmp/210526/merged_study_table.json"

#python Manual_FINNGENs --in_manifest r5_finngen.json --in_EFO finngen_df5_efo_mapping.lastest.nonull.tsv 
#    --in_study_table tmp/version_date/merged_study_table.old.json --outf tmp/version_date/merged_study_table.json

