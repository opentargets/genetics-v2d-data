#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import argparse
from collections import OrderedDict

import numpy as np
import pandas as pd


def main():

    # Parse args
    args = parse_args()
    
    # Only keep required cols
    to_keep = OrderedDict([
        ('code', 'study_id'),
        ('n_total', 'n_total'),
        ('n_cases', 'n_cases')
    ])

    # Load manifest
    manifest = pd.read_csv(args.in_manifest, sep='\t',
                           header=0, dtype=object) \
                           .filter(items=to_keep)

    #
    # Add other columns -------------------------------------------------------
    #

    # Vector to get Neale or SAIGE studies
    is_neale = manifest['study_id'].str.startswith("NEALE2_")
    is_saige = manifest['study_id'].str.startswith("SAIGE_")
    assert (is_neale | is_saige).all()

    # Make columns
    manifest.loc[:, 'pmid'] = ''
    manifest.loc[is_neale, 'pub_date'] = '2018-08-01'
    manifest.loc[is_saige, 'pub_date'] = '2018-10-24'
    manifest.loc[:, 'pub_journal'] = ''
    manifest.loc[:, 'pub_title'] = ''
    manifest.loc[is_neale, 'pub_author'] = 'UKB Neale v2'
    manifest.loc[is_saige, 'pub_author'] = 'UKB SAIGE'
    manifest.loc[:, 'n_initial'] = manifest['n_total'].apply(to_int_safe)
    manifest.loc[:, 'n_cases'] = manifest['n_cases'].apply(to_int_safe)
    manifest.loc[:, 'n_replication'] = 0
    manifest.loc[:, 'ancestry_initial'] = 'European=' + manifest['n_initial'].astype(str)
    manifest.loc[:, 'ancestry_replication'] = ''

    # Ouput required columns
    cols = OrderedDict([
        ('study_id', 'study_id'),
        ('pmid', 'pmid'),
        ('pub_date', 'pub_date'),
        ('pub_journal', 'pub_journal'),
        ('pub_title', 'pub_title'),
        ('pub_author', 'pub_author'),
        ('ancestry_initial', 'ancestry_initial'),
        ('ancestry_replication', 'ancestry_replication'),
        ('n_initial', 'n_initial'),
        ('n_cases', 'n_cases'),
        ('n_replication', 'n_replication')
        ])
    manifest = ( manifest.loc[:, list(cols.keys())]
                         .rename(columns=cols) )

    # Write
    manifest.to_json(args.outf, orient='records', lines=True)

def to_int_safe(i):
    try:
        return int(float(i))
    except ValueError:
        return None

def parse_args():
    """
    Load command line args.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('--in_manifest', metavar="<str>", help=("Input"), type=str, required=True)
    parser.add_argument('--outf', metavar="<str>", help=("Output"), type=str, required=True)
    args = parser.parse_args()

    return args

if __name__ == '__main__':

    main()
