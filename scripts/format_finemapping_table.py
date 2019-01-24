#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
import argparse
import pandas as pd
from pprint import pprint
from collections import OrderedDict
from parquet_writer import write_parquet

def main():

    # Parse args
    args = parse_args()

    # Load
    credset = pd.read_csv(args.inf, sep='\t', header=0)

    # Make study id
    credset.loc[:, 'study_id'] = 'NEALEUKB_' + credset['trait'].astype(str)

    # Filter to remove rows not in a 95% credible set
    credset = credset.loc[credset.is95_credset == 1, :]

    # Decompose variant IDs
    credset[['lead_chrom', 'lead_pos', 'lead_ref', 'lead_alt']] = \
        credset.locus_index_varid.str.split('_', expand=True)
    credset[['tag_chrom', 'tag_pos', 'tag_ref', 'tag_alt']] = \
        credset.varid.str.split('_', expand=True)
    credset['lead_pos'] = credset['lead_pos'].astype(int)
    credset['tag_pos'] = credset['tag_pos'].astype(int)

    # Rename and select columns
    cols = OrderedDict([
        ('study_id', 'study_id'),
        ('lead_chrom', 'lead_chrom'),
        ('lead_pos', 'lead_pos'),
        ('lead_ref', 'lead_ref'),
        ('lead_alt', 'lead_alt'),
        ('tag_chrom', 'tag_chrom'),
        ('tag_pos', 'tag_pos'),
        ('tag_ref', 'tag_ref'),
        ('tag_alt', 'tag_alt'),
        ('logABF', 'log10_ABF'),
        ('postprob', 'posterior_prob')
        ])
    credset = ( credset.loc[:, list(cols.keys())]
                       .rename(columns=cols) )

    # Coerce data types
    dtypes = OrderedDict([
        ('study_id', 'object'),
        ('lead_chrom', 'object'),
        ('lead_pos', 'Int64'),
        ('lead_ref', 'object'),
        ('lead_alt', 'object'),
        ('tag_chrom', 'object'),
        ('tag_pos', 'Int64'),
        ('tag_ref', 'object'),
        ('tag_alt', 'object'),
        ('log10_ABF', 'float64'),
        ('posterior_prob', 'float64')
    ])
    assert(set(dtypes.keys()) == set(credset.columns))
    credset = (
        credset.loc[:, dtypes.keys()]
        .astype(dtype=dtypes)
    )

    # Sort
    credset = credset.sort_values(
        ['study_id', 'lead_chrom', 'lead_pos', 'lead_ref', 'lead_alt',
         'tag_chrom', 'tag_pos', 'tag_ref', 'tag_alt']
    )

    # Save as parquet
    write_parquet(credset,
                  args.outf,
                  compression='snappy',
                  flavor='spark')

    # Save csv
    # credset.to_csv(args.outf, sep='\t', index=None, compression='gzip')

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--inf', metavar="<file>", help=('GWAS Catalog input'), type=str, required=True)
    parser.add_argument('--outf', metavar="<str>", help=("Output"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
