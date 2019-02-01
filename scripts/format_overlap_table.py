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
    df = pd.read_csv(args.inf, sep='\t', header=0)

    # Decompose variant IDs
    df[['A_chrom', 'A_pos', 'A_ref', 'A_alt']] = \
        df.index_variantid_b37_A.str.split('_', expand=True)
    df[['B_chrom', 'B_pos', 'B_ref', 'B_alt']] = \
        df.index_variantid_b37_B.str.split('_', expand=True)
    df['A_pos'] = df['A_pos'].astype(int)
    df['B_pos'] = df['B_pos'].astype(int)

    # Rename and select columns
    cols = OrderedDict([
        ('study_id_A', 'A_study_id'),
        ('A_chrom', 'A_chrom'),
        ('A_pos', 'A_pos'),
        ('A_ref', 'A_ref'),
        ('A_alt', 'A_alt'),
        ('study_id_B', 'B_study_id'),
        ('B_chrom', 'B_chrom'),
        ('B_pos', 'B_pos'),
        ('B_ref', 'B_ref'),
        ('B_alt', 'B_alt'),
        ('set_type', 'set_type'),
        ('distinct_A', 'A_distinct'),
        ('overlap_AB', 'AB_overlap'),
        ('distinct_B', 'B_distinct')
        ])
    df = ( df.loc[:, list(cols.keys())]
                       .rename(columns=cols) )

    # Coerce data types
    dtypes = OrderedDict([
        ('A_study_id', 'object'),
        ('A_chrom', 'object'),
        ('A_pos', 'Int64'),
        ('A_ref', 'object'),
        ('A_alt', 'object'),
        ('B_study_id', 'object'),
        ('B_chrom', 'object'),
        ('B_pos', 'Int64'),
        ('B_ref', 'object'),
        ('B_alt', 'object'),
        ('set_type', 'object'),
        ('A_distinct', 'Int64'),
        ('AB_overlap', 'Int64'),
        ('B_distinct', 'Int64')
        ])
    assert(set(dtypes.keys()) == set(df.columns))
    df = (
        df.loc[:, dtypes.keys()]
        .astype(dtype=dtypes)
    )

    # Sort
    df = df.sort_values(
        ['A_study_id', 'A_chrom', 'A_pos', 'A_ref', 'A_alt', 'B_study_id',
         'B_chrom', 'B_pos', 'B_ref', 'B_alt']
    )

    # Save as parquet
    write_parquet(df,
                  args.outf,
                  compression='snappy',
                  flavor='spark')

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--inf', metavar="<str>", help=('input'), type=str, required=True)
    parser.add_argument('--outf', metavar="<str>", help=("Output"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
