#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#
# Merge top loci table
#

import sys
import argparse
import pandas as pd
from collections import OrderedDict
from parquet_writer import write_parquet
import datetime

def main():

    # Parse args
    args = parse_args()

    # Load
    gwas = pd.read_csv(args.in_gwascat, sep='\t', header=0)
    neale = pd.read_csv(args.in_neale, sep='\t', header=0)

    # Merge together
    merged = pd.concat([gwas, neale], sort=False)

    # Load the number of associated loci per study
    num_loci = ( pd.read_parquet(args.in_toploci, engine='pyarrow')
                   .groupby('study_id')
                   .size().reset_index(name='counts')
                   .rename(columns={'counts': 'num_assoc_loci'}) )
    merged = pd.merge(merged, num_loci, how='left', on='study_id')
    merged['num_assoc_loci'] = merged['num_assoc_loci'].fillna(value=0).astype(int)

    # Coerce data types
    dtypes = OrderedDict([
        ('study_id', 'object'),
        ('pmid', 'object'),
        ('pub_date', 'object'),
        ('pub_journal', 'object'),
        ('pub_title', 'object'),
        ('pub_author', 'object'),
        ('trait_reported', 'object'),
        ('trait_efos', 'object'),
        ('ancestry_initial', 'object'),
        ('ancestry_replication', 'object'),
        ('n_initial', 'Int64'),
        ('n_replication', 'Int64'),
        ('n_cases', 'Int64'),
        ('trait_category', 'object'),
        ('num_assoc_loci', 'Int64')
    ])
    assert(set(dtypes.keys()) == set(merged.columns))
    merged = (
        merged.loc[:, dtypes.keys()]
        .astype(dtype=dtypes)
    )

    # Split array columns
    array_cols = ['trait_efos', 'ancestry_initial', 'ancestry_replication']
    for col in array_cols:
        merged[col] = merged[col].str.split(';')

    # Sort output
    merged = merged.sort_values(['study_id'])

    # Save as parquet
    write_parquet(merged,
                  args.output,
                  str_list_cols=array_cols,
                  compression='snappy',
                  flavor='spark')

    # merged.to_csv(args.output.replace('.parquet', '.test.tsv'),
    #               sep='\t', index=None)


    return 0

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_gwascat', metavar="<str>", type=str, required=True)
    parser.add_argument('--in_neale', metavar="<str>", type=str, required=True)
    parser.add_argument('--in_toploci', metavar="<str>", type=str, required=True)
    parser.add_argument('--output', metavar="<str>", help=("Output merged file"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
