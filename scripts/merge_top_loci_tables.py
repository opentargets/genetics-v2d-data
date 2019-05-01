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

def main():

    # Parse args
    args = parse_args()

    # Load
    gwas = pd.read_csv(args.in_gwascat, sep='\t', header=0)
    sumstat = pd.read_csv(args.in_sumstat, sep='\t', header=0)
    
    #
    # Remove GCST studies from gwas, that appear in sumstat -------------------
    #

    # Get set of GWAS Calaog IDs that have sumstat results
    stid_in_sumstat = set([
        return_original_gwascat_id(idx)
        for idx in sumstat['study_id']
    ])

    # Identify rows to remove from gwas
    to_remove = [
        return_original_gwascat_id(idx) in stid_in_sumstat
        for idx in gwas['study_id']
    ]
    
    # Remove from gwas catalog results
    gwas = gwas.loc[~pd.Series(to_remove), :]

    #
    # Merge -------------------------------------------------------------------
    #

    # Merge
    merged = pd.concat([gwas, sumstat], sort=False)

    # Split variant ID into chrom, pos, ref, alt
    merged[['chrom', 'pos', 'ref', 'alt']] = (
        merged['variant_id_b38'].str.split('_', expand=True)
    )
    merged.pos = merged.pos.astype(int)
    merged.drop(['variant_id_b38', 'rsid'], axis=1, inplace=True)

    # Fix data types
    dtypes = OrderedDict([
        ('study_id', 'object'),
        ('chrom', 'object'),
        ('pos', 'Int64'),
        ('ref', 'object'),
        ('alt', 'object'),
        ('direction', 'object'),
        ('beta', 'float64'),
        ('beta_ci_lower', 'float64'),
        ('beta_ci_upper', 'float64'),
        ('odds_ratio', 'float64'),
        ('oddsr_ci_lower', 'float64'),
        ('oddsr_ci_upper', 'float64'),
        ('pval_mantissa', 'float64'),
        ('pval_exponent', 'Int64')
    ])
    assert(set(dtypes.keys()) == set(merged.columns))
    merged = (
        merged.loc[:, dtypes.keys()]
        .astype(dtype=dtypes)
    )

    # Sort values
    merged = merged.sort_values(
        ['study_id', 'chrom', 'pos', 'ref', 'alt']
    )

    # Save as parquet
    write_parquet(merged,
                  args.output,
                  compression='snappy',
                  flavor='spark')

    return 0

def return_original_gwascat_id(s):
    ''' Strips suffix from gwascatalog IDs
    '''
    if s.startswith('GCST'):
        return s.split('_')[0]
    else:
        return s

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_gwascat', metavar="<str>", type=str, required=True)
    parser.add_argument('--in_sumstat', metavar="<str>", type=str, required=True)
    parser.add_argument('--output', metavar="<str>", help=("Output merged file"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
