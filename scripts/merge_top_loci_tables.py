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
    neale = pd.read_csv(args.in_neale, sep='\t', header=0)

    # Merge
    merged = pd.concat([gwas, neale], sort=False)

    # Split variant ID into chrom, pos, ref, alt
    merged[['chrom', 'pos', 'ref', 'alt']] = \
        merged.variant_id_b37.str.split('_', expand=True)
    merged.pos = merged.pos.astype(int)
    merged.drop('variant_id_b37', axis=1, inplace=True)

    # Fix data types
    dtypes = OrderedDict([
        ('study_id', 'object'),
        ('chrom', 'object'),
        ('pos', 'Int64'),
        ('ref', 'object'),
        ('alt', 'object'),
        # ('rsid', 'object'),
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

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_gwascat', metavar="<str>", type=str, required=True)
    parser.add_argument('--in_neale', metavar="<str>", type=str, required=True)
    parser.add_argument('--output', metavar="<str>", help=("Output merged file"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
