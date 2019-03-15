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
    ld = pd.read_csv(args.inf, sep='\t', header=0)

    # Filter
    ld = ld.loc[ld.R2_overall >= args.min_r2, :]

    # Format
    cols = OrderedDict([
        ('study_id', 'study_id'),
        ('index_variant_id_b37', 'index_variantid_b37'),
        ('tag_variant_id_b37', 'tag_variantid_b37'),
        ('R2_overall', 'overall_r2'),
        ('AFR_prop', 'AFR_1000G_prop'),
        ('AMR_prop', 'AMR_1000G_prop'),
        ('EAS_prop', 'EAS_1000G_prop'),
        ('EUR_prop', 'EUR_1000G_prop'),
        ('SAS_prop', 'SAS_1000G_prop')
        ])
    ld = ( ld.loc[:, list(cols.keys())]
             .rename(columns=cols) )

    # Save
    write_parquet(
        ld,
        args.outf
    )
    # ld.to_csv(args.outf, sep='\t', index=None, compression='gzip')

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--inf', metavar="<file>", help=('GWAS Catalog input'), type=str, required=True)
    parser.add_argument('--outf', metavar="<str>", help=("Output"), type=str, required=True)
    parser.add_argument('--min_r2', metavar="<str>", help=("Minimum R2 to be included"), type=float, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
