#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#
'''
Calculates credible sets using LD data:
- Uses PICS method to adjust p-value for all tag variants
- Convert to approximate Bayes factors using Wakefield et al method
- Credible set analysis
'''

import sys
import os
import pandas as pd
import argparse
import numpy as np

def main():

    pd.set_option('display.max_columns', 500)

    # Parse args
    # args = parse_args()
    args = ArgsPlaceholder()
    args.in_top_loci = '../output/190530/toploci.parquet'
    args.in_ld = '../output/190530/ld.parquet'

    #
    # Load and clean ----------------------------------------------------------
    #

    # Load
    toploci = pd.read_parquet(args.in_top_loci)
    ld = pd.read_parquet(args.in_ld)

    # Drop rows where LD was not available
    ld = ld[ld['ld_available']]

    #
    # Merge P-values onto the LD table ----------------------------------------
    #

    # Calc -log(pval)
    toploci['neglog_p'] = toploci.apply(
        lambda x: -1 * (np.log10(x['pval_mantissa']) + x['pval_exponent']),
        axis=1
    )

    # Drop unneed columns and rename for join
    toploci = (
        toploci
        .loc[:, ['study_id', 'chrom', 'pos',
                 'ref', 'alt', 'neglog_p']]
        .rename(columns={
            'chrom': 'lead_chrom',
            'pos': 'lead_pos',
            'ref': 'lead_ref',
            'alt': 'lead_alt'
        })
    )

    # Join
    ld = pd.merge(
        ld, toploci,
        on=['study_id', 'lead_chrom', 'lead_pos', 'lead_ref', 'lead_alt'],
        how='left'
    )

    # Assert there are no rows without pvals
    problem_studies = (
        ld.loc[
            pd.isnull(ld['neglog_p']),
            ['study_id']
        ].drop_duplicates()
          
    )

    print(problem_studies)
    print(problem_studies.shape)

    # print(ld.loc[pd.isnull(ld['neglog_p']), :])

    # assert (pd.isnull(ld['neglog_p']).sum() == 0)

    # print(ld.head())
    # print(ld.shape)
    


    return 0

class ArgsPlaceholder():
    pass

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_ld', metavar="<str>", help=("Input ld file"), type=str, required=True)
    parser.add_argument('--in_manifest', metavar="<str>", help=("Input manifest file"), type=str, required=True)
    parser.add_argument('--out', metavar="<str>", help=("Output file"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
