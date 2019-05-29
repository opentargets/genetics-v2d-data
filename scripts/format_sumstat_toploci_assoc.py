#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
import argparse
import pandas as pd
from decimal import Decimal
from collections import OrderedDict
import scipy.stats as st
import numpy as np

def main():

    # Parse args
    args = parse_args()

    # Load top loci
    top_loci = pd.read_json(args.inf, orient='records', lines=True)

    # Only keep type == gwas
    top_loci = top_loci.loc[top_loci['type'] == 'gwas', :]

    # Load study information, required for case-control information GCST000964
    study = pd.read_json(args.study_info, orient='records', lines=True)
    study['case_prop'] = study['n_cases'] / study['n_initial']
    study = study.loc[:, ['study_id', 'case_prop']]
    study = study.drop_duplicates()

    # Fix very low p-values
    top_loci.loc[top_loci['pval'] == 0, 'pval'] = (1 / sys.float_info.max)
    
    # Filter p-values
    top_loci = top_loci.loc[top_loci['pval'] <= 5e-8, :]

    # Get P mantissa and exponent
    top_loci.loc[:, 'p_mantissa'] = top_loci['pval'].apply(fman)
    top_loci.loc[:, 'p_exponent'] = top_loci['pval'].apply(fexp)

    # Extract effect size information
    top_loci = pd.merge(top_loci, study, on='study_id', how='left')
    top_loci[['direction',
              'beta',
              'beta_ci_lower',
              'beta_ci_upper',
              'oddsr',
              'oddsr_ci_lower',
              'oddsr_ci_upper']] = top_loci.apply(extract_effect_sizes, axis=1).apply(pd.Series)

    # Make a variant ID
    top_loci['variant_id_b38'] = (
        top_loci.loc[:, ['chrom', 'pos', 'ref', 'alt']].apply(
            lambda row: '_'.join([str(x) for x in row.tolist()]),
            axis=1
        )
    )

    # Extract and rename required columns
    cols = OrderedDict([
        ('study_id', 'study_id'),
        ('variant_id_b38', 'variant_id_b38'),
        # ('locus_index_snp', 'rsid'),
        ('direction', 'direction'),
        ('beta', 'beta'),
        ('beta_ci_lower', 'beta_ci_lower'),
        ('beta_ci_upper', 'beta_ci_upper'),
        ('oddsr', 'odds_ratio'),
        ('oddsr_ci_lower', 'oddsr_ci_lower'),
        ('oddsr_ci_upper', 'oddsr_ci_upper'),
        ('p_mantissa', 'pval_mantissa'),
        ('p_exponent', 'pval_exponent')
    ])
    assert all([col in top_loci.columns for col in cols.keys()])
    top_loci = top_loci[list(cols.keys())].rename(columns=cols)

    # Sort and save
    top_loci = top_loci.sort_values(['study_id', 'pval_exponent', 'pval_mantissa'])
    top_loci.to_csv(args.outf, sep='\t', index=None)

def extract_effect_sizes(row):
    ''' Extract beta, oddsr, lower_ci, upper_ci
        - extract relevant stat
        - impute 95% CI
    Args:
        row (pd.Series)
    Returns:,
        direction, beta, beta_ci_lower, beta_ci_upper, oddsr, oddsr_ci_lower, oddsr_ci_upper
    '''
    # Initiate
    direction, beta, beta_ci_lower, beta_ci_upper, oddsr, oddsr_ci_lower, oddsr_ci_upper = None, None, None, None, None, None, None
    # Check whether quantitative or binary
    is_beta = pd.isnull(row['case_prop'])
    # Estimate z-score
    z = abs(st.norm.ppf((float(row['pval']))/2))
    # Extract relevant stats
    if is_beta:
        beta = row['beta']
        se = row['se']
        beta_ci_lower = beta - 1.96 * se
        beta_ci_upper = beta + 1.96 * se
        direction = '+' if beta >= 0 else '-'
    else:
        # Perform transformation: https://github.com/opentargets/sumstat_data#requirements-when-adding-new-datasets
        log_or = row['beta']
        log_se = row['se']
        oddsr = np.exp(log_or)
        oddsr_ci_lower = np.exp(log_or - 1.96 * log_se)
        oddsr_ci_upper = np.exp(log_or + 1.96 * log_se)
        direction = '+' if oddsr >= 1 else '-'

    return direction, beta, beta_ci_lower, beta_ci_upper, oddsr, oddsr_ci_lower, oddsr_ci_upper

def fexp(number):
    ''' https://stackoverflow.com/a/45359185 '''
    (sign, digits, exponent) = Decimal(number).as_tuple()
    return len(digits) + exponent - 1

def fman(number):
    ''' https://stackoverflow.com/a/45359185 '''
    return Decimal(number).scaleb(-fexp(number)).normalize()

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--inf', metavar="<file>", type=str, required=True)
    parser.add_argument('--study_info', metavar="<file>", type=str, required=True)
    parser.add_argument('--outf', metavar="<str>", type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
