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

    # Load credible set info
    cred = pd.read_csv(args.inf, sep='\t', header=0)

    # Load study information, required for case-control information
    study = pd.read_csv(args.study_info, sep='\t', header=0)
    study['case_prop'] = study['n_cases'] / study['n_initial']
    study = study.loc[:, ['study_id', 'case_prop']]

    # Keep rows where index == var
    top_loci = cred.loc[cred.locus_index_varid == cred.varid, :]

    # Make study id
    top_loci.loc[:, 'study_id'] = 'NEALEUKB_' + top_loci['trait'].astype(str)

    # Fix very low p-values
    top_loci.loc[top_loci.p == 0, 'p'] = (1 / sys.float_info.max)

    # Filter p-values
    top_loci = top_loci.loc[top_loci.p < 1e-5, :]

    # Get P mantissa and exponent
    top_loci.loc[:, 'p_mantissa'] = top_loci.p.apply(fman)
    top_loci.loc[:, 'p_exponent'] = top_loci.p.apply(fexp)

    # Extract effect size information
    top_loci = pd.merge(top_loci, study, how='left')
    top_loci[['direction', 'beta', 'beta_ci_lower', 'beta_ci_upper', 'oddsr', 'oddsr_ci_lower', 'oddsr_ci_upper']] = top_loci.apply(extract_effect_sizes, axis=1).apply(pd.Series)

    # Extract and rename required columns
    cols = OrderedDict([
        ('study_id', 'study_id'),
        ('locus_index_varid', 'variant_id_b37'),
        ('locus_index_snp', 'rsid'),
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
    z = abs(st.norm.ppf((float(row['p']))/2))
    # Extract relevant stats
    if is_beta:
        beta = row['b']
        se = row['se']
        beta_ci_lower = beta - 1.96 * se
        beta_ci_upper = beta + 1.96 * se
        direction = '+' if beta >= 0 else '-'
    else:
        # Perform transformation: https://github.com/opentargets/sumstat_data#requirements-when-adding-new-datasets
        log_or = row['b'] / (row['case_prop'] / (1 - row['case_prop']))
        log_se = row['se'] / (row['case_prop'] / (1 - row['case_prop']))
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
