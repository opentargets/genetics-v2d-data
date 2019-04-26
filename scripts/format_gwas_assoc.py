#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
import re
import argparse
import pandas as pd
from pprint import pprint
from collections import OrderedDict
import logging
import scipy.stats as st
import numpy as np

def main():

    # Parse args
    args = parse_args()

    # Start logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    handler = logging.FileHandler(args.log)
    logger.addHandler(handler)

    #
    # Load data ----------------------------------------------------------------
    #

    # Load assoc data
    gwas = pd.read_csv(args.inf, sep='\t', header=0, dtype={'P-VALUE':str, '95% CI (TEXT)':str})
    logger.info('Total associations: {0}'.format(gwas.shape[0]))

    # Remove Sun et al pQTL study
    gwas = gwas.loc[gwas['STUDY ACCESSION'] != 'GCST005806', :]
    # Remove Huang et al IBD study (GWAS Catalog should not have curated this)
    gwas = gwas.loc[gwas['STUDY ACCESSION'] != 'GCST005837', :]

    #
    # Merge new study IDs ------------------------------------------------------
    #

    # Merge new study IDs
    id_lut = pd.read_csv(args.id_lut, sep='\t', header=0)
    prerowcount = gwas.shape[0]
    gwas = pd.merge(gwas, id_lut,
                    on=['STUDY ACCESSION', 'P-VALUE (TEXT)'],
                    how='left')

    assert(gwas.study_id.isna().sum() == 0) # Assert no missing study IDs
    assert(gwas.shape[0] == prerowcount)    # Assert no duplicated rows

    #
    # Process ------------------------------------------------------------------
    #

    # Drop rows with no variant ID
    gwas = gwas.dropna(subset=['variant_id_b38'])
    logger.info('Total with variant ID: {0}'.format(gwas.shape[0]))

    # Extract pvalue mantissa and exponent from p-value string
    gwas[['pval_mantissa', 'pval_exponent']] = gwas['P-VALUE'].apply(parse_pval_mantissa_exponent).apply(pd.Series)

    # Explode lines
    gwas['variant_id_b38'] = gwas['variant_id_b38'].astype(str).str.split(';')
    gwas['rsid'] = gwas['rsid'].astype(str).str.split(';')
    gwas = explode(gwas, ['variant_id_b38', 'rsid'])
    logger.info('Total associations after exploding compound variant IDs: {0}'.format(gwas.shape[0]))

    # Extract harmonised effect sizes
    gwas[['direction', 'beta', 'beta_ci_lower', 'beta_ci_upper', 'oddsr', 'oddsr_ci_lower', 'oddsr_ci_upper']] = gwas.apply(parse_harmonised_effect, axis=1).apply(pd.Series)

    # Extract and rename required columns
    cols = OrderedDict([
        ('study_id', 'study_id'),
        ('variant_id_b38', 'variant_id_b38'),
        ('rsid', 'rsid'),
        ('direction', 'direction'),
        ('beta', 'beta'),
        ('beta_ci_lower', 'beta_ci_lower'),
        ('beta_ci_upper', 'beta_ci_upper'),
        ('oddsr', 'odds_ratio'),
        ('oddsr_ci_lower', 'oddsr_ci_lower'),
        ('oddsr_ci_upper', 'oddsr_ci_upper'),
        ('pval_mantissa', 'pval_mantissa'),
        ('pval_exponent', 'pval_exponent')
    ])
    gwas = gwas[list(cols.keys())].rename(columns=cols)

    # Drop duplicates by study, variant and sub_phenotype
    gwas = gwas.sort_values(['pval_exponent', 'pval_mantissa'])
    gwas = gwas.drop_duplicates(subset=['study_id', 'variant_id_b38'],
                                keep='first')
    logger.info('Total after duplicates removed: {0}'.format(gwas.shape[0]))

    gwas = gwas.sort_values(['study_id', 'variant_id_b38'])

    gwas.to_csv(args.outf, sep='\t', index=None)

def parse_harmonised_effect(row):
    ''' Takes a row from GWAS Catalog associations.
        - Parses the beta or OR
        - Harmonises to the alt allele
        - Impute SE and CI using p-value
    Args:
        row (pandas.Serires)
    Returns:
        direction, beta, beta_ci_lower, beta_ci_upper, oddsr, oddsr_ci_lower, oddsr_ci_upper
    '''
    # Initiate
    direction, beta, beta_ci_lower, beta_ci_upper, oddsr, oddsr_ci_lower, oddsr_ci_upper = None, None, None, None, None, None, None

    # Extract risk allele, ref and alt
    risk = extract_risk_allele(row['STRONGEST SNP-RISK ALLELE'])
    ref, alt = row['variant_id_b38'].split(';')[0].split('_')[2:4]

    # Stop if no effect is available, or there is no risk allele reported
    if pd.isnull(row['OR or BETA']) or not risk:
        return None, None, None, None, None, None, None

    # Check if palindromic or ambiguous
    is_palin = (ref == revcomp(alt))
    is_ambiguous = len(row['variant_id_b38'].split(';')) > 1

    # Only proceed if not palindromic and not ambiguous and is concordant
    if (is_palin or is_ambiguous):
        return None, None, None, None, None, None, None

    # Check whether harmonisation is requried
    if risk == alt or risk == revcomp(alt):
        needs_harmonising = False
    elif risk == ref or risk == revcomp(ref):
        needs_harmonising = True
    else:
        # Is not concordant with the ref alt alleles
        return None, None, None, None, None, None, None

    # Estimate z-score from p-value
    z = abs(st.norm.ppf((float(row['P-VALUE']))/2))

    # Extract beta or oddsr
    is_beta = ('increase' in str(row['95% CI (TEXT)']) or
               'decrease' in str(row['95% CI (TEXT)']) )
    if is_beta:
        beta = float(row['OR or BETA'])
        # Flip direction based on '95% CI (TEXT)' field
        if 'decrease' in str(row['95% CI (TEXT)']):
            beta = beta * -1
        if needs_harmonising:
            beta = beta * -1
        estimate = beta
        se = estimate / z
        cis = [estimate - 1.96 * se,
               estimate + 1.96 * se]
        direction = '+' if beta >=0 else '-'
        beta_ci_lower = min(cis)
        beta_ci_upper = max(cis)
    else:
        oddsr = float(row['OR or BETA'])
        if needs_harmonising:
            oddsr = oddsr ** -1
        estimate = np.log(oddsr)
        se = estimate / z
        cis = [np.exp(estimate - 1.96 * se),
               np.exp(estimate + 1.96 * se)]
        direction = '+' if oddsr >=1 else '-'
        oddsr_ci_lower = min(cis)
        oddsr_ci_upper = max(cis)

    return direction, beta, beta_ci_lower, beta_ci_upper, oddsr, oddsr_ci_lower, oddsr_ci_upper

def extract_risk_allele(s):
    ''' Takes a string from GWAS Catalog STRONGEST SNP-RISK ALLELE field and
        returns the risk allele, if any.
    Args:
        row (pd.series)
    Returns:
        str
    '''
    mtch = re.match(r'.*-([A|T|G|C]+)', s)
    if mtch:
        allele = mtch.group(1)
    else:
        allele = None
    return allele

def comp(s):
    """ Complementary sequence """
    com_dict = {"A":"T", "T":"A", "G":"C", "C":"G"}
    com_seq = "".join([com_dict.get(nucl, "N") for nucl in s])
    return com_seq

def revcomp(s):
    """ Reverse complement sequence """
    return comp(s[::-1])

def parse_pval_mantissa_exponent(s):
    ''' Takes the p-value str and converts to mantissa and exponent
    '''
    mantissa, exponent = s.split('E')
    return float(mantissa), float(exponent)

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--inf', metavar="<file>", help=('GWAS Catalog input'), type=str, required=True)
    parser.add_argument('--id_lut', metavar="<file>", help=('LUT for new study IDs'), type=str, required=True)
    parser.add_argument('--outf', metavar="<str>", help=("Output"), type=str, required=True)
    parser.add_argument('--log', metavar="<str>", help=("Log output"), type=str, required=True)
    args = parser.parse_args()
    return args

def explode(df, columns):
    ''' Explodes multiple columns
    '''
    idx = np.repeat(df.index, df[columns[0]].str.len())
    a = df.T.reindex(columns).values
    concat = np.concatenate([np.concatenate(a[i]) for i in range(a.shape[0])])
    p = pd.DataFrame(concat.reshape(a.shape[0], -1).T, idx, columns)
    return pd.concat([df.drop(columns, axis=1), p], axis=1).reset_index(drop=True)

if __name__ == '__main__':

    main()
