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
import logging

def main():

    # Parse args
    args = parse_args()

    # Start logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    handler = logging.FileHandler(args.log)
    logger.addHandler(handler)

    # Load assoc data
    gwas = pd.read_csv(args.inf, sep='\t', header=0, dtype={'P-VALUE':str})
    logger.info('Total associations: {0}'.format(gwas.shape[0]))
    # pprint(gwas.columns.tolist())

    # Drop rows with no variant ID
    gwas = gwas.dropna(subset=['variant_id_b37'])
    logger.info('Total with variant ID: {0}'.format(gwas.shape[0]))

    # Extract pvalue mantissa and exponent from p-value string
    gwas[['pval_mantissa', 'pval_exponent']] = gwas['P-VALUE'].apply(parse_pval_mantissa_exponent).apply(pd.Series)

    # Extract and rename required columns
    cols = OrderedDict([
        ('STUDY ACCESSION', 'study_id'),
        ('variant_id_b37', 'variant_id_b37'),
        ('rsid', 'rsid'),
        # ('P-VALUE (TEXT)', 'sub_phenotype'),
        ('pval_mantissa', 'pval_mantissa'),
        ('pval_exponent', 'pval_exponent')
    ])
    gwas = gwas[list(cols.keys())].rename(columns=cols)

    # Clean sub-phenotype column
    # gwas['sub_phenotype'] = gwas['sub_phenotype'].str.strip('()')

    # Drop duplicates by study, variant and sub_phenotype
    gwas = gwas.sort_values(['pval_exponent', 'pval_mantissa'])
    # gwas = gwas.drop_duplicates(subset=['study_id', 'variant_id_b37', 'sub_phenotype'], keep='first')
    gwas = gwas.drop_duplicates(subset=['study_id', 'variant_id_b37'], keep='first')
    logger.info('Total after duplicates removed: {0}'.format(gwas.shape[0]))

    gwas = gwas.sort_values(['study_id', 'variant_id_b37'])

    # Remove Sun et al pQTL study
    gwas = gwas.loc[gwas.study_id != 'GCST005806', :]
    # Remove Huang et al IBD study (GWAS Catalog should not have curated this)
    gwas = gwas.loc[gwas.study_id != 'GCST005837', :]

    gwas.to_csv(args.outf, sep='\t', index=None)

def parse_pval_mantissa_exponent(s):
    ''' Takes the p-value str and converts to mantissa and exponent
    '''
    mantissa, exponent = s.split('E')
    return float(mantissa), float(exponent)

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--inf', metavar="<file>", help=('GWAS Catalog input'), type=str, required=True)
    parser.add_argument('--outf', metavar="<str>", help=("Output"), type=str, required=True)
    parser.add_argument('--log', metavar="<str>", help=("Log output"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
