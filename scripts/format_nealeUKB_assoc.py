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

def main():

    # Parse args
    args = parse_args()

    # Load
    cred = pd.read_csv(args.inf, sep='\t', header=0)

    # Keep rows where index == var
    top_loci = cred.loc[cred.locus_index_varid == cred.varid, :]

    # Make study id
    top_loci.loc[:, 'study_id'] = 'NEALEUKB_' + top_loci['trait'].astype(str)

    # Fix very low p-values
    # top_loci.loc[top_loci.p == 0, 'p'] = sys.float_info.min
    top_loci.loc[top_loci.p == 0, 'p'] = (1 / sys.float_info.max)

    # Filter p-values
    top_loci = top_loci.loc[top_loci.p < 1e-5, :]

    # Get P mantissa and exponent
    top_loci.loc[:, 'p_mantissa'] = top_loci.p.apply(fman)
    top_loci.loc[:, 'p_exponent'] = top_loci.p.apply(fexp)

    # Extract and rename required columns
    cols = OrderedDict([
        ('study_id', 'study_id'),
        ('locus_index_varid', 'variant_id_b37'),
        ('locus_index_snp', 'rsid'),
        ('p_mantissa', 'pval_mantissa'),
        ('p_exponent', 'pval_exponent')
    ])
    top_loci = top_loci[list(cols.keys())].rename(columns=cols)

    # Sort and save
    top_loci = top_loci.sort_values(['study_id', 'pval_exponent', 'pval_mantissa'])
    top_loci.to_csv(args.outf, sep='\t', index=None)


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
    parser.add_argument('--outf', metavar="<str>", type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
