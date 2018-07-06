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

def main():

    # Parse args
    args = parse_args()

    # Load
    credset = pd.read_csv(args.inf, sep='\t', header=0)

    # Make study id
    credset.loc[:, 'study_id'] = 'NEALEUKB_' + credset['trait'].astype(str)

    # Filter
    credset = credset.loc[credset.is95_credset == 1, :]

    # Format
    cols = OrderedDict([
        ('study_id', 'study_id'),
        ('locus_index_varid', 'index_variantid_b37'),
        ('varid', 'tag_variantid_b37'),
        ('logABF', 'log10_ABF'),
        ('postprob', 'posterior_prob')
        ])
    credset = ( credset.loc[:, list(cols.keys())]
                       .rename(columns=cols) )

    # Save
    credset.to_csv(args.outf, sep='\t', index=None, compression='gzip')

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--inf', metavar="<file>", help=('GWAS Catalog input'), type=str, required=True)
    parser.add_argument('--outf', metavar="<str>", help=("Output"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
