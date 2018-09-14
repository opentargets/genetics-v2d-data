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

    # Load manifest
    manifest = pd.read_csv(args.in_manifest, sep='\t', header=0, dtype=object)

    # Load EFOs and merge multiple rows into one
    efos = ( pd.read_csv(args.in_efos, sep='\t', header=0)
               .dropna(how='any') )
    efos_collapsed = (
        efos.groupby('trait_code')
                .agg({'efo_code':combine_rows,
                      'efo_trait':combine_rows})
                .reset_index() )

    # Merge EFO codes
    manifest = pd.merge(manifest, efos_collapsed,
                        how='left',
                        left_on='Field.code', right_on='trait_code')


    # Make columns
    manifest.loc[:, 'study_id'] = 'NEALEUKB_' + manifest['Field.code'].astype(str)
    manifest.loc[:, 'trait_code'] = 'UKB_' + manifest['Field.code'].astype(str)
    manifest.loc[:, 'pmid'] = ''
    manifest.loc[:, 'pub_date'] = '2017-09-15'
    manifest.loc[:, 'pub_journal'] = ''
    manifest.loc[:, 'pub_title'] = ''
    manifest.loc[:, 'pub_author'] = 'Neale et al.'
    manifest.loc[:, 'n_initial'] = manifest['N.non.missing'].astype(int)
    manifest.loc[:, 'n_replication'] = 0
    manifest.loc[:, 'ancestry_initial'] = 'European=' + manifest['n_initial'].astype(str)
    manifest.loc[:, 'ancestry_replication'] = ''


    # Load categories
    categ = load_categories(args.in_categories)
    manifest = pd.merge(manifest, categ,
                        on='trait_code',
                        how='left')

    # Ouput required columns
    cols = OrderedDict([
        ('study_id', 'study_id'),
        ('pmid', 'pmid'),
        ('pub_date', 'pub_date'),
        ('pub_journal', 'pub_journal'),
        ('pub_title', 'pub_title'),
        ('pub_author', 'pub_author'),
        ('Field', 'trait_reported'),
        ('trait_code', 'trait_code'),
        # ('efo_trait', 'trait_mapped'),
        ('efo_code', 'trait_efos'),
        ('category', 'trait_category'),
        ('ancestry_initial', 'ancestry_initial'),
        ('ancestry_replication', 'ancestry_replication'),
        ('n_initial', 'n_initial'),
        ('n_replication', 'n_replication'),
        ('N.cases', 'n_cases')
        ])
    manifest = ( manifest.loc[:, list(cols.keys())]
                         .rename(columns=cols) )
    # print(manifest.columns)

    # Write
    manifest.to_csv(args.outf, sep='\t', index=None)

def combine_rows(items):
    return ';'.join(items)

def load_categories(inf):
    ''' Loads category annotations that will be used in PheWAS plot
    args:
        inf (str)
    returns:
        pd.Df
    '''
    df = (pd.read_csv(inf, sep='\t', header=0)
            .loc[:, ['phenocode', 'category']]
            .rename(columns={'phenocode':'trait_code'})
    )
    df.trait_code = 'UKB_' + df.trait_code
    return df


def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_manifest', metavar="<str>", help=("Input"), type=str, required=True)
    parser.add_argument('--in_efos', metavar="<str>", help=("Input"), type=str, required=True)
    parser.add_argument('--in_categories', metavar="<str>", help=("Input"), type=str, required=True)
    parser.add_argument('--outf', metavar="<str>", help=("Output"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
