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
import hashlib

def main():

    # Parse args
    args = parse_args()


    # Load ancestry info
    anc = parse_ancestry_info(args.in_ancestries)

    # Load study table
    studies = pd.read_csv(args.in_studies, sep='\t', header=0)

    # Merge
    merged = pd.merge(studies, anc,
                      on='STUDY ACCESSION',
                      how='left')

    merged['trait_code'] = [make_hash_trait(val) for val in merged['DISEASE/TRAIT']]

    # Keep required cols
    cols = OrderedDict([
        ['STUDY ACCESSION', 'study_id'],
        ['PUBMEDID', 'pmid'],
        ['DATE', 'pub_date'],
        ['JOURNAL', 'pub_journal'],
        ['STUDY', 'pub_title'],
        ['FIRST AUTHOR', 'pub_author'],
        ['DISEASE/TRAIT', 'trait_reported'],
        ['MAPPED_TRAIT_URI', 'trait_efos'],
        ['trait_code', 'trait_code'],
        # ['MAPPED_TRAIT', 'trait_mapped'],
        ['ancestry_initial', 'ancestry_initial'],
        ['ancestry_replication', 'ancestry_replication'],
        ['n_initial', 'n_initial'],
        ['n_replication', 'n_replication']
    ])
    df = merged.rename(columns=cols).loc[:, cols.values()]
    #
    df['n_cases'] = ''

    # Clean efo codes
    df['trait_efos'] = df['trait_efos'].apply(clean_efo)
    # df['trait_mapped'] = df['trait_mapped'].str.replace(', ', ';')

    # Remove rows where n_initial == nan, these are old studies and their data is
    # inconsistent with the rest of GWAS Catalog (correspondance w Annalisa Buniello)
    df = df.loc[~pd.isnull(df.n_initial), :]

    # Write
    df.to_csv(args.outf, sep='\t', index=None)


def make_hash_trait(traitlabel):
    '''creates an id for each GWAS trait with
    low probability of collision
    '''
    hasher = hashlib.md5(traitlabel.encode('utf-8').lower())
    # truncating it in half gives 64 bit, which for ~10000 traits makes the collision probability < 1e-7
    return 'GT_' + hasher.hexdigest()[0:16]


def clean_efo(efo_str):
    ''' Extracts efo codes from the URI strings provided
    '''
    try:
        codes = []
        for uri in efo_str.split(', '):
            code = uri.split('/')[-1]
            codes.append(code)
        return ';'.join(codes)
    except AttributeError:
        return ''

def parse_ancestry_info(inf):
    ''' Parse the GWAS Catalog ancestry info file to get 1 row per study
    Args:
        inf (str): GWAS Catalog ancestry file
    Returns:
        pandas df
    '''
    anc_dict = {}
    df = pd.read_csv(inf, sep='\t', header=0, index_col=False)
    grouped = df.groupby(['STUDY ACCCESSION', 'STAGE'])
    for name, group in grouped:
        # Add study id to output dict
        study_id, stage = name
        if not study_id in anc_dict:
            anc_dict[study_id] = {}
        # Make N
        key = 'n_{0}'.format(stage)
        value = int(group['NUMBER OF INDIVDUALS'].sum())
        anc_dict[study_id][key] = str(value)
        # Make ancestry string
        key = 'ancestry_{0}'.format(stage)
        values = []
        for cat, n in zip(group['BROAD ANCESTRAL CATEGORY'], group['NUMBER OF INDIVDUALS']):
            try:
                n = int(n)
            except ValueError:
                n = 0
            values.append('{0}={1}'.format(cat, n))
        value = ';'.join(values)
        anc_dict[study_id][key] = value

    # Make df
    anc_df = pd.DataFrame.from_dict(anc_dict, orient='index')
    anc_df.insert(0, 'STUDY ACCESSION', anc_df.index)

    return anc_df

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_studies', metavar="<str>", type=str, required=True)
    parser.add_argument('--in_ancestries', metavar="<str>", type=str, required=True)
    parser.add_argument('--outf', metavar="<str>", help=("Output"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
