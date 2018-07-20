#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
import argparse
import pandas as pd
import numpy as np

def main():

    # Parse args
    args = parse_args()

    # Load top loci
    loci = pd.read_csv(args.in_loci, sep='\t', header=0)

    # Explode lines
    loci['variant_id_b37'] = loci['variant_id_b37'].astype(str).str.split(';')
    loci['rsid'] = loci['rsid'].astype(str).str.split(';')
    loci = explode(loci, ['variant_id_b37', 'rsid'])

    # Load ld
    ld = pd.read_csv(args.in_ld, sep='\t', header=0)

    # Merge
    merged = pd.merge(loci, ld,
                      how='left',
                      left_on='variant_id_b37',
                      right_on='varid_1_b37')

    # Keep required columns
    merged = merged.loc[:, ['study_id', 'variant_id_b37', 'varid_2_b37', 'r2']]
    merged = merged.drop_duplicates()
    merged =merged.rename(columns={'r2': 'overall_r2',
                                   'variant_id_b37': 'index_variantid_b37',
                                   'varid_2_b37': 'tag_variantid_b37'})

    # Drop rows with missing LD information
    merged = merged.dropna()

    # Add columns
    merged['AFR_1000G_prop'] = 0.0
    merged['AMR_1000G_prop'] = 0.0
    merged['EAS_1000G_prop'] = 0.0
    merged['EUR_1000G_prop'] = 1.0
    merged['SAS_1000G_prop'] = 0.0

    # Save
    merged.to_csv(args.outf, sep='\t', index=None, compression='gzip')

def explode(df, columns):
    ''' Explodes multiple columns
    '''
    idx = np.repeat(df.index, df[columns[0]].str.len())
    a = df.T.reindex_axis(columns).values
    concat = np.concatenate([np.concatenate(a[i]) for i in range(a.shape[0])])
    p = pd.DataFrame(concat.reshape(a.shape[0], -1).T, idx, columns)
    return pd.concat([df.drop(columns, axis=1), p], axis=1).reset_index(drop=True)

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_ld', metavar="<file>", type=str, required=True)
    parser.add_argument('--in_loci', metavar="<file>", type=str, required=True)
    parser.add_argument('--outf', metavar="<str>", type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
