#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
import pandas as pd
from pprint import pprint
from collections import OrderedDict
from parquet_writer import write_parquet

def main():

    # Parse args
    inf = 'ld.tsv.gz'
    instudies = '../../output/190219/studies.parquet'
    outf = 'ld.parquet'

    # Load
    df = pd.read_csv(inf, sep='\t', header=0)
    studies = pd.read_parquet(instudies)

    # Map new study ids to old study ids
    studies['old_study_id'] = studies['study_id'].apply(new_to_old)
    studies = studies.loc[:, ['study_id', 'old_study_id']]
    df = (
        df.rename(columns={'study_id': 'old_study_id'})
          .merge(studies, on='old_study_id')
          .drop('old_study_id', axis=1)
    )
    
    # Decompose variant IDs
    df[['lead_chrom', 'lead_pos', 'lead_ref', 'lead_alt']] = \
        df.index_variantid_b37.str.split('_', 3, expand=True)
    df[['tag_chrom', 'tag_pos', 'tag_ref', 'tag_alt']] = \
        df.tag_variantid_b37.str.split('_', 3, expand=True)
    df['lead_pos'] = df['lead_pos'].astype(int)
    df['tag_pos'] = df['tag_pos'].astype(int)

    # Rename and select columns
    cols = OrderedDict([
        ('study_id', 'study_id'),
        ('lead_chrom', 'lead_chrom'),
        ('lead_pos', 'lead_pos'),
        ('lead_ref', 'lead_ref'),
        ('lead_alt', 'lead_alt'),
        ('tag_chrom', 'tag_chrom'),
        ('tag_pos', 'tag_pos'),
        ('tag_ref', 'tag_ref'),
        ('tag_alt', 'tag_alt'),
        ('overall_r2', 'overall_r2'),
        ('AFR_1000G_prop', 'AFR_1000G_prop'),
        ('AMR_1000G_prop', 'AMR_1000G_prop'),
        ('EAS_1000G_prop', 'EAS_1000G_prop'),
        ('EUR_1000G_prop', 'EUR_1000G_prop'),
        ('SAS_1000G_prop', 'SAS_1000G_prop')
        ])
    df = ( df.loc[:, list(cols.keys())]
                       .rename(columns=cols) )

    # Coerce data types
    dtypes = OrderedDict([
        ('study_id', 'object'),
        ('lead_chrom', 'object'),
        ('lead_pos', 'Int64'),
        ('lead_ref', 'object'),
        ('lead_alt', 'object'),
        ('tag_chrom', 'object'),
        ('tag_pos', 'Int64'),
        ('tag_ref', 'object'),
        ('tag_alt', 'object'),
        ('overall_r2', 'float64'),
        ('AFR_1000G_prop', 'float64'),
        ('AMR_1000G_prop', 'float64'),
        ('EAS_1000G_prop', 'float64'),
        ('EUR_1000G_prop', 'float64'),
        ('SAS_1000G_prop', 'float64')
        ])
    assert(set(dtypes.keys()) == set(df.columns))
    df = (
        df.loc[:, dtypes.keys()]
        .astype(dtype=dtypes)
    )

    # Sort
    df = df.sort_values(
        ['study_id', 'lead_chrom', 'lead_pos', 'lead_ref', 'lead_alt',
         'tag_chrom', 'tag_pos', 'tag_ref', 'tag_alt']
    )

    # Save as parquet
    write_parquet(df,
                  outf,
                  compression='snappy',
                  flavor='spark')

    # Save csv
    # df.to_csv(args.outf, sep='\t', index=None, compression='gzip')

def new_to_old(study_id):
    ''' Returns the old study id given the new
    '''
    if study_id.startswith('GCST'):
        return study_id.split('_')[0]
    else:
        return study_id

if __name__ == '__main__':

    main()
