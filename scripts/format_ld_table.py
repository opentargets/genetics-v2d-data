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
from parquet_writer import write_parquet

def main():

    # Parse args
    args = parse_args()

    # Load
    ld = pd.read_csv(args.inf, sep='\t', header=0, low_memory=False)

    # Create table which tells us if LD is available for a given (study, lead)
    ld_available = ld.loc[:, ['study_id', 'index_variant_id']].drop_duplicates()
    ld_available['ld_available'] = True

    # Filter
    ld = ld.loc[ld.R2_overall >= args.min_r2, :]

    # For each lead variant, add a tag variant with R-squared == 1.0
    # Need to do this using the manifest, as leads without LD have already been
    # filtered out by this point
    manifest = pd.read_csv(args.in_manifest, sep='\t', header=0)
    manifest['index_variant_id'] = manifest['variant_id']
    manifest['tag_variant_id'] = manifest['variant_id']
    manifest['R2_overall'] = 1.0
    manifest = manifest.drop(['chrom', 'pos', 'ref', 'alt', 'variant_id'], axis=1)
    ld = pd.concat([ld, manifest], sort=False)

    # Sort by overall R2 and then deduplicate on (study, lead, tag)
    ld = (
        ld.sort_values('R2_overall', ascending=False)
        .drop_duplicates(subset=['study_id', 'index_variant_id', 'tag_variant_id'])
    )

    # Add column showing whether LD was available for each (study, lead)
    ld = pd.merge(ld,
                  ld_available,
                  on=['study_id', 'index_variant_id'],
                  how='left')
    ld['ld_available'] = ld['ld_available'].fillna(False)

    # Decompose variant IDs
    ld[['lead_chrom', 'lead_pos', 'lead_ref', 'lead_alt']] = \
        ld.index_variant_id.str.split('_', 3, expand=True)
    ld[['tag_chrom', 'tag_pos', 'tag_ref', 'tag_alt']] = \
        ld.tag_variant_id.str.split('_', 3, expand=True)
    ld['lead_pos'] = ld['lead_pos'].astype(int)
    ld['tag_pos'] = ld['tag_pos'].astype(int)

    # Format
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
        ('R2_overall', 'overall_r2'),
        ('AFR_prop', 'AFR_1000G_prop'),
        ('AMR_prop', 'AMR_1000G_prop'),
        ('EAS_prop', 'EAS_1000G_prop'),
        ('EUR_prop', 'EUR_1000G_prop'),
        ('SAS_prop', 'SAS_1000G_prop'),
        ('ld_available', 'ld_available')
    ])
    ld = ( ld.loc[:, list(cols.keys())]
             .rename(columns=cols) )
    
    # Sort
    ld = ld.sort_values(
        ['study_id', 'lead_chrom', 'lead_pos', 'lead_ref', 'lead_alt',
         'tag_chrom', 'tag_pos', 'tag_ref', 'tag_alt']
    )

    # Save as parquet
    write_parquet(
        ld,
        args.outf,
        compression='snappy',
        flavor='spark'
    )

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--inf', metavar="<file>", help=('Study weighted LD file'), type=str, required=True)
    parser.add_argument('--in_manifest', metavar="<file>", help=('Input manifest file'), type=str, required=True)
    parser.add_argument('--outf', metavar="<str>", help=("Output"), type=str, required=True)
    parser.add_argument('--min_r2', metavar="<str>", help=("Minimum R2 to be included"), type=float, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
