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
import re
from parquet_writer import write_parquet

def main():

    # Parse args
    args = parse_args()

    # Load top loci table
    loci = pd.read_parquet(args.in_loci,
                           engine='pyarrow',
                           columns=['study_id', 'chrom', 'pos', 'ref', 'alt'])

    # Load study table
    studies = pd.read_parquet(args.in_study,
                           engine='pyarrow',
                           columns=['study_id', 'ancestry_initial', 'ancestry_replication'])

    # Convert numpy array to list
    for col in ['ancestry_initial', 'ancestry_replication']:
        studies[col] = studies[col].apply(numpya_to_list)

    # Merge
    merged = pd.merge(loci, studies, on='study_id', how='inner')

    # Load population map
    pop_map = load_pop_map(args.in_popmap)

    # Get superpopulation proportions
    merged[['AFR_prop', 'AMR_prop', 'EAS_prop', 'EUR_prop', 'SAS_prop']] = (
        merged.apply(to_superpopulation_proportions, pop_map=pop_map, axis=1)
              .apply(pd.Series) )
    
    # Make variant ID
    merged['variant_id'] = (
        merged.loc[:, ['chrom', 'pos', 'ref', 'alt']]
        .apply(lambda row: '_'.join([str(x) for x in row]), axis=1)
    )

    # Write output to parquet
    out_cols = ['study_id', 'variant_id', 'chrom', 'pos', 'ref', 'alt',
                'AFR_prop', 'AMR_prop', 'EAS_prop', 'EUR_prop', 'SAS_prop']
    # write_parquet(merged.loc[:, out_cols],
    #               args.outf,
    #               compression='snappy',
    #               flavor='spark')

    # Write csv
    merged.loc[:, out_cols].to_csv(
        args.outf, sep='\t', index=None, compression='gzip')

def to_superpopulation_proportions(row, pop_map, anc_sep=', '):
    ''' Parses GWAS Catalog ancestries,
        maps to 1000G superpopulations,
        returns the proportion of samples from each superpopulation,
    params:
        anc_sep (str): separator used to split ancesties by GWAS Catalog. Warning, this is due to change.

    '''

    # Parse GWAS cat ancestries
    gwas_anc = {}
    for col in ['ancestry_initial', 'ancestry_replication']:
        for entry in row[col]:
            
            # DEBUG
            try:
                anc, n = entry.split('=')
            except:
                print('col:', col)
                print('row:', row)
                print('entry:', entry)
                print('str(entry):', str(entry))
                print('type(entry):', type(entry))
                print('entry == "":', entry == '')
                if not entry:
                    print('skip')

            # Extract ancestry and n
            anc, n = entry.split('=')

            # Split compounded ancestries into separate parts
            sub_anc_parts = anc.split(', ')
            for sub_anc in sub_anc_parts:

                # Split sample size between them
                sub_n = int(float(n) / len(sub_anc_parts))

                # Add to dict
                try:
                    gwas_anc[sub_anc] += int(sub_n)
                except KeyError:
                    gwas_anc[sub_anc] = int(sub_n)

    # Map to superpopulations
    superpop = {}
    for anc in gwas_anc:
        if anc in pop_map:
            superpop[pop_map[anc]] = gwas_anc[anc]

    # Make output row
    out_row = []
    total_n = sum(superpop.values())
    for pop in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
        if total_n > 0:
            out_row.append(float(superpop.get(pop, 0)) / total_n)
        else:
            out_row.append(np.nan)

    return out_row

def load_pop_map(inf):
    ''' Load dictionary to map from GWAS ancestry to 1000G superpopulation
    Args:
        inf (str): file name
    Returns:
        dict(anc -> super pop)
    '''
    df = pd.read_csv(inf, sep='\t', header=0).dropna()
    pop_map = dict(zip(df['gwascat_population'], df['1000g_superpopulation']))
    return pop_map

def numpya_to_list(numpy_array):
    ''' Converts a numpy array to a list
    '''
    try:
        l = numpy_array.tolist()
        # If list contains only an empty string
        if not l[0]:
            return []      
        # Else return list
        return l
    except AttributeError:
        return []

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_loci', metavar="<str>", type=str, required=True)
    parser.add_argument('--in_study', metavar="<str>", type=str, required=True)
    parser.add_argument('--in_popmap', metavar="<str>", type=str, required=True)
    parser.add_argument('--outf', metavar="<str>", help=("Output"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
