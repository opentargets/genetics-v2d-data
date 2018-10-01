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

def main():

    # Parse args
    args = parse_args()

    # Load top loci table
    loci = pd.read_csv(args.in_loci, sep='\t', header=0,
                       usecols=['study_id', 'variant_id_b37', 'rsid'])

    # Explode lines
    loci['variant_id_b37'] = loci['variant_id_b37'].astype(str).str.split(';')
    loci['rsid'] = loci['rsid'].astype(str).str.split(';')
    loci = explode(loci, ['variant_id_b37', 'rsid'])

    # Split variant ID into parts
    loci[['chrom', 'pos', 'ref', 'alt']] = (
        loci.variant_id_b37.apply(lambda x: x.split('_'))
              .apply(pd.Series) )

    # Load study table
    studies = pd.read_csv(args.in_study, sep='\t', header=0,
                          usecols=['study_id', 'ancestry_initial', 'ancestry_replication'])

    # Merge
    merged = pd.merge(loci, studies, on='study_id', how='left')

    # Load population map
    pop_map = load_pop_map(args.in_popmap)

    # Get superpopulation proportions
    merged[['AFR_prop', 'AMR_prop', 'EAS_prop', 'EUR_prop', 'SAS_prop']] = (
        merged.apply(to_superpopulation_proportions, pop_map=pop_map, axis=1)
              .apply(pd.Series) )

    merged.to_csv(args.outf, sep='\t', index=None, compression='gzip')

def to_superpopulation_proportions(row, pop_map):
    ''' Parses GWAS Catalog ancestries,
        maps to 1000G superpopulations,
        returns the proportion of samples from each superpopulation,
    '''
    # Parse GWAS cat ancestries
    gwas_anc = {}
    for col in ['ancestry_initial', 'ancestry_replication']:
        if not pd.isnull(row[col]):
            for entry in re.split(r',|;', row[col].replace(' ', '')):
                try:
                    anc, n = entry.split('=')
                except ValueError:
                    continue
                try:
                    gwas_anc[anc] += int(n)
                except KeyError:
                    gwas_anc[anc] = int(n)
    # if row['study_id'] == 'GCST006208':
    # # if row['study_id'] == 'GCST004132':
    #     print(row)
    #     print(gwas_anc)
    #     sys.exit()
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
    df.gwascat_population = df.gwascat_population.str.replace(' ', '')
    pop_map = dict(zip(df['gwascat_population'], df['1000g_superpopulation']))
    return pop_map

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
    parser.add_argument('--in_loci', metavar="<str>", type=str, required=True)
    parser.add_argument('--in_study', metavar="<str>", type=str, required=True)
    parser.add_argument('--in_popmap', metavar="<str>", type=str, required=True)
    parser.add_argument('--outf', metavar="<str>", help=("Output"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
