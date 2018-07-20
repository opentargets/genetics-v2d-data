#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
import argparse
import gzip
import pandas as pd

def main():

    # Parse args
    args = parse_args()

    # Load bims
    print('Loading bims...')
    bim_dfs = []
    for bim in args.bims:
        df = pd.read_csv(bim, sep='\t', header=None)
        df.columns = ['chrom_1kg', 'rsid', 'None', 'pos_1kg', 'A1_1kg', 'A2_1kg']
        bim_dfs.append(df)
    bim_df = pd.concat(bim_dfs)
    print(' {0} 1000G variants'.format(bim_df.shape[0]))

    # Make 1000G variant merge key
    bim_df['key_1kg'] = bim_df.loc[:, ['chrom_1kg', 'pos_1kg', 'A1_1kg', 'A2_1kg']].apply(make_variant_key, axis=1)

    # Make variant ID
    bim_df['varid_1kg'] = bim_df.loc[:, ['chrom_1kg', 'pos_1kg', 'A1_1kg', 'A2_1kg']].apply(lambda row: '_'.join([str(x) for x in row]), axis=1)

    # Make set of chrom:pos
    loc_set = set(bim_df.loc[:, ['chrom_1kg', 'pos_1kg']].apply(lambda x: '_'.join([str(x.iloc[0]), str(x.iloc[1])]), axis=1).tolist())

    # Extract required rows from VCF
    print('Loading VCF...')
    vcf = parse_vcf(args.vcf, loc_set)

    # Make 1000G variant merge key
    vcf['key_ens'] = vcf.loc[:, ['chrom_ens', 'pos_ens', 'A1_ens', 'A2_ens']].apply(make_variant_key, axis=1)

    # Make variant ID
    vcf['varid_ens'] = vcf.loc[:, ['chrom_ens', 'pos_ens', 'A1_ens', 'A2_ens']].apply(lambda row: '_'.join([str(x) for x in row]), axis=1)

    # Merge
    print('Merging...')
    merged = pd.merge(bim_df, vcf, left_on='key_1kg', right_on='key_ens', how='inner')

    # Out
    merged.loc[:, ['rsid', 'varid_1kg', 'varid_ens']].to_csv(args.outf, sep='\t', index=None, compression='gzip')

def parse_vcf(inf, loc_set):
    ''' Parses the reference VCF to extract only required rows.
    Args:
        inf (str): filename
        loc_set (set): set of chrom:pos strings to extract
    Returns:
        pd.Df
    '''
    rows = []
    c = 0
    with gzip.open(inf, 'r') as in_h:
        for line in in_h:
            line = line.decode()
            # Skip header
            if line.startswith('#'):
                continue
            # Get info
            parts = line.rstrip().split('\t')
            chrom = parts[0]
            pos = parts[1]
            a1 = parts[3]
            a2_multi = parts[4]
            if '{0}_{1}'.format(chrom, pos) in loc_set:
                # Split multiallelic into separate
                for a2 in a2_multi.split(','):
                    rows.append([chrom, pos, a1, a2])
                    if c % 100000 == 0:
                        print(' {0} variants extracted from vcf'.format(c))
                    c += 1
    # Make df
    df = pd.DataFrame(rows, columns=['chrom_ens', 'pos_ens', 'A1_ens', 'A2_ens'])
    return df

def make_variant_key(row):
    ''' Returns a key to be used to merge on variant. To do this, alleles are
        sorted alphabetically.
    '''
    key = '_'.join(str(x) for x in
        [row.iloc[0], row.iloc[1]] + list(sorted([row.iloc[2], row.iloc[3]]))
        )
    return key

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--bims', metavar="<file>", type=str, nargs='+', required=True)
    parser.add_argument('--vcf', metavar="<file>", type=str, required=True)
    parser.add_argument('--outf', metavar="<str>", type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
