#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
import argparse
import gzip

def main():

    # Parse args
    args = parse_args()

    # Load rsid and chr:pos sets from gwas catalog assocs
    var_set = parse_sets(args.gwas)

    # Parse vcf
    with gzip.open(args.vcf, 'r') as in_h:
        with gzip.open(args.out, 'w') as out_h:
            for line in in_h:
                lined = line.decode()
                # Skip headers
                if lined.startswith('#'):
                    continue
                chrom, pos, rsid, _, _, _, _, _ = lined.rstrip().split('\t')
                chrom_pos = '{0}:{1}'.format(chrom, pos)
                if chrom_pos in var_set or rsid in var_set:
                    out_h.write(line)

    return 0

def parse_sets(in_gwas):
    ''' Parses rsids from SNP_ID_CURRENT and chrom:pos from CHR_ID, CHR_POS
    Args:
        in_gwas (str): GWAS Catalog assoc file name
    Returns:
        set(rsids and chrom:pos)
    '''
    var_set = set([])
    with open(in_gwas, 'r') as in_h:
        header = in_h.readline().rstrip().split('\t')
        for line in in_h:

            parts = line.rstrip().split('\t')
            # print(parts)
            # print()

            # Add SNP_ID_CURRENT to set
            value = parts[header.index('SNP_ID_CURRENT')]
            if value:
                try:
                    rsid = 'rs{0}'.format(int(value))
                    var_set.add(rsid)
                except ValueError:
                    pass

            # Add SNPS to set
            value = parts[header.index('SNPS')]
            if value:
                for rsid in value.split('; '):
                    if rsid.startswith('rs'):
                        var_set.add(rsid)
                    elif rsid.startswith('chr'):
                        chrom, pos = rsid.replace('_', ':').replace('.', ':').replace('-', ':').split(':')[:2]
                        chrom = chrom.replace('chr', '')
                        pos = int(pos)
                        try:
                            chrom_pos = '{chrom}:{pos}'.format(chrom=chrom, pos=int(pos))
                            var_set.add(chrom_pos)
                        except ValueError:
                            pass

            # Add chrom:pos to set
            # chrom_value = parts[header.index('CHR_ID')]
            # pos_value = parts[header.index('CHR_POS')]
            # if chrom_value and pos_value:
            #     chroms = chrom_value.split(';')
            #     positions = pos_value.split(';')
            #     assert(len(chroms) == len(positions))
            #     for chrom, pos in zip(chroms, positions):
            #         try:
            #             chrom_pos = '{chrom}:{pos}'.format(chrom=chrom, pos=int(pos))
            #             var_set.add(chrom_pos)
            #         except ValueError:
            #             pass

    return var_set

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
    parser.add_argument('--gwas', metavar="<file>", help=('GWAS Catalog input'), type=str, required=True)
    parser.add_argument('--vcf', metavar="<file>", help=("VCF input"), type=str, required=True)
    parser.add_argument('--out', metavar="<str>", help=("Output"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
