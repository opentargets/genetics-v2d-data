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
    allowed_chroms = set(
        [str(x) for x  in range(1, 23)] + ['X', 'Y', 'MT']
    )

    # Load rsid and chr:pos sets from gwas catalog assocs
    rsid_set, chrom_pos_b37_set, chrom_pos_b38_set = parse_sets(args.gwas)

    # Parse vcf
    with gzip.open(args.vcf, 'r') as in_h:
        with gzip.open(args.out, 'w') as out_h:
            in_h.readline() # Skip header
            for line in in_h:
                
                parts = line.decode().rstrip().split('\t')
                
                # Extract b37 chrom:pos
                chrom_b37 = parts[2]
                pos_b37 = parts[3]
                chrom_pos_b37 = '{0}:{1}'.format(chrom_b37, pos_b37)
                # Extract b38 chrom:pos
                chrom_b38 = parts[4]
                pos_b38 = parts[5]
                chrom_pos_b38 = '{0}:{1}'.format(chrom_b38, pos_b38)
                # Extract rsid
                rsid = parts[8]

                # Skip invalid chroms
                if ((chrom_b37 not in allowed_chroms) or
                    (chrom_b38 not in allowed_chroms)):
                    continue
                
                # Write line if variant is a gwas association
                if ((chrom_pos_b37 in chrom_pos_b37_set) or 
                    (chrom_pos_b38 in chrom_pos_b38_set) or
                    (rsid in rsid_set)):

                    out_h.write(line)

    return 0

def parse_sets(in_gwas):
    ''' Parses 
        1. rsids from SNP_ID_CURRENT, add rs prefix
        2. rsids from SNPS if startswith "rs"
        3. build37 chrom:pos from SNPS if startswith "chr"
        4. build38 chrom:pos from CHR_ID, CHR_POS
    Args:
        in_gwas (str): GWAS Catalog assoc file name
    Returns:
        set(rsids), set(chrom_pos_b37), set(chrom_pos_b38)
    '''
    rsid_set = set([])
    chrom_pos_b37_set = set([])
    chrom_pos_b38_set = set([])

    with open(in_gwas, 'r') as in_h:
        header = in_h.readline().rstrip().split('\t')
        for line in in_h:
            
            parts = line.rstrip().split('\t')

            # Add rsid from SNP_ID_CURRENT to rsid set
            value = parts[header.index('SNP_ID_CURRENT')]
            if value:
                try:
                    rsid = 'rs{0}'.format(int(value))
                    rsid_set.add(rsid)
                except ValueError:
                    pass

            # Add info from SNPS column
            value = parts[header.index('SNPS')]
            if value:
                for item in value.split('; '):
                    # Add rsid to rsid_set
                    if item.startswith('rs'):
                        rsid_set.add(item)
                    # Add chrom_pos to chrom_pos_b37_set
                    elif item.startswith('chr'):
                        chrom, pos = item.replace('_', ':').replace('.', ':').replace('-', ':').split(':')[:2]
                        chrom = chrom.replace('chr', '')
                        pos = int(pos)
                        try:
                            chrom_pos = '{chrom}:{pos}'.format(chrom=chrom, pos=int(pos))
                            chrom_pos_b37_set.add(chrom_pos)
                        except ValueError:
                            pass

            # Add chrom:pos from CHR_ID, CHR_POS to chrom_pos_b38_set
            chrom_value = parts[header.index('CHR_ID')]
            pos_value = parts[header.index('CHR_POS')]
            if chrom_value and pos_value:
                chroms = chrom_value.split(';')
                positions = pos_value.split(';')
                assert(len(chroms) == len(positions))
                for chrom, pos in zip(chroms, positions):
                    try:
                        chrom_pos = '{chrom}:{pos}'.format(chrom=chrom, pos=int(pos))
                        chrom_pos_b38_set.add(chrom_pos)
                    except ValueError:
                        pass

    return rsid_set, chrom_pos_b37_set, chrom_pos_b38_set

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
