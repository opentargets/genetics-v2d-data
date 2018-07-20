#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
import argparse
from pprint import pprint
import pandas as pd
import gzip

def main():

    # Parse args
    args = parse_args()

    # Load postgap r2
    print('Loading postgap...')
    ld = pd.read_csv(args.in_postgap,
                     sep='\t',
                     usecols=['gwas_snp', 'ld_snp_rsID', 'r2'])
    ld = ld.loc[:, ['gwas_snp', 'ld_snp_rsID', 'r2']]
    ld.columns = ['rsid1', 'rsid2', 'r2']
    ld = ld.drop_duplicates()

    # Load rsid to variant id map
    print('Loading var map...')
    rsids = set(ld.rsid1.tolist() + ld.rsid2.tolist())
    var_map = load_map(args.in_lut, rsids)
    print('Mapping variant ids...')
    ld['varid_1_b37'] = ld.rsid1.apply(lambda x: var_map.get(x, None))
    ld['varid_2_b37'] = ld.rsid2.apply(lambda x: var_map.get(x, None))

    # Output
    print('Saving...')
    ld = ld.loc[:, ['varid_1_b37', 'rsid1', 'varid_2_b37', 'rsid2', 'r2']]
    ld.to_csv(args.outf, sep='\t', index=None, compression='gzip')

def load_map(inf, rsids):
    ''' Loads rsid to variant id map
    Args:
        inf (str): File containing lut
        rsids (set): rsids to extract
    Returns:
        dict(rsid -> variantid)
    '''
    d = {}
    with gzip.open(inf, 'r') as in_h:
        for line in in_h:
            rsid, varid_1kg, varid_ens = line.decode().rstrip().split('\t')
            if rsid in rsids:
                d[rsid] = varid_ens
    return d

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_postgap', metavar="<file>", type=str, required=True)
    parser.add_argument('--in_lut', metavar="<file>", type=str, required=True)
    parser.add_argument('--outf', metavar="<str>", type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
