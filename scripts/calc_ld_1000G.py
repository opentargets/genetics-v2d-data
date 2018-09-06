#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#
# Calculates LD in 1000 Genomes using plink
#

import sys
import os
import argparse
import pandas as pd
import subprocess as sp
from functools import reduce

def main():

    # Parse args
    args = parse_args()

    # Calculate R for each population
    results = []
    for pop in args.pops:

        # Make command variables
        plink_varid = args.varid.replace('_', ':')
        inbfile = args.bfile.replace('POPULATION', pop).replace('CHROM', args.varid.split('_')[0])
        outtemp = args.outf.replace('.ld.gz', '.plink')

        # Calc LD
        res = calc_ld(plink_varid,
                      inbfile,
                      pop,
                      args.ld_window,
                      outtemp)
        results.append(res)

    # Merge results together
    merged = reduce(lambda left, right: pd.merge(left, right), results)
    merged = merged.loc[:, ['index_variant_id_b37', 'tag_variant_id_b37', 'R_AFR', 'R_AMR', 'R_EAS', 'R_EUR', 'R_SAS']]

    # DEBUG, i'm expected this to fail but don't know what the error will be
    # assert(merged.shape[0] > 0)

    # Save
    merged.to_csv(args.outf, sep='\t', index=None, compression='gzip')


def calc_ld(varid, bfile, pop, ld_window, outf):
    ''' Uses plink to calc LD for a single variant
    Args:
        varid (str): variant ID as it appears in the plink bim
        bfile (str): plink file prefix
        pop (str): name of population
        ld_window (int): window are variant to calc LD for
        outf (file): location to save temp plink output
    Returns:
        pd.DataFrame
    '''
    # Make command
    cmd = [
        'plink',
        '--bfile', bfile,
        '--ld-snp', varid,
        '--ld-window-kb', ld_window,
        # '--ld-window-r2', 0,
        '--r',
        # '--r2',
        '--memory', 1000,
        '--threads', 1,
        '--out', outf
    ]
    cmd_str = ' '.join([str(x) for x in cmd])
    print(cmd_str)

    # Run command
    sp.call(cmd_str, shell=True)

    # Load result to df
    try:
        res_file = outf + '.ld'
        res = pd.read_table(res_file, header=0, sep=r"\s+", engine='python')
        res = res.loc[:, ['SNP_A', 'SNP_B', 'R']]
        res.columns = ['index_variant_id_b37', 'tag_variant_id_b37', 'R_{}'.format(pop)]
        # Replace : with _ in variant fields
        res.index_variant_id_b37 = res.index_variant_id_b37.str.replace(':', '_')
        res.tag_variant_id_b37 = res.tag_variant_id_b37.str.replace(':', '_')
    except FileNotFoundError:
        res = pd.DataFrame(columns=['index_variant_id_b37', 'tag_variant_id_b37', 'R_{}'.format(pop)])

    return res

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--varid', metavar="<str>", help=("Input variant"), type=str, required=True)
    parser.add_argument('--bfile', metavar="<str>", help=("Input plink file pattern"), type=str, required=True)
    parser.add_argument('--pops', metavar="<str>", help=("Populations"), nargs='+', type=str, required=True)
    parser.add_argument('--ld_window', metavar="<int>", help=("Window to calc LD in (kb)"), type=str, required=True)
    parser.add_argument('--outf', metavar="<str>", help=("Output"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
