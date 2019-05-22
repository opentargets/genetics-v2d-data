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
from concurrent.futures import ProcessPoolExecutor

def main():

    # Parse args
    global args
    args = parse_args()

    # Make output dir
    os.makedirs(args.outdir, exist_ok=True)

    # Load a list of variant IDs
    varids = set([])
    with open(args.varfile, 'r') as in_h:
        for line in in_h:
            varids.add(line.rstrip())

    # Run all in parallel
    with ProcessPoolExecutor(max_workers=args.max_cores) as executor:
        executor.map(run_single_variant, varids)


def run_single_variant(varid):
    ''' Runs the LD pipeline for a single variant
    '''

    # Calculate R for each population
    results = []
    for pop in args.pops:

        # Make command variables
        inbfile = args.bfile.replace('POPULATION', pop).replace(
            'CHROM', varid.split(':')[0])
        outtemp = os.path.join(
            args.outdir,
            varid.replace(':', '_') + '.plink'
        )

        # Calc LD
        res = calc_ld(varid,
                      inbfile,
                      pop,
                      args.ld_window,
                      outtemp)
        results.append(res)

    # Merge results together
    merged = reduce(lambda left, right: pd.merge(
        left, right, how='outer'), results)
    merged = merged.loc[:, ['index_variant_id', 'tag_variant_id',
                            'R_AFR', 'R_AMR', 'R_EAS', 'R_EUR', 'R_SAS']]

    # Caluclate R2 and remove rows where the max R2 < min_r2
    # This is needed to save storage space in the output files
    r2 = merged.loc[:, ['R_AFR', 'R_AMR', 'R_EAS', 'R_EUR', 'R_SAS']] ** 2
    to_keep = (r2.max(axis=1) >= args.min_r2)
    merged = merged.loc[to_keep, :]

    # DEBUG, i'm expected this to fail but don't know what the error will be
    # assert(merged.shape[0] > 0)

    # Save
    outf = os.path.join(
        args.outdir,
        varid.replace(':', '_') + '.ld.tsv.gz'
    )
    merged.to_csv(outf, sep='\t', index=None, compression='gzip')

    return 0

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
        '--ld-window', 99999999,
        '--r', 'gz',
        '--memory', 1000,
        '--threads', 1,
        '--allow-extra-chr',
        '--out', outf
    ]
    cmd_str = ' '.join([str(x) for x in cmd])
    # print(cmd_str)

    # Run command
    FNULL = open(os.devnull, 'w')
    sp.call(cmd_str, shell=True, stdout=FNULL, stderr=sp.STDOUT)

    # Load result to df
    try:
        res_file = outf + '.ld.gz'
        res = pd.read_csv(res_file, header=0, sep=r"\s+", engine='python')
        res = res.loc[:, ['SNP_A', 'SNP_B', 'R']]
        res.columns = ['index_variant_id', 'tag_variant_id', 'R_{}'.format(pop)]
    except FileNotFoundError:
        '''
        TODO. I don't want this to pick up ANY missing. It should only create a
        file if the following error is in the plink log file:
            Error: No valid variants specified by --ld-snp/--ld-snps/--ld-snp-list.
        '''
        res = pd.DataFrame(columns=['index_variant_id', 'tag_variant_id', 'R_{}'.format(pop)])
    
    # Delete temp
    if args.delete_temp:
        for tempfile in [res_file, outf + '.nosex']:
            if os.path.exists(tempfile):
                os.remove(tempfile)

    return res

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--varfile', metavar="<str>", help=("Input file specifying a variant on each line"), type=str, required=True)
    parser.add_argument('--bfile', metavar="<str>", help=("Input plink file pattern"), type=str, required=True)
    parser.add_argument('--pops', metavar="<str>", help=("Populations"), nargs='+', type=str, required=True)
    parser.add_argument('--ld_window', metavar="<int>", help=("Window to calc LD in (kb)"), type=int, required=True)
    parser.add_argument('--min_r2', metavar="<float>", help=("Minimum R2 to be kept"), type=float, required=True)
    parser.add_argument('--outdir', metavar="<str>", help=("Output directory"), type=str, required=True)
    parser.add_argument('--max_cores', metavar="<int>", help=("Maximum cores to use"), type=int, default=os.cpu_count())
    parser.add_argument('--delete_temp', help=("Remove temporary files"), action='store_true')
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
