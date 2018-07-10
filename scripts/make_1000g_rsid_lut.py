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

    with gzip.open(args.outf, 'w') as out_h:
        for inf in args.bims:
            with open(inf, 'r') as in_h:
                for line in in_h:
                    chrom, rsid, _, pos, a1, a2 = line.rstrip().split('\t')
                    varid = '{0}_{1}_{2}_{3}'.format(chrom, pos, a1, a2)
                    out_h.write(('\t'.join([rsid, varid]) + '\n').encode())

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--bims', metavar="<file>", type=str, nargs='+', required=True)
    parser.add_argument('--outf', metavar="<str>", type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
