#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#
# Merge all LD files
#

import sys
import os
import argparse
import gzip

def main():

    # Parse args
    args = parse_args()
    header_done = False

    # Open output file
    with gzip.open(args.output, 'wt') as out_h:

        # Process each input file
        for inf in args.infiles:
            with gzip.open(inf, 'rt') as in_h:

                # Process header
                header = in_h.readline()
                if not header_done:
                    out_h.write(header)
                    header_done = True

                # Process remaining lines
                for line in in_h:
                    out_h.write(line)

    return 0

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--infiles', metavar="<str>", help=("Input files to merge"), nargs='+', type=str, required=True)
    parser.add_argument('--output', metavar="<str>", help=("Output merged file"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
