#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
import argparse
import pandas as pd

def main():

    # Parse args
    args = parse_args()

    # Load
    df = pd.read_csv(args.inf, sep='\t', header=0)

    # Parse ancestries
    anc_list = df.ancestry_initial.dropna().tolist() + df.ancestry_replication.dropna().tolist()
    anc_clean = []
    for anc_str in anc_list:
        for entry in anc_str.split(';'):
            anc = entry.split('=')[0]
            anc_clean.append(anc)

    # Output
    out = pd.DataFrame(list(set(anc_clean)), columns=['gwascat_population'])
    out.to_csv(args.outf, sep='\t', index=None)


def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--inf', metavar="<str>", help=("Input"), type=str, required=True)
    parser.add_argument('--outf', metavar="<str>", help=("Output"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
