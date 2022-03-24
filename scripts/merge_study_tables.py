"""
Merges studies from different sources into one
"""
#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import argparse

import pandas as pd


def main(in_gwascat: str, in_ukb: str, in_finngen: str, output: str) -> None:

    # Load input files
    gwas = pd.read_json(in_gwascat, orient='records', lines=True)
    ukb = pd.read_json(in_ukb, orient='records', lines=True)
    finngen = pd.read_json(in_finngen, orient='records', lines=True)

    # Merge
    merged = pd.concat([gwas, ukb, finngen], sort=False)

    # Write
    merged.to_json(output, orient='records', lines=True)



def parse_args():
    """
    Load command line args.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_gwascat', metavar="<str>", type=str, required=True)
    parser.add_argument('--in_ukb', metavar="<str>", type=str, required=True)
    parser.add_argument('--in_finngen', metavar="<str>", type=str, required=True)
    parser.add_argument('--output', metavar="<str>", help=("Output merged file"), type=str, required=True)
    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_args()
    main(in_gwascat=args.in_gwascat, in_ukb=args.in_ukb, in_finngen=args.in_finngen, output=args.output)
