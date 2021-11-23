#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#
# Merge top loci table
#

import sys
import argparse
import pandas as pd
from collections import OrderedDict
from parquet_writer import write_parquet
import datetime

def main():

    # Parse args
    args = parse_args()

    # Load
    gwas = pd.read_json(args.in_gwascat, orient='records', lines=True)
    ukb = pd.read_json(args.in_ukb, orient='records', lines=True)

    # Merge
    merged = pd.concat([gwas, ukb], sort=False)

    # Write
    merged.to_json(args.output, orient='records', lines=True)

    return 0

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_gwascat', metavar="<str>", type=str, required=True)
    parser.add_argument('--in_ukb', metavar="<str>", type=str, required=True)
    parser.add_argument('--output', metavar="<str>", help=("Output merged file"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
