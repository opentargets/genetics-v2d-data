"""
Merges studies from different sources into one
"""
#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import argparse
import logging

import pandas as pd


def main(in_gwascat: str, in_ukb: str, in_finngen: str, output_path: str) -> None:

    logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.INFO)

    # Load input files
    gwas = pd.read_json(in_gwascat, orient='records', lines=True).drop_duplicates()
    ukb = pd.read_json(in_ukb, orient='records', lines=True).drop_duplicates()
    finngen = pd.read_json(in_finngen, orient='records', lines=True).drop_duplicates()

    logging.info(f"{len(gwas)} studies from GWAS Catalog have been loaded. Formatting...")
    logging.info(f"{len(ukb)} studies from UK Biobank have been loaded. Formatting...")
    logging.info(f"{len(finngen)} studies from Finngen have been loaded. Formatting...")
    

    # Merge
    merged = pd.concat([gwas, ukb, finngen], sort=False).drop_duplicates()

    assert gwas.shape[0] + ukb.shape[0] + finngen.shape[0] == merged.shape[0], 'Merged table has different number of rows'

    # Write
    merged.to_json(output_path, orient='records', lines=True)
    logging.info(f"{len(merged)} studies have been saved in {output_path}. Exiting.")



def parse_args():
    """
    Load command line args.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_gwascat', metavar="<str>", type=str, required=True)
    parser.add_argument('--in_ukb', metavar="<str>", type=str, required=True)
    parser.add_argument('--in_finngen', metavar="<str>", type=str, required=True)
    parser.add_argument('--output', metavar="<str>", help=("Output merged file in parquet format"), type=str, required=True)
    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_args()
    main(in_gwascat=args.in_gwascat, in_ukb=args.in_ukb, in_finngen=args.in_finngen, output=args.output)
