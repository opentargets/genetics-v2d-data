#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import re
import sys
import os
import argparse
<<<<<<< HEAD
import logging
from collections import OrderedDict

import pandas as pd

=======
import pandas as pd
from pprint import pprint
from collections import OrderedDict
from operator import itemgetter
import json
>>>>>>> d58d0d4c6c7e6822fe566ee129890ba82045f2b7

def main(input_path: str, output_path: str) -> None:

    logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.INFO)

    # Only keep required cols
    to_keep = OrderedDict(
        [("code", "study_id"), ("n_total", "n_total"), ("n_cases", "n_cases")]
    )

    # Load manifest
    manifest = (
        pd.read_csv(input_path, sep="\t", header=0, dtype=object)
        .filter(items=to_keep)
        .rename(columns=to_keep)
    )

    logging.info(f"{input_path} has been loaded. Formatting...")

    #
    # Add other columns -------------------------------------------------------
    #

    # Vector to get Neale or SAIGE studies
    is_neale = manifest["study_id"].str.startswith("NEALE2_")
    is_saige = manifest["study_id"].str.startswith("SAIGE_")
    assert (is_neale | is_saige).all()

    # Make columns
<<<<<<< HEAD
    manifest.loc[:, "pmid"] = ""
    manifest.loc[is_neale, "pub_date"] = "2018-08-01"
    manifest.loc[is_saige, "pub_date"] = "2018-10-24"
    manifest.loc[:, "pub_journal"] = ""
    manifest.loc[:, "pub_title"] = ""
    manifest.loc[is_neale, "pub_author"] = "UKB Neale v2"
    manifest.loc[is_saige, "pub_author"] = "UKB SAIGE"
    manifest.loc[:, "n_initial"] = manifest["n_total"].apply(to_int_safe)
    manifest.loc[:, "n_cases"] = manifest["n_cases"].apply(to_int_safe)
    manifest.loc[:, "n_replication"] = 0
    manifest.loc[:, "ancestry_initial"] = "European=" + manifest["n_initial"].astype(
        str
    )
    manifest.loc[:, "ancestry_replication"] = ""

    # Ouput required columns
    cols = OrderedDict(
        [
            ("study_id", "study_id"),
            ("pmid", "pmid"),
            ("pub_date", "pub_date"),
            ("pub_journal", "pub_journal"),
            ("pub_title", "pub_title"),
            ("pub_author", "pub_author"),
            ("ancestry_initial", "ancestry_initial"),
            ("ancestry_replication", "ancestry_replication"),
            ("n_initial", "n_initial"),
            ("n_cases", "n_cases"),
            ("n_replication", "n_replication"),
        ]
    )
    manifest = manifest.loc[:, list(cols.keys())].rename(columns=cols)
=======
    manifest.loc[:, 'pmid'] = ''
    manifest.loc[is_neale, 'pub_date'] = '2018-08-01'
    manifest.loc[is_saige, 'pub_date'] = '2018-10-24'
    manifest.loc[:, 'pub_journal'] = ''
    manifest.loc[:, 'pub_title'] = ''
    manifest.loc[is_neale, 'pub_author'] = 'UKB Neale v2'
    manifest.loc[is_saige, 'pub_author'] = 'UKB SAIGE'
    manifest.loc[:, 'n_initial'] = manifest['n_total'].apply(to_int_safe)
    manifest.loc[:, 'n_cases'] = manifest['n_cases'].apply(to_int_safe)
    manifest.loc[:, 'n_replication'] = 0
    manifest.loc[:, 'ancestry_initial'] = 'European=' + manifest['n_initial'].astype(str)
    manifest.loc[:, 'ancestry_replication'] = ''

    # Load efos annotations
    efo_mapper = {}
    with open(args.in_efos, 'r') as in_h:
        for line in in_h:
            parts = json.loads(line)
            efo_mapper[parts['study_id']] = parts['efos']
    
    # Map efos
    manifest['trait_efos'] = manifest['study_id'].apply(
        lambda stid: efo_mapper.get(stid, None)
    )

    # Ouput required columns
    cols = OrderedDict([
        ('study_id', 'study_id'),
        ('pmid', 'pmid'),
        ('pub_date', 'pub_date'),
        ('pub_journal', 'pub_journal'),
        ('pub_title', 'pub_title'),
        ('pub_author', 'pub_author'),
        ('trait_reported', 'trait_reported'),
        ('trait_efos', 'trait_efos'),
        # ('trait_category', 'trait_category'),
        ('ancestry_initial', 'ancestry_initial'),
        ('ancestry_replication', 'ancestry_replication'),
        ('n_initial', 'n_initial'),
        ('n_cases', 'n_cases'),
        ('n_replication', 'n_replication')
        ])
    manifest = ( manifest.loc[:, list(cols.keys())]
                         .rename(columns=cols) )
>>>>>>> d58d0d4c6c7e6822fe566ee129890ba82045f2b7

    # Write
    manifest.to_json(args.output, orient="records", lines=True)
    logging.info(f"{len(manifest)} studies have been saved in {output_path}. Exiting.")


def to_int_safe(i):
    try:
        return int(float(i))
    except ValueError:
        return None

<<<<<<< HEAD
=======
def make_trait_reported_string(s_raw):
    ''' Takes the raw trait name and outputs trnasformed name
    '''

    # Replace any double spaces with single
    s_raw = re.sub(r' +', r' ', s_raw)

    # Assert no "|" in trait name
    assert "|" not in s_raw

    # Split prefix
    parts = s_raw.split(': ', 1)

    # Move prefix to end if exists
    if len(parts) == 2:
        trait = " | ".join([parts[1], parts[0]])
    else:
        trait = s_raw

    # Capitalise the frist letter
    trait = trait.capitalize()
    
    return trait

def count_prefixes(l, sep=': '):
    ''' Counts the occurence of prefixes based on sep
    '''
    # Extract prefixes
    prefixes = []
    for entry in l:
        parts = entry.split(sep)
        if len(parts) > 1:
            prefixes.append(parts[0])
    # Count
    counts = {}
    for prefix in prefixes:
        try:
            counts[prefix] += 1
        except KeyError:
            counts[prefix] = 1
    
    return counts

def combine_rows(items):
    return ';'.join(items)
>>>>>>> d58d0d4c6c7e6822fe566ee129890ba82045f2b7

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
<<<<<<< HEAD
    parser.add_argument(
        "--input",
        metavar="<str>",
        help=("TSV file with UK Biobank (Neale V2 and SAIGE) manifest"),
        type=str,
        required=True,
    )
    parser.add_argument(
        "--output",
        metavar="<str>",
        help=("JSON file with all Saige and Neale studies and their metadata."),
        type=str,
        required=True,
    )
=======
    parser.add_argument('--in_manifest', metavar="<str>", help=("Input"), type=str, required=True)
    parser.add_argument('--in_efos', metavar="<str>", help=("EFO mapping file"), type=str, required=True)
    parser.add_argument('--prefix_counts', metavar="<str>", help=("File to output prefix counts to"), type=str, required=True)
    parser.add_argument('--outf', metavar="<str>", help=("Output"), type=str, required=True)
>>>>>>> d58d0d4c6c7e6822fe566ee129890ba82045f2b7
    args = parser.parse_args()
    return args


if __name__ == "__main__":

    # Parse args
    args = parse_args()

    main(
        input_path=args.input,
        output_path=args.output,
    )
