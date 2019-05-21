#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import re
import sys
import os
import argparse
import pandas as pd
from pprint import pprint
from collections import OrderedDict
from operator import itemgetter
import json

def main():

    # Parse args
    args = parse_args()

    # Load manifest
    manifest = pd.read_csv(args.in_manifest, sep='\t',
                           header=0, dtype=object)
    
    # Only keep required cols
    to_keep = OrderedDict([
        ('code', 'study_id'),
        ('trait', 'trait_raw'),
        ('n_total', 'n_total'),
        ('n_cases', 'n_cases')
    ])
    manifest = (
        manifest
        .loc[:, to_keep.keys()]
        .rename(columns=to_keep)
    )

    #
    # Clean trait names -------------------------------------------------------
    #

    # Get list of phenotype prefixes with counts
    counts = count_prefixes(manifest['trait_raw'].tolist())
    with open(args.prefix_counts, 'w') as out_h:
        for prefix, count in sorted(counts.items(),
                                    key=itemgetter(1),
                                    reverse=True):
            out_h.write('{}\t{}\n'.format(prefix, count))
    
    # Move prefixes to suffix (helps when viewed in the UI)
    manifest['trait_reported'] = (
        manifest['trait_raw'].apply(make_trait_reported_string)
    )   

    #
    # Add other columns -------------------------------------------------------
    #

    # Vector to get Neale or SAIGE studies
    is_neale = manifest['study_id'].str.startswith("NEALE2_")
    is_saige = manifest['study_id'].str.startswith("SAIGE_")
    assert (is_neale | is_saige).all()

    # Make columns
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

    # Write
    manifest.to_json(args.outf, orient='records', lines=True)

def to_int_safe(i):
    try:
        return int(float(i))
    except ValueError:
        return None

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
        return " | ".join([parts[1], parts[0]])
    else:
        return s_raw

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

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_manifest', metavar="<str>", help=("Input"), type=str, required=True)
    parser.add_argument('--in_efos', metavar="<str>", help=("EFO mapping file"), type=str, required=True)
    parser.add_argument('--prefix_counts', metavar="<str>", help=("File to output prefix counts to"), type=str, required=True)
    parser.add_argument('--outf', metavar="<str>", help=("Output"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
