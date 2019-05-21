#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#
# Merge top loci table
#

import os
import sys
import argparse
import pandas as pd
from collections import OrderedDict
from parquet_writer import write_parquet
import datetime
import json

def main():

    # Parse args
    args = parse_args()

    # Load
    merged = pd.read_json(args.in_study_table, orient='records', lines=True)

    # Load the number of associated loci per study
    top_loci = pd.read_parquet(args.in_toploci, engine='pyarrow')
    num_loci = ( top_loci
                   .groupby('study_id')
                   .size().reset_index(name='counts')
                   .rename(columns={'counts': 'num_assoc_loci'}) )
    merged = pd.merge(merged, num_loci, how='left', on='study_id')
    merged['num_assoc_loci'] = merged['num_assoc_loci'].fillna(value=0).astype(int)

    #
    # Annotate whether summary stats are available ----------------------------
    #

    # Load list of study IDs that have sumstats
    from_sumstats = set([])
    with open(args.sumstat_studies, 'r') as in_h:
        for line in in_h:
            # Get study_id
            line = line.rstrip().rstrip('/')
            stid = os.path.basename(line).replace('.parquet', '')
            # If GCST id, add "_1" suffix
            if stid.startswith('GCST'):
                from_sumstats.add(stid + '_1')
            else:
                from_sumstats.add(stid)
    
    # Annotate study table with field showing if there are sumstats
    merged['has_sumstats'] = [
        stid in from_sumstats for stid in merged['study_id']
    ]

    #
    # Annotate EFOs with therapeutic area (trait category) --------------------
    #

    # Make sure that the trait_efos column is a list
    merged['trait_efos'] = merged['trait_efos'].apply(clean_efo)

    # Load efo to category mappings
    efo_anno = {}
    with open(args.efo_categories, 'r') as in_h:
        for line in in_h:
            parts = json.loads(line)
            efo_anno[parts['efo_term']] = parts['therapeutic_areas']
    
    # Load theraputic areas to sort results against
    ta_dict = load_therapeutic_area_labels(args.in_ta)

    # Annotate efos with therapeutic areas
    merged['trait_category_list'] = merged['trait_efos'].apply(
        annotate_efos,
        efo_anno_dict=efo_anno,
        sort_list=list(ta_dict.values())
    )

    # Select first
    merged['trait_category'] = merged['trait_category_list'].apply(
        select_best_category,
        unknown_label=args.unknown_label
    )

    # Remove category_list
    merged = merged.drop('trait_category_list', axis=1)

    #
    # Format output parquet ---------------------------------------------------
    #


    # # Debug
    # merged.to_csv('tmp_study_table.tsv', sep='\t', index=None)
    # sys.exit()

    # Coerce data types
    dtypes = OrderedDict([
        ('study_id', 'object'),
        ('pmid', 'object'),
        ('pub_date', 'object'),
        ('pub_journal', 'object'),
        ('pub_title', 'object'),
        ('pub_author', 'object'),
        ('trait_reported', 'object'),
        ('trait_efos', 'object'),
        ('ancestry_initial', 'object'),
        ('ancestry_replication', 'object'),
        ('n_initial', 'Int64'),
        ('n_replication', 'Int64'),
        ('n_cases', 'Int64'),
        ('trait_category', 'object'),
        ('num_assoc_loci', 'Int64'),
        ('has_sumstats', 'bool')
        
    ])
    assert(set(dtypes.keys()) == set(merged.columns))
    merged = (
        merged.loc[:, dtypes.keys()]
        .astype(dtype=dtypes)
    )

    # Split array columns
    split_cols = ['ancestry_initial', 'ancestry_replication']
    for col in split_cols:
        merged[col] = merged[col].str.split(';')

    # Sort output
    merged = merged.sort_values(['study_id'])

    # DEBUG output study table
    merged.to_csv('tmp/study_table.tsv', sep='\t', index=None)

    # Save as parquet
    array_cols = ['trait_efos', 'ancestry_initial', 'ancestry_replication']
    write_parquet(merged,
                  args.output,
                  str_list_cols=array_cols,
                  compression='snappy',
                  flavor='spark')

    return 0


def clean_efo(value):
    ''' Always returns a list of EFOs or empty list
    '''
    if isinstance(value, list):
        return value
    elif pd.isnull(value):
        return []
    elif isinstance(value, str):
        return [value]
    else:
        assert True, 'Error: unrecognised EFO column type {} {}'.format(type(value), value)

def select_best_category(cat_list, unknown_label):
    ''' Selects the best (first) category for use in the portal
    '''
    try:
        return cat_list[0].capitalize()
    except IndexError:
        return unknown_label.capitalize()


def annotate_efos(efo_list, efo_anno_dict, sort_list=None):
    ''' Returns therapeutic area annotations for one or more efos
    Params:
        efo_list (list): list of EFO terms
        efo_anno_dict (dict): dictionary of efo to TA mappings
        sort_list (list): list against which to sort terms
    Return:
        list of therapeutic ares
    '''
    # Make list
    ta_list = []
    for efo in efo_list:
        tas = efo_anno_dict.get(efo, [])
        ta_list = ta_list + tas

    # Sort list
    if sort_list:
        ta_list = sorted(
            ta_list, key=lambda x: sort_list.index(x)
        )

    return ta_list

def load_therapeutic_area_labels(inf):
    ''' Loads therapeutic labels and display labels
    '''
    d = OrderedDict()
    with open(inf, 'r') as in_h:
        in_h.readline()  # skip header
        for line in in_h:
            
            if line.startswith('#'):
                continue

            try:
                category, term_label, display_label = line.rstrip().split('\t')
            except ValueError:
                sys.exit('Error for in {}: {}'.format(inf, line))

            if term_label in d:
                sys.exit(
                    'Error: duplicate term_label in therapuetic'
                    ' area list: {}'.format(term_label)
                )
            d[term_label] = display_label
    return d

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_study_table', metavar="<str>", type=str, required=True)
    parser.add_argument('--in_toploci', metavar="<str>", type=str, required=True)
    parser.add_argument('--sumstat_studies', metavar="<str>", help=("List of studies with sumstats"), type=str, required=True)
    parser.add_argument('--efo_categories', metavar="<str>", help=("Efo category mappings"), type=str, required=True)
    parser.add_argument('--in_ta', metavar="<str>", help=("Therapuetic are list"), type=str, required=True)
    parser.add_argument('--unknown_label', metavar="<str>", help=("Category label when unkown"), type=str, required=True)
    parser.add_argument('--output', metavar="<str>", help=("Output merged file"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
