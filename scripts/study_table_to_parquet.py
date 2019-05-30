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
            from_sumstats.add(stid)
    
    # Annotate study table with field showing if there are sumstats
    merged['has_sumstats'] = merged['study_id'].isin(from_sumstats)

    #
    # Fail if there are any sumstat studies not in table ----------------------
    #

    # Get list of study IDs that are in the sumstats but not the study table
    from_sumstats_series = pd.Series(list(from_sumstats))
    missing_ids = from_sumstats_series[
        ~from_sumstats_series.isin(merged['study_id'])
    ]

    # Add missing studies manually (we shouldn't have to do this but
    # occasionally GWAS Catalog will remove studies without warning due to
    # a curation error being detected)
    manual_additions = pd.DataFrame([
        {'ancestry_initial': 'European=360063',
        'ancestry_replication': None,
        'n_cases': 108473,
        'n_initial': 360063,
        'n_replication': None,
        'pmid': None,
        'pub_author': 'UKB Neale v2',
        'pub_date': '2018-08-01',
        'pub_journal': '',
        'pub_title': '',
        'study_id': 'NEALE2_6160_1',
        'trait_efos': None,
         'trait_reported': 'Leisure/social activities: Sports club or gym',
        'has_sumstats': True},
        {'ancestry_initial': 'European=360088',
        'ancestry_replication': None,
        'n_cases': 326854,
        'n_initial': 360088,
        'n_replication': None,
        'pmid': None,
        'pub_author': 'UKB Neale v2',
        'pub_date': '2018-08-01',
        'pub_journal': '',
        'pub_title': '',
         'study_id': 'NEALE2_670_1',
        'trait_efos': None,
        'trait_reported': 'Type of accommodation lived in: A house or bungalow',
        'has_sumstats': True},
        {'ancestry_initial': 'European=359931',
        'ancestry_replication': None,
        'n_cases': 2818,
        'n_initial': 359931,
        'n_replication': None,
        'pmid': None,
        'pub_author': 'UKB Neale v2',
        'pub_date': '2018-08-01',
        'pub_journal': '',
        'pub_title': '',
        'study_id': 'NEALE2_6142_7',
        'trait_efos': None,
         'trait_reported': 'Current employment status: Full or part-time student',
        'has_sumstats': True},
        {'ancestry_initial': 'European=20883',
         'ancestry_replication': 'European=41309;Greater Middle Eastern (Middle Eastern / North African / Persian)=493;South Asian=1174;East Asian=5409',
         'n_cases': 5956,
         'n_initial': 20883,
         'n_replication': 48385,
         'pmid': 'PMID:26192919',
         'pub_author': 'Liu JZ',
         'pub_date': '2015-07-20',
         'pub_journal': 'Nat Genet',
         'pub_title': 'Association analyses identify 38 susceptibility loci for inflammatory bowel disease and highlight shared genetic risk across populations.',
         'study_id': 'GCST003044',
         'trait_efos': ['EFO_0000384'],
         'trait_reported': "Crohn's disease"}
    ])

    # Only keep manual additions if they are missing (when GWAS Catalog re-add
    # them at a later I won't have to do anything)
    manual_additions = manual_additions.loc[
        manual_additions['study_id'].isin(missing_ids)
    ]
    
    # Merge manual additions
    merged = pd.concat([merged, manual_additions], axis=0, sort=True)

    # Assert that there are no longer any missing studies
    missing_ids = from_sumstats_series[
        ~from_sumstats_series.isin(merged['study_id'])
    ]
    if len(missing_ids) > 0:
        sys.exit(
            'ERROR: the {0} following study IDs are in the summary statistics '
            'but are missing from the study table, they will need adding '
            'manually:\n{1}'.format(len(missing_ids), missing_ids)
        )

    #
    # Annotate number of associated loci from top loci table ------------------
    #

    # Load the number of associated loci per study
    top_loci = pd.read_parquet(args.in_toploci, engine='pyarrow')
    num_loci = (top_loci
                .groupby('study_id')
                .size().reset_index(name='counts')
                .rename(columns={'counts': 'num_assoc_loci'}))
    merged = pd.merge(merged, num_loci, how='left', on='study_id')
    merged['num_assoc_loci'] = merged['num_assoc_loci'].fillna(
        value=0).astype(int)

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
