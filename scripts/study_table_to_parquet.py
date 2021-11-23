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

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_study_table', metavar="<str>", type=str, required=True)
    parser.add_argument('--in_toploci', metavar="<str>", type=str, required=True)
    parser.add_argument('--sumstat_studies', metavar="<str>", help=("List of studies with sumstats"), type=str, required=True)
    parser.add_argument('--output', metavar="<str>", help=("Output merged file"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
