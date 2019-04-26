#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
import argparse
import pandas as pd
from pprint import pprint
import numpy as np
import requests
import json

def main():

    # Parse args
    args = parse_args()

    # Load ICD10 curation
    icd10 = pd.read_csv(args.in_icd10, sep=',', header=0)
    icd10 = ( icd10.loc[icd10['Confidence (High/Medium/Low/None)'].isin(['High']), :]
                   .loc[:, ['Field.code', 'Manual curation']]
                   .rename(columns={'Field.code': 'trait_code',
                                    'Manual curation': 'efo_code'}) )
    # Split
    icd10['efo_code'] = icd10['efo_code'].astype(str).str.split('\|\|')
    icd10 = explode(icd10, ['efo_code'])

    # Load self-reported curation
    selfrep = pd.read_csv(args.in_self, sep=',', header=0)
    selfrep = ( selfrep.loc[selfrep['Confidence (High/Medium/Low/None)'].isin(['High', 'Medium']), :]
                   .loc[:, ['Field.code', 'Manual curation']]
                   .rename(columns={'Field.code': 'trait_code',
                                    'Manual curation': 'efo_code'}) )
    # Split
    selfrep['efo_code'] = selfrep['efo_code'].astype(str).str.split('\|\|')
    selfrep = explode(selfrep, ['efo_code'])

    # Combine
    merge = pd.concat([icd10, selfrep])
    merge['efo_code'] = merge['efo_code'].str.replace(':', '_')

    # Get labels for EFO codes
    print('Getting EFO labels from OLS API...')
    merge['efo_trait'] = merge['efo_code'].apply(get_efo_label)

    # Write
    merge.to_csv(args.outf, sep='\t', index=None)

def get_efo_label(code):
    ''' Gets the mapped trait name from the efo code from the OLS API
    '''
    # Get from API
    url = 'https://www.ebi.ac.uk/ols/api/search?q={code}&queryFields=short_form'.format(code=code)
    resp = requests.get(url)
    # Check the response was ok
    if not (resp.ok):
        resp.raise_for_status()
    # Load the response data into a dict variable
    data = json.loads(resp.content.decode('utf-8'))

    # Extract label
    label = None
    for entry in data['response']['docs']:
        if entry['short_form'] == code:
            label = entry['label']
            break

    return label

def combine_rows(items):
    return ';'.join(items)

def explode(df, columns):
    ''' Explodes multiple columns
    '''
    idx = np.repeat(df.index, df[columns[0]].str.len())
    a = df.T.reindex(columns).values
    concat = np.concatenate([np.concatenate(a[i]) for i in range(a.shape[0])])
    p = pd.DataFrame(concat.reshape(a.shape[0], -1).T, idx, columns)
    return pd.concat([df.drop(columns, axis=1), p], axis=1).reset_index(drop=True)

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_icd10', metavar="<str>", help=("Input"), type=str, required=True)
    parser.add_argument('--in_self', metavar="<str>", help=("Input"), type=str, required=True)
    parser.add_argument('--outf', metavar="<str>", help=("Output"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
