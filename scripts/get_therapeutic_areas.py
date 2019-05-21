#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

'''
http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0002508

'''

import sys
import os
import pandas as pd
import requests
import json
from pprint import pprint
import numpy as np
from collections import OrderedDict
import argparse

def main():

    pd.set_option('display.max_columns', 500)

    # Args
    args = parse_args()

    # Load therapeutic areas dict
    ta_dict = load_therapeutic_area_labels(args.in_ta)

    # Load study info and explode efo column
    std = pd.read_json(args.in_study, orient='records', lines=True)
    std['trait_efos'] = std['trait_efos'].apply(
        lambda x: x if isinstance(x, list) else [None])
    std = explode(std, ['trait_efos'])
    
    # # DEBUG
    # std = std.head(100)

    # Get efo therapeutic areas
    efo_res = get_efo_therapeutic_areas_multi(
        efo_list=std['trait_efos'],
        ta_dict=ta_dict,
        order=True
    )

    # Write as json
    with open(args.output, 'w') as out_h:
        for efo, ta in efo_res.items():
            out_h.write(
                json.dumps(
                    {'efo_term': efo,
                     'therapeutic_areas': ta}
                ) + '\n'
            )

    return 0

def get_efo_therapeutic_areas_multi(efo_list, ta_dict, order=False):
    ''' For a list of efo terms, get a list of therapeutic areas
    params:
        efo (str): EFO short form code
        ta_dict (dict): Dictionary of therapeutic area EFO labels -> display labels
        order (bool): whether to order efo's by order in input file
    returns:
        dict of efo -> list of therapeutic areas
    '''
    d = {}
    efo_set = set([x for x in efo_list if isinstance(x, str)])
    # efo_set = set([x for x in efo_list])
    for i, efo in enumerate(set(efo_set)):
        if i % 10 == 0:
            print('Processed {} of {}...'.format(i, len(efo_set)))
        d[efo] = get_efo_therapeutic_areas(efo, ta_dict, order)
    return d

def get_efo_therapeutic_areas(efo, ta_dict, order=False):
    ''' For a single efo term, get a list of therapeutic areas
    params:
        efo (str): EFO short form code
        ta_dict (dict): Dictionary of therapeutic area EFO labels -> display labels
        order (bool): whether to order efo's by order in input file
    returns list of therapeutic area display labels
    '''
    ta_set = set([])

    # Get set of labels from ancestors
    for anc in get_efo_ancestor_terms(efo):
        if anc['label'] in ta_dict:
            ta_set.add(ta_dict[anc['label']])
    
    # Check if the efo term is a therapeutic area itself
    # This is only likely if we get to a root (phenotype, disease, measurement)
    if len(ta_set) < 2:
        efo_label = get_efo_label(efo)
        if efo_label in ta_dict:
            ta_set.add(ta_dict[efo_label])
    
    # Order labels same as input file
    if order:
        ta_list = sorted(ta_set, key=lambda x: list(ta_dict.values()).index(x))
    # Random order list
    else:
        ta_list = list(ta_set)
    
    return ta_list

def load_therapeutic_area_labels(inf):
    ''' Loads therapeutic labels and display labels
    '''
    d = OrderedDict()
    with open(inf, 'r') as in_h:
        in_h.readline() # skip header
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

def get_efo_ancestor_terms(efo):
    ''' Uses OLS API to get all ancestors for an efo term.
    Params:
        efo (str): efo short form code
    Returns:
        yields dicts of json respons for each ancestor term
    '''

    # Query OLS
    url = ("http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/"
           "http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252F{efo}/ancestors")
    url = url.format(efo=efo.replace(':', '_'))
    # print(url)

    # Get first page
    page = query_rest(url)
    
    # Process data
    while True:
        
        # Stop if no terms
        try:
            page['_embedded']['terms']
        except KeyError:
            break

        # Iterate over ancestors
        for anc in page['_embedded']['terms']:
            yield anc

        # Get next page
        if 'next' in page['_links']:
            page = query_rest(page['_links']['next']['href'])
        else:
            break

def get_efo_label(code):
    ''' Gets the mapped trait name from the efo code from the OLS API
    '''
    # Get from API
    url = ('https://www.ebi.ac.uk/ols/api/'
           'search?q={code}&queryFields=short_form')
    url = url.format(code=code)
    
    # Make query
    data = query_rest(url)

    # Extract label
    label = None
    for entry in data['response']['docs']:
        if entry['short_form'] == code:
            label = entry['label']
            break

    return label

def query_rest(url):
    ''' Queries rest api, checks response is ok and parases json
    '''
    try:
        resp = requests.get(url)
    except:
        sys.exit('Error fetching: {}'.format(url))

    # Check the response was ok
    if not (resp.ok):
        resp.raise_for_status()
    
    # Load the response data into a dict variable
    data = json.loads(resp.content.decode('utf-8'))

    return data

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
    parser.add_argument('--in_study', metavar="<str>", type=str, required=True)
    parser.add_argument('--in_ta', metavar="<str>", type=str, required=True)
    parser.add_argument('--output', metavar="<str>", help=("Output"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
