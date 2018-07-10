#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
import argparse
import pandas as pd
import requests
import json
from pprint import pprint

def main():

    # Parse args
    args = parse_args()

    # Get each study from the api
    data = []
    for i, study in enumerate(yield_studies()):

        if i % 20 == 0:
            print('Processing study {0}...'.format(i))

        row = {}
        row['study_id'] = study['accessionId']
        # Get publication info
        row['pmid'] = str(study['publicationInfo']['pubmedId'])
        row['pub_date'] = study['publicationInfo']['publicationDate']
        row['pub_journal'] = study['publicationInfo']['publication']
        row['pub_title'] = study['publicationInfo']['title']
        row['pub_author'] = study['publicationInfo']['author']['fullname']
        # Get trait information
        row['trait_reported'] = study['diseaseTrait']['trait']
        mapped_traits, efos = parse_efo_traits(study['_links']['efoTraits']['href'])
        row['trait_mapped'] = ';'.join(mapped_traits)
        row['trait_efos'] = ';'.join(efos)
        # Get ancestry information
        initial, repl = parse_ancestries(study['ancestries'])
        row['ancestry_initial'] = ';'.join(['='.join([str(x) for x in pair]) for pair in initial])
        row['ancestry_replication'] = ';'.join(['='.join([str(x) for x in pair]) for pair in repl])
        # Get sample size info
        row['n_initial'] = sum([entry[1] for entry in initial])
        row['n_replication'] = sum([entry[1] for entry in repl])

        # Add to output
        data.append(row)

        # DEBUG
        # if i == 5:
        #     break

    # Make into df
    df = pd.DataFrame(data)
    # Sort columns
    df = df.loc[:, [
            'study_id',
            'pmid',
            'pub_date',
            'pub_journal',
            'pub_title',
            'pub_author',
            'trait_reported',
            'trait_mapped',
            'trait_efos',
            'ancestry_initial',
            'ancestry_replication',
            'n_initial',
            'n_replication'
    ]]
    df['n_cases'] = ''

    # Remove rows where n_initial == 0, these are old studies and their data is
    # inconsistent with the rest of GWAS Catalog (correspondance w Annalisa Buniello)
    df = df.loc[df.n_initial != 0, :]

    # Write
    df.to_csv(args.outf, sep='\t', index=None)

def parse_ancestries(data):
    ''' Parse the ancestry field
    '''
    out = {'initial': [],
           'replication': []}
    for entry in data:
        try:
            anc_group = entry['ancestralGroups'][0]['ancestralGroup']
        except IndexError:
            anc_group = 'NR'
        n = entry['numberOfIndividuals']
        if n:
            out[entry['type']].append((anc_group, n))
    return out['initial'], out['replication']

def parse_efo_traits(url):
    ''' Parse EFO traits from GWAS Catalog API
    '''
    studyResponse = requests.get(url)
    # Check the response was ok
    if not (studyResponse.ok):
        studyResponse.raise_for_status()
    # Load the response data into a dict variable
    data = json.loads(studyResponse.content.decode('utf-8'))

    traits = []
    efos = []
    for entry in data['_embedded']['efoTraits']:
        traits.append(entry['trait'])
        efos.append(entry['shortForm'])

    return traits, efos

def yield_studies():
    """ Get all studies from GWAS Catalog API
    Returns:
        study_data json dict for a single study
    """

    # Search gwas cat by pubmed ID
    url = 'https://www.ebi.ac.uk/gwas/labs/rest/api/studies'
    studyResponse = requests.get(url)

    # Check the response was ok
    if not (studyResponse.ok):
        studyResponse.raise_for_status()

    # Load the response data into a dict variable
    data = json.loads(studyResponse.content.decode('utf-8'))

    while True:
        for study in data['_embedded']['studies']:
            yield study
        # Check if other pages exist
        if 'next' in data['_links']:
            # Get next page
            url = data['_links']['next']['href']
            studyResponse = requests.get(url)
            # Check the response was ok
            if not (studyResponse.ok):
                studyResponse.raise_for_status()
            # Extract json
            data = json.loads(studyResponse.content.decode('utf-8'))
        else:
            break

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--outf', metavar="<str>", help=("Output"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
