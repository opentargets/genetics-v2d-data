#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
import pandas as pd

def main():

    # Args
    in_custom = 'input_data/ld.custom.tsv.gz'
    in_postgap = 'input_data/ld.postgap.tsv.gz'
    out_merge = 'output/ld.merged.tsv.gz'

    # Load
    datasets = {}
    datasets['custom'] = pd.read_csv(in_custom, sep='\t', header=0)
    datasets['custom']['in_custom'] = True
    datasets['postgap'] = pd.read_csv(in_postgap, sep='\t', header=0)
    datasets['postgap']['in_postgap'] = True
    print(datasets['custom'].columns)

    # Unique studies
    print('Number of unique studies')
    for name, ds in datasets.items():
        print(' {0}: {1}'.format(name, len(ds.study_id.unique())))

    # Unique index
    print('Number of unique index variants')
    for name, ds in datasets.items():
        print(' {0}: {1}'.format(name, len(ds.index_variantid_b37.unique())))

    # Unique tag
    print('Number of unique tag variants')
    for name, ds in datasets.items():
        print(' {0}: {1}'.format(name, len(ds.tag_variantid_b37.unique())))

    # Tag per index
    print('Number of tag variants per index')
    for name, ds in datasets.items():
        print(' {0}: {1}'.format(name, len(ds.tag_variantid_b37.unique()) / len(ds.index_variantid_b37.unique()) ))

    # Merge custom and postgap
    print('Merging...')
    merged = pd.merge(datasets['custom'],
                      datasets['postgap'],
                      on=['study_id', 'index_variantid_b37', 'tag_variantid_b37'],
                      how='outer')
    merged = merged.loc[merged.index_variantid_b37.str.startswith('21_'), :]
    merged = merged.sort_values(['study_id', 'index_variantid_b37', 'tag_variantid_b37'])
    merged.to_csv(out_merge, sep='\t', index=None, compression='gzip')


    return 0

if __name__ == '__main__':

    main()
