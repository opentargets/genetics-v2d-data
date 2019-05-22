#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#
# Calculates overall R2 weighted across populations and outputs final file
#

import sys
import os
import argparse
import pandas as pd
import numpy as np

def main():

    # Parse args
    args = parse_args()
    round_1_to = 0.9999995

    # Load
    manifest = pd.read_csv(args.in_manifest, sep='\t', header=0)
    ld = pd.read_csv(args.in_ld, sep='\t', header=0)

    # Replace ":" with "_" in ld to match manifest
    for col in ['index_variant_id', 'tag_variant_id']:
        ld[col] = ld[col].str.replace(':', '_')

    # Merge ld to manifest
    man_ld = pd.merge(manifest, ld,
                      left_on='variant_id', right_on='index_variant_id',
                      how='inner')
    # print(man_ld.head())
    # print(man_ld.shape)
    # print(man_ld.columns)

    # Convert correlations and weights to numpy arrays
    r = man_ld.loc[:, ['R_AFR', 'R_AMR', 'R_EAS', 'R_EUR', 'R_SAS']].values
    w = man_ld.loc[:, ['AFR_prop', 'AMR_prop', 'EAS_prop', 'EUR_prop', 'SAS_prop']].values

    # Have to convert r = 1 to r = 0.9999995 or we get divide by 0 error
    r[r == 1.0] = round_1_to

    # Fisher transform correlations to z-scores
    r_z = np.arctanh(r)

    # Compute weighted average
    r_z_masked = np.ma.masked_array(r_z, np.isnan(r_z))
    r_z_weighted = np.average(r_z_masked, axis=1, weights=w)

    # Inverse Fisher transformation z-scores back to correlations
    r_overall = np.tanh(r_z_weighted)

    # Round back down to 6dp
    r_overall_round = np.around(r_overall, 6)

    # Add as column in df and calculate R2
    man_ld['R_overall'] = r_overall_round.filled(np.nan)
    man_ld['R2_overall'] = man_ld['R_overall'].pow(2)

    # Write
    man_ld.to_csv(args.out, sep='\t', index=None, compression='gzip')


    return 0

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_ld', metavar="<str>", help=("Input ld file"), type=str, required=True)
    parser.add_argument('--in_manifest', metavar="<str>", help=("Input manifest file"), type=str, required=True)
    parser.add_argument('--out', metavar="<str>", help=("Output file"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
