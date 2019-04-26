#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
import argparse
import pandas as pd
import logging

def main():

    # Parse args
    args = parse_args()

    # Start logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    handler = logging.FileHandler(args.log)
    logger.addHandler(handler)

    #
    # Load data ----------------------------------------------------------------
    #

    # Load assoc data
    loci = pd.read_csv(args.inf, sep='\t', header=0)
    logger.info('Total inital associations: {0}'.format(loci.shape[0]))

    # Filter pbased on p-value
    pval = loci.pval_mantissa * (10 ** loci.pval_exponent)
    loci = loci.loc[pval <= args.min_p, :]
    logger.info('N associations after P-value filter: {0}'.format(loci.shape[0]))

    # Get chrom and pos from variant ID
    loci[['chrom', 'pos', 'alleles']] = (
        loci.variant_id_b38
            .str.split('_', 2, expand=True) )
    loci['pos'] = loci['pos'].astype(int)

    #
    # Initial clustering -------------------------------------------------------
    # Produce statistics on how many loci there are pre/post clustering
    #

    # WARNING, this must be done on a df deduplicated on (study_id, chrom, pos)
    # as the process of mapping GWAS Catalog RSIDs to variant IDs can be one
    # to many, meaning that false independent loci will have been introduced
    loci_dedup = loci.drop_duplicates(subset=['study_id', 'chrom', 'pos'])
    logger.info('N associations after deduplication on study, chrom, pos: {0}'.format(loci_dedup.shape[0]))

    # Perform pre-clustering
    loci_initial_cluster  = (
        loci_dedup.groupby('study_id')
                  .apply(distance_clumping, dist=args.cluster_dist_kb)
                  .reset_index()
    )
    logger.info('N associations after initial clustering: {0}'.format(loci_initial_cluster.shape[0]))

    # Calc number of loci per study
    pre_num = ( loci_dedup.groupby('study_id')
                          .pos
                          .count()
                          .reset_index()
                          .rename(columns={'pos': 'pre_clustering'}) )
    post_num = ( loci_initial_cluster.groupby('study_id')
                        .pos
                        .count()
                        .reset_index()
                        .rename(columns={'pos': 'post_clustering'}) )
    cluster_stats = pd.merge(pre_num, post_num)

    # Calculate stats on whether they will be clustered
    cluster_stats['proportion_multi'] = (
        ((cluster_stats.pre_clustering - cluster_stats.post_clustering) /
          cluster_stats.pre_clustering) )
    cluster_stats['to_cluster'] = (
        (cluster_stats.proportion_multi > args.cluster_multi_prop) &
        (cluster_stats.pre_clustering >= args.cluster_min_loci) )
    # cluster_stats.to_csv('cluster_stats.tsv', sep='\t', index=None)

    # Get list of studies for which clustering should be applied
    studies_to_cluster = cluster_stats.loc[cluster_stats.to_cluster, 'study_id']

    #
    # Main clustering ----------------------------------------------------------
    #

    # Split studies that require clustering from those that don't
    loci_no_cluster = loci.loc[~loci.study_id.isin(studies_to_cluster), :]
    loci_to_cluster = loci.loc[loci.study_id.isin(studies_to_cluster), :]

    logger.info('Clustering will be applied to N studies: {0}'.format(loci_to_cluster.study_id.nunique()))
    logger.info('Clustering will not be applied to N studies: {0}'.format(loci_no_cluster.study_id.nunique()))
    logger.info('Total N studies: {0}'.format(loci.study_id.nunique()))
    logger.info('Clustering will be applied to the following studies: {0}'.format(list(studies_to_cluster)))

    # Apply main distance clustering
    loci_to_clustered = (
        loci_to_cluster.groupby('study_id')
                       .apply(distance_clumping, dist=args.cluster_dist_kb)
                       .reset_index()
                       .drop('level_1', axis=1)
    )

    # Concatenate clustered and none clustered
    final_loci = pd.concat([loci_no_cluster, loci_to_clustered])
    logger.info('Final num associations: {0}'.format(final_loci.shape[0]))

    # Checks
    assert final_loci.study_id.nunique() == loci.study_id.nunique()
    assert final_loci.shape[1] == loci.shape[1]

    # Write output
    (
        final_loci.drop(['chrom', 'pos', 'alleles'], axis=1)
                  .to_csv(args.outf, sep='\t', index=None)
    )

def distance_clumping(df, dist=500):
    """ Does distance based clumping.
    Args:
        df:   pandas df in standard format
        dist: (kb) distance around index SNP to clump
    Returns:
        pandas df with additional columns showing cluster number
    """

    # Sort by pval and deduplicate
    df = (
        df.sort_values(['pval_exponent', 'pval_mantissa'])
          .drop_duplicates(subset='pos')
    )

    # Initiate clustering
    df["cluster"] = None
    clusnum = 1
    unclustered = pd.isnull(df["cluster"])

    # Continue clumping whilst there are unclustered SNPs
    while unclustered.any():

        # Get index row
        index_row = df.loc[unclustered, :].iloc[0]

        # Find other rows within set distance
        in_cluster = ( (df["chrom"] == index_row["chrom"]) &
                       (df["pos"] >= index_row["pos"] - dist * 1000) &
                       (df["pos"] <= index_row["pos"] + dist * 1000) &
                       unclustered )
        df.loc[in_cluster, "cluster"] = clusnum

        # Increase cluster number
        clusnum += 1
        unclustered = pd.isnull(df["cluster"])

    # Deduplicate on cluster num
    df = (
        df.drop_duplicates(subset='cluster', keep='first')
          .drop(['study_id', 'cluster'], axis=1)
    )

    return df

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--inf', metavar="<str>", help=('Input'), type=str, required=True)
    parser.add_argument('--outf', metavar="<str>", help=("Output"), type=str, required=True)
    parser.add_argument('--log', metavar="<str>", help=("Log output"), type=str, required=True)
    parser.add_argument('--min_p', metavar="<float>", help=("Minimum p-value to be included"), type=float, required=True)
    parser.add_argument('--cluster_dist_kb', metavar="<int>", help=("± Distance in Kb"), type=int, required=True)
    parser.add_argument('--cluster_min_loci', metavar="<int>", help=("Minimum number of reported loci for that study to be included in the clustering analysis"), type=int, required=True)
    parser.add_argument('--cluster_multi_prop', metavar="<float>", help=("For a given study, if more than this proportion of loci are multi-signals (>1 signal within gwas_cat_cluster_dist_kb), the study will be clustered"), type=float, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
