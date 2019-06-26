#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#
'''
- Reads LD information output by plink
- Calculates study/population weighted LD
- Calulcates PICS probabilities, and credible sets

# Set SPARK_HOME and PYTHONPATH to use 2.4.0
export PYSPARK_SUBMIT_ARGS="--driver-memory 8g pyspark-shell"
export SPARK_HOME=/Users/em21/software/spark-2.4.0-bin-hadoop2.7
export PYTHONPATH=$SPARK_HOME/python:$SPARK_HOME/python/lib/py4j-2.4.0-src.zip:$PYTHONPATH
'''

import os
import sys
import numpy as np
import argparse
import pyspark.sql
from pyspark.sql.types import *
from pyspark.sql.functions import *
from pyspark.sql.window import Window
from scipy.stats import norm

def main():

    # Args
    args = parse_args()
    # args.in_ld_pattern = 'input_data/ld_each_variant/*.ld.tsv.gz'
    # args.in_manifest = 'input_data/190625/ld_analysis_input.tsv'
    # args.in_top_loci = 'input_data/190625/toploci.parquet'
    # args.out = 'output/ld_w_crediblesets.parquet'
    # args.min_r2 = 0.5

    # Make spark session
    global spark
    spark = (
        pyspark.sql.SparkSession.builder
        .config("spark.master", "local[*]")
        .getOrCreate()
    )
    print('Spark version: ', spark.version)
    
    #
    # Load data ---------------------------------------------------------------
    #

    # Load LD
    ld = (
        load_ld(args.in_ld_pattern)
        .withColumn('index_variant_id', regexp_replace(col('index_variant_id'), ':', '_'))
        .withColumn('tag_variant_id', regexp_replace(col('tag_variant_id'), ':', '_'))
        # .limit(10000) # Debug
    )
    
    # Load manifest
    manifest = (
        load_manifest(args.in_manifest)
        .withColumnRenamed('variant_id', 'index_variant_id')
    )
    
    #
    # Weight correlations by study population ---------------------------------
    #

    # Join LD to manifest
    data = manifest.join(
        ld,
        on='index_variant_id'
    )

    # Replace all R values == 1 with 0.9999995, otherwise we get error
    # This is reverted later by rounding to 6 dp
    for coln in ['R_AFR', 'R_AMR', 'R_EAS', 'R_EUR', 'R_SAS']:
        data = data.withColumn(
            coln, 
            when(col(coln) == 1, 0.9999995).otherwise(col(coln))
        )

    # Fisher transform correlations to z-scores
    for coln in ['R_AFR', 'R_AMR', 'R_EAS', 'R_EUR', 'R_SAS']:
        data = data.withColumn(
            coln.replace('R_', 'Z_'),
            arctanh(col(coln))
        )
    
    # Compute weighted average across populations
    data = data.withColumn('Z_overall',
        (
            (col('AFR_prop') * col('Z_AFR')) +
            (col('AMR_prop') * col('Z_AMR')) +
            (col('EAS_prop') * col('Z_EAS')) +
            (col('EUR_prop') * col('Z_EUR')) +
            (col('SAS_prop') * col('Z_SAS'))
        )
    )

    # Inverse Fisher transform weigthed z-score back to correlation
    data = data.withColumn('R_overall', tanh(col('Z_overall')))

    # Round R_overall to 6 dp
    data = data.withColumn('R_overall', round6dp(col('R_overall')))


    # Convert R to R2
    data = data.withColumn('R2_overall',
        pow(col('R_overall'), 2)
    )

    # Drop rows where R2 is null
    data = data.filter(col('R2_overall').isNotNull())

    # Filter based on overall R2
    data = data.filter(col('R2_overall') >= args.min_r2)

    # Drop unneeded columns
    data = data.drop(*['Z_overall','R_overall', 'R_AFR',
                       'R_AMR', 'R_EAS', 'R_EUR', 'R_SAS', 'Z_AFR',
                       'Z_AMR', 'Z_EAS', 'Z_EUR', 'Z_SAS', 'index_variant_id'])
    
    # Denormalise variant IDs
    data = (
        data
        .withColumnRenamed('chrom', 'lead_chrom')
        .withColumnRenamed('pos', 'lead_pos')
        .withColumnRenamed('ref', 'lead_ref')
        .withColumnRenamed('alt', 'lead_alt')
        .withColumn('tag_split', split(col('tag_variant_id'), '_'))
        .withColumn('tag_chrom', col('tag_split').getItem(0))
        .withColumn('tag_pos', col('tag_split').getItem(1).cast('int'))
        .withColumn('tag_ref', col('tag_split').getItem(2))
        .withColumn('tag_alt', col('tag_split').getItem(3))
        .drop('tag_split', 'tag_variant_id')
    )   
    
    #
    # Conduct credible set analysis using PICS adjustment ---------------------
    #

    ''' Probabilistic Identification of Causal SNPs (PICS) from Farh (2014):
            https://www.nature.com/articles/nature13835

        Adjusts the p-values for tag SNPs based on the p-value of the lead SNP
        and it's LD.
    '''

    # Empiric constant that can be adjusted to fit the curve, 6.4 recommended.
    k = 6.4

    # Load toploci
    toploci = spark.read.parquet(args.in_top_loci)

    # Join negative log pvalue from toploci onto data
    toploci = (
        toploci
        .withColumn('neglog_p', 
                    -1 * (log10(col('pval_mantissa')) + col('pval_exponent')))
        .withColumnRenamed('chrom', 'lead_chrom')
        .withColumnRenamed('pos', 'lead_pos')
        .withColumnRenamed('ref', 'lead_ref')
        .withColumnRenamed('alt', 'lead_alt')
        .select('study_id', 'lead_chrom', 'lead_pos', 'lead_ref',
                'lead_alt', 'neglog_p')
    )
    data = data.join(
        toploci,
        on=['study_id', 'lead_chrom', 'lead_pos',
            'lead_ref', 'lead_alt']
    )

    # Calculate PICS statistics
    data = (
        data
        .withColumn('pics_mu',
            col('R2_overall') * col('neglog_p')
        )
        .withColumn('pics_std',
            sqrt(1 - sqrt(col('R2_overall'))**k) *
            sqrt(col('neglog_p')) / 2
        )
        .withColumn('pics_relative_prob',
            when(col('pics_std') == 0, 1.0).otherwise(
                norm_sf(col('pics_mu'), col('pics_std'), col('neglog_p'))
            )
        )
    )

    # Calculate the sum of the posterior probabilities at each locus
    pics_prob_sums = (
        data.groupby('study_id', 'lead_chrom', 'lead_pos',
                     'lead_ref', 'lead_alt')
        .agg(sum('pics_relative_prob').alias('pics_relative_prob_sum'))
    )

    # Merge back onto data
    data = data.join(
        pics_prob_sums,
        on=['study_id', 'lead_chrom', 'lead_pos',
            'lead_ref', 'lead_alt']
    )

    # Calculate posterior probability at each locus
    data = (
        data
        .withColumn(
            'pics_postprob',
            col('pics_relative_prob') / col('pics_relative_prob_sum')
        )
        .drop('pics_relative_prob_sum', 'neglog_p')
    )

    # Calculate cumulative sum per locus
    window_spec = (
        Window.partitionBy('study_id', 'lead_chrom', 'lead_pos',
                           'lead_ref', 'lead_alt')
        .orderBy(desc('pics_postprob'))
        .rowsBetween(Window.unboundedPreceding, Window.currentRow)
    )
    data = (
        data.withColumn(
            'pics_postprob_cumsum',
            sum('pics_postprob').over(window_spec)
        )
    )

    # Label whether each row is in the 95 and 99% credible sets
    window_spec = (
        Window.partitionBy('study_id', 'lead_chrom', 'lead_pos',
                           'lead_ref', 'lead_alt')
        .orderBy('pics_postprob_cumsum')
    )
    data = (
        data
        .withColumn(
            'pics_95perc_credset',
            when(
                lag('pics_postprob_cumsum', 1).over(window_spec) >= 0.95,
                False).otherwise(True)
        )
        .withColumn(
            'pics_99perc_credset',
            when(
                lag('pics_postprob_cumsum', 1).over(window_spec) >= 0.99,
                False).otherwise(True)
        )
    )

    #
    # Write output ------------------------------------------------------------
    #

    # Rename columns and format
    data = (
        data
        .withColumnRenamed('AFR_prop', 'AFR_1000G_prop')
        .withColumnRenamed('AMR_prop', 'AMR_1000G_prop')
        .withColumnRenamed('EAS_prop', 'EAS_1000G_prop')
        .withColumnRenamed('EUR_prop', 'EUR_1000G_prop')
        .withColumnRenamed('SAS_prop', 'SAS_1000G_prop')
        .withColumnRenamed('R2_overall', 'overall_r2')
        .select('study_id',
                'lead_chrom',
                'lead_pos',
                'lead_ref',
                'lead_alt',
                'tag_chrom',
                'tag_pos',
                'tag_ref',
                'tag_alt',
                'overall_r2',
                'pics_mu',
                'pics_postprob',
                'pics_95perc_credset',
                'pics_99perc_credset',
                'AFR_1000G_prop',
                'AMR_1000G_prop',
                'EAS_1000G_prop',
                'EUR_1000G_prop',
                'SAS_1000G_prop')
    )

    # Save output
    (
        data
        .repartitionByRange('study_id', 'lead_chrom', 'lead_pos')
        .write
        .parquet(
            args.out,
            mode='overwrite'
        )
    )
    
    return 0

@udf(DoubleType())
def norm_sf(mu, std, neglog_p):
    ''' Wrap norm(mean, std).sf(neglog_p) * 2 in a UDF'''
    return float(norm(mu, std).sf(neglog_p) * 2)

@udf(DoubleType())
def round6dp(x):
    ''' Wrap np.around(x, 6) in a UDF'''
    try:
        return float(np.around(x, 6))
    except TypeError:
        return None

@udf(DoubleType())
def arctanh(x):
    ''' Wrap np.arctanh in a UDF'''
    try:
        return float(np.arctanh(x))
    except AttributeError:
        return None

@udf(DoubleType())
def tanh(x):
    ''' Wrap np.tanh in a UDF'''
    try:
        return float(np.tanh(x))
    except AttributeError:
        return None

def load_manifest(inf):
    ''' Loads manifest file
    '''
    # Specify schema
    import_schema = (
        StructType()
        .add('study_id', StringType())
        .add('variant_id', StringType())
        .add('chrom', StringType())
        .add('pos', IntegerType())
        .add('ref', StringType())
        .add('alt', StringType())
        .add('AFR_prop', DoubleType())
        .add('AMR_prop', DoubleType())
        .add('EAS_prop', DoubleType())
        .add('EUR_prop', DoubleType())
        .add('SAS_prop', DoubleType())
    )
    # Load
    df = (
        spark.read.csv(inf,
                       sep='\t',
                       schema=import_schema,
                       enforceSchema=True,
                       header=True)
    )

    return df

def load_ld(in_pattern):
    ''' Loads all LD information from individual files
    '''
    # Specify schema
    import_schema = (
        StructType()
        .add('index_variant_id', StringType())
        .add('tag_variant_id', StringType())
        .add('R_AFR', DoubleType())
        .add('R_AMR', DoubleType())
        .add('R_EAS', DoubleType())
        .add('R_EUR', DoubleType())
        .add('R_SAS', DoubleType())
    )
    # Load
    df = (
        spark.read.csv(in_pattern,
                       sep='\t',
                       schema=import_schema,
                       enforceSchema=True,
                       header=True)
    )

    return df

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_ld_pattern', metavar="<str>", help=("Pattern to match all individual variant LD files"), type=str, required=True)
    parser.add_argument('--in_manifest', metavar="<str>", help=("Input manifest file"), type=str, required=True)
    parser.add_argument('--in_top_loci', metavar="<str>", help=("Input top loci table"), type=str, required=True)
    parser.add_argument('--min_r2', metavar="<float>", help=("Minimum R2"), type=float, required=True)
    parser.add_argument('--out', metavar="<str>", help=("Output file"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
