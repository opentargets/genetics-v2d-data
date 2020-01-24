#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import os
import sys
import argparse
import pyspark.sql
from pyspark.sql.types import *
from pyspark.sql.functions import *
from glob import glob
import gzip

def main():

    # Args
    args = parse_args()

    # Make spark session
    global spark
    spark = (
        pyspark.sql.SparkSession.builder
        .config("spark.master", "local[*]")
        .getOrCreate()
    )
    print('Spark version: ', spark.version)
    
    # Load data
    credset = (
        spark.read.json(arg.inf)
        .filter(col('type') == 'gwas')
        .filter(col('is95_credset'))
    )

    # Rename columns and cast types
    credset = (
        credset
        .select(
            col('study_id').cast('str').alias('study_id'),
            col('lead_chrom').cast('str').alias('lead_chrom'),
            col('lead_pos').cast('int').alias('lead_pos'),
            col('lead_ref').cast('str').alias('lead_ref'),
            col('lead_alt').cast('str').alias('lead_alt'),
            col('tag_chrom').cast('str').alias('tag_chrom'),
            col('tag_pos').cast('int').alias('tag_pos'),
            col('tag_ref').cast('str').alias('tag_ref'),
            col('tag_alt').cast('str').alias('tag_alt'),
            col('logABF').cast('float').alias('log10_ABF'),
            col('postprob').cast('float').alias('posterior_prob'),
        )
    )

    # Repartition and save
    (
        credset
        .repartitionByRane('lead_chrom', 'lead_pos')
        .write.parquet(
            args.outf,
            mode='overwrite'
        )
    )

    return 0

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--inf', metavar="<file>", help=('Credible set json'), type=str, required=True)
    parser.add_argument('--outf', metavar="<str>", help=("Output"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
