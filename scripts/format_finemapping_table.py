#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import argparse

from pyspark.sql.types import *
from pyspark.sql.functions import *

from common.utils import initialize_sparksession

def main():

    # Args
    args = parse_args()
    
    # Load data
    credset = (
        spark.read.json(args.inf)
        .filter(col('type') == 'gwas')
        .filter(col('is95_credset'))
    )

    # Rename columns and cast types
    credset = (
        credset
        .select(
            col('study_id').cast(StringType()).alias('study_id'),
            col('lead_chrom').cast(StringType()).alias('lead_chrom'),
            col('lead_pos').cast(IntegerType()).alias('lead_pos'),
            col('lead_ref').cast(StringType()).alias('lead_ref'),
            col('lead_alt').cast(StringType()).alias('lead_alt'),
            col('tag_chrom').cast(StringType()).alias('tag_chrom'),
            col('tag_pos').cast(IntegerType()).alias('tag_pos'),
            col('tag_ref').cast(StringType()).alias('tag_ref'),
            col('tag_alt').cast(StringType()).alias('tag_alt'),
            col('logABF').cast(DoubleType()).alias('log10_ABF'),
            col('postprob').cast(DoubleType()).alias('posterior_prob'),
        )
    )

    # Repartition and save
    (
        credset
        .repartitionByRange('lead_chrom', 'lead_pos')
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

    global spark
    spark = initialize_sparksession()

    main()
