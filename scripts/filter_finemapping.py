#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#
'''
Filters finemapping credible sets to only keep those of type gwas
'''

import os
import sys
import argparse
import pyspark.sql
from pyspark.sql.types import *
from pyspark.sql.functions import *
from shutil import copyfile
from glob import glob

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
    
    # Write intermediate file
    (
        spark.read.json(args.in_path)
        .filter(col('type') == 'gwas')
        .coalesce(1)
        .write.json(
            args.intermediate_path,
            mode='overwrite',
            compression='gzip'
        )
    )

    # Copy intermediate file to final output
    src_files = glob(os.path.join(args.intermediate_path), '*.json.gz')
    assert len(src_files) == 1
    src_file = src_files[0]
    copyfile(src_file, args.out_path)

    
    return 0


def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_path', metavar="<str>", help=("Input file"), type=str, required=True)
    parser.add_argument('--intermediate_path', metavar="<str>", help=("Intermediate path"), type=str, required=True)
    parser.add_argument('--out_path', metavar="<str>", help=("Output file"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
