#!/usr/bin/env snakemake
'''
Makes:
    1. Locus overlap table
'''

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
import pandas as pd

from datetime import date

# Load configuration
configfile: "configs/config.yaml"
tmpdir = config['temp_dir']
KEEP_LOCAL = False
UPLOAD = True

if 'version' not in config:
    config['version'] = date.today().strftime("%y%m%d")

targets = []

# Make targets for locus overlap table (pseudo-coloc)
targets.append(
    'output/{version}/locus_overlap.parquet'.format(version=config['version']) )
if UPLOAD:
    # targets.append(GSRemoteProvider().remote(
    # '{gs_dir}/{version}/locus_overlap.tsv.gz'.format(gs_dir=config['gs_dir'],
    #     version=config['version']) ))
    targets.append(GSRemoteProvider().remote(
    '{gs_dir}/{version}/locus_overlap.parquet'.format(gs_dir=config['gs_dir'],
        version=config['version']) ))

# Trigger making of targets
rule all:
    input:
        targets

rule parquet_to_tsv:
    input:
        '{filename}.parquet'
    output:
        temp('{filename}.tsv.gz')
        # '{filename}.tsv.gz'
    run:
        print(input)
        (
            pd.read_parquet(input[0], engine='pyarrow')
              .to_csv(output[0], sep='\t', index=None)
        )

rule calculate_overlaps:
    ''' Calcs overlap between trait associated loci
    '''
    input:
        top_loci='output/{version}/toploci.tsv.gz',
        ld='output/{version}/ld.tsv.gz',
        finemap='output/{version}/finemapping.tsv.gz'
    output:
        tmpdir + '/{version}/locus_overlap.tsv.gz'
    shell:
        'pypy3 scripts/calculate_locus_set_overlaps.py '
        '--top_loci {input.top_loci} '
        '--ld {input.ld} '
        '--finemap {input.finemap} '
        '--outf {output}'

rule format_overlap:
    ''' Formats the overlap table to parquets
    '''
    input:
        tmpdir + '/{version}/locus_overlap.tsv.gz'
    output:
        'output/{version}/locus_overlap.parquet'
    shell:
        'python scripts/format_overlap_table.py '
        '--inf {input} '
        '--outf {output}'

rule overlap_to_GCS:
    ''' Copy to GCS
    '''
    input:
        'output/{version}/locus_overlap.parquet'
    output:
        GSRemoteProvider().remote(
            '{gs_dir}/{{version}}/locus_overlap.parquet'.format(gs_dir=config['gs_dir'])
            )
    shell:
        'cp {input} {output}'
