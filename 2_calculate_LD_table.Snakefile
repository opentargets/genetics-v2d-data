#!/usr/bin/env snakemake
'''
Makes:
  1. LD table
'''

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

from datetime import date
import sys
import pandas as pd

# Load configuration
configfile: "configs/config.yaml"
tmpdir = config['temp_dir']
KEEP_LOCAL = False
if 'version' not in config:
    config['version'] = date.today().strftime("%y%m%d")

#
# Load LD manifest and variant ID LUT
#

in_manifest = 'output/{version}/ld_analysis_input.tsv'.format(version=config['version'])

# Load manifest
manifest = pd.read_csv(
    in_manifest,
    sep='\t',
    header=0,
    dtype='object'
)

# Ony keep chromosomes in 1000G
valid_chroms = [str(x) for x in range(1, 23)] + ['X'] # DEBUG
# valid_chroms = ['22']
manifest.chrom = manifest.chrom.astype(str)
manifest = manifest.loc[manifest.chrom.isin(valid_chroms), :]

# manifest = manifest.iloc[:20, :] # DEBUG

# Make variant id list
varid_list = manifest.variant_id.unique().tolist()

#
# Make main LD target and workflow
#

targets = []

# # Make targets for ld table
targets.append(
    'output/{version}/ld.parquet'.format(version=config['version']) )

# Trigger making of targets
rule all:
    input:
        targets

# Add workflows
include: 'snakefiles/ld_table_2.Snakefile'
