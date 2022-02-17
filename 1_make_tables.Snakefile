#!/usr/bin/env snakemake
'''
Makes:
  1. Top loci table
  2. Study table
  3. Finemapping table
  4. Input manifest for LD calculation table
'''

from datetime import date

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

from src.utils import *

# Load configuration
configfile: "configs/config.yaml"
tmpdir = config['temp_dir']
KEEP_LOCAL = True

if 'version' not in config:
    config['version'] = date.today().strftime("%y%m%d")

targets = []

# Make targets for top loci table
targets.append(
    'output/{version}/toploci.parquet'.format(version=config['version']) )

# Make targets for study table
targets.append(
    'output/{version}/studies.parquet'.format(version=config['version']) )

# Make targets for study table
targets.append(
    'output/{version}/trait_efo.parquet'.format(version=config['version']) )

# Make targets for finemapping table
targets.append(
    'output/{version}/finemapping.parquet'.format(version=config['version']) )

# Make targets for LD input table
targets.append(
    'output/{version}/ld_analysis_input.tsv'.format(version=config['version']) )

# Make lut table for Irene
# Removing this currently as mappings may be generated as a separate process
#targets.append(
#    'output/{version}/trait_efo.parquet'.format(version=config['version']) )

# Trigger making of targets
rule all:
    input:
        targets
    log: f"logs/{config['version']}/1_make_tables.log"

# Add workflows
include: 'snakefiles/study_and_top_loci_tables.Snakefile'
include: 'snakefiles/finemapping_table.Snakefile'
include: 'snakefiles/ld_table_1.Snakefile'
