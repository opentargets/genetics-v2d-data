#!/usr/bin/env snakemake
'''
Makes:
  1. Top loci table
  2. Study table
  3. Trait to EFO look up table
  4. Finemapping table
  5. Input manifest for LD calculation table
'''

from datetime import date

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

# Load configuration
configfile: "configs/config.yaml"
tmpdir = config['temp_dir']
KEEP_LOCAL = False

# Extract run version from configuration if exists:
version = config['version'] if 'version' not in config else date.today().strftime("%y%m%d")

# Trigger making of targets
rule all:
    input:  
        f'output/{version}/toploci.parquet'  # Make targets for top loci table
        f'output/{version}/studies.parquet'  # Make targets for study table
        f'output/{version}/trait_efo.parquet'  # Make targets for study table
        f'output/{version}/finemapping.parquet'  # Make targets for finemapping table
        f'output/{version}/ld_analysis_input.tsv'  # Make targets for LD input table
        f'output/{version}/trait_efo.parquet'  # Make trait to EFO lut table
    log: f"logs/{version}/1_make_tables.log"

# Add workflows
include: 'snakefiles/study_and_top_loci_tables.Snakefile'
include: 'snakefiles/finemapping_table.Snakefile'
include: 'snakefiles/ld_table_1.Snakefile'