#!/usr/bin/env snakemake
'''
Makes:
  1. Top loci table
  2. Study table
  3. Finemapping table
  4. Input manifest for LD table
'''

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

from datetime import date

# Load configuration
configfile: "configs/config.yaml"
tmpdir = config['temp_dir']
KEEP_LOCAL = True
UPLOAD = False # DEBUG

if 'version' not in config:
    config['version'] = date.today().strftime("%y%m%d")

targets = []

# Make targets for top loci table
targets.append(
    'output/{version}/toploci.parquet'.format(version=config['version']) )

# Make targets for study table
## targets.append(
##     'output/{version}/studies.tsv'.format(version=config['version']) )
# targets.append(
#     'output/{version}/studies.json'.format(version=config['version']) )
targets.append(
    'output/{version}/studies.parquet'.format(version=config['version']) )


#
# # Make targets for finemapping table
# targets.append(
#     'output/{version}/finemapping.tsv.gz'.format(version=config['version']) )
#
# # Make targets for LD input table
# targets.append(
#     'output/{version}/ld_analysis_input.tsv.gz'.format(version=config['version']) )
#
# if UPLOAD:
#     targets.append(GSRemoteProvider().remote(
#     '{gs_dir}/{version}/toploci.tsv'.format(gs_dir=config['gs_dir'],
#         version=config['version']) ))
#     targets.append(GSRemoteProvider().remote(
#     '{gs_dir}/{version}/studies.tsv'.format(gs_dir=config['gs_dir'],
#         version=config['version']) ))
#     targets.append(GSRemoteProvider().remote(
#     '{gs_dir}/{version}/studies.json'.format(gs_dir=config['gs_dir'],
#         version=config['version']) ))
#     targets.append(GSRemoteProvider().remote(
#     '{gs_dir}/{version}/finemapping.tsv.gz'.format(gs_dir=config['gs_dir'],
#         version=config['version']) ))
#     targets.append(GSRemoteProvider().remote(
#     '{gs_dir}/{version}/extras/ld_analysis_input.tsv.gz'.format(gs_dir=config['gs_dir'],
#         version=config['version']) ))

# Trigger making of targets
rule all:
    input:
        targets

# Add workflows
include: 'scripts/top_loci_table.Snakefile'
include: 'scripts/study_table.Snakefile'
include: 'scripts/finemapping_table.Snakefile'
include: 'scripts/ld_table_1.Snakefile'
