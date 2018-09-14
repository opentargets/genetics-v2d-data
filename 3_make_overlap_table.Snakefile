#!/usr/bin/env snakemake
'''
Makes:
    1. Locus overlap table
'''

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

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
    'output/{version}/locus_overlap.tsv.gz'.format(version=config['version']) )
if UPLOAD:
    targets.append(GSRemoteProvider().remote(
    '{gs_dir}/{version}/locus_overlap.tsv.gz'.format(gs_dir=config['gs_dir'],
        version=config['version']) ))

# Trigger making of targets
rule all:
    input:
        targets

# Add workflows
include: 'scripts/locus_overlap_table.Snakefile'
