#!/usr/bin/env snakemake
'''
Ed Mountjoy (June 2018)

'''

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

# Load configuration
configfile: "configs/config.yaml"
tmpdir = config['temp_dir']
keep_local = True

targets = []

# Make targets for top loci table
outfile = 'output/ot_genetics_toploci.{version}.tsv'.format(version=config['version'])
targets.append(outfile)

# Trigger making of targets
rule all:
    input:
        targets
        # tmpdir + '/nealeUKB-associations_ot-format.{version}.tsv'.format(version=config['version'])

# Add workflows
include: 'scripts/top_loci_table.Snakefile'
