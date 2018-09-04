#!/usr/bin/env snakemake
'''
Ed Mountjoy (June 2018)

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

config['version'] = date.today().strftime("%y%m%d")

targets = []

# Make targets for top loci table
# targets.append(
#     'output/{version}/toploci.tsv'.format(version=config['version']) )
# if UPLOAD:
#     targets.append(GSRemoteProvider().remote(
#     '{gs_dir}/{version}/toploci.tsv'.format(gs_dir=config['gs_dir'],
#                                                                         version=config['version']) ))

# Make targets for study table
targets.append(
    'output/{version}/studies.tsv'.format(version=config['version']) )
# if UPLOAD:
#     targets.append(GSRemoteProvider().remote(
#     '{gs_dir}/{version}/studies.tsv'.format(gs_dir=config['gs_dir'],
#                                                                         version=config['version']) ))

# Make targets for finemapping table
# targets.append(
#     'output/{version}/finemapping.tsv.gz'.format(version=config['version']) )
# if UPLOAD:
#     targets.append(GSRemoteProvider().remote(
#     '{gs_dir}/{version}/finemapping.tsv.gz'.format(gs_dir=config['gs_dir'],
#                                                                             version=config['version']) ))
#
# # Make targets for ld table
# targets.append(
#     'output/{version}/ld.tsv.gz'.format(version=config['version']) )
# if UPLOAD:
#     targets.append(GSRemoteProvider().remote(
#     '{gs_dir}/{version}/ld.tsv.gz'.format(gs_dir=config['gs_dir'],
#                                                                       version=config['version']) ))
#
# # Make targets for ld table query variant inputs with population information
# targets.append(
#     'output/{version}/ld_analysis_input.tsv.gz'.format(version=config['version']) )
# if UPLOAD:
#     targets.append(GSRemoteProvider().remote(
#     '{gs_dir}/{version}/extras/ld_analysis_input.tsv.gz'.format(gs_dir=config['gs_dir'],
#                                                                       version=config['version']) ))
#
# # Make targets for locus overlap table (pseudo-coloc)
# targets.append(
#     'output/{version}/locus_overlap.tsv.gz'.format(version=config['version']) )
# if UPLOAD:
#     targets.append(GSRemoteProvider().remote(
#     '{gs_dir}/{version}/locus_overlap.tsv.gz'.format(gs_dir=config['gs_dir'],
#                                                                       version=config['version']) ))

# Trigger making of targets
rule all:
    input:
        targets
        # tmpdir + '/postgap_ld.temp.tsv.gz'

# Add workflows
include: 'scripts/top_loci_table.Snakefile'
include: 'scripts/study_table.Snakefile'
include: 'scripts/finemapping_table.Snakefile'
include: 'scripts/ld_table.Snakefile'
include: 'scripts/locus_overlap_table.Snakefile'
