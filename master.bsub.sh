#!/bin/sh
#BSUB -J ot_ld
#BSUB -q long
#BSUB -n 16
#BSUB -R "select[mem>32000] rusage[mem=32000] span[hosts=1]" -M32000
#BSUB -o output.%J.%I # %J=jobid; %I=array index
#BSUB -e errorfile.%J.%I
#BSUB -E <script> # Execute this script on host before main script

# Run interactive:   bsub -q normal -J interactive -n 1 -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" -M8000 -Is bash

version_date=`date +%y%m%d`
cores=16
snakemake -s 1_make_tables.Snakefile --config version=$version_date --cores $cores
snakemake -s 2_calculate_LD_table.Snakefile --config version=$version_date --cores $cores
snakemake -s 3_make_overlap_table.Snakefile --config version=$version_date --cores $cores
