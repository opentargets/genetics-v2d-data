#!/bin/sh

# set -euo pipefail

# Set args
cores=32
version_date=`date +%y%m%d`

# Load environment
source activate v2d_data

# Run pipelines
# snakemake -s 1_make_tables.Snakefile --config version=$version_date --cores 1 --rerun-incomplete
snakemake -s 2_calculate_LD_table.Snakefile --config version=$version_date --cores $cores --rerun-incomplete
# snakemake -s 3_make_overlap_table.Snakefile --config version=$version_date --cores $cores --rerun-incomplete

# Shutdown instance
gcloud compute instances stop "em-ld" --zone="europe-west1-d"

echo COMPLETE
