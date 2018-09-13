#!/bin/sh

set -euo pipefail

# Set args
cores=64
version_date=`date +%y%m%d`

# Load environment
source activate v2d_data

# Run pipelines
snakemake -s 1_make_tables.Snakefile --config version=$version_date --cores 1 --rerun-incomplete
# snakemake -s 2_calculate_LD_table.Snakefile --unlock --config version=$version_date
snakemake -s 2_calculate_LD_table.Snakefile --config version=$version_date --cores $cores --rerun-incomplete
snakemake -s 3_make_overlap_table.Snakefile --config version=$version_date --cores $cores --rerun-incomplete

# Copy output to gcs
gsutil -m -o GSUtil:parallel_composite_upload_threshold=150M rsync -r -x ".*DS_Store$" output gs://genetics-portal-staging/v2d

# Shutdown instance
gcloud compute instances stop em-ld --zone="europe-west1-d"

echo COMPLETE
