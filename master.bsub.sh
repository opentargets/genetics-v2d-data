#!/bin/sh
#BSUB -J ot_ld
#BSUB -q long
#BSUB -n 16
#BSUB -R "select[mem>32000] rusage[mem=32000] span[hosts=1]" -M32000
#BSUB -o output.%J
#BSUB -e errorfile.%J

# Run interactive:   bsub -q normal -J interactive -n 1 -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" -M8000 -Is bash

set -euo pipefail

# Set args
cores=16
version_date=`date +%y%m%d`

# Load environment
source activate v2d_data

# Run pipelines
snakemake -s 1_make_tables.Snakefile --config version=$version_date --cores 1 --rerun-incomplete
# snakemake -s 2_calculate_LD_table.Snakefile --unlock --config version=$version_date
snakemake -s 2_calculate_LD_table.Snakefile --config version=$version_date --cores $cores --rerun-incomplete
snakemake -s 3_make_overlap_table.Snakefile --config version=$version_date --cores $cores --rerun-incomplete

echo COMPLETE
