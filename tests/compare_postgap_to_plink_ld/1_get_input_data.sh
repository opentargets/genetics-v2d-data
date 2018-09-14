#!/usr/bin/env bash
#

set -euo pipefail

mkdir -p input_data
gsutil cp -n gs://genetics-portal-staging/v2d/180913/ld.tsv.gz input_data/ld.custom.tsv.gz
gsutil cp -n gs://genetics-portal-staging/v2d/180904/ld.tsv.gz input_data/ld.postgap.tsv.gz


echo COMPLETE
