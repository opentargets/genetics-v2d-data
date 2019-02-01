#!/usr/bin/env bash
#

set -euo pipefail

gsutil cp gs://genetics-portal-data/v2d/ld.tsv.gz .

echo COMPLETE
