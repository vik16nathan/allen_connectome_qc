#!/bin/bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 <excluded_exps_file> <list_rgn_numbers_file> <suffix>" >&2
    exit 1
fi

source ../.venv/bin/activate
excluded_exps_file=$1
list_rgn_numbers_file=$2
suffix=$3

python build_homogeneous_model_new_excluded.py "${excluded_exps_file}" "${list_rgn_numbers_file}" "${suffix}"
