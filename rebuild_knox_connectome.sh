#!/bin/bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <excluded_exps_file> <suffix>" >&2
    exit 1
fi

source .venv/bin/activate
excluded_exps_file=$1
suffix=$2
python3 build_model_new_excluded.py "${excluded_exps_file}" "${suffix}"

