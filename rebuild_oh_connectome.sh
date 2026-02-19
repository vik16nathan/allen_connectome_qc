#!/bin/bash

source .venv/bin/activate
excluded_exps_file=$1
list_rgn_numbers_file=$2
suffix=$3

python3 build_homogeneous_model_new_excluded.py ${excluded_exps_file} ${list_rgn_numbers_file} ${suffix}
