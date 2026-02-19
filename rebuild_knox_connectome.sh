#!/bin/bash

source .venv/bin/activate
excluded_exps_file=$1
suffix=$2
python3 build_model_new_excluded.py ${excluded_exps_file} ${suffix}

