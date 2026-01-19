#!/bin/bash

cd ./mouse_connectivity_models/
source .venv/bin/activate
python3 ./paper/figures/model_comparison/build_homogeneous_model_new_excluded.py
#python3 ./paper/figures/model_comparison/build_homogeneous_model.py
