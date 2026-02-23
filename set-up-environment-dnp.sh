#!/bin/bash
source load-all-modules.sh

##############SET UP DIRECTORY STRUCTURE#####################

mkdir "../derivatives"
mkdir "../preprocessed"

#############################################################
set -euo pipefail
#########DOWNLOAD PyQC BEFORE MANUAL QC###################
git clone https://github.com/CoBrALab/PyQC
cd PyQC
uv sync
cd ..

############SETUP VIRTUAL ENVIRONMENT#####################
##download mouse_connectivity_models to rebuild connectomes
# Create virtual environment with Python 3.9
uv venv --python 3.9
source .venv/bin/activate

###download mouse_connectivity_models###


# Install using pyproject.toml
git clone https://github.com/AllenInstitute/mouse_connectivity_models
uv pip install six setuptools==80.9.0 numpy==1.23.5 allensdk scikit-learn==1.0.2
uv pip install --no-build-isolation -r pyproject.toml
uv pip install --no-build-isolation .

#############DOWNLOAD REQUIRED PACKAGES####################
##download ggslicer for visualizations
git clone https://github.com/yohanyee/ggslicer
###download Brain Connectivity Toolbox for graph theory
set -e
FILE_ID="1DmMvRnferBfGe057O-sZwB5jL4j8w1Hu"
ZIP_NAME="BCT.zip"
# install gdown if needed
command -v gdown >/dev/null 2>&1
# download
gdown --id "$FILE_ID" -O "$ZIP_NAME"
# unzip into current directory
unzip -o "$ZIP_NAME"