#!/bin/bash                                                                                           
module load anaconda R
module load minc-toolkit-v2
module load RMINC
module load MATLAB
module load ANTs

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

# Clone the repository
git clone https://github.com/AllenInstitute/mouse_connectivity_models
cd mouse_connectivity_models

# Install using pyproject.toml
uv pip install .

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

