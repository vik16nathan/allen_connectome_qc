#!/bin/bash
# Run all steps before manual QC
# This includes: downloading data, binarizing, automated QC, and making QC images

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

echo "Starting pre-manual QC pipeline..."

# Load required modules
module load anaconda R
module load minc-toolkit-v2
module load RMINC

# Download PyQC
if [[ ! -d "PyQC" ]]; then
    echo "Downloading PyQC..."
    git clone https://github.com/CoBrALab/PyQC
    cp ./PyQC/PyQC.py ./
else
    echo "PyQC already downloaded"
fi

# Download and binarize inj/proj data at different thresholds
echo "Downloading Knox connectome data..."
./download_knox_conn_data.sh

if [[ -d "./allen_template_inputs" ]]; then
    mkdir -p ../preprocessed/
    mv ./allen_template_inputs ../preprocessed/
fi

# Optional: see how label/mask files were downloaded and transformed
echo "Processing Allen template inputs..."
./download_process_allen_template_inputs.sh

# Binarize/threshold inj/proj files at different thresholds
echo "Applying binary thresholds..."
./multiply_threshold_inj_proj.sh 0.5 0.1
./multiply_threshold_inj_proj.sh 0.01 0.01

# Create directory for voxel count results (including OOB/vent voxels)
mkdir -p tables

echo "Running automated QC analysis..."
Rscript create_automated_qc_csv.R 0.5 0.1
Rscript create_automated_qc_csv.R 0.01 0.01

# Make images for manual QC
echo "Creating QC images..."
./make_qc_images_bin_threshold.sh 0.5 0.1
./make_qc_images_bin_threshold.sh 0.01 0.01

echo "Pre-manual QC pipeline completed successfully"
echo "Next steps: Perform manual QC using PyQC, then run run-all-after-manual-QC.sh"