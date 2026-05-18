#!/bin/bash
source load-all-modules.sh
set -euo pipefail
cd before_manual_qc
##download and binarize inj/proj data @ diff thresholds
###all command-line arguments are in the order inj_thresh proj_thresh
#mv ../allen_template_inputs ../../preprocessed/ ###load allen template inputs without needing to download

#./patch_api_timeout.sh
#./download_process_allen_template_inputs.sh

##download allen template inputs
#./download_knox_conn_data.sh

###binarize/threshold inj/proj files
./multiply_threshold_inj_proj.sh 0.01 0.01
./multiply_threshold_inj_proj.sh 0.5 0.1
./multiply_threshold_inj_proj.sh 0.4 0.05
./multiply_threshold_inj_proj.sh 0.6 0.2

##make directory for voxel count results (including OOB/vent voxels)
mkdir -p ../tables
Rscript create_automated_qc_csv.R 0.01 0.01
Rscript create_automated_qc_csv.R 0.5 0.1
Rscript create_automated_qc_csv.R 0.4 0.05
Rscript create_automated_qc_csv.R 0.6 0.2


##make images for manual QC
./make_qc_images_bin_threshold.sh 0.01 0.01
./make_qc_images_bin_threshold.sh 0.5 0.1