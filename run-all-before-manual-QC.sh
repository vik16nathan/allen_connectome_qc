#!/bin/bash                                                                                           
module load anaconda R
module load minc-toolkit-v2
module load RMINC


##download PyQC
git clone https://github.com/CoBrALab/PyQC
cp ./PyQC/PyQC.py ./

##download and binarize inj/proj data @ diff thresholds
###all command-line arguments are in the order inj_thresh proj_thresh
./download_knox_conn_data.sh
mv ./allen_template_inputs ../preprocessed/

##optional: to see how the label/mask files were downloaded and transformed
./download_process_allen_template_inputs.sh

###binarize/threshold inj/proj files
./multiply_threshold_inj_proj.sh 0.5 0.1
./multiply_threshold_inj_proj.sh 0.01 0.01

##make directory for voxel count results (including OOB/vent voxels)
mkdir -p tables
Rscript create_automated_qc_csv.R 0.5 0.1
Rscript create_automated_qc_csv.R 0.01 0.01

##make images for manual QC
./make_qc_images_bin_threshold.sh 0.5 0.1
./make_qc_images_bin_threshold.sh 0.01 0.01