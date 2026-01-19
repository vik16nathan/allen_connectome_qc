#!/bin/bash
module load minc-toolkit-v2
module load ANTs
# Define file paths
###NOTE: can change threshold as needed
csv_file="knox_experiments_included.csv"
inj_dir="/data/chamal/projects/natvik/knox_qc_full_06232025/preprocessed/knox_connectome_tracers/"
proj_dir=${inj_dir}
qc_jpg_dir_inj="../derivatives/knox_inj/"
qc_jpg_dir_proj="../derivatives/knox_proj/"
#qc_jpg_dir_inj="../derivatives/knox_inj/bin0.1/"
#qc_jpg_dir_proj="../derivatives/knox_proj/bin0.1/"
allen_50um_template_path="/data/chamal/projects/yohan/common/allenbrain/data/mouse_atlas/ccfv3/average_template_50.mnc"

mkdir -p ${qc_jpg_dir_inj}
mkdir -p ${qc_jpg_dir_proj}

while IFS=, read -r tracer; do 
     ##change threshold from 0.01 as needed
     inj_file_thresh="${inj_dir}/${tracer}/injection_dens_times_frac_bin0.5.mnc"
     proj_file_thresh="${proj_dir}/${tracer}/projection_density_bin0.1.mnc"

     ./make_slice_images.sh \
     --label-overlay-opacity 0.5 \
     --label-overlay $inj_file_thresh \
     "${allen_50um_template_path}" \
     "${qc_jpg_dir_inj}/${tracer}.jpg"

    ./make_slice_images.sh \
     --label-overlay-opacity 0.5 \
     --label-overlay $proj_file_thresh \
     "${allen_50um_template_path}" \
     "${qc_jpg_dir_proj}/${tracer}.jpg"

done < <(tail -n +2 "$csv_file")

###list of tracers to create cropped images (no injection seen in original QC image)
allen_50um_template_path="/data/chamal/projects/yohan/common/allenbrain/data/mouse_atlas/ccfv3/average_template_50.mnc"
inj_dir="/data/chamal/projects/natvik/knox_qc_full_06232025/preprocessed/knox_connectome_tracers/"
proj_dir=${inj_dir}
tracers=("180568155" "180708524" "146856593" "120494729" "180404418" "272970039" "112424813" "272825299" "100148443" "147790181" "120493315" "114400640" "112460257" "273055501" "277618054" "148197327" "114045733" "100148554" "180523704" "146012184" "158738894" "100141598" "126351299" "141601779" "126352037")
qc_jpg_dir_inj="../derivatives/knox_inj/bin0.1/cropped/"
tracers=("100148443" "100148554" "112460257" "114045733" "114400640" "125801739" "126115436" "126352037" "127085717" "127711098" "141601779" "146012184" "158738894" "180568155" "272970039")
qc_jpg_dir_proj="../derivatives/knox_proj/bin0.1/cropped/"
for tracer in "${tracers[@]}"; do
     
     #replace inj/proj as needed
     inj_file_thresh="${inj_dir}/${tracer}/injection_dens_times_frac_bin0.5.mnc"

      ./make_slice_images.sh \
     --label-overlay-opacity 0.5 \
     --crop-file $inj_file_thresh \
     --label-overlay $inj_file_thresh \
     "${allen_50um_template_path}" \
     "${qc_jpg_dir_inj}/${tracer}.jpg"

     proj_file_thresh="${proj_dir}/${tracer}/projection_density_bin0.1.mnc"

      ./make_slice_images.sh \
     --label-overlay-opacity 0.5 \
     --crop-file $proj_file_thresh \
     --label-overlay $proj_file_thresh \
     "${allen_50um_template_path}" \
     "${qc_jpg_dir_proj}/${tracer}.jpg"
done