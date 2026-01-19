#!/bin/bash

#GOAL: get "true" injection density in annotated area by multiplying injection fraction with injection density
#use a binary threshold to keep only inj/proj files over the threshold
inj_thresh="0.5"
proj_thresh="0.1"
###note: can change threshold as needed
#csv_file="knox_experiments_included.csv" ###add 20 tracers with no injections 
# Read CSV into an array, skipping the header
#mapfile -t csv_data < <(tail -n +2 "$csv_file")

output_dir="/data/chamal/projects/natvik/knox_qc_full_06232025/preprocessed/knox_connectome_tracers/"
csv_file="../derivatives/knox_excluded_original/knox_excluded_tracers.csv"
# Read CSV into an array, skipping the header
mapfile -t csv_data < <(tail -n +2 "$csv_file")

output_dir="/data/chamal/projects/natvik/knox_qc_full_06232025/preprocessed/knox_connectome_tracers/orig_removed/"
tracer_dir="${output_dir}/orig_removed"
#for tracer in "${csv_data[@]}"; do
for tracer in 127909584 ; do
    if [ "$tracer" -eq 310207648 ]; then
        continue
    fi
    tracer_dir="${output_dir}/${tracer}/"
    mincmath -mult "${tracer_dir}injection_fraction.mnc" "${tracer_dir}injection_density.mnc" "${tracer_dir}injection_dens_times_frac.mnc"

    ##binary thresholds
    mincmath -2 -seg -const2 ${inj_thresh} 999999 "${tracer_dir}injection_dens_times_frac.mnc" "${tracer_dir}injection_dens_times_frac_bin${inj_thresh}.mnc" 
    mincmath -2 -seg -const2 ${proj_thresh} 999999 "${tracer_dir}projection_density.mnc" "${tracer_dir}projection_density_bin${proj_thresh}.mnc" 

    ##one percent thresholds

done
