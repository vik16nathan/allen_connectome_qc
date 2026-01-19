#!/bin/bash

#csv_file="missing_tracers_proj.csv"
#csv_file="knox_experiments_included.csv" ###add 20 tracers with no injections
csv_file="../derivatives/knox_excluded_original/knox_excluded_tracers.csv"
# Read CSV into an array, skipping the header
mapfile -t csv_data < <(tail -n +2 "$csv_file")

output_dir="/data/chamal/projects/natvik/knox_qc_full_06232025/preprocessed/knox_connectome_tracers/orig_removed/"
for tracer in "${csv_data[@]}"; do  

    ##download injection fraction/injection density/projection density from Allen API
    mkdir -p "${output_dir}/${tracer}"

    wget -O "${output_dir}/${tracer}/projection_density_raw.nrrd"  "http://api.brain-map.org/grid_data/download_file/${tracer}?image=projection_density&resolution=50"
    sleep 5 ##can make less time if needed
    python3 transform_space.py "${output_dir}/$tracer/projection_density_raw.nrrd" "${output_dir}/${tracer}/projection_density_raw.mnc" -v=RAS -w=MICe -x=1
    mincresample -2 \
    "${output_dir}/$tracer/projection_density_raw.mnc"  \
    "${output_dir}/$tracer/projection_density.mnc" \
     -like "${output_dir}/$tracer/projection_density_raw.mnc" 

    ##repeat for injection fraction (if needed)
    wget -O "${output_dir}/${tracer}/injection_fraction_raw.nrrd"  "http://api.brain-map.org/grid_data/download_file/${tracer}?image=injection_fraction&resolution=50" 
    sleep 5 ##can make less time if needed
    python3 transform_space.py "${output_dir}/$tracer/injection_fraction_raw.nrrd" "${output_dir}/${tracer}/injection_fraction_raw.mnc" -v=RAS -w=MICe -x=1
    mincresample -2 \
    "${output_dir}/$tracer/injection_fraction_raw.mnc"  \
    "${output_dir}/$tracer/injection_fraction.mnc" \
    -like "${output_dir}/$tracer/injection_fraction_raw.mnc" 

    ##repeat for injection density (if needed)
    wget -O "${output_dir}/${tracer}/injection_density_raw.nrrd"  "http://api.brain-map.org/grid_data/download_file/${tracer}?image=injection_density&resolution=50" 
    sleep 5 ##can make less time if needed
    python3 transform_space.py "${output_dir}/$tracer/injection_density_raw.nrrd" "${output_dir}/${tracer}/injection_density_raw.mnc" -v=RAS -w=MICe -x=1
    mincresample -2 \
    "${output_dir}/$tracer/injection_density_raw.mnc"  \
    "${output_dir}/$tracer/injection_density.mnc" \
    -like "${output_dir}/$tracer/injection_density_raw.mnc" 
done

###download data masks
#output_dir="/data/chamal/projects/natvik/knox_qc_full_06232025/preprocessed/knox_connectome_tracers/"
#for tracer in "${csv_data[@]}"; do 
#    wget -O "${output_dir}/${tracer}/data_mask_raw.nrrd"  "http://api.brain-map.org/grid_data/download_file/${tracer}?image=data_mask&resolution=50" 
#    sleep 5 ##can make less time if needed
#    python3 transform_space.py "${output_dir}/$tracer/data_mask_raw.nrrd" "${output_dir}/${tracer}/data_mask_raw.mnc" -v=RAS -w=MICe -x=1
#    mincresample -2 \
#    "${output_dir}/$tracer/data_mask_raw.mnc"  \
#    "${output_dir}/$tracer/data_mask.mnc" \
#    -like "${output_dir}/$tracer/data_mask_raw.mnc" 
#done
