#!/bin/bash
set -euo pipefail
cd ..

#############################################################################
########STEP 2: Download Experiment Tracer Files#############################
###have these files in the knox_experiments_csvs/ directory (pulled from GitHub)

csv_file_main_exps="./knox_experiment_csvs/knox_experiments_included.csv"
csv_file_removed_exps="./knox_experiment_csvs/knox_excluded_tracers.csv"

##experiments that are in Allen API and downloaded when rebuilding connectomes,
##but that aren't within the Knox connectome
csv_file_not_in_knox_exps="./knox_experiment_csvs/wt_experiments_allen_api_not_knox.csv"
###
main_exps_output_dir="../preprocessed/knox_connectome_tracers/"
removed_exps_output_dir="../preprocessed/knox_connectome_tracers/orig_removed/"
not_in_knox_output_dir="../preprocessed/knox_connectome_tracers/allen_api_not_in_knox/"

# Read CSV into an array, skipping the header
mapfile -t csv_data_main < <(tail -n +2 "$csv_file_main_exps")
mapfile -t csv_data_removed < <(tail -n +2 "$csv_file_removed_exps")
mapfile -t csv_file_allen_api_not_knox < <(tail -n +2 "$csv_file_not_in_knox_exps")


output_dir="../preprocessed/knox_connectome_tracers/"

##make any directories that aren't there if needed
mkdir -p "${main_exps_output_dir}"
mkdir -p "${removed_exps_output_dir}"
mkdir -p "${not_in_knox_output_dir}"

declare -A file_to_outdir=(
    ["$csv_file_not_in_knox_exps"]="$not_in_knox_output_dir"
    ["$csv_file_main_exps"]="$main_exps_output_dir"
    ["$csv_file_removed_exps"]="$removed_exps_output_dir"
)

for csv_file in "${!file_to_outdir[@]}"; do
    output_dir="${file_to_outdir[$csv_file]}"

    mapfile -t csv_data < <(tail -n +2 "$csv_file")
    for tracer in "${csv_data[@]}"; do  
        ###get rid of tracer 310207648
        if [ "$tracer" -eq "310207648" ]; then
            continue
        fi
        ##download injection fraction/injection density/projection density from Allen API
        mkdir -p "${output_dir}/${tracer}"

        wget -O "${output_dir}/${tracer}/projection_density_raw.nrrd"  "http://api.brain-map.org/grid_data/download_file/${tracer}?image=projection_density&resolution=50"
        sleep 5 ##can make less time if needed
        python transform_space.py "${output_dir}/$tracer/projection_density_raw.nrrd" "${output_dir}/${tracer}/projection_density_raw.mnc" -v=RAS -w=MICe -x=1
        mincresample -2 \
        "${output_dir}/$tracer/projection_density_raw.mnc"  \
        "${output_dir}/$tracer/projection_density.mnc" \
        -like "${output_dir}/$tracer/projection_density_raw.mnc" 

        ##repeat for injection fraction (if needed)
        wget -O "${output_dir}/${tracer}/injection_fraction_raw.nrrd"  "http://api.brain-map.org/grid_data/download_file/${tracer}?image=injection_fraction&resolution=50" 
        sleep 5 ##can make less time if needed
        python transform_space.py "${output_dir}/$tracer/injection_fraction_raw.nrrd" "${output_dir}/${tracer}/injection_fraction_raw.mnc" -v=RAS -w=MICe -x=1
        mincresample -2 \
        "${output_dir}/$tracer/injection_fraction_raw.mnc"  \
        "${output_dir}/$tracer/injection_fraction.mnc" \
        -like "${output_dir}/$tracer/injection_fraction_raw.mnc" 

        ##repeat for injection density (if needed)
        wget -O "${output_dir}/${tracer}/injection_density_raw.nrrd"  "http://api.brain-map.org/grid_data/download_file/${tracer}?image=injection_density&resolution=50" 
        sleep 5 ##can make less time if needed
        python transform_space.py "${output_dir}/$tracer/injection_density_raw.nrrd" "${output_dir}/${tracer}/injection_density_raw.mnc" -v=RAS -w=MICe -x=1
        mincresample -2 \
        "${output_dir}/$tracer/injection_density_raw.mnc"  \
        "${output_dir}/$tracer/injection_density.mnc" \
        -like "${output_dir}/$tracer/injection_density_raw.mnc"
        done
done

