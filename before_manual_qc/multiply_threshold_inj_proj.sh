#!/bin/bash
set -euo pipefail
cd ..
#GOAL: get "true" injection density in annotated area by multiplying injection fraction with injection density
#use a binary threshold to keep only inj/proj files over the threshold
if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <inj_thresh> <proj_thresh>" >&2
    exit 1
fi

inj_thresh=$1
proj_thresh=$2
###note: can change threshold as needed
#csv_file="knox_experiments_included.csv" ###add 20 tracers with no injections 
# Read CSV into an array, skipping the header
#mapfile -t csv_data < <(tail -n +2 "$csv_file")

csv_file_main_exps="./knox_experiment_csvs/knox_experiments_included.csv"
csv_file_removed_exps="./knox_experiment_csvs/knox_excluded_tracers.csv"
csv_file_not_in_knox_exps="./knox_experiment_csvs/wt_experiments_allen_api_not_knox.csv"

###
main_exps_output_dir="../preprocessed/knox_connectome_tracers/"
removed_exps_output_dir="../preprocessed/knox_connectome_tracers/orig_removed/"
not_in_knox_output_dir="../preprocessed/knox_connectome_tracers/allen_api_not_in_knox/"
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

        if [ "$tracer" -eq 310207648 ]; then
        ##missing data
            continue
        fi
        tracer_dir="${output_dir}/${tracer}/"
        mincmath -mult "${tracer_dir}injection_fraction.mnc" "${tracer_dir}injection_density.mnc" "${tracer_dir}injection_dens_times_frac.mnc"

        ##binary thresholds
        mincmath -2 -seg -const2 ${inj_thresh} 999999 "${tracer_dir}injection_dens_times_frac.mnc" "${tracer_dir}injection_dens_times_frac_bin${inj_thresh}.mnc" 
        mincmath -2 -seg -const2 ${proj_thresh} 999999 "${tracer_dir}projection_density.mnc" "${tracer_dir}projection_density_bin${proj_thresh}.mnc" 

        ##one percent thresholds

    done
done
