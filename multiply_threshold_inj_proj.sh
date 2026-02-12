#!/bin/bash
# Multiply injection fraction and density, then apply binary thresholds
# Usage: ./multiply_threshold_inj_proj.sh <inj_thresh> <proj_thresh>
# Example: ./multiply_threshold_inj_proj.sh 0.5 0.1

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

module load minc-toolkit-v2

# Validate input arguments
if [[ $# -ne 2 ]]; then
    echo "Error: Incorrect number of arguments" >&2
    echo "Usage: $0 <injection_threshold> <projection_threshold>" >&2
    echo "Example: $0 0.5 0.1" >&2
    exit 1
fi

inj_thresh=$1
proj_thresh=$2

# Validate thresholds are numeric
if ! [[ "$inj_thresh" =~ ^[0-9]+\.?[0-9]*$ ]]; then
    echo "Error: Injection threshold must be numeric" >&2
    exit 1
fi

if ! [[ "$proj_thresh" =~ ^[0-9]+\.?[0-9]*$ ]]; then
    echo "Error: Projection threshold must be numeric" >&2
    exit 1
fi

echo "Using thresholds: injection=${inj_thresh}, projection=${proj_thresh}"

csv_file_main_exps="./knox_experiment_csvs/knox_experiments_included.csv"
csv_file_removed_exps="./knox_experiment_csvs/knox_excluded_tracers.csv"

# Validate CSV files exist
if [[ ! -f "$csv_file_main_exps" ]]; then
    echo "Error: CSV file not found: $csv_file_main_exps" >&2
    exit 1
fi

if [[ ! -f "$csv_file_removed_exps" ]]; then
    echo "Error: CSV file not found: $csv_file_removed_exps" >&2
    exit 1
fi

main_exps_output_dir="../preprocessed/knox_connectome_tracers/"
removed_exps_output_dir="../preprocessed/knox_connectome_tracers/orig_removed/"

declare -A file_to_outdir=(
    ["$csv_file_main_exps"]="$main_exps_output_dir"
    ["$csv_file_removed_exps"]="$removed_exps_output_dir"
)

for csv_file in "${!file_to_outdir[@]}"; do
    output_dir="${file_to_outdir[$csv_file]}"
    
    echo "Processing experiments from: $csv_file"
    mapfile -t csv_data < <(tail -n +2 "$csv_file")
    
    for tracer in "${csv_data[@]}"; do
        # Skip tracer 310207648 (missing data)
        if [[ "$tracer" -eq 310207648 ]]; then
            echo "  Skipping tracer $tracer (missing data)"
            continue
        fi
        
        echo "  Processing tracer: $tracer"
        tracer_dir="${output_dir}/${tracer}/"
        
        # Validate input files exist
        if [[ ! -f "${tracer_dir}injection_fraction.mnc" ]]; then
            echo "    Warning: injection_fraction.mnc not found for tracer $tracer" >&2
            continue
        fi
        
        if [[ ! -f "${tracer_dir}injection_density.mnc" ]]; then
            echo "    Warning: injection_density.mnc not found for tracer $tracer" >&2
            continue
        fi
        
        if [[ ! -f "${tracer_dir}projection_density.mnc" ]]; then
            echo "    Warning: projection_density.mnc not found for tracer $tracer" >&2
            continue
        fi
        
        # Multiply injection fraction and density
        mincmath -mult "${tracer_dir}injection_fraction.mnc" \
                       "${tracer_dir}injection_density.mnc" \
                       "${tracer_dir}injection_dens_times_frac.mnc"

        # Apply binary thresholds
        mincmath -2 -seg -const2 "${inj_thresh}" 999999 \
                 "${tracer_dir}injection_dens_times_frac.mnc" \
                 "${tracer_dir}injection_dens_times_frac_bin${inj_thresh}.mnc"
        
        mincmath -2 -seg -const2 "${proj_thresh}" 999999 \
                 "${tracer_dir}projection_density.mnc" \
                 "${tracer_dir}projection_density_bin${proj_thresh}.mnc"
    done
done

echo "Processing completed successfully"
