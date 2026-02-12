#!/bin/bash
# Download Allen Institute connectivity data and transform to MINC format
# This script downloads template files and experiment tracer data from the Allen Brain API

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

# Configuration
API_DELAY=5  # Delay between API requests in seconds

module load minc-toolkit-v2


#######STEP 1: Download Template and Label Files############################
allen_template_dir="../preprocessed/allen_template_inputs/"
mkdir -p "${allen_template_dir}"

echo "Downloading 50 micron template..."
if wget -O "${allen_template_dir}/average_template_50.nrrd" "https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/average_template/average_template_50.nrrd"; then
    echo "Template downloaded successfully"
else
    echo "Error: Failed to download template" >&2
    exit 1
fi

python transform_space.py "${allen_template_dir}/average_template_50.nrrd" "${allen_template_dir}/average_template_50.mnc" -v RAS -w MICe -x 1 --volume_type uint

echo "Downloading label file..."
if wget -O "${allen_template_dir}/annotation_25.nrrd" "https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/annotation_25.nrrd"; then
    echo "Label file downloaded successfully"
else
    echo "Error: Failed to download label file" >&2
    exit 1
fi

python transform_space.py "${allen_template_dir}/annotation_25.nrrd" "${allen_template_dir}/annotation_25.mnc" -v RAS -w MICe -x 1 --volume_type uint

#############################################################################
########STEP 2: Download Experiment Tracer Files#############################

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

# Create output directories
mkdir -p "${main_exps_output_dir}"
mkdir -p "${removed_exps_output_dir}"

declare -A file_to_outdir=(
    ["$csv_file_main_exps"]="$main_exps_output_dir"
    ["$csv_file_removed_exps"]="$removed_exps_output_dir"
)

for csv_file in "${!file_to_outdir[@]}"; do
    output_dir="${file_to_outdir[$csv_file]}"
    
    echo "Processing experiments from: $csv_file"
    mapfile -t csv_data < <(tail -n +2 "$csv_file")
    
    for tracer in "${csv_data[@]}"; do
        echo "Downloading data for tracer: $tracer"
        mkdir -p "${output_dir}/${tracer}"

        # Download projection density
        if wget -O "${output_dir}/${tracer}/projection_density_raw.nrrd" "http://api.brain-map.org/grid_data/download_file/${tracer}?image=projection_density&resolution=50"; then
            echo "  Projection density downloaded"
        else
            echo "  Warning: Failed to download projection density for tracer $tracer" >&2
        fi
        sleep "$API_DELAY"
        
        python transform_space.py "${output_dir}/$tracer/projection_density_raw.nrrd" "${output_dir}/${tracer}/projection_density_raw.mnc" -v=RAS -w=MICe -x=1
        mincresample -2 \
            "${output_dir}/$tracer/projection_density_raw.mnc" \
            "${output_dir}/$tracer/projection_density.mnc" \
            -like "${output_dir}/$tracer/projection_density_raw.mnc"

        # Download injection fraction
        if wget -O "${output_dir}/${tracer}/injection_fraction_raw.nrrd" "http://api.brain-map.org/grid_data/download_file/${tracer}?image=injection_fraction&resolution=50"; then
            echo "  Injection fraction downloaded"
        else
            echo "  Warning: Failed to download injection fraction for tracer $tracer" >&2
        fi
        sleep "$API_DELAY"
        
        python transform_space.py "${output_dir}/$tracer/injection_fraction_raw.nrrd" "${output_dir}/${tracer}/injection_fraction_raw.mnc" -v=RAS -w=MICe -x=1
        mincresample -2 \
            "${output_dir}/$tracer/injection_fraction_raw.mnc" \
            "${output_dir}/$tracer/injection_fraction.mnc" \
            -like "${output_dir}/$tracer/injection_fraction_raw.mnc"

        # Download injection density
        if wget -O "${output_dir}/${tracer}/injection_density_raw.nrrd" "http://api.brain-map.org/grid_data/download_file/${tracer}?image=injection_density&resolution=50"; then
            echo "  Injection density downloaded"
        else
            echo "  Warning: Failed to download injection density for tracer $tracer" >&2
        fi
        sleep "$API_DELAY"
        
        python transform_space.py "${output_dir}/$tracer/injection_density_raw.nrrd" "${output_dir}/${tracer}/injection_density_raw.mnc" -v=RAS -w=MICe -x=1
        mincresample -2 \
            "${output_dir}/$tracer/injection_density_raw.mnc" \
            "${output_dir}/$tracer/injection_density.mnc" \
            -like "${output_dir}/$tracer/injection_density_raw.mnc"
    done
done

echo "All downloads completed successfully"

