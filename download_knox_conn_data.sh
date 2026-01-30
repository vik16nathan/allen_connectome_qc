#!/bin/bash
module load minc-toolkit-v2


#######STEP 1: Download Template and Label Files############################
allen_template_dir="../preprocessed/allen_template_inputs/"
mkdir -p "${allen_template_dir}"

##download 50 micron template
wget -O "${allen_template_dir}/average_template_50.nrrd" "https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/average_template/average_template_50.nrrd"
python transform_space.py "${allen_template_dir}/average_template_50.nrrd" "${allen_template_dir}/average_template_50.mnc" -v RAS -w MICe -x 1 --volume_type uint

##download label file
wget -O "${allen_template_dir}/annotation_25.nrrd"  "https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/annotation_25.nrrd"
python transform_space.py "${allen_template_dir}/annotation_25.nrrd" "${allen_template_dir}/annotation_25.mnc" -v RAS -w MICe -x 1 --volume_type uint

#####TODO: fill-in conversion to relabelled label file from Yohan######
######and mapping from relabelled label file --> mask file########

#############################################################################
########STEP 2: Download Experiment Tracer Files#############################
###have these files in the knox_experiments_csvs/ directory (pulled from GitHub)

csv_file_main_exps="./knox_experiment_csvs/knox_experiments_included.csv"
csv_file_removed_exps="./knox_experiment_csvs/knox_excluded_tracers.csv"
###
main_exps_output_dir="../preprocessed/knox_connectome_tracers/"
removed_exps_output_dir="../preprocessed/knox_connectome_tracers/orig_removed/"

# Read CSV into an array, skipping the header
mapfile -t csv_data_main < <(tail -n +2 "$csv_file_main_exps")
mapfile -t csv_data_removed < <(tail -n +2 "$csv_file_removed_exps")


output_dir="../preprocessed/knox_connectome_tracers/"

##make any directories that aren't there if needed
mkdir -p "${main_exps_output_dir}"
mkdir -p "${removed_exps_output_dir}"

declare -A file_to_outdir=(
    ["$csv_file_main_exps"]="$main_exps_output_dir"
    ["$csv_file_removed_exps"]="$removed_exps_output_dir"
)

for csv_file in "${!file_to_outdir[@]}"; do
    output_dir="${file_to_outdir[$csv_file]}"

    mapfile -t csv_data < <(tail -n +2 "$csv_file")
    for tracer in "${csv_data[@]}"; do  
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

