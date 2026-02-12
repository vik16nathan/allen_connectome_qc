#!/bin/bash
# Run all steps after manual QC
# This includes: connectome rebuilding, graph theory analyses, and visualizations

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

echo "Starting post-manual QC pipeline..."

# Load required modules
module load anaconda R
module load minc-toolkit-v2
module load RMINC
module load MATLAB

############DOWNLOAD REQUIRED PACKAGES#####################
echo "Setting up Python environment..."
if [[ ! -d ".venv" ]]; then
    uv venv --python 3.9
fi
source .venv/bin/activate

echo "Installing Python dependencies..."
uv pip install -r requirements.txt

# Download mouse_connectivity_models to rebuild connectomes
if [[ ! -d "mouse_connectivity_models" ]]; then
    echo "Downloading mouse_connectivity_models..."
    git clone https://github.com/AllenInstitute/mouse_connectivity_models
    cd mouse_connectivity_models
    uv pip install --no-build-isolation .
    cd ..
else
    echo "mouse_connectivity_models already downloaded"
fi

# Download ggslicer for visualizations
if [[ ! -d "ggslicer" ]]; then
    echo "Downloading ggslicer..."
    git clone https://github.com/yohanyee/ggslicer
else
    echo "ggslicer already downloaded"
fi

# Download Brain Connectivity Toolbox for graph theory
if [[ ! -f "BCT.zip" ]]; then
    echo "Downloading Brain Connectivity Toolbox..."
    FILE_ID="1DmMvRnferBfGe057O-sZwB5jL4j8w1Hu"
    ZIP_NAME="BCT.zip"
    
    # Check if gdown is available
    if ! command -v gdown &> /dev/null; then
        echo "Installing gdown..."
        uv pip install gdown
    fi
    
    # Download and extract
    gdown --id "$FILE_ID" -O "$ZIP_NAME"
    unzip -o "$ZIP_NAME"
else
    echo "Brain Connectivity Toolbox already downloaded"
fi

##################################################################################################
echo "Running overall QC exclusion analysis..."
Rscript overall_qc_exclusion.R

# Apply slight modifications to mouse_connectivity_models (legacy python compatibility)
echo "Patching mouse_connectivity_models..."
cp scorers_updated.py ./mouse_connectivity_models/paper/figures/model_comparison/helpers/scorers.py
./patch_api_timeout.sh  # Increase time limits for file loads

oh_rgn_list="../preprocessed/allen_template_inputs/oh_connectome_rgn_numbers_ccfv3.txt"
knox_region_list="../preprocessed/allen_template_inputs/knox_connectome_rgn_numbers_ccfv3.txt"

# Validate region list files exist
if [[ ! -f "$oh_rgn_list" ]]; then
    echo "Error: Oh region list not found: $oh_rgn_list" >&2
    exit 1
fi

if [[ ! -f "$knox_region_list" ]]; then
    echo "Error: Knox region list not found: $knox_region_list" >&2
    exit 1
fi

# Rebuild connectomes with original list of experiments to exclude from Knox et al., 2018
echo "Rebuilding original connectomes..."
source .venv/bin/activate
./rebuild_oh_connectome.sh "experiments_exclude.json" "${oh_rgn_list}" "original"
./rebuild_oh_connectome.sh "experiments_exclude.json" "${knox_region_list}" "original_291"
./rebuild_knox_connectome.sh "experiments_exclude.json" "original"
./rebuild_knox_connectome.sh "experiments_exclude.json" "original_oh_211_regions"
./rebuild_knox_connectome.sh "experiments_exclude_updated.json" "rebuilt_oh_211_regions"

# Rebuild connectomes with increased list of experiments to exclude post-QC
echo "Rebuilding updated connectomes..."
./rebuild_oh_connectome.sh "experiments_exclude_updated.json" "${oh_rgn_list}" "rebuilt"
./rebuild_oh_connectome.sh "experiments_exclude_updated.json" "${knox_region_list}" "rebuilt_291"
./rebuild_knox_connectome.sh "experiments_exclude_updated.json" "rebuilt"
./rebuild_knox_connectome.sh "experiments_exclude_updated.json" "rebuilt_oh_211_regions"

#################################################
#############GRAPH THEORY ANALYSES###############
echo "Running graph theory analyses..."
mkdir -p ../derivatives/regionalized_connectomes/
Rscript process_flip_regionalized_knox_connectomes.R

mkdir -p ../derivatives/rich_club/
matlab -batch "dbstop if error; rich_club_connectome_bin"

mkdir -p ../derivatives/community_louvain/
matlab -batch "dbstop if error; community_connectome_bin"

###########VISUALIZATIONS##########################
echo "Generating figures..."
mkdir -p figures
mkdir -p figures/oh/

Rscript figure_1_workflow.R
Rscript figure_2.R

mkdir -p "../derivatives/excluded_tracer_aggregate_volumes/"
Rscript figure_3.R

# Visualize each combination of thresholds and connectome models
echo "Generating supplementary figures..."
Rscript supp_fig_major_div_connectomes.R "knox" 0.2
Rscript supp_fig_major_div_connectomes.R "knox" 0.05
Rscript supp_fig_major_div_connectomes.R "oh" 0.2
Rscript supp_fig_major_div_connectomes.R "oh" 0.05

Rscript figure_4.R "knox" 0.2
Rscript figure_4.R "oh" 0.2
Rscript figure_4.R "knox" 0.05
Rscript figure_4.R "oh" 0.05

Rscript figure_5.R

echo "Post-manual QC pipeline completed successfully"