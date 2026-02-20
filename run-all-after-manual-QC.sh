#!/bin/bash                                                                                           
##################################################################################################

#### NOTE: All scripts here run with injection thresholds of 0.5 and projection thresholds of 0.1###
Rscript compare_manual_qc_both_raters.R
Rscript overall_qc_exclusion.R

###very slight modifications needed to mouse_connectivity_models to get things running (legacy python)
###single change to line 5 of scorers.py
cp scorers_updated.py ./mouse_connectivity_models/paper/figures/model_comparison/helpers/scorers.py
./patch_api_timeout.sh ###increase time limits for file loads to allow download of all conn. files

oh_rgn_list="../preprocessed/allen_template_inputs/oh_connectome_rgn_numbers_ccfv3.txt"
knox_region_list="../preprocessed/allen_template_inputs/knox_connectome_rgn_numbers_ccfv3.txt"

###rebuild connectomes with the original list of experiments to exclude from Knox et al., 2018
source .venv/bin/activate
./rebuild_oh_connectome.sh "experiments_exclude.json" ${oh_rgn_list} "original"
./rebuild_oh_connectome.sh "experiments_exclude.json" ${knox_region_list} "original_291"
./rebuild_knox_connectome.sh "experiments_exclude.json" "original"
./rebuild_knox_connectome.sh "experiments_exclude.json" "original_oh_211_regions"
./rebuild_knox_connectome.sh "experiments_exclude_updated.json" "rebuilt_oh_211_regions" 

###rebuild connectomes with increased list of experiments to exclude post-QC
./rebuild_oh_connectome.sh "experiments_exclude_updated.json" ${oh_rgn_list} "rebuilt"
./rebuild_oh_connectome.sh "experiments_exclude_updated.json" ${knox_region_list} "rebuilt_291"
./rebuild_knox_connectome.sh "experiments_exclude_updated.json" "rebuilt" ##automatically writes out 291 regions
./rebuild_knox_connectome.sh "experiments_exclude_updated.json" "rebuilt_oh_211_regions" ##automatically writes out 211 rgns from Oh et al.

#################################################
#############graph theory analyses###############
mkdir -p ../derivatives/regionalized_connectomes/
Rscript process_flip_regionalized_connectomes.R

mkdir -p ../derivatives/rich_club/
matlab -batch "dbstop if error; rich_club_connectome_bin"

mkdir -p ../derivatives/community_louvain/
matlab -batch "dbstop if error; community_connectome_bin"

###########visualizations##########################
mkdir -p figures
mkdir -p figures/oh/

###misc - manually fill in experiments in query.csv with missing injection regions (available on website)
Rscript fill_in_missing_knox_tracer_regions.R

Rscript figure_1_workflow.R
Rscript figure_2.R

mkdir -p "../derivatives/excluded_tracer_aggregate_volumes/"
Rscript figure_3.R

##visualize each combination of thresholds and connectome models
Rscript supp_fig_major_div_connectomes.R "knox" 0.2
Rscript supp_fig_major_div_connectomes.R "knox" 0.05
Rscript supp_fig_major_div_connectomes.R "oh" 0.2
Rscript supp_fig_major_div_connectomes.R "oh" 0.05


Rscript figure_4.R "knox" 0.2
Rscript figure_4.R "oh" 0.2

Rscript figure_4.R "knox" 0.05
Rscript figure_4.R "oh" 0.05

Rscript figure_5.R