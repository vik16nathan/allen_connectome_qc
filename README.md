# Allen Connectome QC

## About

This is the code to reproduce all results in *Experimental Quality Control Induces Changes in Allen Mouse Brain Connectomes* (Nathan et al., 2025). 

## Dependencies

All analyses were run on a Linux OS with the following packages/scripts downloaded:
* `minc-toolkit-v2`
* `python 3.9` (to run legacy software in mouse_connectivity_models)
* `mouse_connectivity_models`: https://github.com/AllenInstitute/mouse_connectivity_models -  follow download instructions to install within working directory  
* python packages: see virtual environment file; install packages using `uv` package manager
* `make_slices_images.sh` (https://github.com/CoBrALab/make_slice_images) and PyQC for manual QC images (https://github.com/CoBrALab/PyQC)
* Brain Connectivity Toolbox: https://sites.google.com/site/bctnet/ - download within working directory
* `rstudio` version 2022.02.3+492: for visualizations
    * `install.packages("pacman")` to allow further package management/package downloads in subsequent R scripts
* `ggslicer` for visualizations, within working directory (courtesy of Yohan Yee https://github.com/yohanyee/ggslicer) 


## Quickstart to Reproduce Results
Note that there are two stages: before manual QC (downloading/thresholding/automated QC/making manual QC images) and after manual QC (rebuilding connectomes/visualizing changes/graph theory analyses). To reproduce our results, you can use our manual QC ratings provided within a subdirectory (stored as a .csv file). 

Within working directory:
1. `mkdir preprocessed derivatives`
2. `git clone https://github.com/vik16nathan/allen_connectome_qc/`
3. `cd allen_connectome_qc`

Before running: make sure PyQC.py, mouse_connectivity_models, ggslicer, and Brain Connectivity Toolbox are all downloaded within the `allen_connectome_qc` directory (see above). 

4. `./run-all-before-manual-QC.sh` (steps 1-4 of full analysis)

Manual QC: These require manual adjustment to redo the QC. The final results from this QC are stored within the harmonized_ratings directory within this repo; to reproduce results, skip this step and go to `./run-all-after-manual-QC.sh`.

5. Use `python PyQC.py` and select folders `../derivatives/knox_inj/bin0.5` and `../derivatives/knox_proj/bin0.1` (thresholds used for main results in the paper). Store results in the same directory as the images being QCed.
6. `compare_manual_QC_both_raters.R`: Takes the QC .csv file from another rater and determines disagreements, before deciding on consensus ratings.
7. `write_harmonized_qc_ratings.R`: Writes the consensus QC ratings to a file with the same format as the original QC files (compatible with PyQC).

After Manual QC: 

7. `./run-all-after-manual-QC.sh` (steps 7+ of full analysis)

## Full Analysis

### Setup 
1. `download_knox_conn_data.sh` : Downloads all 3D volumes for injection_fraction, injection_density, and projection_density from the Allen API (For definitions of these quantities, see https://community.brain-map.org/t/api-allen-brain-connectivity/2988), to the directory specified at the top of the file. Uses `transform_space.py` (from Yohan Yee) to convert from Allen PIR --> RAS standard imaging coordinates (see manuscript for more info). Downloads WT tracer experiments that were included and excluded in the connectome from Knox et al. into respective subdirectories. Also downloads template/label files. 
2. `multiply_threshold_inj_proj.sh inj_thresh proj_thresh`: multiplies injection fraction and injection density together (as was done in the code by Knox et al.). Then, uses separate the numerical thresholds (between 0 and 1) provided as command-line arguments to binarize the injection and projection files before generating QC images for PyQC.

### Automated QC
3. `create_automated_qc_csv.R`: calculates number of out-of-brain/ventricular projection voxels and overall injection and projection voxel counts for each experiment, based on the binarized files above. Note: overlaps with all ventricles except the third ventricle and cerebral aqueduct are considered, since these regions are too thin and prone to partial volume effects. Outputs results in "tables" subdirectory. 

### Manual QC
4. `make_qc_images_bin_threshold.sh`: makes all manual QC images using `make_slice_images.sh` (modified from https://github.com/CoBrALab/make_slice_images with _arg_row_width="5760" and _arg_row_height="1080" for higher-resolution images. Be sure to change binary thresholds to match whichever thresholds were chosen in (2); we use 0.5 for injection fraction and 0.1 for projection density by default

5. Use PyQC: `python PyQC.py` to load the GUI; navigate to (1) the injection slice images (2) the projection slice images and separately QC using the numerical ratings outlined in the *Manual Quality Control* section of Nathan et al. Save the injection and projection QC results as separate .csv files.

6. `compare_manual_qc_both_raters.R`: come up with a table of experiments with discrepancies across both raters, and each rater's rating. This facilitates easier discussion over consensus ratings for these images.

7. `overall_qc_exclusion.R`: Creates a binary matrix of all experiments to remove (tracers_to_remove.csv), with rows representing experiment IDs and columns representing all possible failure modes. An entry of 1 in row i and column j indicates that experiment i failed in mode j. Use the experiment IDs to update experiments_to_exclude.json. Uses outlier detection for right-skewed count data for automated QC (https://rpubs.com/dario-dellevedove/601843). 

### Rebuilding Connectomes
8. `rebuild_knox_connectome.sh / rebuild_oh_connectome.sh`: Before running, create a new experiments_exclude.json file using the overall experiments excluded in manual/automated QC in addition to all experiments originally within `experiments_exclude.json` (https://github.com/AllenInstitute/mouse_connectivity_models/tree/master/paper). The python scripts build_model_new_excluded.py and build_homogeneous_model_new_excluded.py are almost identical to build_model.py and build_homogeneous_model.py within https://github.com/AllenInstitute/mouse_connectivity_models/tree/master/paper/figures/model_comparison, except the new file replacing `experiments_exclude.json`.

* Note: I also modified `get_regional_data()` in `model_data.py` (https://github.com/AllenInstitute/mouse_connectivity_models/tree/master/paper/figures/model_comparison/helpers) to include the 211 Oh connectome regions (minus SUBd and SUBv) that are included in the Allen CCFv3 labels, rather than the 291 regions defined in Knox et al. This choice enables direct comparison/downstream application of the rebuilt homogeneous model in place of the Oh connectome.

### Graph Theory Analyses
Code for (10) and (11) provided by Lizette Herrera-Portillo from her recent manuscript (https://doi.org/10.1101/2025.05.10.653278) 

9. `process_flip_regionalized_knox_connectomes.R`: Since connectivity is only defined for right-hemisphere injections, but we need a full n x n matrix for all graph theoretical analyses, we assume that connectivity is symmetric across the midline to create a square connectivity matrix. We assume connectivity to the contralateral regions is the same as connectivity to the ipsilateral regions, but flipped across the midline (ex: if we know that the right DG is connected to the left CA1, we assume the left DG is connected to the right CA1. This assumption has been verified in the literature).
10. `rich_club_connectome_bin.m`: Calculates rich club coefficient and topological rich club using the top 20% of region x region connections above. See *Rich Club Algorithm* in Supplementary Methods for more info.
11. `community_connectome_bin.m`: Calculates Louvain community identity for the top 20% of region x region connections above; uses `call_randind.py` to calculate the mutual information between community assignments at various "resolution parameters" (gammas) to identify a stable community assignment. See *Community Detection Algorithm* in Supplementary Methods for more info.

### Visualizations
12. `figure_1_workflow.R`: Makes Methods workflow (Figure 1)
13. `figure_2.R`: Makes individual figure plots for Figure 2
14. `figure_3.R`: Makes figure plots for Figure 3, Supplementary Figure 1, and Supplementary Figure 2.
15. `supp_fig_major_div_connectomes.R`: Makes Supplementary Figures 3, 5, 6, and 8 based on the choice of connectome/binarization threshold.
16. `figure_4.R`: Makes figure plots visualizing the changes in the full regionalized voxel model (Figure 4) or homogeneous model (Supplementary Figure 4) based on the choice of Oh vs. Knox connectome. Also makes Supplementary Figure 7 at a 5% binarization threshold.
17. `figure_5.R`: Visualizes graph theory results in Figure 5. Can also make Supplementary Figure 9 if using normalized connection density metric instead of normalized connection strength

