# Allen Connectome QC

## About

This is the code to reproduce all results in *Experimental Quality Control Induces Changes in Allen Mouse Brain Connectomes* (Nathan et al., 2025). 

## Dependencies

All analyses were run on a Linux OS with the following packages/scripts downloaded:
* minc-toolkit-v2
* python 3.9 (to run legacy software in mouse_connectivity_models)
* mouse_connectivity_models: https://github.com/AllenInstitute/mouse_connectivity_models -  follow download instructions to install within working directory  
* python packages: see virtual environment file; install packages using uv package manager
* make_slices_images.sh (https://github.com/CoBrALab/make_slice_images) and PyQC for manual QC images (https://github.com/CoBrALab/PyQC)
* Brain Connectivity Toolbox: https://sites.google.com/site/bctnet/ - download within working directory
* rstudio version 2022.02.3+492: for visualizations

## Analysis

### Setup 
1. download_knox_conn_data.sh : downloads all 3D volumes for injection_fraction, injection_density, and projection_density from the Allen API (For definitions of these quantities, see https://community.brain-map.org/t/api-allen-brain-connectivity/2988), to the directory specified at the top of the file. Uses transform_space.py (from Yohan Yee) to convert from Allen PIR --> RAS standard imaging coordinates (see manuscript for more info)
2. multiply_threshold_inj_proj.sh: multiplies injection fraction and injection density together (as was done in the code by Knox et al.); applies a binary threshold to the new multiplied injection density, and applies a binary threshold to the existing projection density to make these into "label files" for PyQC

### Automated QC
3. create_automated_qc_csv.R: calculates number of out-of-brain/ventricular projection voxels and overall injection and projection voxel counts for each experiment, based on the binarized files above. Note: overlaps with all ventricles except the third ventricle and cerebral aqueduct are considered, since these regions are too thin and prone to partial volume effects. Outputs results in "tables" subdirectory. 

### Manual QC
4. make_qc_images_bin_threshold.sh: makes all manual QC images using make_slice_images.sh (modified from https://github.com/CoBrALab/make_slice_images with _arg_row_width="5760" and _arg_row_height="1080" for higher-resolution images. Be sure to change binary thresholds to match whichever thresholds were chosen in (2); we use 0.5 for injection fraction and 0.1 for projection density by default

5. Use PyQC: python PyQC.py to load the GUI; navigate to (1) the injection slice images (2) the projection slice images and separately QC using the numerical ratings outlined in the *Manual Quality Control* section of Nathan et al. Save the injection and projection QC results as separate .csv files.

6. compare_manual_qc_both_raters.R: come up with a table of experiments with discrepancies across both raters, and each rater's rating. This facilitates easier discussion over consensus ratings for these images. 

### Rebuilding Connectomes
7. rebuild_knox_connectome.sh / rebuild_oh_connectome.sh: Before running, create a new experiments_exclude.json file using the overall experiments excluded in manual/automated QC in addition to all experiments originally within experiments_exclude.json (https://github.com/AllenInstitute/mouse_connectivity_models/tree/master/paper). The python scripts build_model_new_excluded.py and build_homogeneous_model_new_excluded.py are almost identical to build_model.py and build_homogeneous_model.py within https://github.com/AllenInstitute/mouse_connectivity_models/tree/master/paper/figures/model_comparison, except the NEW files replacing experiments_exclude.json.

* Note: I also modified get_regional_data() in model_data.py to include the 211 Oh connectome regions, rather than the 291 regions defined in Knox et al. This choice enables direct comparison/downstream application of the rebuilt homogeneous model in place of the Oh connectome.

### Graph Theory Analyses
Code provided by Lizette Herrera-Portillo from her recent manuscript.

8. rich_club_vikram.m
9. community_vikram_updated.m

### Visualizations
7. figure_2.R:
8. figure_3.R:
9. figure_4.R:
10. supplementary_figure_4.R:
11. figure_5.R:

