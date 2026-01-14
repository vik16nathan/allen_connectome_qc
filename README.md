# Allen Connectome QC

## About

This is the code to reproduce all results in *Experimental Quality Control Induces Changes in Allen Mouse Brain Connectomes* (Nathan et al., 2025). 

## Requirements

All analyses were run on a Linux OS with:
* minc-toolkit-v2
* python 3.9 (to run legacy software in mouse_connectivity_models)
* python packages: see virtual environment file; install packages using uv package manager
* rstudio version 2022.02.3+492: for visualizations

## Analysis
1. ./download_knox_conn_data.sh : downloads all 3D volumes for injection_fraction, injection_density, and projection_density from the Allen API, to the directory specified at the top of the file.
2. ./multiply_threshold_inj_proj.sh: multiplies injection fraction and injection density together (as was done in the code by Knox et al.); applies a binary threshold to the new multiplied injection density, and applies a binary threshold to the existing projection density to make it a "label file" to be used with PyQC
3. 
