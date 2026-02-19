#!/bin/bash

output_dir="../preprocessed/allen_template_inputs/"

cd ${output_dir}
##download template file (50 microns)
wget -O "average_template_50.nrrd" "https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/average_template/average_template_50.nrrd"
python ../../analysis/transform_space.py "average_template_50.nrrd" "average_template_50.mnc"  -v RAS -w MICe -x 1

###download label file (25 microns)
wget -O "annotations_25.nrrd" "https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/annotation_25.nrrd"
python ../../analysis/transform_space.py "annotations_25.nrrd" "annotations_25.mnc" -v RAS -w MICe -x 1 --volume_type uint

mincresample -2 -like "average_template_50.mnc" -nearest_neighbour "annotations_25.mnc" "AMBA_25um_resampled_50um_int.mnc"
##resample to 50 um minc volume

##produce brain mask
mincmath -2 -seg -const2 0.5 9999999999 "annotations_25.mnc" "mask_25um.mnc"
mincresample -2 -like "average_template_50.mnc" -nearest_neighbour -keep_real_range "mask_25um.mnc" "mask_50um.mnc"