#install.packages("pacman")
library("pacman")
pacman::p_load(dplyr, readr, readxl, stringr, foreach, doParallel, RMINC)
setwd(".")

####CHANGE THESE IF NEEDED####
###set thresholds for the entire script - can change later as needed##
args <- commandArgs(trailingOnly = TRUE)
inj_thresh <- args[1]
proj_thresh <- args[2]
inj_manual_qc_dir <- "../derivatives/knox_inj/bin0.5/"###should contain QC .csv file and .jpg images @ inj_thresh
proj_manual_qc_dir <- "../derivatives/knox_proj/bin0.1/"###should contain QC .csv file and .jpg images @ proj_thresh

#############################
####################SETUP##########################################

allen_input_dir <- "../preprocessed/allen_template_inputs/"

##load primary/secondary injection sites for each experiment(query.csv)
tracer_sites <- as.data.frame(read_csv(paste0(allen_input_dir, "query.csv"), col_names=TRUE))

##load mask - see if all ABC injection sites are actually in the brain 
mask <- paste0(allen_input_dir,"mask_50um.mnc")
maskVol <- mincGetVolume(mask)

##load in excel file with region name/ontology before mapping regions from our label file to GM/WM/CSF
# Read the Excel file
aba_region_filepath <- paste0(allen_input_dir, "allen_ccfv3_tree_wang_2020_s2.xlsx")

##50 microns
aba_label_filepath <- paste0(allen_input_dir, "AMBA_25um_resampled_50um_int.mnc")
aba_label_file <- round(mincGetVolume(aba_label_filepath))

##both
df <- as.data.frame(read_excel(aba_region_filepath))  # Replace with your actual file path
colnames(df) <- df[1,]
df <- df[2:nrow(df),]

# Convert to numeric if needed
df$second_number <- str_extract(df$structure_id_path, "^/\\d+/\\d+")  # Extract first two numbers
df$second_number <- as.numeric(str_extract(df$second_number, "\\d+$"))  # Extract the second number
df <- df %>%
  mutate(gm_wm_vent = case_when(
    second_number == 8 ~ "GM",
    second_number == 1009 ~ "WM",
    second_number == 73 ~ "ventricles",
    TRUE ~ NA_character_  # For any other values
  ))

#use this conversion to keep only GM regions (screen every voxel)
wm_labels <- df[which(df$gm_wm_vent == "WM") , "structure ID"]
wm_voxels <- which(aba_label_file %in% wm_labels)
csf_labels <- df[which(df$gm_wm_vent == "ventricles") , "structure ID"]
csf_voxels <- which(aba_label_file %in% csf_labels)
gm_labels <- df[which(df$gm_wm_vent == "GM") , "structure ID"]
gm_voxels <- which(aba_label_file %in% gm_labels)

###write csf voxels as a .mnc file for future analysis
vent_volume <- rep(0, length(aba_label_file))
vent_volume[csf_voxels] <- 1
mincWriteVolume(vent_volume, "../preprocessed/allen_template_inputs/vent_voxels_label_file_50um.mnc",
                like=aba_label_filepath)

##output a version of the ventricles that DOES NOT include the third ventricle
##or the cerebral aqueduct (too thin, close to midline)
csf_labels_filt <- csf_labels[which(! csf_labels %in% c("129","140"))]
csf_voxels_filt <- which(aba_label_file %in% csf_labels_filt)
vent_volume <- rep(0, length(aba_label_file))
vent_volume[csf_voxels_filt] <- 1
mincWriteVolume(vent_volume, "../preprocessed/allen_template_inputs/vent_voxels_label_file_50um_filt.mnc",
                like=aba_label_filepath)
##rationale: because the medial ventricle is so thin,it's more susceptible
##to partial volume effects


##load list of Knox connectome tracers
knox_experiments_included <- as.data.frame(read.csv("knox_experiment_csvs/knox_experiments_included.csv"))
###remove missing experiment
knox_experiments_included <- knox_experiments_included[which(knox_experiments_included$id != 310207648),]
tracer_list <- unlist(knox_experiments_included)

##load Knox experiments originally excluded
knox_experiments_orig_excluded <- as.data.frame(read.csv("knox_experiment_csvs/knox_excluded_tracers.csv"))
tracer_list_orig_excluded <- unlist(knox_experiments_orig_excluded)
###################################################

#######################AUTOMATED QC#############################################
csf_voxels <- csf_voxels_filt ###exclude ventricles that are too clse to midline
##STEP 1: see if injection densities are misaligned to manually annotated ROI
tracer_dir <- "../preprocessed/knox_connectome_tracers/"
tracer_dir_orig_removed <- "../preprocessed/knox_connectome_tracers/orig_removed/"

####STEP 2: Screen for OOB/vent voxels in injection AND projection files #####
calculate_size_oob_vent_voxels_thresh <- function(tracer_dir, tracer, 
                                                  inj_thresh, proj_thresh, csf_voxels) {
  proj_file <- tryCatch({
    # Code that may produce an error
    mincGetVolume(paste0(tracer_dir, tracer,"/projection_density_bin",proj_thresh,".mnc"))
  }, error = function(e) {
    # Handle the error
    cat("Could not read file:", conditionMessage(e), "\n")
    return(c(tracer, rep(NA, 8)))
  })
  inj_file <- tryCatch({
    # Code that may produce an error
    mincGetVolume(paste0(tracer_dir, tracer,"/injection_dens_times_frac_bin", inj_thresh,".mnc"))
  }, error = function(e) {
    # Handle the error
    cat("Could not read file:", conditionMessage(e), "\n")
    return(c(tracer, rep(NA, 8)))
  })
  proj_file_vox_list <- which(proj_file > 0)
  proj_file_voxels <- length(proj_file_vox_list)
  ###figure out number of nonzero elements in proj file x mask 
  proj_file_voxels_in_brain <- sum(proj_file*maskVol > 0)
  ###see the difference between these two
  proj_file_voxels_outside_brain <- proj_file_voxels - proj_file_voxels_in_brain
  
  ###figure out number of voxels in ventricles (based on allen region names)
  proj_file_voxels_in_vent <- length(intersect(proj_file_vox_list, csf_voxels))
  
  ##repeat for injection 
  inj_file_vox_list <- which(inj_file > 0)
  inj_file_voxels <- length(inj_file_vox_list)
  ###figure out number of nonzero elements in inj file x mask 
  inj_file_voxels_in_brain <- sum(inj_file*maskVol > 0)
  ###see the difference between these two
  inj_file_voxels_outside_brain <- inj_file_voxels - inj_file_voxels_in_brain
  
  ###figure out number of voxels in ventricles (based on allen region names)
  inj_file_voxels_in_vent <- length(intersect(inj_file_vox_list, csf_voxels))
  
  ##output results
  c(tracer, proj_file_voxels, proj_file_voxels_outside_brain,
    proj_file_voxels_in_vent, 
    inj_file_voxels, inj_file_voxels_outside_brain, 
    inj_file_voxels_in_vent)
}



#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-3) #not to overload your computer
registerDoParallel(cl)
tracer_oob_vent_prop_df <-  foreach(
  tracer = tracer_list, 
  .combine = rbind, 
  .export = c("calculate_size_oob_vent_voxels_thresh", "inj_thresh", "proj_thresh", "tracer_dir", "csf_voxels", "mincGetVolume", "maskVol")
  
) %dopar% {
  tempMatrix <- calculate_size_oob_vent_voxels_thresh(tracer_dir, tracer, inj_thresh, proj_thresh, csf_voxels)
  tempMatrix
}

tracer_oob_vent_prop_df_orig_removed <-  foreach(
  tracer = tracer_list_orig_excluded, 
  .combine = rbind, 
  .export = c("calculate_size_oob_vent_voxels_thresh", "inj_thresh", "proj_thresh", "tracer_dir_orig_removed", "csf_voxels", "mincGetVolume", "maskVol")
  
) %dopar% {
  tempMatrix <- calculate_size_oob_vent_voxels_thresh(tracer_dir_orig_removed, tracer, inj_thresh, proj_thresh, csf_voxels)
  tempMatrix
}
#stop cluster
stopCluster(cl)

colnames(tracer_oob_vent_prop_df) <- c("tracer", "proj_num_vox", "proj_vox_oob", "proj_vox_vent", 
                                       "inj_num_vox", "inj_vox_oob", "inj_vox_vent")

colnames(tracer_oob_vent_prop_df_orig_removed) <- c("tracer", "proj_num_vox", "proj_vox_oob", "proj_vox_vent", 
                                       "inj_num_vox", "inj_vox_oob", "inj_vox_vent")
write.csv(tracer_oob_vent_prop_df, paste0("tables/knox_oob_vent_df_inj",inj_thresh,"_proj",proj_thresh,".csv"))
write.csv(tracer_oob_vent_prop_df_orig_removed, paste0("tables/knox_oob_vent_df_inj",inj_thresh,"_proj",proj_thresh,"_orig_removed.csv"))
