library("dplyr")
library("readr")
library("readxl")
library("stringr")
setwd("/data/chamal/projects/natvik/knox_qc_full_06232025/analysis/")

####CHANGE THESE IF NEEDED####
###set thresholds for the entire script - can change later as needed##
harmonized_qc_dir <- "/data/chamal/projects/natvik/knox_qc_full_06232025/analysis/qc_comparison_with_steph/"


##########combine all the exclusion criteria#################################
##binary matrix documenting tracers + failure reasons 
tracer_removal_df <- data.frame(matrix(ncol=12, nrow=0)) 
colnames(tracer_removal_df) <- c("Tracer", "Auto OOB Proj.", "Auto Vent Proj.", "Auto Large Inj.", "Auto Small Inj.", "Auto Large Proj.", "Auto Small Proj.",
                                 "Manual Nonspecific Inj.", "Manual Cortical Leaking Inj.", 
                                 "Manual Nonspecific Proj.", "Manual Small Proj.", "Manual Misaligned Proj.")
###read in data
zscore_thresh <- 4
tracer_num_vox_oob_vent_df <- as.data.frame(read.csv("./tables/knox_oob_vent_df_inj0.5_proj0.1_1018.csv"))
tracer_num_vox_oob_vent_df <- tracer_num_vox_oob_vent_df[which(tracer_num_vox_oob_vent_df$tracer != 310207648),]

#################automated OOB/ventricular voxels##########################
auto_oob_proj_tracers <- tracer_num_vox_oob_vent_df[which(scale(tracer_num_vox_oob_vent_df$proj_vox_oob) > zscore_thresh),"tracer"]
for(tracer in auto_oob_proj_tracers) {
  tracer_removal_df[nrow(tracer_removal_df)+1, "Tracer"] <- tracer
  tracer_removal_df[nrow(tracer_removal_df), "Auto OOB Proj."] <- 1
}

##things might start overlapping at this point: see if tracer row already exists
auto_vent_proj_tracers <- tracer_num_vox_oob_vent_df[which(scale(tracer_num_vox_oob_vent_df$proj_vox_vent) > zscore_thresh),"tracer"]
for(tracer in auto_vent_proj_tracers) {
  if(tracer %in% tracer_removal_df$Tracer) {
    row_index <- which(tracer_removal_df$Tracer == tracer)
    tracer_removal_df[row_index,"Auto Vent Proj."] <- 1
  } else {
    tracer_removal_df[nrow(tracer_removal_df)+1, "Tracer"] <- tracer
    tracer_removal_df[nrow(tracer_removal_df), "Auto Vent Proj."] <- 1
  }
}



################large injections and projections#######################
auto_large_inj_list <- tracer_num_vox_oob_vent_df[which(scale(tracer_num_vox_oob_vent_df$inj_num_vox) > zscore_thresh),"tracer"]
for(tracer in auto_large_inj_list) {
  if(tracer %in% tracer_removal_df$Tracer) {
    row_index <- which(tracer_removal_df$Tracer == tracer)
    tracer_removal_df[row_index,"Auto Large Inj."] <- 1
  }
  else {
    tracer_removal_df[nrow(tracer_removal_df)+1, "Tracer"] <- tracer
    tracer_removal_df[nrow(tracer_removal_df), "Auto Large Inj."] <- 1
  }
}

auto_large_proj_list <- tracer_num_vox_oob_vent_df[which(scale(tracer_num_vox_oob_vent_df$proj_num_vox) > zscore_thresh),"tracer"]
for(tracer in auto_large_proj_list) {
  if(tracer %in% tracer_removal_df$Tracer) {
    row_index <- which(tracer_removal_df$Tracer == tracer)
    tracer_removal_df[row_index,"Auto Large Proj."] <- 1
  }
  else {
    tracer_removal_df[nrow(tracer_removal_df)+1, "Tracer"] <- tracer
    tracer_removal_df[nrow(tracer_removal_df), "Auto Large Proj."] <- 1
  }
}

#############small injections and projections#########################
####get cutoff values from figure 2###
sorted_df <- tracer_num_vox_oob_vent_df %>%
  arrange(inj_num_vox) %>%
  slice(1:50) %>%
  mutate(Rank = 1:50)

# Add group coloring
sorted_df <- sorted_df %>%
  mutate(ColorGroup = ifelse(Rank <= 4, "Bottom 4", "Others"))

# Get cutoff value
inj_cutoff_val <- max(sorted_df$inj_num_vox[sorted_df$ColorGroup == "Bottom 4"])

##repeat for projection
sorted_df <- tracer_num_vox_oob_vent_df %>%
  arrange(proj_num_vox) %>%
  slice(1:50) %>%
  mutate(Rank = 1:50)

# Add group coloring
sorted_df <- sorted_df %>%
  mutate(ColorGroup = ifelse(Rank <= 2, "Bottom 2", "Others"))

# Get cutoff value
proj_cutoff_val <- max(sorted_df$proj_num_vox[sorted_df$ColorGroup == "Bottom 2"])


auto_small_inj_list <- tracer_num_vox_oob_vent_df[which(tracer_num_vox_oob_vent_df$inj_num_vox <= inj_cutoff_val),"tracer"]
for(tracer in auto_small_inj_list) {
  if(tracer %in% tracer_removal_df$Tracer) {
    row_index <- which(tracer_removal_df$Tracer == tracer)
    tracer_removal_df[row_index,"Auto Small Inj."] <- 1
  }
  else {
    tracer_removal_df[nrow(tracer_removal_df)+1, "Tracer"] <- tracer
    tracer_removal_df[nrow(tracer_removal_df), "Auto Small Inj."] <- 1
  }
}

auto_small_proj_list <- tracer_num_vox_oob_vent_df[which(tracer_num_vox_oob_vent_df$proj_num_vox <= proj_cutoff_val),"tracer"]
for(tracer in auto_small_proj_list) {
  if(tracer %in% tracer_removal_df$Tracer) {
    row_index <- which(tracer_removal_df$Tracer == tracer)
    tracer_removal_df[row_index,"Auto Small Proj."] <- 1
  }
  else {
    tracer_removal_df[nrow(tracer_removal_df)+1, "Tracer"] <- tracer
    tracer_removal_df[nrow(tracer_removal_df), "Auto Small Proj."] <- 1
  }
}



############################MANUAL QC#################################
################load consensus QC ratings from Steph##################

injection_manual_qc_file <- as.data.frame(read.csv(paste0(harmonized_qc_dir,"knox_inj_prod_bin0.5_qc_harmonized.csv"), header=FALSE))

colnames(injection_manual_qc_file)<- c("file", "tracer", "rating", "extra")
inj_fail_large_subcort  <- injection_manual_qc_file[which(injection_manual_qc_file$rating == 1), "tracer"]
inj_fail_cort_leak <- injection_manual_qc_file[which(injection_manual_qc_file$rating == 2), "tracer"]

##filter injections that are poorly defined
tracers_with_no_inj_at_thresh <- injection_manual_qc_file[which(injection_manual_qc_file$rating == 9), "tracer"]

for(tracer in inj_fail_large_subcort) {
  if(tracer %in% tracer_removal_df$Tracer) {
    row_index <- which(tracer_removal_df$Tracer == tracer)
    tracer_removal_df[row_index,"Manual Nonspecific Inj."] <- 1
  }
  else {
    tracer_removal_df[nrow(tracer_removal_df)+1, "Tracer"] <- tracer
    tracer_removal_df[nrow(tracer_removal_df), "Manual Nonspecific Inj."] <- 1
  }
}

for(tracer in inj_fail_cort_leak) {
  if(tracer %in% tracer_removal_df$Tracer) {
    row_index <- which(tracer_removal_df$Tracer == tracer)
    tracer_removal_df[row_index,"Manual Cortical Leaking Inj."] <- 1
  }
  else {
    tracer_removal_df[nrow(tracer_removal_df)+1, "Tracer"] <- tracer
    tracer_removal_df[nrow(tracer_removal_df), "Manual Cortical Leaking Inj."] <- 1
  }
}

######################PROJECTION############################################
projection_manual_qc_file <- as.data.frame(read.csv(paste0(harmonized_qc_dir, "knox_proj_bin0.1_qc_harmonized.csv"), header=FALSE))
colnames(projection_manual_qc_file)<- c("file", "tracer", "rating", "extra")


proj_fail_large_subcort <- projection_manual_qc_file[which(projection_manual_qc_file$rating == 1), "tracer"]
proj_fail_no_spread <- projection_manual_qc_file[which(projection_manual_qc_file$rating == 2), "tracer"]
proj_fail_misaligned <- projection_manual_qc_file[which(projection_manual_qc_file$rating == 3), "tracer"]


for(tracer in proj_fail_large_subcort) {
  if(tracer %in% tracer_removal_df$Tracer) {
    row_index <- which(tracer_removal_df$Tracer == tracer)
    tracer_removal_df[row_index,"Manual Nonspecific Proj."] <- 1
  }
  else {
    tracer_removal_df[nrow(tracer_removal_df)+1, "Tracer"] <- tracer
    tracer_removal_df[nrow(tracer_removal_df), "Manual Nonspecific Proj."] <- 1
  }
}

for(tracer in proj_fail_no_spread) {
  if(tracer %in% tracer_removal_df$Tracer) {
    row_index <- which(tracer_removal_df$Tracer == tracer)
    tracer_removal_df[row_index,"Manual Small Proj."] <- 1
  }
  else {
    tracer_removal_df[nrow(tracer_removal_df)+1, "Tracer"] <- tracer
    tracer_removal_df[nrow(tracer_removal_df), "Manual Small Proj."] <- 1
  }
}

for(tracer in proj_fail_misaligned) {
  if(tracer %in% tracer_removal_df$Tracer) {
    row_index <- which(tracer_removal_df$Tracer == tracer)
    tracer_removal_df[row_index,"Manual Misaligned Proj."] <- 1
  }
  else {
    tracer_removal_df[nrow(tracer_removal_df)+1, "Tracer"] <- tracer
    tracer_removal_df[nrow(tracer_removal_df), "Manual Misaligned Proj."] <- 1
  }
}

# 1) Replace NA with 0
df0 <- tracer_removal_df %>%
  mutate(across(-Tracer, ~ ifelse(is.na(.x), 0, .x)))
write_csv(df0, "tracers_to_remove.csv", col_names=TRUE)

