library("pacman")
pacman::p_load(tidyverse, dplyr, readr, readxl, patchwork)
setwd(".")
allen_input_dir <- "../preprocessed/allen_template_inputs/"

process_regionalized_conn_contra_ipsi <- function(regionalized_conn_strength_path) {
  knox_conn_old <- as.data.frame(read_csv(regionalized_conn_strength_path))
  knox_conn_region_numbers <- knox_conn_old[2:nrow(knox_conn_old),1]
  ###filter into contralateral/ipsilateral
  knox_conn_ipsi <- knox_conn_old[,which(substr(colnames(knox_conn_old),0,4) == "ipsi")]
  colnames(knox_conn_ipsi) <- as.integer(knox_conn_ipsi[1,])
  knox_conn_ipsi <- knox_conn_ipsi[2:nrow(knox_conn_ipsi),]
  rownames(knox_conn_ipsi) <- knox_conn_region_numbers
  
  knox_conn_contra <- knox_conn_old[,which(substr(colnames(knox_conn_old),0,6) == "contra")]
  colnames(knox_conn_contra) <- as.integer(knox_conn_contra[1,])
  knox_conn_contra <- knox_conn_contra[2:nrow(knox_conn_contra),]
  rownames(knox_conn_contra) <- knox_conn_region_numbers
  list(knox_conn_contra, knox_conn_ipsi)
}

binarize_regionalized_conn_contra_ipsi <- function(knox_conn_contra, knox_conn_ipsi, thresh_prop=0.2) {
  knox_conn_full <- cbind(knox_conn_ipsi, knox_conn_contra)
  cutoff <- quantile(unlist(knox_conn_full), 1-thresh_prop) ###keep the TOP THRESH_PROP conns
  
  knox_conn_contra_binary <- as.data.frame(
    lapply(knox_conn_contra, function(col) ifelse(col >= cutoff, 1, 0)))
  rownames(knox_conn_contra_binary) <- rownames(knox_conn_contra)
  colnames(knox_conn_contra_binary) <- colnames(knox_conn_contra)
  
  knox_conn_ipsi_binary <- as.data.frame(
    lapply(knox_conn_ipsi, function(col) ifelse(col >= cutoff, 1, 0)))
  rownames(knox_conn_ipsi_binary) <- rownames(knox_conn_ipsi)
  colnames(knox_conn_ipsi_binary) <- colnames(knox_conn_ipsi)
  
  list(knox_conn_contra_binary, knox_conn_ipsi_binary)
}
####load parcellated Knox connectomes (old vs. rerun)

for(suffix in c("strength", "density")) {
  

  knox_conn_old <- as.data.frame(read_csv(paste0("mouse_connectivity_models/paper/connectivity/voxel-standard-model/normalized_connection_",suffix,"_original.csv")))
  knox_conn_region_numbers <- knox_conn_old[2:nrow(knox_conn_old),1]
  
  knox_conn_contra_ipsi <- process_regionalized_conn_contra_ipsi(paste0("mouse_connectivity_models/paper/connectivity/voxel-standard-model/normalized_connection_",suffix,"_original.csv"))
  knox_conn_contra_old <- knox_conn_contra_ipsi[[1]]
  knox_conn_ipsi_old <- knox_conn_contra_ipsi[[2]]
  
  knox_conn_contra_ipsi <- process_regionalized_conn_contra_ipsi(paste0("mouse_connectivity_models/paper/connectivity/voxel-standard-model/normalized_connection_",suffix,"_rebuilt.csv"))
  knox_conn_contra_new <- knox_conn_contra_ipsi[[1]]
  knox_conn_ipsi_new <- knox_conn_contra_ipsi[[2]]
  
  #######STACK CONNECTOMES TO REPRESENT SYMMETRY ACROSS MIDLINE (ASSUMPTION)#####
  ###fill in missing columns of original contra. matrix with 0's (turn 291 x 286 --> 291 x 291 for stacking)
  knox_conn_contra_old_full <- matrix(nrow=nrow(knox_conn_contra_old), ncol=nrow(knox_conn_contra_old))
  for(j in c(1:nrow(knox_conn_contra_old))) {
    if(rownames(knox_conn_contra_old)[j] %in% colnames(knox_conn_contra_old)) {
      knox_conn_contra_old_full[,j] <- knox_conn_contra_old[,which(colnames(knox_conn_contra_old) == rownames(knox_conn_contra_old)[j])]
    }
    else {
      knox_conn_contra_old_full[,j] <- 0
    }
  }
  colnames(knox_conn_contra_old_full) <- colnames(knox_conn_ipsi_old)
  ###output matrix:
  ###top left: ipsi-ipsi
  ###top right: ipsi-contra
  ###bottom left: contra-ipsi
  ###bottom right: contra-contra
  
  ###build the top of the output matrix
  knox_conn_ipsi_contra_old_full <- cbind(knox_conn_ipsi_old, knox_conn_contra_old_full)
  
  ###build the bottom of the output matrix
  ###the contra-contra connection should be the same as the ipsi-ipsi connection above
  ###and the contra-ipsi connection should be the same as the ipsi-contra conn. above
  knox_conn_contra_ipsi_old_full <- cbind(knox_conn_contra_old_full, knox_conn_ipsi_old)
  
  ###stack reflected regionalized connectomes
  knox_conn_old_bilateral <- rbind(knox_conn_ipsi_contra_old_full, knox_conn_contra_ipsi_old_full)
  write.csv(knox_conn_old_bilateral, paste0("../derivatives/regionalized_connectomes/knox_conn_",suffix,"_old_bilateral.csv"))
  
  
  ###repeat for new
  knox_conn_contra_new_full <- matrix(nrow=nrow(knox_conn_contra_new), ncol=nrow(knox_conn_contra_new))
  for(j in c(1:nrow(knox_conn_contra_new))) {
    if(rownames(knox_conn_contra_new)[j] %in% colnames(knox_conn_contra_new)) {
      knox_conn_contra_new_full[,j] <- knox_conn_contra_new[,which(colnames(knox_conn_contra_new) == rownames(knox_conn_contra_new)[j])]
    }
    else {
      knox_conn_contra_new_full[,j] <- 0
    }
  }
  colnames(knox_conn_contra_new_full) <- colnames(knox_conn_ipsi_new)
  
  knox_conn_ipsi_contra_new_full <- cbind(knox_conn_ipsi_new, knox_conn_contra_new_full)
  knox_conn_contra_ipsi_new_full <- cbind(knox_conn_contra_new_full, knox_conn_ipsi_new)
  
  knox_conn_new_bilateral <- rbind(knox_conn_ipsi_contra_new_full, knox_conn_contra_ipsi_new_full)
  write.csv(knox_conn_new_bilateral, paste0("../derivatives/regionalized_connectomes/knox_conn_",suffix,"_new_bilateral.csv"))
}

