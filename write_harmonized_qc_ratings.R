library("pacman")
pacman::p_load(dplyr, readr, tidyr)


############LOAD HARMONIZED QC TO CREATE FINAL QC CSV########################
setwd("./harmonized_ratings/")
harmonized_qc_inj <- read_csv("vikram_steph_inj_qc_differences_harmonized.csv",col_names=FALSE)
harmonized_qc_proj <- read_csv("vikram_steph_proj_qc_differences_harmonized.csv",col_names=FALSE)
colnames(harmonized_qc_inj) <- c("file", "tracer", "vikram_rating", "steph_rating", "harmonized", "comments")
colnames(harmonized_qc_proj) <- c("file", "tracer", "vikram_rating", "steph_rating", "harmonized", "comments")

vikram_inj_qc <- as.data.frame(read_csv(paste0(knox_inj_dir, "knox_inj_prod_bin0.5_qc.csv"), col_names=FALSE))
vikram_proj_qc <- as.data.frame(read_csv(paste0(knox_proj_dir,"knox_proj_bin0.1_qc.csv"), col_names=FALSE))

colnames(vikram_inj_qc) <- c("file", "tracer", "rating", "extra")
colnames(vikram_proj_qc) <-  c("file", "tracer", "rating", "extra")

vikram_inj_qc_final <- vikram_inj_qc
vikram_inj_qc_final[which(vikram_inj_qc_final$tracer %in% harmonized_qc_inj$tracer), "rating"] <- harmonized_qc_inj$harmonized

vikram_proj_qc_final <- vikram_proj_qc
vikram_proj_qc_final[which(vikram_proj_qc_final$tracer %in% harmonized_qc_proj$tracer), "rating"] <- harmonized_qc_proj$harmonized

write.table(vikram_inj_qc_final, "knox_inj_prod_bin0.5_qc_harmonized.csv", sep=",",row.names=FALSE, col.names=FALSE, quote=FALSE, na="")
write.table(vikram_proj_qc_final, "knox_proj_bin0.1_qc_harmonized.csv", sep=",",row.names=FALSE, col.names=FALSE, quote=FALSE, na="")

