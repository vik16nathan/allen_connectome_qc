library("dplyr")
library("readr")
library("tidyr")
setwd("..")
##load Steph's qc csv
steph_qc_csv <- as.data.frame(read_delim("./rater2/vikram-qc-full.csv", delim="/"))
steph_qc_csv <- steph_qc_csv[,c(1,2,4,5)]
steph_inj_qc <- as.data.frame(cbind(as.numeric(substr(steph_qc_csv[,1], 1,9)), as.numeric(steph_qc_csv[,2])))
colnames(steph_inj_qc) <- c("tracer","rating")
steph_proj_qc <- as.data.frame(cbind(as.numeric(substr(steph_qc_csv[,3], 1,9)), as.numeric(steph_qc_csv[,4])))
colnames(steph_proj_qc) <- c("tracer","rating")

###load Vikram's qc csv's (one for inj, one for proj)
knox_inj_dir <- "../derivatives/knox_inj/bin0.5/"
knox_proj_dir <- "../derivatives/knox_proj/bin0.1/"
vikram_inj_qc <- as.data.frame(read_csv(paste0(knox_inj_dir, "knox_inj_prod_bin0.5_qc.csv"), col_names=FALSE))
vikram_proj_qc <- as.data.frame(read_csv(paste0(knox_proj_dir,"knox_proj_bin0.1_qc.csv"), col_names=FALSE))
vikram_inj_qc_cropped <- as.data.frame(read_csv(paste0(knox_inj_dir,"cropped/knox_inj_prod_bin0.5_qc_cropped.csv"), col_names=FALSE))
vikram_proj_qc_cropped <- as.data.frame(read_csv(paste0(knox_proj_dir,"cropped/knox_proj_bin0.1_qc_cropped.csv"), col_names=FALSE))
vikram_inj_qc_not_in_knox <- as.data.frame(read_csv(paste0(knox_inj_dir,"allen_api_not_in_knox/knox_inj_prod_bin0.5_qc_not_in_knox.csv"), col_names=FALSE))
vikram_proj_qc_not_in_knox <- as.data.frame(read_csv(paste0(knox_proj_dir,"allen_api_not_in_knox/knox_proj_bin0.1_qc_not_in_knox.csv"), col_names=FALSE))

colnames(vikram_inj_qc) <- c("file", "tracer", "rating", "extra")
colnames(vikram_proj_qc) <-  c("file", "tracer", "rating", "extra")
colnames(vikram_inj_qc_cropped) <- c("file", "tracer", "rating", "extra")
colnames(vikram_proj_qc_cropped) <-  c("file", "tracer", "rating", "extra")
colnames(vikram_inj_qc_not_in_knox) <- c("file", "tracer", "rating", "extra")
colnames(vikram_proj_qc_not_in_knox) <-  c("file", "tracer", "rating", "extra")

###replace the 9's (cropped rating) with actual rating
vikram_inj_qc <- rbind(vikram_inj_qc, vikram_inj_qc_not_in_knox)
colnames(vikram_inj_qc) <- c("file", "tracer", "rating", "extra") 
vikram_inj_qc <- left_join(vikram_inj_qc, vikram_inj_qc_cropped, by='tracer')
vikram_inj_qc[which(vikram_inj_qc$rating.x == 9),"rating.x"] <- vikram_inj_qc[which(vikram_inj_qc$rating.x == 9),"rating.y"]

vikram_proj_qc <- rbind(vikram_proj_qc, vikram_proj_qc_not_in_knox)
colnames(vikram_proj_qc) <- c("file", "tracer", "rating", "extra") 
vikram_proj_qc <- left_join(vikram_proj_qc, vikram_proj_qc_cropped, by='tracer')
vikram_proj_qc[which(vikram_proj_qc$rating.x == 2),"rating.x"] <- vikram_proj_qc[which(vikram_proj_qc$rating.x == 2),"rating.y"]

vikram_inj_qc <- vikram_inj_qc[,c("file.x","tracer","rating.x")]
colnames(vikram_inj_qc) <- c("file","tracer", "rating")

vikram_proj_qc <- vikram_proj_qc[,c("file.x","tracer","rating.x")]
colnames(vikram_proj_qc) <- c("file", "tracer", "rating")
##########INJECTION#########################
vikram_steph_inj_qc <- inner_join(vikram_inj_qc, steph_inj_qc, by="tracer")
colnames(vikram_steph_inj_qc) <- c("file", "tracer", "vikram_rating", "steph_rating")
inj_tracers_to_flag <- vikram_steph_inj_qc[which(vikram_steph_inj_qc$vikram_rating != vikram_steph_inj_qc$steph_rating),]

###########PROJECTION########################
vikram_steph_proj_qc <- inner_join(vikram_proj_qc, steph_proj_qc, by="tracer")
colnames(vikram_steph_proj_qc) <- c("file","tracer", "vikram_rating", "steph_rating")
proj_tracers_to_flag <- vikram_steph_proj_qc[which(vikram_steph_proj_qc$vikram_rating != vikram_steph_proj_qc$steph_rating),]

######write out comparison files############
write.table(inj_tracers_to_flag, "harmonized_ratings/vikram_steph_inj_qc_differences.csv", sep=",",row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(proj_tracers_to_flag, "harmonized_ratings/vikram_steph_proj_qc_differences.csv", sep=",",row.names=FALSE, col.names=FALSE, quote=FALSE)


############LOAD HARMONIZED QC TO CREATE FINAL QC CSV########################
harmonized_qc_inj <- read_csv("harmonized_ratings/vikram_steph_inj_qc_differences_harmonized.csv",col_names=FALSE)
harmonized_qc_proj <- read_csv("harmonized_ratings/vikram_steph_proj_qc_differences_harmonized.csv",col_names=FALSE)
colnames(harmonized_qc_inj) <- c("file", "tracer", "vikram_rating", "steph_rating", "harmonized", "comments")
colnames(harmonized_qc_proj) <- c("file", "tracer", "vikram_rating", "steph_rating", "harmonized", "comments")

vikram_inj_qc_final <- vikram_inj_qc
vikram_inj_qc_final[which(vikram_inj_qc_final$tracer %in% harmonized_qc_inj$tracer), "rating"] <- harmonized_qc_inj$harmonized

vikram_proj_qc_final <- vikram_proj_qc
vikram_proj_qc_final[which(vikram_proj_qc_final$tracer %in% harmonized_qc_proj$tracer), "rating"] <- harmonized_qc_proj$harmonized

write.table(vikram_inj_qc_final, "harmonized_ratings/knox_inj_prod_bin0.5_qc_harmonized.csv", sep=",",row.names=FALSE, col.names=FALSE, quote=FALSE, na="")
write.table(vikram_proj_qc_final, "harmonized_ratings/knox_proj_bin0.1_qc_harmonized.csv", sep=",",row.names=FALSE, col.names=FALSE, quote=FALSE, na="")

