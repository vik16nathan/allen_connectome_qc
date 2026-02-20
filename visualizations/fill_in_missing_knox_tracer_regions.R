library("dplyr")
library("readr")
library("readxl")

setwd("..")
allen_input_dir <- "../preprocessed/allen_template_inputs/"
knox_experiments_included <- as.data.frame(read.csv("./knox_experiment_csvs/knox_experiments_included.csv"))
allen_api_exps <- as.data.frame(read.csv("./knox_experiment_csvs/wt_experiments_allen_api_not_knox.csv"))
knox_experiments_included_full <- unlist(c(knox_experiments_included, allen_api_exps))
knox_experiments_included_full <- data.frame(knox_experiments_included_full)
colnames(knox_experiments_included_full) <- c("id")
tracer_sites <- as.data.frame(read_csv(paste0(allen_input_dir, "query.csv"), col_names=TRUE))

aba_region_filepath <- paste0(allen_input_dir, "allen_ccfv3_tree_wang_2020_s2.xlsx")
aba_region_labels <- as.data.frame(read_excel(aba_region_filepath))    
colnames(aba_region_labels) <- aba_region_labels[1,]
aba_region_labels <- aba_region_labels[2:nrow(aba_region_labels),]


###repeat the same process for the included experiments
knox_joined_full <- knox_experiments_included_full %>%
  left_join(tracer_sites %>% select(id, `structure-id`),
            by = c("id" = "id")) 

missing_ids <- knox_joined_full[which(is.na(knox_joined_full$`structure-id`)),"id"]
print(knox_joined_full[which(knox_joined_full$id == 638314843),])
print(missing_ids)
###NOTE: 310207648 is truly missing

###fill in experiments manually from ABC web viewer
knox_joined_full[which(knox_joined_full$id == 112670853), "structure-id"] <- aba_region_labels[which(aba_region_labels$abbreviation == "MOs"),"structure ID"]
knox_joined_full[which(knox_joined_full$id == 100142354), "structure-id"] <- aba_region_labels[which(aba_region_labels$abbreviation == "POST"),"structure ID"]
knox_joined_full[which(knox_joined_full$id == 638314843), "structure-id"] <- aba_region_labels[which(aba_region_labels$abbreviation == "VISp"),"structure ID"]

write_csv(knox_joined_full,paste0(allen_input_dir, "knox_tracer_ids_inj_regions_full.csv"))
