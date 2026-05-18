library(tidyverse)
setwd("..")

# --- Load tables ---
new_table <- "mouse_connectivity_models/paper/figures/model_comparison/output/cv_results_voxel-standard_homogeneous-standard.csv"

#####downloaded from https://github.com/AllenInstitute/mouse_connectivity_models/blob/master/paper/figures/model_comparison/output/cv_results_voxel-standard_homogeneous-standard.csv
#####renamed to cv_results_voxel-standard_homogeneous-standard-original.csv
old_table <- "mouse_connectivity_models/paper/figures/model_comparison/output/cv_results_voxel-standard_homogeneous-standard-original.csv"

new_df <- read.csv(new_table, header = TRUE)
old_df <- read.csv(old_table, header = TRUE)

# --- Drop row 1 ---
new_df <- new_df[-1, ]
old_df <- old_df[-1, ]

# --- Identify numeric columns (all except the first two: Major division, Model) ---
numeric_cols <- colnames(new_df)[3:ncol(new_df)]

# --- Build the output table ---
out_df <- new_df

for (col in numeric_cols) {
  new_vals <- as.numeric(new_df[[col]]) * 100
  old_vals <- as.numeric(old_df[[col]]) * 100
  diff_vals <- new_vals - old_vals

  out_df[[col]] <- ifelse(
    !is.na(new_vals),
    sprintf("%d%% (%.1f)", round(new_vals), diff_vals),
    ""
  )
}

# --- Rename columns ---
colnames(out_df) <- c("Major division", "Model",
                      "Voxel MSErel", "Voxel MSErel, Training",
                      "Region MSErel", "Region MSErel, Training",
                      "Region PTP", "Region PTP, Training")

# --- Relabel Model names ---
out_df[["Model"]] <- recode(out_df[["Model"]],
  "voxel-standard"       = "Voxel",
  "homogeneous-standard" = "Homogeneous"

)
# --- Write output ---
write.csv(out_df,
          file = "tables/cv_results_comparison.csv",
          row.names = FALSE)

print(out_df)