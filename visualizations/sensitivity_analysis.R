library("pacman")
pacman::p_load(dplyr, readr, readxl, stringr, ggrepel, tidyr, ggplot2, ggtext)

####################ARGUMENTS##################################################
args <- commandArgs(trailingOnly = TRUE)
knox_or_oh <- args[1]  # expects "knox" or "oh"
stopifnot(knox_or_oh %in% c("knox", "oh"))

setwd("/data/chamal/projects/natvik/copy_of_allen_qc_reproducible_rerun_cv_04122026/allen_connectome_qc/")
output_supp_fig_dir_knox <- "figures/"
output_supp_fig_dir_oh   <- "figures/oh/"
supp_fig_dir <- if (knox_or_oh == "oh") output_supp_fig_dir_oh else output_supp_fig_dir_knox

#############HELPER FUNCTION: PROCESS CONNECTOMES#############################
process_regionalized_conn_contra_ipsi <- function(regionalized_conn_strength_path) {
  knox_conn_old <- as.data.frame(read_csv(regionalized_conn_strength_path))
  knox_conn_region_numbers <- knox_conn_old[2:nrow(knox_conn_old),1]

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

#############HELPER FUNCTION: REORDER CONNECTOME BY MAJOR DIVISION############
reorder_conn_by_major_division <- function(conn_ipsi, conn_contra, conn_region_numbers,
                                           row_major_divisions, ipsi_col_major_divisions,
                                           contra_col_major_divisions, desired_order) {
  row_index        <- order(match(row_major_divisions,       desired_order))
  col_index        <- order(match(ipsi_col_major_divisions,  desired_order))
  contra_col_index <- order(match(contra_col_major_divisions, desired_order))

  list(
    conn_ipsi                  = conn_ipsi[row_index, col_index],
    conn_contra                = conn_contra[row_index, contra_col_index],
    conn_region_numbers        = conn_region_numbers[row_index],
    row_major_divisions        = row_major_divisions[row_index],
    ipsi_col_major_divisions   = ipsi_col_major_divisions[col_index],
    contra_col_major_divisions = contra_col_major_divisions[contra_col_index]
  )
}

#############HELPER FUNCTION: LOAD + FIND MAJOR DIVISIONS#####################
load_major_divisions_for_conn <- function(conn_ipsi, conn_contra, conn_region_numbers,
                                          major_division_dict, aba_region_labels) {
  row_major_divisions    <- sapply(conn_region_numbers,  find_major_division,
                                   major_division_dict, aba_region_labels)
  ipsi_col_major_divisions  <- sapply(colnames(conn_ipsi),  find_major_division,
                                      major_division_dict, aba_region_labels)
  contra_col_major_divisions <- sapply(colnames(conn_contra), find_major_division,
                                       major_division_dict, aba_region_labels)
  list(row = row_major_divisions,
       ipsi = ipsi_col_major_divisions,
       contra = contra_col_major_divisions)
}

#############HELPER FUNCTION: MAKE HEATMAP####################################
# conn_comparison_list: list of one or more connectomes to correlate against
# conn_full_removal. When length > 1, cell values are averaged across comparisons.
make_corr_major_div_heatmap_generic <- function(conn_full_removal, conn_comparison_list,
                                                row_major_divisions, col_major_divisions,
                                                plot_title = "Spearman Correlation of Connection Strengths") {
  div_order  <- c("Isocortex","OLF","HPF","CTXsp","STR","PAL",
                  "Thal","Hypothal","Midbrain","Pons","Medulla","CB")
  div_levels <- div_order[div_order %in% unique(row_major_divisions)]
  n          <- length(div_levels)

  corr_matrix_sum <- matrix(0,  nrow = n, ncol = n, dimnames = list(div_levels, div_levels))
  corr_matrix_cnt <- matrix(0L, nrow = n, ncol = n, dimnames = list(div_levels, div_levels))

  for (conn_comparison in conn_comparison_list) {
    for (row_div in div_levels) {
      for (col_div in div_levels) {
        row_idx   <- which(row_major_divisions == row_div)
        col_idx   <- which(col_major_divisions == col_div)
        vals_full <- unlist(conn_full_removal[row_idx, col_idx])
        vals_comp <- unlist(conn_comparison[row_idx,   col_idx])
        complete  <- complete.cases(vals_full, vals_comp)
        if (sum(complete) >= 3) {
          r <- cor(vals_comp[complete], vals_full[complete], method = "spearman")
          corr_matrix_sum[row_div, col_div] <- corr_matrix_sum[row_div, col_div] + r
          corr_matrix_cnt[row_div, col_div] <- corr_matrix_cnt[row_div, col_div] + 1L
        }
      }
    }
  }

  corr_matrix_avg <- ifelse(corr_matrix_cnt > 0,
                            corr_matrix_sum / corr_matrix_cnt,
                            NA_real_)

  df_heat           <- as.data.frame(as.table(corr_matrix_avg))
  colnames(df_heat) <- c("RowDiv", "ColDiv", "Value")
  df_heat$RowDiv    <- factor(df_heat$RowDiv, levels = div_levels)
  df_heat$ColDiv    <- factor(df_heat$ColDiv, levels = div_levels)

  ggplot(df_heat, aes(x = ColDiv, y = RowDiv, fill = Value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(is.na(Value), "NA", sprintf("%.2f", Value))), size = 8) +
    scale_fill_gradient(low = "cyan", high = "red", limits = c(0, 1), name = "Spearman r") +
    scale_y_discrete(
      limits = rev(div_levels),
      labels = function(x) paste0("<span style='color:", div_cols[x], "'>", x, "</span>")
    ) +
    scale_x_discrete(
      labels = function(x) paste0("<span style='color:", div_cols[x], "'>", x, "</span>")
    ) +
    coord_fixed() +
    labs(title = plot_title, x = "Target Major Division", y = "Source Major Division") +
    theme_minimal(base_size = 30) +
    theme(
      axis.text.x  = element_markdown(angle = 45, hjust = 1, vjust = 1),
      axis.text.y  = element_markdown(hjust = 1),
      axis.title.x = element_text(margin = margin(t = 15)),
      axis.title.y = element_text(margin = margin(r = 15)),
      plot.title   = element_text(hjust = 0.5)
    )
}

#############LOAD MAJOR DIVISION LABELS########################################
allen_input_dir  <- "../preprocessed/allen_template_inputs/"
allen_tracer_dir <- "../preprocessed/knox_connectome_tracers/"
tracer_sites     <- as.data.frame(read_csv(paste0(allen_input_dir, "knox_tracer_ids_inj_regions_full.csv")))

aba_region_filepath <- paste0(allen_input_dir, "allen_ccfv3_tree_wang_2020_s2.xlsx")
aba_region_labels   <- as.data.frame(read_excel(aba_region_filepath))
colnames(aba_region_labels) <- aba_region_labels[1,]
aba_region_labels   <- aba_region_labels[2:nrow(aba_region_labels),]

major_division_dict <- data.frame(Isocortex=315, OLF=698, HPF=1089, CTXsp=703, STR=477, PAL=803,
                                  Thal=549, Hypothal=1097, Midbrain=313, Pons=771,
                                  Medulla=354, CB=512)

find_major_division <- function(label, major_division_dict, aba_region_labels) {
  sample_string        <- aba_region_labels[which(aba_region_labels[,"structure ID"] == label), "structure_id_path"]
  region_hierarchy_list <- as.numeric(unlist(strsplit(sample_string, "/")))
  colnames(major_division_dict)[which(major_division_dict %in% region_hierarchy_list)]
}

iwant_hex <- c("#a83537","#4fc79c","#63348a","#73c161","#6280d6","#d3a046",
               "#ca78cd","#4c792a","#bb467a","#aab248","#a96126","#ff846b")
div_lvls  <- colnames(major_division_dict)
div_cols  <- setNames(iwant_hex[seq_along(div_lvls)], div_lvls)
desired_order <- colnames(major_division_dict)

###############A: BINARIZATION THRESHOLD SENSITIVITY##########################
inj_thresholds  <- c(0.4, 0.5, 0.6)
proj_thresholds <- c(0.05, 0.1, 0.2)

for (i in seq_along(inj_thresholds)) {
  inj_thresh  <- inj_thresholds[i]
  proj_thresh <- proj_thresholds[i]

  tracer_num_vox_oob_vent_df <- as.data.frame(read.csv(paste0(
    "tables/knox_oob_vent_df_inj", inj_thresh, "_proj", proj_thresh, "_437exps.csv")))
  tracer_num_vox_oob_vent_df <- tracer_num_vox_oob_vent_df[
    which(tracer_num_vox_oob_vent_df$tracer != 310207648),]

  qc_table_removal_full <- as.data.frame(read_csv(paste0(
    "tables/tracers_to_remove_adjboxStats_skew_outliers_auto_inj", inj_thresh, "_proj", proj_thresh, ".csv")))
  oob_outliers  <- qc_table_removal_full[which(qc_table_removal_full[,"Auto OOB Proj."]  == 1), "Tracer"]
  vent_outliers <- qc_table_removal_full[which(qc_table_removal_full[,"Auto Vent Proj."] == 1), "Tracer"]

  long_df <- tracer_num_vox_oob_vent_df %>%
    select(tracer, proj_vox_oob, proj_vox_vent) %>%
    pivot_longer(cols = -tracer, names_to = "metric", values_to = "value") %>%
    mutate(
      facet_label = case_when(
        metric == "proj_vox_oob"  ~ paste0("OOB Voxels, Projection > ",        proj_thresh),
        metric == "proj_vox_vent" ~ paste0("Ventricular Voxels, Projection > ", proj_thresh)
      ),
      type = ifelse(grepl("^inj", metric), "Injection", "Projection")
    )

  long_df <- long_df %>%
    group_by(facet_label) %>%
    mutate(
      bin = cut(value, breaks = seq(min(value), max(value), length.out = 51), include.lowest = TRUE),
      upper_tail = case_when(
        facet_label == paste0("OOB Voxels, Projection > ",        proj_thresh) & tracer %in% oob_outliers  ~ TRUE,
        facet_label == paste0("Ventricular Voxels, Projection > ", proj_thresh) & tracer %in% vent_outliers ~ TRUE,
        TRUE ~ FALSE
      ),
      tracer_label = ifelse(upper_tail, tracer, NA)
    )

  p <- ggplot(long_df, aes(x = value, fill = type)) +
    geom_histogram(bins = 50, color = "black", alpha = 0.7) +
    geom_text_repel(
      data = filter(long_df, upper_tail),
      aes(label = tracer_label), stat = "identity",
      y = 0, nudge_y = 100, size = 8, segment.color = "black"
    ) +
    facet_wrap(~ facet_label, scales = "free", ncol = 1) +
    scale_fill_manual(values = c("Injection" = "#B0E57C", "Projection" = "#F4C2C2")) +
    labs(x = "Voxel Count", y = "Count", fill = "Type") +
    theme_minimal(base_size = 40)

  ggsave(paste0("figures/figure_oob_vent_hist_inj", inj_thresh, "_proj", proj_thresh, ".png"),
         p, width = 24, height = 16, dpi = 300)

  ##########INJ/PROJ VOXEL COUNT SCATTERPLOT###################################
  inj_size_outliers  <- qc_table_removal_full[which(qc_table_removal_full[,"Auto Large Inj."]  == 1), "Tracer"]
  proj_size_outliers <- qc_table_removal_full[which(qc_table_removal_full[,"Auto Large Proj."] == 1), "Tracer"]

  tracer_num_vox_oob_vent_df <- tracer_num_vox_oob_vent_df %>%
    mutate(
      sd_outlier   = tracer %in% inj_size_outliers | tracer %in% proj_size_outliers,
      outlier_type = case_when(sd_outlier ~ "Large Inj/Proj.", TRUE ~ "None"),
      label        = ifelse(outlier_type != "None", tracer, NA)
    )

  p <- ggplot(tracer_num_vox_oob_vent_df, aes(x = inj_num_vox, y = proj_num_vox)) +
    geom_point(aes(color = outlier_type), size = 6) +
    scale_color_manual(values = c("None" = "black", "Large Inj/Proj." = "red")) +
    geom_text_repel(aes(label = label), size = 15, na.rm = TRUE) +
    labs(
      x     = paste0("Inj > ", inj_thresh, " Voxel Count"),
      y     = paste0("Proj > ", proj_thresh, " Voxel Count"),
      color = "Outlier Type"
    ) +
    theme_minimal(base_size = 45)

  ggsave(paste0("figures/figure_large_inj_proj_scatterplot_inj", inj_thresh, "_proj", proj_thresh, ".png"),
         p, width = 30, height = 16, dpi = 300)
}

###############LOAD CONNECTOMES FOR PARTS B AND C#############################
connectome_dir <- "mouse_connectivity_models/paper/figures/model_comparison/output/"

knox_paths <- list(
  full = "mouse_connectivity_models/paper/connectivity/voxel-standard-model/normalized_connection_strength_rebuilt.csv",
  auto = "mouse_connectivity_models/paper/connectivity/voxel-standard-model/normalized_connection_strength_rebuilt_auto.csv",
  lo   = c("mouse_connectivity_models/paper/connectivity/voxel-standard-model/normalized_connection_strength_rebuilt_lo_6.csv",
           "mouse_connectivity_models/paper/connectivity/voxel-standard-model/normalized_connection_strength_rebuilt_lo_8.csv")
)

oh_paths <- list(
  full = paste0(connectome_dir, "homogeneous-standard-model_rebuilt.csv"),
  auto = paste0(connectome_dir, "homogeneous-standard-model_rebuilt_auto.csv"),
  lo   = c(paste0(connectome_dir, "homogeneous-standard-model_rebuilt_lo_6.csv"),
           paste0(connectome_dir, "homogeneous-standard-model_rebuilt_lo_8.csv"))
)

paths <- if (knox_or_oh == "oh") oh_paths else knox_paths

conn_full <- process_regionalized_conn_contra_ipsi(paths$full)
conn_auto <- process_regionalized_conn_contra_ipsi(paths$auto)
conn_lo   <- lapply(paths$lo, process_regionalized_conn_contra_ipsi)

## Extract region numbers from full-removal file
conn_raw         <- as.data.frame(read_csv(paths$full))
conn_region_numbers <- conn_raw[2:nrow(conn_raw), 1]

## Major divisions (computed once, reused for B and C)
divs <- load_major_divisions_for_conn(conn_full[[2]], conn_full[[1]], conn_region_numbers,
                                      major_division_dict, aba_region_labels)

## Reorder full-removal connectome
ord_full <- reorder_conn_by_major_division(conn_full[[2]], conn_full[[1]], conn_region_numbers,
                                           divs$row, divs$ipsi, divs$contra, desired_order)

reorder_one <- function(m, row_divs, col_divs) {
  m[order(match(row_divs, desired_order)), order(match(col_divs, desired_order))]
}

###############B: AUTO-ONLY VS. FULL-REMOVAL##################################
ord_auto_ipsi   <- reorder_one(conn_auto[[2]], divs$row, divs$ipsi)
ord_auto_contra <- reorder_one(conn_auto[[1]], divs$row, divs$contra)

p_b_ipsi <- make_corr_major_div_heatmap_generic(
  ord_full$conn_ipsi, list(ord_auto_ipsi),
  ord_full$row_major_divisions, ord_full$ipsi_col_major_divisions,
  "Spearman Correlation of Connection Strengths\n(Auto-Only vs. Full-Removal) by Major Division Pair"
)
ggsave(paste0(supp_fig_dir, "cor_heatmap_ipsi_auto_full_10b.png"),
       p_b_ipsi, width = 20, height = 16, dpi = 300)

p_b_contra <- make_corr_major_div_heatmap_generic(
  ord_full$conn_contra, list(ord_auto_contra),
  ord_full$row_major_divisions, ord_full$contra_col_major_divisions,
  "Spearman Correlation of Connection Strengths\n(Auto-Only vs. Full-Removal) by Major Division Pair"
)
ggsave(paste0(supp_fig_dir, "cor_heatmap_contra_auto_full_10b.png"),
       p_b_contra, width = 20, height = 16, dpi = 300)

###############C: LOWER OUTLIERS SENSITIVITY (AVERAGED)#######################
lo_ipsi_list   <- lapply(conn_lo, function(ci) reorder_one(ci[[2]], divs$row, divs$ipsi))
lo_contra_list <- lapply(conn_lo, function(ci) reorder_one(ci[[1]], divs$row, divs$contra))

p_c_ipsi <- make_corr_major_div_heatmap_generic(
  ord_full$conn_ipsi, lo_ipsi_list,
  ord_full$row_major_divisions, ord_full$ipsi_col_major_divisions,
  "Avg. Spearman Correlation of Connection Strengths\n(6/8 Lower Outliers vs. 4) by Major Division Pair"
)
ggsave(paste0(supp_fig_dir, "cor_heatmap_ipsi_lo_full_10c.png"),
       p_c_ipsi, width = 20, height = 16, dpi = 300)

p_c_contra <- make_corr_major_div_heatmap_generic(
  ord_full$conn_contra, lo_contra_list,
  ord_full$row_major_divisions, ord_full$contra_col_major_divisions,
  "Avg. Spearman Correlation of Connection Strengths\n(6/8 Lower Outliers vs. 4) by Major Division Pair"
)
ggsave(paste0(supp_fig_dir, "cor_heatmap_contra_lo_full_10c.png"),
       p_c_contra, width = 20, height = 16, dpi = 300)


###########D: MUTUAL INFORMATION IN BINARY CONN. LOSSES ACROSS THRESHOLDS#########
binarize_regionalized_conn_contra_ipsi <- function(knox_conn_contra, knox_conn_ipsi, thresh_prop=0.2) {
  knox_conn_full <- cbind(knox_conn_ipsi, knox_conn_contra)
  unlisted_conn <- unlist(knox_conn_full)
  nonzero_unlisted_conn <- unlisted_conn[which(unlisted_conn > 0)]
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
###########SETUP: TAKEN FROM FIGURE_4.R#############

knox_conn_old <- as.data.frame(read_csv("mouse_connectivity_models/paper/connectivity/voxel-standard-model/normalized_connection_density_original.csv"))
knox_conn_region_numbers <- knox_conn_old[2:nrow(knox_conn_old),1]
####load parcellated Knox/Oh connectomes (old vs. rerun)
knox_conn_contra_old <- read.csv("../derivatives/regionalized_connectomes/knox_conn_contra_old.csv", check.names=FALSE, row.names=1)
knox_conn_ipsi_old   <- read.csv("../derivatives/regionalized_connectomes/knox_conn_ipsi_old.csv", check.names=FALSE, row.names=1)

knox_conn_contra_new <- read.csv("../derivatives/regionalized_connectomes/knox_conn_contra_new.csv", check.names=FALSE, row.names=1)
knox_conn_ipsi_new   <- read.csv("../derivatives/regionalized_connectomes/knox_conn_ipsi_new.csv", check.names=FALSE, row.names=1)

oh_conn_contra_old   <- read.csv("../derivatives/regionalized_connectomes/oh_conn_contra_old_211.csv", check.names=FALSE, row.names=1)
oh_conn_ipsi_old     <- read.csv("../derivatives/regionalized_connectomes/oh_conn_ipsi_old_211.csv", check.names=FALSE, row.names=1)

oh_conn_contra_new   <- read.csv("../derivatives/regionalized_connectomes/oh_conn_contra_new_211.csv", check.names=FALSE, row.names=1)
oh_conn_ipsi_new     <- read.csv("../derivatives/regionalized_connectomes/oh_conn_ipsi_new_211.csv", check.names=FALSE, row.names=1)
knox_conn_ipsi_regions_in_oh <- rownames(knox_conn_ipsi_old)[which(rownames(knox_conn_ipsi_old) %in% rownames(oh_conn_ipsi_old))]
knox_conn_contra_regions_in_oh <- colnames(knox_conn_contra_old)[which(colnames(knox_conn_contra_old) %in% colnames(oh_conn_contra_old))]

connectome_dir <- "mouse_connectivity_models/paper/figures/model_comparison/output/"
oh_conn_old <- as.data.frame(read_csv(paste0(connectome_dir,"homogeneous-standard-model_original.csv")))
oh_conn_region_numbers <- oh_conn_old[2:nrow(oh_conn_old),1]

if(knox_or_oh == "oh") {
  fig_dir <- output_supp_fig_dir_oh
  conn_region_numbers <- oh_conn_region_numbers
  conn_ipsi_old <- oh_conn_ipsi_old
  conn_contra_old <- oh_conn_contra_old
  conn_ipsi_new <- oh_conn_ipsi_new
  conn_contra_new <- oh_conn_contra_new
}

if(knox_or_oh == "knox") {
  fig_dir <- output_supp_fig_dir_knox
  conn_region_numbers <- knox_conn_region_numbers
  conn_ipsi_old <- knox_conn_ipsi_old
  conn_contra_old <- knox_conn_contra_old
  conn_ipsi_new <- knox_conn_ipsi_new
  conn_contra_new <- knox_conn_contra_new
}
###################################################

#####calculate binarization over a range of thresholds#######################
thresh_list <- c(0.15, 0.2, 0.25)

####calculate gains and losses in connectivity at each of the thresholds####

conn_ipsi_loss_list <- list()
conn_contra_loss_list <- list()

i <- 1
for(binary_pct_threshold in thresh_list) {
  ##binarize connectomes
    bin_list <- binarize_regionalized_conn_contra_ipsi(conn_contra_old, conn_ipsi_old, binary_pct_threshold)
    conn_contra_old_bin <- bin_list[[1]]
    conn_ipsi_old_bin <- bin_list[[2]]

    bin_list <- binarize_regionalized_conn_contra_ipsi(conn_contra_new, conn_ipsi_new, binary_pct_threshold)
    conn_contra_new_bin <- bin_list[[1]]
    conn_ipsi_new_bin <- bin_list[[2]]

    ###########ipsilateral connectivity losses##########
    mat_old  <- as.matrix(conn_ipsi_old_bin)
    mat_new  <- as.matrix(conn_ipsi_new_bin)
    mat_diff <- mat_new - mat_old  # +1 gained, -1 lost

    ## Build composite code matrix: 0,1,2,3 (see legend below)
    code <- matrix(NA_integer_, nrow(mat_old), ncol(mat_old))
    code[mat_diff == 0 & mat_old == 0] <- 0   
    code[mat_diff == 0 & mat_old == 1] <- 0   
    code[mat_diff ==  1]               <- 0 
    code[mat_diff == -1]               <- 1   # lost       -> blue
    rownames(code) <- rownames(mat_old)
    colnames(code) <- colnames(mat_old)
    conn_ipsi_loss_list[[i]] <- code
    

    ############contralateral connectivity losses####################
    mat_old  <- as.matrix(conn_contra_old_bin)
    mat_new  <- as.matrix(conn_contra_new_bin)
    mat_diff <- mat_new - mat_old  # +1 gained, -1 lost

    ## Build composite code matrix: 0,1,2,3 (see legend below)
    code <- matrix(NA_integer_, nrow(mat_old), ncol(mat_old))
    code[mat_diff == 0 & mat_old == 0] <- 0   # unchanged 0 -> white
    code[mat_diff == 0 & mat_old == 1] <- 0   # unchanged 1 -> black
    code[mat_diff ==  1]               <- 0   # gained     -> red
    code[mat_diff == -1]               <- 1   # lost       -> blue
    rownames(code) <- rownames(mat_old)
    colnames(code) <- colnames(mat_old)
    conn_contra_loss_list[[i]] <- code

    i <- i + 1

}


#############HELPER FUNCTION: MUTUAL INFORMATION###################
# Returns mutual information (arbitrary units)
# Returns NA if either vector is constant (entropy = 0)
compute_nmi <- function(x, y, eps = 1e-10) {
  stopifnot(length(x) == length(y))
  
  entropy <- function(v) {
    p <- table(v) / length(v)
    p <- p[p > 0]
    -sum(p * log(p))
  }
  
  
  ## Joint distribution from 2x2 contingency table
  joint <- table(x, y) / length(x)
  px    <- rowSums(joint)
  py    <- colSums(joint)
  
  mi <- 0
  for (xi in rownames(joint)) {
    for (yi in colnames(joint)) {
      pxy <- joint[xi, yi]
      if (pxy > eps) mi <- mi + pxy * log(pxy / (px[xi] * py[yi]))
    }
  }
  
  #mi / sqrt(hx * hy) ##optional: normalize MI (less interpretable)
  mi
}

#############HELPER FUNCTION: AVG PAIRWISE NMI HEATMAP#######################
# conn_loss_list: list of binary loss matrices (one per threshold),
#                 already reordered by major division
make_nmi_major_div_heatmap <- function(conn_loss_list,
                                       row_major_divisions, col_major_divisions,
                                       plot_title = "Avg. Pairwise NMI of Connectivity Losses by Major Division Pair") {
  div_order  <- c("Isocortex","OLF","HPF","CTXsp","STR","PAL",
                  "Thal","Hypothal","Midbrain","Pons","Medulla","CB")
  div_levels <- div_order[div_order %in% unique(row_major_divisions)]
  n          <- length(div_levels)
  n_thresh   <- length(conn_loss_list)
  
  ## All pairwise threshold index combinations: (1,2), (1,3), (2,3), ...
  thresh_pairs <- combn(n_thresh, 2, simplify = FALSE)
  
  nmi_matrix_sum <- matrix(0,  nrow = n, ncol = n, dimnames = list(div_levels, div_levels))
  nmi_matrix_cnt <- matrix(0L, nrow = n, ncol = n, dimnames = list(div_levels, div_levels))
  
  for (pair in thresh_pairs) {
    mat_a <- as.matrix(conn_loss_list[[pair[1]]])
    mat_b <- as.matrix(conn_loss_list[[pair[2]]])
    
    for (row_div in div_levels) {
      for (col_div in div_levels) {
        row_idx <- which(row_major_divisions == row_div)
        col_idx <- which(col_major_divisions == col_div)
        
        vec_a <- as.integer(mat_a[row_idx, col_idx])
        vec_b <- as.integer(mat_b[row_idx, col_idx])
        
        ## Drop positions where either matrix has NA
        complete <- complete.cases(vec_a, vec_b)
        if (sum(complete) < 3) next
        
        nmi <- compute_nmi(vec_a[complete], vec_b[complete])
        if (!is.na(nmi)) {
          nmi_matrix_sum[row_div, col_div] <- nmi_matrix_sum[row_div, col_div] + nmi
          nmi_matrix_cnt[row_div, col_div] <- nmi_matrix_cnt[row_div, col_div] + 1L
        }
      }
    }
  }
  
  ## Element-wise mean (NA where no valid pairs)
  nmi_matrix_avg <- ifelse(nmi_matrix_cnt > 0,
                           nmi_matrix_sum / nmi_matrix_cnt,
                           NA_real_)
  
  df_heat           <- as.data.frame(as.table(nmi_matrix_avg))
  colnames(df_heat) <- c("RowDiv", "ColDiv", "Value")
  df_heat$RowDiv    <- factor(df_heat$RowDiv, levels = div_levels)
  df_heat$ColDiv    <- factor(df_heat$ColDiv, levels = div_levels)
  
  ggplot(df_heat, aes(x = ColDiv, y = RowDiv, fill = Value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(is.na(Value), "NA", sprintf("%.2f", Value))), size = 8) +
    scale_fill_gradient(low = "cyan", high = "red", name = "Avg. MI") +
    scale_y_discrete(
      limits = rev(div_levels),
      labels = function(x) paste0("<span style='color:", div_cols[x], "'>", x, "</span>")
    ) +
    scale_x_discrete(
      labels = function(x) paste0("<span style='color:", div_cols[x], "'>", x, "</span>")
    ) +
    coord_fixed() +
    labs(title = plot_title, x = "Target Major Division", y = "Source Major Division") +
    theme_minimal(base_size = 30) +
    theme(
      axis.text.x  = element_markdown(angle = 45, hjust = 1, vjust = 1),
      axis.text.y  = element_markdown(hjust = 1),
      axis.title.x = element_text(margin = margin(t = 15)),
      axis.title.y = element_text(margin = margin(r = 15)),
      plot.title   = element_text(hjust = 0.5)
    )
}

###############D: NMI OF CONNECTIVITY LOSSES ACROSS THRESHOLDS################

## Reorder each loss matrix to match the major-division ordering from ord_full
reorder_loss <- function(loss_mat, row_divs, col_divs) {
  loss_mat[order(match(row_divs, desired_order)), order(match(col_divs, desired_order))]
}

ipsi_loss_list_ord   <- lapply(conn_ipsi_loss_list,   reorder_loss,
                                row_divs = divs$row, col_divs = divs$ipsi)
contra_loss_list_ord <- lapply(conn_contra_loss_list, reorder_loss,
                                row_divs = divs$row, col_divs = divs$contra)

p_d_ipsi <- make_nmi_major_div_heatmap(
  ipsi_loss_list_ord,
  ord_full$row_major_divisions, ord_full$ipsi_col_major_divisions,
  "Avg. Pairwise NMI of Ipsilateral Connectivity Losses\nAcross Binarization Thresholds by Major Division Pair"
)
ggsave(paste0(supp_fig_dir, "nmi_heatmap_ipsi_losses_10d.png"),
       p_d_ipsi, width = 20, height = 16, dpi = 300)

p_d_contra <- make_nmi_major_div_heatmap(
  contra_loss_list_ord,
  ord_full$row_major_divisions, ord_full$contra_col_major_divisions,
  "Avg. Pairwise NMI of Contralateral Connectivity Losses\nAcross Binarization Thresholds by Major Division Pair"
)
ggsave(paste0(supp_fig_dir, "nmi_heatmap_contra_losses_10d.png"),
       p_d_contra, width = 20, height = 16, dpi = 300)
