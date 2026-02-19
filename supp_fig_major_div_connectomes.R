library("pacman")
pacman::p_load(tidyverse, pheatmap, dplyr, readr, readxl, patchwork, ggrepel, scales)


###SET DIRECTORY PATHS AND HYPERPARAMETERS
setwd(".")
allen_input_dir <- "../preprocessed/allen_template_inputs/"
source("./ggslicer/plotting_functions/plotting_functions.R")
source("./ggslicer/plotting_functions/plotting_functions_labels.R")
output_supp_fig_dir_knox <- "figures/"
output_supp_fig_dir_oh <- "figures/oh/"

args <- commandArgs(trailingOnly=TRUE)
knox_or_oh <- args[1]
binary_pct_threshold <- as.numeric(args[2])

###load major division colors
# install.packages("ggtext") # if needed
library(ggtext)

##palette##
iwant_hex <- c("#a83537","#4fc79c","#63348a","#73c161","#6280d6","#d3a046",
               "#ca78cd","#4c792a","#bb467a","#aab248","#a96126","#ff846b")

# --- Build color map for major_division ---
# If you also have fmd/fmd_col, you can include them in the union as you showed.
# If not, we just use the levels present in tracer_joined:
# --- Fixed legend/order from major_division_dict ---
# Preserve the column order exactly as written in your dict
major_division_dict <- data.frame(Isocortex=315, OLF=698, HPF=1089, CTXsp=703, STR=477, PAL=803,
                                  Thal=549, Hypothal=1097, Midbrain=313, Pons=771,
                                  Medulla=354, CB=512)

div_lvls <- colnames(major_division_dict)

# Map colors to levels in dict order
div_cols <- setNames(iwant_hex[seq_along(div_lvls)], div_lvls)

###organize region numbers --> names --> major subdivisions
#load the dictionary from label numbers <-- --> region acronyms
aba_region_filepath <- paste0(allen_input_dir, "allen_ccfv3_tree_wang_2020_s2.xlsx")
aba_region_labels <- as.data.frame(read_excel(aba_region_filepath))  
colnames(aba_region_labels) <- aba_region_labels[1,]
aba_region_labels <- aba_region_labels[2:nrow(aba_region_labels),]

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
knox_conn_old <- as.data.frame(read_csv("mouse_connectivity_models/paper/connectivity/voxel-standard-model/normalized_connection_strength_original.csv"))
knox_conn_region_numbers <- knox_conn_old[2:nrow(knox_conn_old),1]

knox_conn_contra_ipsi <- process_regionalized_conn_contra_ipsi("mouse_connectivity_models/paper/connectivity/voxel-standard-model/normalized_connection_strength_original.csv")
knox_conn_contra_old <- knox_conn_contra_ipsi[[1]]
knox_conn_ipsi_old <- knox_conn_contra_ipsi[[2]]

knox_conn_contra_ipsi <- process_regionalized_conn_contra_ipsi("mouse_connectivity_models/paper/connectivity/voxel-standard-model/normalized_connection_strength_rebuilt.csv")
knox_conn_contra_new <- knox_conn_contra_ipsi[[1]]
knox_conn_ipsi_new <- knox_conn_contra_ipsi[[2]]

###repeat for Oh connectomes#############################
####load parcellated Oh connectomes (old vs. rerun)
##NOTE: changed directories to rerun w/ 211 regions for OHBM abstract (more similar to original 213 regions)
##not the same as taking the 291 mesoscale regions --> looking at original OH regions
connectome_dir <- "mouse_connectivity_models/paper/figures/model_comparison/output/"
oh_conn_old <- as.data.frame(read_csv(paste0(connectome_dir,"homogeneous-standard-model_original.csv")))
oh_conn_region_numbers <- oh_conn_old[2:nrow(oh_conn_old),1]

oh_conn_contra_ipsi <- process_regionalized_conn_contra_ipsi(paste0(connectome_dir,"homogeneous-standard-model_original.csv"))
oh_conn_contra_old <- oh_conn_contra_ipsi[[1]]
oh_conn_ipsi_old <- oh_conn_contra_ipsi[[2]]

oh_conn_contra_ipsi <- process_regionalized_conn_contra_ipsi(paste0(connectome_dir,"homogeneous-standard-model_rebuilt.csv"))
oh_conn_contra_new <- oh_conn_contra_ipsi[[1]] 
oh_conn_ipsi_new <- oh_conn_contra_ipsi[[2]]

#################TOGGLE BETWEEN KNOX/OH CONNECTOMES HERE#######################
##original: don't exclude any of the regions

#knox_conn_ipsi_regions_in_oh <- rownames(knox_conn_ipsi_old)[which(rownames(knox_conn_ipsi_old) %in% rownames(oh_conn_ipsi_old))]
#knox_conn_contra_regions_in_oh <- colnames(knox_conn_contra_old)[which(colnames(knox_conn_contra_old) %in% colnames(oh_conn_contra_old))]

if(knox_or_oh == "oh") {
  supp_fig_dir <- output_supp_fig_dir_oh
  conn_region_numbers <- oh_conn_region_numbers
  conn_ipsi_old <- oh_conn_ipsi_old
  conn_contra_old <- oh_conn_contra_old
  conn_ipsi_new <- oh_conn_ipsi_new
  conn_contra_new <- oh_conn_contra_new
}

if(knox_or_oh == "knox") {
  supp_fig_dir <- output_supp_fig_dir_knox
  conn_region_numbers <- knox_conn_region_numbers
  conn_ipsi_old <- knox_conn_ipsi_old
  conn_contra_old <- knox_conn_contra_old
  conn_ipsi_new <- knox_conn_ipsi_new
  conn_contra_new <- knox_conn_contra_new
}

######################START HERE################################################
###########reorganize all functions here to run on Oh vs. Knox connectomes######
###organize region numbers --> names --> major subdivisions
#load the dictionary from label numbers <-- --> region acronyms
aba_region_filepath <- paste0(allen_input_dir, "allen_ccfv3_tree_wang_2020_s2.xlsx")
aba_region_labels <- as.data.frame(read_excel(aba_region_filepath))  # Replace with your actual file path
colnames(aba_region_labels) <- aba_region_labels[1,]
aba_region_labels <- aba_region_labels[2:nrow(aba_region_labels),]

###major division dictionary (depth 4/5 in Allen region tree)
major_division_dict <- data.frame(Isocortex=315, OLF=698, HPF=1089, CTXsp=703, STR=477, PAL=803,
                                  Thal=549, Hypothal=1097, Midbrain=313, Pons=771,
                                  Medulla=354, CB=512)

find_major_division <- function(label, major_division_dict, aba_region_labels) {
  sample_string <- aba_region_labels[which(aba_region_labels[,"structure ID"] == label),"structure_id_path"]
  region_hierarchy_list <- as.numeric(unlist(strsplit(sample_string, "/")))
  colnames(major_division_dict)[which(major_division_dict %in% region_hierarchy_list)]
}


#############FIND MAJOR DIVISIONS FOR HEATMAP#########################
row_major_divisions <- c()
for(label in conn_region_numbers) {
  row_major_divisions <- c(row_major_divisions, find_major_division(label, major_division_dict, aba_region_labels))
}

ipsi_col_major_divisions <- c()
for(label in colnames(conn_ipsi_new)) {
  ipsi_col_major_divisions <- c(ipsi_col_major_divisions, find_major_division(label, major_division_dict, aba_region_labels))
}

contra_col_major_divisions <- c()
for(label in colnames(conn_contra_new)) {
  contra_col_major_divisions <- c(contra_col_major_divisions, find_major_division(label, major_division_dict, aba_region_labels))
}

###enforce SAME major division/region ordering across connectomes############
desired_order <- colnames(major_division_dict)
row_order <- match(row_major_divisions, desired_order)
col_order <- match(ipsi_col_major_divisions, desired_order)
row_index <- order(row_order)
row_major_divisions <- row_major_divisions[row_index]

col_index <- order(col_order)
ipsi_col_major_divisions <- ipsi_col_major_divisions[col_index]

conn_ipsi_old <- conn_ipsi_old[row_index, col_index]
conn_ipsi_new <- conn_ipsi_new[row_index, col_index]
conn_region_numbers <- conn_region_numbers[row_index]

contra_col_order <- match(contra_col_major_divisions, desired_order)
contra_col_index <- order(contra_col_order)
contra_col_major_divisions <- contra_col_major_divisions[contra_col_index]

conn_contra_old <- conn_contra_old[row_index, contra_col_index]
conn_contra_new <- conn_contra_new[row_index, contra_col_index]

##binarize connectomes
bin_list <- binarize_regionalized_conn_contra_ipsi(conn_contra_old, conn_ipsi_old, binary_pct_threshold)
conn_contra_old_bin <- bin_list[[1]]
conn_ipsi_old_bin <- bin_list[[2]]

bin_list <- binarize_regionalized_conn_contra_ipsi(conn_contra_new, conn_ipsi_new, binary_pct_threshold)
conn_contra_new_bin <- bin_list[[1]]
conn_ipsi_new_bin <- bin_list[[2]]

################################################################################
##########Supplementary supp_figure: Distribution of Connection Strengths############

########change these lines to work with Knox vs. Oh connectome#############
# Function to replace outliers beyond 3*SD with max non-outlier
replace_upper_outliers <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sigma <- sd(x, na.rm = TRUE)
  upper <- mu + 3 * sigma
  lower <- mu - 3 * sigma
  
  # Max non-outlier value
  max_non_outlier <- max(x[x <= upper & x >= lower], na.rm = TRUE)
  min_non_outlier <- min(x[x <= upper & x >= lower], na.rm = TRUE)
  
  # Replace outliers
  print(paste("Number of upper outliers:", sum(x > upper, na.rm=TRUE)))
  print(paste("Max non outlier:", max_non_outlier))
  x[x > upper] <- max_non_outlier
  
  return(x)
}


replace_upper_lower_outliers <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sigma <- sd(x, na.rm = TRUE)
  upper <- mu + 3 * sigma
  lower <- mu - 3 * sigma
  
  # Max non-outlier value
  max_non_outlier <- max(x[x <= upper & x >= lower], na.rm = TRUE)
  min_non_outlier <- min(x[x <= upper & x >= lower], na.rm = TRUE)
  
  # Replace outliers
  print(paste("Number of upper outliers:", sum(x > upper, na.rm=TRUE)))
  print(paste("Max non outlier:", max_non_outlier))
  x[x > upper] <- max_non_outlier
  
  print(paste("Number of lower outliers:", sum(x < lower, na.rm=TRUE)))
  print(paste("Min non outlier:", min_non_outlier))
  x[x < lower] <- min_non_outlier
  
  return(x)
}


# ---- Build a tidy dataframe ----
vals_old <- as.numeric(unlist(c(conn_ipsi_old,  conn_contra_old)))
vals_old <- vals_old[which(vals_old > 0)]
vals_old <- replace_upper_lower_outliers(vals_old)
vals_old <- vals_old/max(vals_old)

vals_new <- as.numeric(unlist(c(conn_ipsi_new,  conn_contra_new)))
vals_new <- vals_new[which(vals_new > 0)]
vals_new <- replace_upper_lower_outliers(vals_new)
vals_new <- vals_new/max(vals_new)

plot_hist_conn_strength <- function(vals_old, vals_new, output_dir) {
  ##get rid of zeros

  df <- data.frame(
    value = c(vals_old, vals_new),
    group = factor(rep(c("old","new"),
                       times = c(length(vals_old), length(vals_new))),
                   levels = c("old","new"))
  )
  
  # --- log-transform values ---
  df$value_log <- log(df$value + 1e-20)
  
  # --- compute medians for each group ---
  medians <- df %>%
    group_by(group) %>%
    summarise(median_val = median(value_log, na.rm = TRUE)) %>%
    ungroup()
  
  # --- compute histogram counts per group to position labels ---
  bins <- 50
  hist_data <- df %>%
    group_by(group) %>%
    do({
      h <- hist(.$value_log, breaks = bins, plot = FALSE)
      data.frame(mid = h$mids, count = h$counts)
    }) %>%
    ungroup()
  
  # --- get approximate y position for labels (above the nearest histogram bin) ---
  label_df <- medians %>%
    rowwise() %>%
    mutate(
      y = hist_data %>%
        filter(group == group) %>%
        slice(which.min(abs(mid - median_val))) %>%
        pull(count),
      y = y + max(hist_data$count) * 0.03  # add small vertical offset
    ) %>%
    ungroup()
  
  # --- plot ---
  rng <- range(df$value_log, na.rm = TRUE)
  
  # Shared color palette
  cols <- c(old = "#1f77b4", new = "#ff7f0e")
  
  p <- ggplot(df, aes(x = value_log)) +
    # Histogram: fill color by group, no border outlines
    geom_histogram(aes(fill = group),
                   position = "identity", alpha = 0.45, bins = bins, color = NA) +
    # Median lines: color by group (uses the same palette)
    geom_vline(data = medians,
               aes(xintercept = median_val, color = group),
               linetype = "dashed", size = 1.2, show.legend = FALSE) +
    # Median labels: also use same group colors
    geom_text_repel(
      data = label_df,
      aes(x = median_val, y = y,
          label = paste0("Median = ", round(median_val, 3)),
          color = group),
      size = 10,
      fontface = "bold",
      direction = "y",
      segment.color = "grey50",
      show.legend = FALSE
    ) +
    # Apply consistent color/fill scales
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    coord_cartesian(xlim = rng) +
    labs(
      title = "Distribution of Logged, Regionalized Connections in Old vs. New Connectome",
      x = "ln(Conn. Strength/max(Conn. Strength) + 1e-20)",
      y = "Number of Regionalized Connections",
      fill = NULL
    ) +
    theme_minimal(base_size = 30) +
    theme(legend.position = "top")
  
  ggsave(paste0(output_dir, "supp_supp_fig_hist_conn_strength_median.png"), p,
         width = 32, height = 16, dpi = 300)
  
}

plot_hist_conn_strength(vals_old, vals_new, supp_fig_dir)

################################################################################
##########supp_figURE 4A: Diff. of New vs. Old Connection Strengths##################
make_md_table_mean <- function(row_levels, col_levels,
                               row_div_map, col_div_map, vals,
                               na.rm = TRUE, fill = NA_real_) {
  vals_mat <- as.matrix(vals)
  R <- nrow(vals_mat)
  C <- ncol(vals_mat)
  rf <- factor(as.character(row_div_map), levels = row_levels)   # length R
  cf <- factor(as.character(col_div_map), levels = col_levels)   # length C
  
  # Expand to element-wise factors that align with as.vector() (column-major)
  rf_vec <- rep(rf, times = C)   # repeats each row label for each column
  cf_vec <- rep(cf, each = R)    # repeats each column label down all rows
  vals_vec <- as.numeric(as.vector(vals_mat))
  
  # use interaction of the expanded factors as a single INDEX (length matches vals_vec)
  idx <- interaction(rf_vec, cf_vec, drop = TRUE)
  mat_raw <- tapply(vals_vec, INDEX = idx, FUN = function(x) mean(x, na.rm = na.rm))
  
  # convert back to matrix of size length(row_levels) x length(col_levels)
  mat <- matrix(as.numeric(mat_raw),
                nrow = length(row_levels),
                ncol = length(col_levels),
                dimnames = list(row_levels, col_levels))
  
  if (!is.na(fill)) mat[is.na(mat)] <- fill
  storage.mode(mat) <- "double"
  mat
}


####main function
make_diff_major_div_heatmap <- function(diff_no_outliers, md, md_col) {
  
  ## Preserve existing major_division order
  fmd <- if (is.factor(md)) factor(md, levels = levels(md)) else factor(md, levels = unique(md))
  fmd_col <- if (is.factor(md_col)) factor(md_col, levels = levels(md_col)) else factor(md_col, levels = unique(md_col))
  
  # Get region -> major-division mappings (character vectors)
  row_div_map <- fmd
  col_div_map <- fmd_col
  
  # Preserve order of major divisions as in your factors (fmd, fmd_col)
  md_row_levels <- as.character(levels(fmd))
  md_col_levels <- as.character(levels(fmd_col))
  
  
  # Build matrices (major-division × major-division)
  zdiff_mean_major_div <- make_md_table_mean(md_row_levels, md_col_levels, 
                                             row_div_map, col_div_map, as.matrix(diff_no_outliers))
  
  
  df_heat <- as.data.frame(as.table(zdiff_mean_major_div))
  colnames(df_heat) <- c("RowDiv", "ColDiv", "Value")
  df_heat$Value <- 100*df_heat$Value ##convert proportion to percent
  # Plot heatmap
  p <- ggplot(df_heat, aes(x = ColDiv, y = RowDiv, fill = Value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", Value)), size = 8) +
    scale_fill_gradient2(
      low = "cyan", mid = "white", high = "red", midpoint = 0,
      name = "Mean Difference in Percentile"
    ) +
    scale_y_discrete(
      limits = rev(levels(df_heat$RowDiv)),
      labels = function(x)
        paste0("<span style='color:", div_cols[x], "'>", x, "</span>")
    ) +
    scale_x_discrete(
      labels = function(x)
        paste0("<span style='color:", div_cols[x], "'>", x, "</span>")
    ) +
    coord_fixed() +
    labs(
      title = "Mean Percentile Difference by Major Division",
      x = "Target Major Division",
      y = "Source Major Division"
    ) +
    theme_minimal(base_size = 30) +
    theme(
      axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_markdown(hjust = 1),
      axis.title.x = element_text(margin = margin(t = 15)),
      axis.title.y = element_text(margin = margin(r = 15)),
      plot.title = element_text(hjust = 0.5)
    )

  
  return(p)
}

# --- compute percentiles on combined matrices (no outlier removal) ---

# combined matrices for each condition
combined_new <- as.matrix(cbind(conn_ipsi_new, conn_contra_new))
combined_old <- as.matrix(cbind(conn_ipsi_old, conn_contra_old))

# function to compute percentile (ecdf) values in matrix form
matrix_percentile <- function(mat) {
  # ecdf returns a function; apply to each entry of the matrix
  f <- ecdf(as.vector(mat))
  pct_vec <- f(as.vector(mat))   # values in [0,1]
  matrix(pct_vec, nrow = nrow(mat), ncol = ncol(mat), byrow = FALSE,
         dimnames = dimnames(mat))
}

pct_new_all   <- matrix_percentile(combined_new)
pct_old_all   <- matrix_percentile(combined_old)

# split back to ipsi / contra parts
n_ipsi <- ncol(conn_ipsi_new)   # assume same # cols in new/old counterparts
# indices
ipsi_idx  <- seq_len(n_ipsi)
contra_idx <- seq(n_ipsi + 1, n_ipsi + ncol(conn_contra_new))

pct_ipsi_new   <- pct_new_all[, ipsi_idx, drop = FALSE]
pct_ipsi_old   <- pct_old_all[, ipsi_idx, drop = FALSE]

pct_contra_new <- pct_new_all[, contra_idx, drop = FALSE]
pct_contra_old <- pct_old_all[, contra_idx, drop = FALSE]

# compute percentile differences (new - old)
diff_ipsi_pct   <- pct_ipsi_new - pct_ipsi_old
diff_contra_pct <- pct_contra_new - pct_contra_old

# force difference to 0 where both original values are exactly zero (optional safeguard)
both_zero_ipsi <- (as.matrix(conn_ipsi_new) == 0) & (as.matrix(conn_ipsi_old) == 0)
both_zero_contra <- (as.matrix(conn_contra_new) == 0) & (as.matrix(conn_contra_old) == 0)

diff_ipsi_pct[both_zero_ipsi] <- 0
diff_contra_pct[both_zero_contra] <- 0

# --- visualize: use your existing heatmap function(s) ---
md  <- row_major_divisions
md_col <- ipsi_col_major_divisions
md_col_contra <- contra_col_major_divisions
p_ipsi <- make_diff_major_div_heatmap(diff_ipsi_pct, md, md_col)
ggsave(paste0(supp_fig_dir, "supp_figure_4a_diff_pct_ipsi.png"), p_ipsi, width = 20, height = 16, dpi = 300)

p_contra <- make_diff_major_div_heatmap(diff_contra_pct, md, md_col_contra)
ggsave(paste0(supp_fig_dir, "supp_figure_4a_diff_pct_contra.png"), p_contra, width = 20, height = 16, dpi = 300)


###########repeat for contra##################

################################################################################
###########supp_figURE 4B: CONDENSED HEATMAPS WITH CHANGES IN BINARIZED CONNECTOMES#########
##################HELPER FUNCTIONS#######################################
# Helper to build major div. x major div. table from index pairs
make_md_table <- function(i_idx, j_idx, row_levels, col_levels,
                          row_div_map, col_div_map) {
  if (length(i_idx) == 0) {
    # return zero matrix with proper dimnames
    m <- matrix(0L, nrow = length(row_levels), ncol = length(col_levels),
                dimnames = list(row_levels, col_levels))
    return(m)
  }
  df <- data.frame(row_div = factor(row_div_map[i_idx], levels = row_levels),
                   col_div = factor(col_div_map[j_idx], levels = col_levels))
  tbl <- as.matrix(with(df, table(row_div, col_div)))
  # ensure integer
  storage.mode(tbl) <- "integer"
  return(tbl)
}

# ---- Prepare a long-format dataframe for plotting ----
make_long_df <- function(mat, metric_name) {
  as.data.frame(as.table(mat), stringsAsFactors = FALSE) %>%
    setNames(c("major_row", "major_col", "value")) %>%
    mutate(metric = metric_name)
}

# helper: map percent -> hex color from white -> color
pct_to_hex <- function(pct, color = c("indianred1","deepskyblue2")) {
  color <- match.arg(color)
  # clip pct to [0,100], replace NA with 0 for mapping (we'll handle NA display later)
  pct_clipped <- ifelse(is.na(pct), 0, pmin(pmax(pct, 0), 100))/100
  # get RGB values from white->color ramp
  ramp_fun <- grDevices::colorRamp(c("white", color))
  cols_rgb <- ramp_fun(matrix(pct_clipped, ncol=1))
  hex <- grDevices::rgb(cols_rgb[,1], cols_rgb[,2], cols_rgb[,3], maxColorValue = 255)
  # for NA percentages return NA so geom will be transparent
  hex[is.na(pct)] <- NA
  return(hex)
}

####main function
make_condensed_major_div_heatmap <- function(mat_old, mat_new, mat_diff, md, md_col) {
  
  ## Preserve existing major_division order
  fmd <- if (is.factor(md)) factor(md, levels = levels(md)) else factor(md, levels = unique(md))
  fmd_col <- if (is.factor(md_col)) factor(md_col, levels = levels(md_col)) else factor(md_col, levels = unique(md_col))
  
  # Get region -> major-division mappings (character vectors)
  row_div_map <- fmd
  col_div_map <- fmd_col
  
  # Preserve order of major divisions as in your factors (fmd, fmd_col)
  md_row_levels <- as.character(levels(fmd))
  md_col_levels <- as.character(levels(fmd_col))
  
  # Indices for each condition
  gain_idx  <- which(mat_diff ==  1, arr.ind = TRUE)
  lost_idx  <- which(mat_diff == -1, arr.ind = TRUE)
  old_idx   <- which(mat_old != 0, arr.ind = TRUE)
  
  # Build matrices (major-division × major-division)
  mat_gained <- make_md_table(gain_idx[,1], gain_idx[,2], md_row_levels, md_col_levels,
                              row_div_map, col_div_map)
  mat_lost   <- make_md_table(lost_idx[,1], lost_idx[,2], md_row_levels, md_col_levels,
                              row_div_map, col_div_map)
  mat_total  <- make_md_table(old_idx[,1], old_idx[,2], md_row_levels, md_col_levels,
                              row_div_map, col_div_map)
  
  # Add informative dimnames (if not already)
  rownames(mat_gained) <- rownames(mat_lost) <- rownames(mat_total) <- md_row_levels
  colnames(mat_gained) <- colnames(mat_lost) <- colnames(mat_total) <- md_col_levels
  
  df_long <- bind_rows(
    make_long_df(mat_gained, "gained"),
    make_long_df(mat_lost,   "lost"),
    make_long_df(mat_total,  "total_old")
  )
  
  # mark metrics as ordered factor for nice plotting order
  df_long$metric <- factor(df_long$metric, levels = c("gained", "lost", "total_old"))
  
  # --- Build a data.frame with percent values per cell ---
  row_names <- rownames(mat_total)
  col_names <- colnames(mat_total)
  
  df <- expand.grid(major_col = col_names, major_row = row_names, stringsAsFactors = FALSE) %>%
    arrange(factor(major_row, levels = row_names), factor(major_col, levels = col_names))
  
  # fill numeric values (use t() so that expand.grid ordering matches columns-major)
  df$gain_count  <- as.integer(as.vector(t(mat_gained)))
  df$loss_count  <- as.integer(as.vector(t(mat_lost)))
  df$total_count <- as.integer(as.vector(t(mat_total)))
  
  # percent (0-100). If total==0 we set NA
  df <- df %>%
    mutate(gain_pct = ifelse(total_count > 0, 100 * gain_count / total_count, NA_real_),
           loss_pct = ifelse(total_count > 0, 100 * loss_count / total_count, NA_real_),
           x = as.numeric(factor(major_col, levels = col_names)),
           y = as.numeric(factor(major_row, levels = row_names)))
  
  # compute fill hex colors per triangle
  df$fill_gain_hex <- pct_to_hex(df$gain_pct, color = "indianred1")
  df$fill_loss_hex <- pct_to_hex(df$loss_pct, color = "deepskyblue2")
  
  # --- Build polygon rows for two triangles per cell ---
  # square half-size
  h <- 0.45
  
  # For each cell we will create two triangles split along top-right -> bottom-left diagonal:
  # vertices of the square:
  # tl = (x-h, y+h)
  # tr = (x+h, y+h)
  # br = (x+h, y-h)
  # bl = (x-h, y-h)
  # Triangles:
  #  - upper-left (blue/loss): tl, tr, bl
  #  - lower-right (red/gain): tr, br, bl
  
  tri_list <- lapply(seq_len(nrow(df)), function(i) {
    r <- df[i, ]
    tl <- c(r$x - h, r$y + h)
    tr <- c(r$x + h, r$y + h)
    br <- c(r$x + h, r$y - h)
    bl <- c(r$x - h, r$y - h)
    # upper-left triangle (loss / blue)
    loss_poly <- data.frame(
      major_row = r$major_row, major_col = r$major_col,
      x = c(tl[1], tr[1], bl[1]),
      y = c(tl[2], tr[2], bl[2]),
      type = "loss",
      pct = r$loss_pct,
      fill = r$fill_loss_hex,
      stringsAsFactors = FALSE
    )
    # lower-right triangle (gain / red)
    gain_poly <- data.frame(
      major_row = r$major_row, major_col = r$major_col,
      x = c(tr[1], br[1], bl[1]),
      y = c(tr[2], br[2], bl[2]),
      type = "gain",
      pct = r$gain_pct,
      fill = r$fill_gain_hex,
      stringsAsFactors = FALSE
    )
    rbind(loss_poly, gain_poly)
  })
  
  poly_df <- do.call(rbind, tri_list)
  
  # Mark rows with total_count == 0 -> give them a light grey background tile and no colored triangles
  empty_cells <- df %>% filter(total_count == 0) %>% mutate(x = x, y = y)
  
    # --- Plot ---
  p <- ggplot() +
    geom_tile(data = df, aes(x = x, y = y), width = 0.98, height = 0.98,
              fill = "white", color = "grey90") +
    geom_tile(data = empty_cells, aes(x = x, y = y), width = 0.98, height = 0.98,
              fill = "grey97", color = "grey85") +
    geom_polygon(data = poly_df %>% filter(type == "loss"),
                aes(x = x, y = y, group = interaction(major_row, major_col, type)),
                fill = poly_df$fill[poly_df$type == "loss"], color = NA) +
    geom_polygon(data = poly_df %>% filter(type == "gain"),
                aes(x = x, y = y, group = interaction(major_row, major_col, type)),
                fill = poly_df$fill[poly_df$type == "gain"], color = NA) +
    geom_text(data = df,
              aes(x = x - 0.18, y = y + 0.12,
                  label = ifelse(is.na(loss_pct), "", sprintf("%.0f%%", loss_pct))),
              size = 5) +
    geom_text(data = df,
              aes(x = x + 0.18, y = y - 0.12,
                  label = ifelse(is.na(gain_pct), "", sprintf("%.0f%%", gain_pct))),
              size = 5) +
    scale_x_continuous(
      breaks = seq_along(col_names),
      labels = function(x) {
        labs <- col_names[x]
        paste0("<span style='color:", div_cols[labs], "'>", labs, "</span>")
      },
      expand = c(0,0)
    ) +
    scale_y_reverse(
      breaks = seq_along(row_names),
      labels = function(y) {
        labs <- row_names[y]
        paste0("<span style='color:", div_cols[labs], "'>", labs, "</span>")
      },
      expand = c(0,0)
    ) +
    coord_fixed() +
    labs(
      x = "Target division",
      y = "Source division (right hemisphere)",
      title = "Percent of Connections Gained and Lost after Rebuilding"
    ) +
    theme_minimal(base_size = 24) +
    theme(
      axis.text.x = element_markdown(angle = 45, hjust = 1),
      axis.text.y = element_markdown(),
      panel.grid = element_blank()
    )
  return(p)
}

####load inputs for visualization
mat_old  <- as.matrix(conn_ipsi_old_bin)
mat_new  <- as.matrix(conn_ipsi_new_bin)
mat_diff <- mat_new - mat_old  # +1 gained, -1 lost
md  <- row_major_divisions
md_col <- ipsi_col_major_divisions

p <- make_condensed_major_div_heatmap(mat_old, mat_new, mat_diff, md, md_col)
ggsave(paste0(supp_fig_dir,"supp_figure_4a_ipsi_bin",binary_pct_threshold,".png"), p,
       width = 16, height = 16, dpi = 300)

##repeat for contra
mat_old  <- as.matrix(conn_contra_old_bin)
mat_new  <- as.matrix(conn_contra_new_bin)
mat_diff <- mat_new - mat_old  # +1 gained, -1 lost
md  <- row_major_divisions
md_col <- contra_col_major_divisions
p <- make_condensed_major_div_heatmap(mat_old, mat_new, mat_diff, md, md_col)
ggsave(paste0(supp_fig_dir,"supp_figure_4a_contra_bin",binary_pct_threshold,".png"), p,
       width = 16, height = 16, dpi = 300)
