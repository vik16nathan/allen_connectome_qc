library("pacman")
pacman::p_load(tidyverse, pheatmap, dplyr, readr, readxl, patchwork, ggrepel, RColorBrewer)
setwd(".")

###LOAD VISUALIZATION SCRIPTS, ARGUMENTS; SET PATHS AND INPUTS#############
source("./ggslicer/plotting_functions/plotting_functions.R")
source("./ggslicer/plotting_functions/plotting_functions_labels.R")
allen_input_dir <- "../preprocessed/allen_template_inputs/"

output_fig_dir_knox <- "./figures/"
output_fig_dir_oh <- "./figures/oh/"
args <- commandArgs(trailingOnly=TRUE)
knox_or_oh <- args[1]
binary_pct_threshold <- as.numeric(args[2])

##########################################################################################
##########################HELPER FUNCTION################################################
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
################################################################################################

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

#################TOGGLE BETWEEN KNOX/OH CONNECTOMES HERE#######################
##original: don't exclude any of the regions
knox_conn_ipsi_regions_in_oh <- rownames(knox_conn_ipsi_old)[which(rownames(knox_conn_ipsi_old) %in% rownames(oh_conn_ipsi_old))]
knox_conn_contra_regions_in_oh <- colnames(knox_conn_contra_old)[which(colnames(knox_conn_contra_old) %in% colnames(oh_conn_contra_old))]

connectome_dir <- "mouse_connectivity_models/paper/figures/model_comparison/output/"
oh_conn_old <- as.data.frame(read_csv(paste0(connectome_dir,"homogeneous-standard-model_original.csv")))
oh_conn_region_numbers <- oh_conn_old[2:nrow(oh_conn_old),1]

if(knox_or_oh == "oh") {
  fig_dir <- output_fig_dir_oh
  conn_region_numbers <- oh_conn_region_numbers
  conn_ipsi_old <- oh_conn_ipsi_old
  conn_contra_old <- oh_conn_contra_old
  conn_ipsi_new <- oh_conn_ipsi_new
  conn_contra_new <- oh_conn_contra_new
}

if(knox_or_oh == "knox") {
  fig_dir <- output_fig_dir_knox
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
###########FIGURE 4A: HEATMAPS WITH CHANGES IN BINARIZED CONNECTOMES#########
mat_old  <- as.matrix(conn_ipsi_old_bin)
mat_new  <- as.matrix(conn_ipsi_new_bin)
mat_diff <- mat_new - mat_old  # +1 gained, -1 lost

## Build composite code matrix: 0,1,2,3 (see legend below)
code <- matrix(NA_integer_, nrow(mat_old), ncol(mat_old))
code[mat_diff == 0 & mat_old == 0] <- 0   # unchanged 0 -> white
code[mat_diff == 0 & mat_old == 1] <- 1   # unchanged 1 -> black
code[mat_diff ==  1]               <- 2   # gained     -> red
code[mat_diff == -1]               <- 3   # lost       -> blue
rownames(code) <- rownames(mat_old)
colnames(code) <- colnames(mat_old)

## Preserve existing major_division order
md  <- row_major_divisions
fmd <- if (is.factor(md)) factor(md, levels = levels(md)) else factor(md, levels = unique(md))

md_col_ipsi <- ipsi_col_major_divisions
md <- md_col_ipsi
fmd_col <- if (is.factor(md)) factor(md, levels = levels(md)) else factor(md, levels = unique(md))

###Row/column name annotations
annotation_row <- data.frame(major_division = fmd)
annotation_col <- data.frame(major_division = fmd_col)
rownames(annotation_row) <- rownames(code)
rownames(annotation_col) <- colnames(code)

##Diverging categorical colors from iWantHue (replace with your own hexes if you like)
# Example: 12 distinct hues from iWantHue-style palettes
iwant_hex <- c("#a83537",
               "#4fc79c",
               "#63348a",
               "#73c161",
               "#6280d6",
               "#d3a046",
               "#ca78cd",
               "#4c792a",
               "#bb467a",
               "#aab248",
               "#a96126",
               "#ff846b")

## Map colors to the UNION of row/col division levels (keeps one legend)
div_lvls <- unique(c(as.character(fmd), as.character(fmd_col)))
div_cols <- setNames(iwant_hex[seq_along(div_lvls)], div_lvls)

annotation_colors <- list(
  major_division = div_cols
)

## 4) Discrete cell palette for your 4 states (unchanged 0/1, gained, lost)
cols   <- c("#F2F2F2", "black", "red", "blue")
breaks <- c(-0.5, 0.5, 1.5, 2.5, 3.5)

## 5) Optional gaps between groups (rows & cols)
gaps     <- cumsum(table(fmd))[ -length(levels(fmd)) ]
gaps_col <- cumsum(table(fmd_col))[ -length(levels(fmd_col)) ]

## 6) Plot
p <- pheatmap(code,
              cluster_rows = FALSE, cluster_cols = FALSE,
              color = cols, breaks = breaks,
              legend_breaks = 0:3,
              legend_labels = c("old=0 (unchanged)",
                                "old=1 (unchanged)",
                                "gained (+1)",
                                "lost (-1)"),
              title="Top 20% of Regionalized Knox Connections, Before and After QC",
              annotation_row = annotation_row,
              annotation_col = annotation_col,
              annotation_colors = annotation_colors,
              annotation_names_row = FALSE,   # hides the label "major_division" beside the row bar
              annotation_names_col = FALSE,   # hides the label "major_division" above the col bar
              show_rownames = FALSE, show_colnames = FALSE,
              gaps_row = gaps, gaps_col = gaps_col,
              border_color = NA,
              fontsize=30)

ggsave(paste0(fig_dir,"figure_4b_ipsi_heatmap_bin", binary_pct_threshold,".png"), p, width = 24, height = 16, dpi = 300)  # adjust width/height as needed


#############REPEAT FOR CONTRA############################
mat_old  <- as.matrix(conn_contra_old_bin)
mat_new  <- as.matrix(conn_contra_new_bin)
mat_diff <- mat_new - mat_old  # +1 gained, -1 lost

## Build composite code matrix: 0,1,2,3 (see legend below)
code <- matrix(NA_integer_, nrow(mat_old), ncol(mat_old))
code[mat_diff == 0 & mat_old == 0] <- 0   # unchanged 0 -> white
code[mat_diff == 0 & mat_old == 1] <- 1   # unchanged 1 -> black
code[mat_diff ==  1]               <- 2   # gained     -> red
code[mat_diff == -1]               <- 3   # lost       -> blue
rownames(code) <- rownames(mat_old)
colnames(code) <- colnames(mat_old)

## Preserve existing major_division order
md  <- row_major_divisions
fmd <- if (is.factor(md)) factor(md, levels = levels(md)) else factor(md, levels = unique(md))

md_col_contra <- contra_col_major_divisions
md <- md_col_contra
fmd_col <- if (is.factor(md)) factor(md, levels = levels(md)) else factor(md, levels = unique(md))

###Row/column name annotations
annotation_row <- data.frame(major_division = fmd)
annotation_col <- data.frame(major_division = fmd_col)
rownames(annotation_row) <- rownames(code)
rownames(annotation_col) <- colnames(code)

## Map colors to the UNION of row/col division levels (keeps one legend)
div_lvls <- unique(c(as.character(fmd), as.character(fmd_col)))
div_cols <- setNames(iwant_hex[seq_along(div_lvls)], div_lvls)

annotation_colors <- list(
  major_division = div_cols
)

## 4) Discrete cell palette for your 4 states (unchanged 0/1, gained, lost)
cols   <- c("#F2F2F2", "black", "red", "blue")
breaks <- c(-0.5, 0.5, 1.5, 2.5, 3.5)

## 5) Optional gaps between groups (rows & cols)
gaps     <- cumsum(table(fmd))[ -length(levels(fmd)) ]
gaps_col <- cumsum(table(fmd_col))[ -length(levels(fmd_col)) ]

## 6) Plot
p <- pheatmap(code,
              cluster_rows = FALSE, cluster_cols = FALSE,
              color = cols, breaks = breaks,
              legend_breaks = 0:3,
              legend_labels = c("old=0 (unchanged)",
                                "old=1 (unchanged)",
                                "gained (+1)",
                                "lost (-1)"),
              title="Top 20% of Regionalized Knox Connections, Before and After QC",
              annotation_row = annotation_row,
              annotation_col = annotation_col,
              annotation_colors = annotation_colors,
              annotation_names_row = FALSE,   # hides the label "major_division" beside the row bar
              annotation_names_col = FALSE,   # hides the label "major_division" above the col bar
              gaps_row = gaps, gaps_col = gaps_col,
              show_rownames = FALSE, show_colnames = FALSE,
              border_color = NA,
              fontsize=30)

ggsave(paste0(fig_dir,"figure_4b_contra_heatmap_bin",binary_pct_threshold, ".png"), p, width = 24, height = 16, dpi = 300)  # adjust width/height as needed


#############REPEAT FOR CONNECTION STRENGTHS/RATIOS###########################
## Preserve existing major_division order
md  <- row_major_divisions
fmd <- if (is.factor(md)) factor(md, levels = levels(md)) else factor(md, levels = unique(md))

md_col_ipsi <- ipsi_col_major_divisions
md <- md_col_ipsi
fmd_col <- if (is.factor(md)) factor(md, levels = levels(md)) else factor(md, levels = unique(md))

major_order <- names(major_division_dict)

fmd     <- factor(row_major_divisions, levels = major_order)
fmd_col <- factor(ipsi_col_major_divisions, levels = major_order)

###Row/column name annotations
annotation_row <- data.frame(major_division = fmd)
annotation_col <- data.frame(major_division = fmd_col)
rownames(annotation_row) <- rownames(conn_ipsi_old)
rownames(annotation_col) <- colnames(conn_ipsi_old)

## Map colors to the UNION of row/col division levels (keeps one legend)
div_lvls <- unique(c(as.character(fmd), as.character(fmd_col)))
div_cols <- setNames(iwant_hex[seq_along(div_lvls)], div_lvls)

annotation_colors <- list(
  major_division = div_cols
)

##zero out diagonals (connection strength not meaningful)
for(j in colnames(conn_ipsi_old)){
  conn_ipsi_old[j,j] <- NA
  conn_ipsi_new[j,j] <- NA
}

for(j in colnames(conn_contra_old)){
  conn_contra_old[j,j] <- NA
  conn_contra_new[j,j] <- NA
}


##replace outliers with max non-outlier value
# Compute the ratio data frame
z_diff <- scale(conn_ipsi_new) - scale(conn_ipsi_old)

# Function to replace outliers beyond 3*SD with max non-outlier
replace_outliers <- function(x) {
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

# Apply column-wise
z_diff_no_outliers <- as.data.frame(replace_outliers(as.matrix(z_diff)))

# Define a diverging color palette
n_colors <- 100
my_colors <- colorRampPalette(c("blue", "white", "red"))(n_colors)

# Define breaks from -1 to 1
#my_breaks <- seq(-1, 1, length.out = n_colors + 1)
my_breaks <- seq(-1, 1, length.out = n_colors + 1)

## 5) Optional gaps between groups (rows & cols)
gaps     <- cumsum(table(fmd))[ -length(levels(fmd)) ]
gaps_col <- cumsum(table(fmd_col))[ -length(levels(fmd_col)) ]

# Plot heatmap with forced color scale
p <- pheatmap(
  z_diff_no_outliers,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  title = "Connection Strength Z Post - Pre-QC",
  annotation_row = annotation_row,
  annotation_col = annotation_col,
  gaps_row=gaps, gaps_col = gaps_col,
  annotation_colors = annotation_colors,
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  fontsize = 30,
  color = my_colors,
  breaks = my_breaks
)

p
ggsave(paste0(fig_dir,"fig_4a_ipsi_conn_strength_diff_z_heatmap.png"), p, width = 20, height = 16, dpi = 300)  # adjust width/height as needed

##########################repeat for contra###########################
md  <- row_major_divisions
fmd <- if (is.factor(md)) factor(md, levels = levels(md)) else factor(md, levels = unique(md))

md_col_contra <- contra_col_major_divisions
md <- md_col_contra
fmd_col <- if (is.factor(md)) factor(md, levels = levels(md)) else factor(md, levels = unique(md))

###Row/column name annotations
annotation_row <- data.frame(major_division = fmd)
annotation_col <- data.frame(major_division = fmd_col)
rownames(annotation_row) <- rownames(conn_contra_old)
rownames(annotation_col) <- colnames(conn_contra_old)

## 5) Optional gaps between groups (rows & cols)
gaps     <- cumsum(table(fmd))[ -length(levels(fmd)) ]
gaps_col <- cumsum(table(fmd_col))[ -length(levels(fmd_col)) ]


z_diff <- scale(conn_contra_new) - scale(conn_contra_old)

# Apply column-wise
z_diff_no_outliers <- as.data.frame(replace_outliers(as.matrix(z_diff)))

# Plot heatmap with forced color scale
p <- pheatmap(
  z_diff_no_outliers,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  title = "Connection Strength Z-Score Post - Pre-QC",
  annotation_row = annotation_row,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  gaps_row = gaps, gaps_col = gaps_col,
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  fontsize = 30,
  color = my_colors,
  breaks = my_breaks
)
p
ggsave(paste0(fig_dir,"fig_4a_contra_conn_strength_diff_z_heatmap.png"), p, width = 20, height = 16, dpi = 300)  # adjust width/height as needed

