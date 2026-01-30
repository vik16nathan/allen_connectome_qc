########################Vikram Nathan, 08/13/2025##########################
###prerequisites: run create_automated_qc_csv.R with threshold of interest##
library("pacman")
library("RMINC")
pacman::p_load(dplyr, readr, readxl, stringr, foreach, doParallel, ggplot2, ggrepel, tidyr)
setwd(".")
source("./ggslicer/plotting_functions/plotting_functions.R")
source("./ggslicer/plotting_functions/plotting_functions_labels.R")

inj_thresh <- 0.5
proj_thresh <- 0.1
allen_50um_template_path <- "../preprocessed/allen_template_inputs/average_template_50.mnc"
allen_mask_path_50um <- "../preprocessed/allen_template_inputs/mask_50um.mnc"
vent_path <- "../preprocessed/allen_template_inputs/vent_voxels_label_file_50um_filt.mnc"
##directory containing .mnc files
tracer_dir <- "../preprocessed/knox_connectome_tracers/"

##################HELPER FUNCTIONS#############################################
#####function to plot tracer + template at slice coordinate using Yohan's vis#####
plot_tracer_proj_slice <- function(tracer, tracer_dir, slice_coords, proj_thresh, template_path, mask_path) {
  allen_template_df <- prepare_masked_anatomy(template_path, mask_path, "y", seq(slice_coords, slice_coords, 0))[[2]]
  allen_template_df <- allen_template_df %>% filter(mask_value == 1)
  tracer_proj_df <- prepare_anatomy(paste0(tracer_dir, tracer,"/projection_density_bin",proj_thresh,".mnc"), "y", seq(slice_coords, slice_coords, 0))[[2]]
  tracer_proj_df <- tracer_proj_df %>% filter(intensity > 0)
  
  # Combine "slice_world" with "y" for unique facet labels
  tracer_proj_df <- tracer_proj_df %>% 
    mutate(facet_label = paste0("y = ", y))
  
  ggplot(data = allen_template_df, mapping = aes(x = x, y = z)) +
    geom_raster(mapping = aes(fill = intensity), interpolate = TRUE) +
    scale_fill_gradient(
      low = 'black',
      high = 'white',
      oob = scales::squish,
      guide = 'none'
    ) +
    ggnewscale::new_scale_fill() + 
    geom_tile(data = tracer_proj_df, mapping = aes(fill = intensity), alpha=0.5) +
    scale_fill_gradient2(
      name = paste0("Projection > ",proj_thresh),
      low = "red",
      mid = "red",
      high = "red",
      midpoint = 1,
      space = "Lab",
      na.value = "grey50",
      transform = "identity",
      guide = "colourbar",
      aesthetics = "fill"
    ) + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_wrap(~ slice_world, ncol = 1, labeller = labeller(slice_world = function(x) paste0("Slice: ", x, " mm"))) +  
    coord_fixed(ratio = 1) +
    theme_void() +
    theme(
      panel.spacing = unit(0, "npc"),                     # Keep spacing minimal
      strip.text = element_text(size = 60, face = "bold"), # Increase facet label text size
      plot.title = element_text(size = 70, face = "bold", hjust = 0.5), # Make title bigger and center it
    ) +
    labs(
      title = tracer # Add title
    )
  
}

plot_tracer_proj_slice_vent <- function(
    tracer, tracer_dir, slice_coords, proj_thresh,
    template_path, mask_path, vent_path
) {
  allen_template_df <- prepare_masked_anatomy(template_path, mask_path, "y",
                                              seq(slice_coords, slice_coords, 0))[[2]] %>%
    filter(mask_value == 1)
  
  tracer_proj_df <- prepare_anatomy(
    paste0(tracer_dir, tracer, "/projection_density_bin", proj_thresh, ".mnc"),
    "y", seq(slice_coords, slice_coords, 0)
  )[[2]] %>% filter(intensity > 0)
  
  vent_vol_df <- prepare_masked_anatomy(
    vent_path, mask_path, "y", seq(slice_coords, slice_coords, 0)
  )[[2]] %>% filter(mask_value == 1, intensity > 0)
  
  # ---- NEW: Identify overlap (voxels with proj > 0 AND vent > 0) ----
  overlap_df <- inner_join(
    tracer_proj_df %>% select(x, y, z),
    vent_vol_df %>% select(x, y, z),
    by = c("x", "y", "z")
  )
  
  # Add facet labels
  tracer_proj_df <- tracer_proj_df %>%
    mutate(facet_label = paste0("y = ", y))
  
  ggplot(data = allen_template_df, aes(x = x, y = z)) +
    geom_raster(aes(fill = intensity), interpolate = TRUE) +
    scale_fill_gradient(low = "black", high = "white", guide = "none") +
    
    ggnewscale::new_scale_fill() +
    geom_tile(data = tracer_proj_df, aes(fill = intensity), alpha = 0.5) +
    scale_fill_gradient2(
      name = paste0("Projection > ", proj_thresh),
      low = "red", mid = "red", high = "red",
      midpoint = 1, guide = "colourbar"
    ) +
    
    ggnewscale::new_scale_fill() +
    geom_tile(data = vent_vol_df, aes(fill = intensity), alpha = 0.5) +
    scale_fill_gradient2(
      name = "Ventricles",
      low = "cyan", mid = "cyan", high = "cyan",
      midpoint = 1, guide = "colourbar"
    ) +
    
    # ---- NEW OVERLAP LAYER ----
  # bright red highlight: full opacity
  geom_tile(data = overlap_df, aes(x = x, y = z), fill = "orange", alpha = 1) +
    
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_wrap(~ slice_world, ncol = 1,
               labeller = labeller(slice_world = function(x) paste0("Slice: ", x, " mm"))) +
    coord_fixed() +
    theme_void() +
    theme(
      panel.spacing = unit(0, "npc"),
      strip.text = element_text(size = 60, face = "bold"),
      plot.title = element_text(size = 70, face = "bold", hjust = 0.5)
    ) +
    labs(title = tracer)
}


plot_tracer_inj_slice <- function(tracer, tracer_dir, slice_coords, inj_thresh, template_path, mask_path) {
  allen_template_df <- prepare_masked_anatomy(template_path, mask_path, "y", seq(slice_coords, slice_coords, 0))[[2]]
  allen_template_df <- allen_template_df %>% filter(mask_value == 1)
  tracer_inj_df <- prepare_anatomy(paste0(tracer_dir, tracer,"/injection_dens_times_frac_bin",inj_thresh,".mnc"), "y", seq(slice_coords, slice_coords, 0))[[2]]
  tracer_inj_df <- tracer_inj_df %>% filter(intensity > 0)
  
  # Combine "slice_world" with "y" for unique facet labels
  tracer_inj_df <- tracer_inj_df %>% 
    mutate(facet_label = paste0("y = ", y))
  
  ggplot(data = allen_template_df, mapping = aes(x = x, y = z)) +
    geom_raster(mapping = aes(fill = intensity), interpolate = TRUE) +
    scale_fill_gradient(
      low = 'black',
      high = 'white',
      oob = scales::squish,
      guide = 'none'
    ) +
    ggnewscale::new_scale_fill() + 
    geom_tile(data = tracer_inj_df, mapping = aes(fill = intensity), alpha=0.5) +
    scale_fill_gradient2(
      name = paste0("Injection > ",inj_thresh),
      low = "purple",
      mid = "purple",
      high = "purple",
      midpoint = 1,
      space = "Lab",
      na.value = "grey50",
      transform = "identity",
      guide = "colourbar",
      aesthetics = "fill"
    ) + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_wrap(~ slice_world, ncol = 1, labeller = labeller(slice_world = function(x) paste0("Slice: ", x, " mm"))) +  
    coord_fixed(ratio = 1) +
    theme_void() +
    theme(
      panel.spacing = unit(0, "npc"),                     # Keep spacing minimal
      strip.text = element_text(size = 60, face = "bold"), # Increase facet label text size
      plot.title = element_text(size = 70, face = "bold", hjust = 0.5), # Make title bigger and center it
    ) +
    labs(
      title = tracer # Add title
    )
  
}


###############HISTOGRAMS OF OOB and VENT VOXELS: FIGURE 2A##############
###read in data
tracer_num_vox_oob_vent_df <- as.data.frame(read.csv(paste0("tables/knox_oob_vent_df_inj",inj_thresh,"_proj",proj_thresh,".csv")))
tracer_num_vox_oob_vent_df <- tracer_num_vox_oob_vent_df[which(tracer_num_vox_oob_vent_df$tracer != 310207648),]
qc_table_removal_full <- as.data.frame(read_csv("tables/tracers_to_remove_adjboxStats_skew_outliers.csv"))
oob_outliers <- qc_table_removal_full[which(qc_table_removal_full[,"Auto OOB Proj."] == 1), "Tracer"]
vent_outliers <- qc_table_removal_full[which(qc_table_removal_full[,"Auto Vent Proj."] == 1), "Tracer"]

# Step 1: Gather all 4 vectors into a long dataframe with labels
long_df <- tracer_num_vox_oob_vent_df %>%
  select(tracer, proj_vox_oob, proj_vox_vent) %>%
  pivot_longer(
    cols = -tracer,
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    facet_label = case_when(
      metric == "proj_vox_oob" ~ "OOB Voxels, Projection > 0.1",
      metric == "proj_vox_vent" ~ "Ventricular Voxels, Projection > 0.1"
    ),
    type = ifelse(grepl("^inj", metric), "Injection", "Projection")
  )

# Step 2: Define 50 bins per facet
bin_breaks <- function(x) {
  seq(min(x), max(x), length.out = 51)
}

# Apply binning and tag upper-tail outliers (identified previously using robustbase)
long_df <- long_df %>%
  group_by(facet_label) %>%
  mutate(
    bin = cut(value, breaks = bin_breaks(value), include.lowest = TRUE),

    upper_tail = case_when(
      facet_label == "OOB Voxels, Projection > 0.1" &
        tracer %in% oob_outliers ~ TRUE,

      facet_label == "Ventricular Voxels, Projection > 0.1" &
        tracer %in% vent_outliers ~ TRUE,

      TRUE ~ FALSE
    ),

    tracer_label = ifelse(upper_tail, tracer, NA)
  )


p <- ggplot(long_df, aes(x = value, fill = type)) +
  geom_histogram(bins = 50, color = "black", alpha = 0.7) +
  geom_text_repel(
    data = filter(long_df, upper_tail),
    aes(label = tracer_label),
    stat = "identity",
    y = 0,
    nudge_y = 100,
    size = 8,
    segment.color = "black"
  ) +
  facet_wrap(~ facet_label, scales = "free", ncol = 1) +  # vertical facets
  scale_fill_manual(values = c("Injection" = "#B0E57C",   # pastel blue
                               "Projection" = "#F4C2C2")) +  # pastel pink
  labs(x = "Voxel Count", y = "Count", fill = "Type") +
  theme_minimal(base_size = 40)

###save plot
ggsave("figures/figure_2a_hist.png", p, width = 32, height = 16, dpi = 300)  # adjust width/height as needed

#####plot representative slices of outlier images
tracer <- 112672974
p <- plot_tracer_proj_slice(tracer, tracer_dir, -4, 0.1, allen_50um_template_path, allen_mask_path_50um)
ggsave("figures/figure_2a_112672974.png", p, width = 16, height = 16, dpi = 300)

tracer <- 113036264
p <- plot_tracer_proj_slice(tracer, tracer_dir, -1, 0.1, allen_50um_template_path, allen_mask_path_50um)
ggsave("figures/figure_2a_113036264.png", p, width = 16, height = 16, dpi = 300)

tracer <- 112670853
p <- plot_tracer_proj_slice(tracer, tracer_dir, 0.9, 0.1, allen_50um_template_path, allen_mask_path_50um)
ggsave("figures/figure_2a_112670853.png", p, width = 16, height = 16, dpi = 300)


#ventricle plots
tracer <- 158257355
p <- plot_tracer_proj_slice_vent(tracer, tracer_dir, -0.122, 0.1, allen_50um_template_path, allen_mask_path_50um, vent_path)
ggsave("figures/figure_2a_158257355.png", p, width = 16, height = 16, dpi = 300)

tracer <- 180435652
p <- plot_tracer_proj_slice_vent(tracer, tracer_dir, -0.1, 0.1, allen_50um_template_path, allen_mask_path_50um, vent_path)
ggsave("figures/figure_2a_180435652.png", p, width = 16, height = 16, dpi = 300)

tracer <- 180296424
p <- plot_tracer_proj_slice_vent(tracer, tracer_dir, -1.2, 0.1, allen_50um_template_path, allen_mask_path_50um, vent_path)
ggsave("figures/figure_2a_180296424.png", p, width = 16, height = 16, dpi = 300)


##################FIGURE 2B###############################################
##load automated QC for previously removed experiments
tracer_num_vox_oob_vent_df_prev <- as.data.frame(read.csv(paste0("tables/knox_oob_vent_df_inj",inj_thresh,"_proj",proj_thresh,"_orig_removed.csv")))
# Step 1: Identify outliers (robust) after automated QC: overall_qc_exclusion.R
qc_table_removal_full <- as.data.frame(read_csv("tables/tracers_to_remove_adjboxStats_skew_outliers.csv"))
inj_size_outliers <- qc_table_removal_full[which(qc_table_removal_full[,"Auto Large Inj."] == 1), "Tracer"]
proj_size_outliers <- qc_table_removal_full[which(qc_table_removal_full[,"Auto Large Proj."] == 1), "Tracer"]

zscore_thresh <- 4
tracer_num_vox_oob_vent_df <- tracer_num_vox_oob_vent_df %>%
  mutate(
    sd_outlier = tracer %in% inj_size_outliers | tracer %in% proj_size_outliers,
    # Categorize outlier type (including overlapping cases)
    outlier_type = case_when(
      sd_outlier ~ "Large Inj/Proj.",
      TRUE ~ "None"
    ),
    
    # Label only the outliers
    label = ifelse(outlier_type != "None", tracer, NA)
  )

# Step 4: Plot with different colors
p <- ggplot(tracer_num_vox_oob_vent_df, aes(x = inj_num_vox, y = proj_num_vox)) +
  geom_point(aes(color = outlier_type),size=6) +
  scale_color_manual(values = c(
    "None" = "black",
    "Large Inj/Proj." = "red"
  )) +
  geom_text_repel(aes(label = label), size = 15, na.rm = TRUE) +  # increase label size
  labs(
    x = "Injection Voxel Count",
    y = "Projection Voxel Count",
    color = "Outlier Type"
  ) +
  theme_minimal(base_size = 45)  # increase base text size for the whole plot

ggsave("figures/figure_2b_scatterplot.png", p, width=20, height=16, dpi=300)

# Step 4: Plot with updated legend and label
p <- ggplot(tracer_num_vox_oob_vent_df, aes(x = inj_num_vox, y = proj_num_vox)) +
  # Main dataset (current)
  geom_point(aes(color = outlier_type), size = 6) +
  
  # Overlay previously removed experiments in darkorange (with legend entry)
  geom_point(
    data = tracer_num_vox_oob_vent_df_prev,
    aes(x = inj_num_vox, y = proj_num_vox, color = "Removed Knox et al."),
    size = 6,
    shape = 17
  ) +
  
  # Define all color categories
  scale_color_manual(
    name = "Removal Type",
    values = c(
      "None" = "black",
      "Large Inj/Proj." = "red",
      "Removed Knox et al." = "darkorange"
    )
  ) +
  
  geom_text_repel(aes(label = label), size = 15, na.rm = TRUE) +  # increase label size
  labs(
    x = "Injection Voxel Count",
    y = "Projection Voxel Count"
  ) +
  theme_minimal(base_size = 45)

# Save updated plot
ggsave("figures/figure_2b_scatterplot_with_prev.png", p, width = 30, height = 16, dpi = 300)

###########PLOT REPRESENTATIVE TRACERS########################
tracer <- 174957972
p <- plot_tracer_inj_slice(tracer, tracer_dir, -0.8, inj_thresh, allen_50um_template_path, allen_mask_path_50um)
ggsave("figures/figure_2b_174957972_inj.png", p, width=16, height=16, dpi=300)

p <- plot_tracer_proj_slice(tracer, tracer_dir, -0.8, proj_thresh, allen_50um_template_path, allen_mask_path_50um)
ggsave("figures/figure_2b_174957972_proj.png", p, width=16, height=16, dpi=300)

tracer <- 180436360
p <- plot_tracer_inj_slice(tracer, tracer_dir, 0, inj_thresh, allen_50um_template_path, allen_mask_path_50um)
ggsave("figures/figure_2b_180436360_inj.png", p, width=16, height=16, dpi=300)

p <- plot_tracer_proj_slice(tracer, tracer_dir, 0, proj_thresh, allen_50um_template_path, allen_mask_path_50um)
ggsave("figures/figure_2b_180436360_proj.png", p, width=16, height=16, dpi=300)

tracer <- 181057754
p <- plot_tracer_inj_slice(tracer, tracer_dir, -0.86, inj_thresh, allen_50um_template_path, allen_mask_path_50um)
ggsave("figures/figure_2b_181057754_inj.png", p, width=16, height=16, dpi=300)

p <- plot_tracer_proj_slice(tracer, tracer_dir, -0.86, proj_thresh, allen_50um_template_path, allen_mask_path_50um)
ggsave("figures/figure_2b_181057754_proj.png", p, width=16, height=16, dpi=300)

# add source column and combine
tracer_num_vox_oob_vent_df <- tracer_num_vox_oob_vent_df %>% mutate(source = "current")
tracer_num_vox_oob_vent_df_prev <- tracer_num_vox_oob_vent_df_prev %>% mutate(source = "prev")
combined_df <- bind_rows(tracer_num_vox_oob_vent_df, tracer_num_vox_oob_vent_df_prev)

##############2C: SMALL INJ DENSITIES##############################
sorted_df <- combined_df %>%
  arrange(inj_num_vox) %>%
  slice(1:50) %>%
  mutate(Rank = 1:50)

# bottom 9 (now include source info)
sorted_df <- sorted_df %>%
  mutate(ColorGroup = case_when(
    Rank <= 9 & source == "prev" ~ "Prev",
    Rank <= 9 & source == "current" ~ "Bottom 9",
    TRUE ~ "Others"
  ))

# cutoff based only on current tracers
inj_cutoff_val <- max(sorted_df$inj_num_vox[sorted_df$ColorGroup == "Bottom 9"], na.rm = TRUE)

p <- ggplot(sorted_df, aes(x = Rank, y = inj_num_vox, color = ColorGroup, shape = ColorGroup)) +
  geom_point(size = 8) +
  scale_color_manual(values = c("Bottom 9" = "red", "Prev" = "darkorange", "Others" = "blue")) +
  scale_shape_manual(values = c("Bottom 9" = 16, "Prev" = 17, "Others" = 16)) +
  geom_hline(yintercept = inj_cutoff_val, linetype = "dashed", color = "black") +
  geom_text_repel(
    data = subset(sorted_df, ColorGroup %in% c("Bottom 9", "Prev")),
    aes(label = tracer, color = ColorGroup),
    nudge_y = 0.2 * max(sorted_df$inj_num_vox),
    size = 15,
    show.legend = FALSE
  ) +
  labs(
    x = "Rank",
    y = "Number of Voxels",
    title = "Smallest Injection Densities (> 0.5)"
  ) +
  theme_minimal(base_size = 40)

ggsave("figures/figure_2c_smallest_inj.png", p, width = 20, height = 16, dpi = 300)


##############2C: SMALL PROJ DENSITIES##############################
sorted_df <- combined_df %>%
  arrange(proj_num_vox) %>%
  slice(1:50) %>%
  mutate(Rank = 1:50)

# bottom 4 (now include source info)
sorted_df <- sorted_df %>%
  mutate(ColorGroup = case_when(
    Rank <= 4 & source == "prev" ~ "Prev",
    Rank <= 4 & source == "current" ~ "Bottom 4",
    TRUE ~ "Others"
  ))

# cutoff based only on current tracers
proj_cutoff_val <- max(sorted_df$proj_num_vox[sorted_df$ColorGroup == "Bottom 4"], na.rm = TRUE)

p <- ggplot(sorted_df, aes(x = Rank, y = proj_num_vox, color = ColorGroup, shape = ColorGroup)) +
  geom_point(size = 8) +
  scale_color_manual(values = c("Bottom 4" = "red", "Prev" = "darkorange", "Others" = "blue")) +
  scale_shape_manual(values = c("Bottom 4" = 16, "Prev" = 17, "Others" = 16)) +
  geom_hline(yintercept = proj_cutoff_val, linetype = "dashed", color = "black") +
  geom_text_repel(
    data = subset(sorted_df, ColorGroup %in% c("Bottom 4", "Prev")),
    aes(label = tracer, color = ColorGroup),
    nudge_y = 0.2 * max(sorted_df$proj_num_vox),
    size = 15,
    show.legend = FALSE
  ) +
  labs(
    x = "Rank",
    y = "Number of Voxels",
    title = "Smallest Projection Densities (> 0.1)"
  ) +
  theme_minimal(base_size = 40)

ggsave("figures/figure_2c_smallest_proj.png", p, width = 20, height = 16, dpi = 300)




############PLOT REPRESENTATIVE SLICES###################################
tracer <- 277618054
p <- plot_tracer_inj_slice(tracer, tracer_dir, -5.7, inj_thresh, allen_50um_template_path, allen_mask_path_50um)
ggsave("figures/figure_2c_277618054_inj.png", p, width=16, height=16, dpi=300)

p <- plot_tracer_proj_slice(tracer, tracer_dir, -5.7, proj_thresh, allen_50um_template_path, allen_mask_path_50um)
ggsave("figures/figure_2c_277618054_proj.png", p, width=16, height=16, dpi=300)

tracer <- 180568155
p <- plot_tracer_inj_slice(tracer, tracer_dir, -0.44, inj_thresh, allen_50um_template_path, allen_mask_path_50um)
ggsave("figures/figure_2c_180568155_inj.png", p, width=16, height=16, dpi=300)

p <- plot_tracer_proj_slice(tracer, tracer_dir, -0.44, proj_thresh, allen_50um_template_path, allen_mask_path_50um)
ggsave("figures/figure_2c_180568155_proj.png", p, width=16, height=16, dpi=300)
