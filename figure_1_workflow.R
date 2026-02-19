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
allen_50um_mask_path <- "../preprocessed/allen_template_inputs/mask_50um.mnc"
##directory containing .mnc files
tracer_dir <- "../preprocessed/knox_connectome_tracers/"

###########1B scatterplot schematic (removing largest/smallest inj/proj)#####
##################FIGURE 2B###############################################
tracer_num_vox_oob_vent_df <- as.data.frame(read.csv(paste0("tables/knox_oob_vent_df_inj",inj_thresh,"_proj",proj_thresh,".csv")))
tracer_num_vox_oob_vent_df <- tracer_num_vox_oob_vent_df[which(tracer_num_vox_oob_vent_df$tracer != 310207648),]
# Step 1: Identify SD outliers (unchanged)
zscore_thresh <- 4
tracer_num_vox_oob_vent_df <- tracer_num_vox_oob_vent_df %>%
  mutate(
    sd_outlier = abs(scale(inj_num_vox)) > zscore_thresh | abs(scale(proj_num_vox)) > zscore_thresh
  )

# Step 2: Flag "bottom-corner" points
tracer_num_vox_oob_vent_df <- tracer_num_vox_oob_vent_df %>%
  mutate(
    inj_rank  = dplyr::min_rank(inj_num_vox),
    proj_rank = dplyr::min_rank(proj_num_vox),
    bottom_corner = replace_na(inj_rank <= 4 | proj_rank <= 2, FALSE)
  )

# Step 3: Plot — draw blue first, then red on top
p <- ggplot(tracer_num_vox_oob_vent_df, aes(x = inj_num_vox, y = proj_num_vox)) +
  # base layer: all "normal" tracers in blue
  geom_point(data = subset(tracer_num_vox_oob_vent_df, !sd_outlier & !bottom_corner),
             color = "blue", size = 8) +
  # overlay: red points (both >4 SD and bottom corner)
  geom_point(data = subset(tracer_num_vox_oob_vent_df, sd_outlier | bottom_corner),
             color = "red", size = 8) +
  labs(
    x = "Injection Voxel Count",
    y = "Projection Voxel Count"
  ) +
  theme_minimal(base_size = 60) +
  theme(legend.position = "none")  # remove legend completely
ggsave("figures/figure_1_scatterplot_schematic.png", p, width = 20, height = 16, dpi = 300)

##############PLOT BRAIN MASK/VENTRICLE MASK#######################

###brain mask
slice_coords <-  seq(-1, -1, 1.5)
template_path <- allen_50um_template_path
allen_template_df <- prepare_masked_anatomy(template_path, allen_50um_mask_path, "y", slice_coords)[[2]]
allen_template_df <- allen_template_df %>% filter(mask_value == 1)
mask_df <- prepare_anatomy(allen_50um_mask_path, "y", slice_coords)[[2]]
mask_df <- mask_df %>% filter(intensity > 0)

ventricle_df <- prepare_anatomy("../preprocessed/allen_template_inputs/vent_voxels_label_file_50um_filt.mnc", "y", slice_coords)[[2]]
ventricle_df <- ventricle_df %>% filter(intensity > 0)

# Combine "slice_world" with "y" for unique facet labels
mask_df <- mask_df %>% 
  mutate(facet_label = paste0("y = ", y))

p <- ggplot(data = allen_template_df, mapping = aes(x = x, y = z)) +
  geom_raster(mapping = aes(fill = intensity), interpolate = TRUE) +
  scale_fill_gradient(
    low = 'black',
    high = 'white',
    oob = scales::squish,
    guide = 'none'
  ) +
  ggnewscale::new_scale_fill() + 
  geom_tile(data = mask_df, mapping = aes(fill = intensity), alpha=0.5) +
  scale_fill_gradient2(

    low = "brown",
    mid = "brown",
    high = "brown",
    space = "Lab",
    na.value = "grey50",
    transform = "identity",
    guide = "colourbar",
    aesthetics = "fill"
  ) + 
  ggnewscale::new_scale_fill() + 
  geom_tile(data = ventricle_df, mapping = aes(fill = intensity)) +
  scale_fill_gradient2(
    
    low = "cyan",
    mid = "cyan",
    high = "cyan",
    space = "Lab",
    na.value = "grey50",
    transform = "identity",
    guide = "colourbar",
    aesthetics = "fill"
  ) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~ slice_world, ncol = 4, labeller = labeller(slice_world = function(x) paste0("Slice: ", x, " mm"))) +  
  coord_fixed(ratio = 1) +
  theme_void() +
  theme(
    panel.spacing = unit(0, "npc"),                     # Keep spacing minimal
    strip.text = element_text(size = 48, face = "bold"), # Increase facet label text size
    plot.title = element_text(size = 54, face = "bold", hjust = 0.5), # Make title bigger and center it
  ) 

ggsave("figures/fig1c_mask_vent.png", p, width=16, height=16, dpi=300)

#####functions to plot tracer + template at slice coordinate using Yohan's vis#####
####################################################################################
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
      strip.text = element_text(size = 48, face = "bold"), # Increase facet label text size
      plot.title = element_text(size = 54, face = "bold", hjust = 0.5), # Make title bigger and center it
    ) +
    labs(
      title = tracer # Add title
    )
  
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
      strip.text = element_text(size = 48, face = "bold"), # Increase facet label text size
      plot.title = element_text(size = 54, face = "bold", hjust = 0.5), # Make title bigger and center it
    ) +
    labs(
      title = tracer # Add title
    )
  
}



#####plot representative slices for each QC failure mode 
tracer <- 158914182
p <- plot_tracer_inj_slice(tracer, tracer_dir, -3.04, 0.5, allen_50um_template_path, allen_50um_mask_path)
ggsave("figures/figure_1e_158914182.png", p, width = 16, height = 16, dpi = 300)


tracer <- 100148443
p <- plot_tracer_proj_slice(tracer, tracer_dir, -3.6, 0.1, allen_50um_template_path, allen_50um_mask_path)
ggsave("figures/figure_1g_100148443.png", p, width = 16, height = 16, dpi = 300)

tracer <- 112670853
p <- plot_tracer_proj_slice(tracer, tracer_dir, 0.9, 0.1, allen_50um_template_path, allen_50um_mask_path)
ggsave("figures/figure_1h_112670853.png", p, width = 16, height = 16, dpi = 300)

tracer <- 158257355
p <- plot_tracer_proj_slice(tracer, tracer_dir, -2.75, 0.1, allen_50um_template_path, allen_50um_mask_path)
ggsave("figures/figure_1f_158257355.png", p, width = 16, height = 16, dpi = 300)

tracer <- 174957972
p <- plot_tracer_proj_slice(tracer, tracer_dir, -1, 0.1, allen_50um_template_path, allen_50um_mask_path)
ggsave("figures/figure_1d_174957972.png", p, width = 16, height = 16, dpi = 300)