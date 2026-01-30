library("pacman")
library("RMINC")
pacman::p_load(dplyr, readr, readxl, stringr, ggrepel, tidyr, ggplot2, ggalluvial, tidyr)
setwd(".")
source("./ggslicer/plotting_functions/plotting_functions.R")
source("./ggslicer/plotting_functions/plotting_functions_labels.R")

###find the intended injection region for each tracer (and major division)

allen_input_dir <- "../preprocessed/allen_template_inputs/"
allen_tracer_dir <- "../preprocessed/knox_connectome_tracers/"
##50 microns
aba_label_filepath <- paste0(allen_input_dir, "AMBA_25um_resampled_50um_int.mnc")
aba_label_file <- round(mincGetVolume(aba_label_filepath))

##load primary injection sites for each experiment + filled in two missing experiments
#modified version of query.csv: see fill_in_missing_knox_tracer_regions.R
tracer_sites <- as.data.frame(read_csv(paste0(allen_input_dir,"knox_tracer_ids_inj_regions_full.csv")))

###organize region numbers --> names --> major subdivisions
#load the dictionary from label numbers <-- --> region acronyms
aba_region_filepath <- paste0(allen_input_dir, "allen_ccfv3_tree_wang_2020_s2.xlsx")
aba_region_labels <- as.data.frame(read_excel(aba_region_filepath))  
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

#############Compare manual vs. automated QC tracers (with regions)###########
##original Knox conn tracer list before removal
##load list of Knox connectome tracers
knox_experiments_included <- as.data.frame(read.csv("knox_experiment_csvs/knox_experiments_included.csv"))
###remove missing experiment
knox_experiments_included <- as.data.frame(knox_experiments_included[which(knox_experiments_included$id != 310207648),])
colnames(knox_experiments_included) <- "id"
##experiments to remove
tracer_removal_df <- as.data.frame(read_csv("tables/tracers_to_remove_adjboxStats_skew_outliers.csv", col_names=TRUE))

###add major_division column for each tracer
# join tracer_removal_df with tracer_sites by Tracer ↔ id
# If you’ll color axis text:
# install.packages("ggtext") # if needed
library(ggtext)

# --- 0) Your palette ---
iwant_hex <- c("#a83537","#4fc79c","#63348a","#73c161","#6280d6","#d3a046",
               "#ca78cd","#4c792a","#bb467a","#aab248","#a96126","#ff846b")

# --- 1) Join to get structure-id, then add major_division ---
tracer_joined <- tracer_removal_df %>%
  left_join(tracer_sites %>% select(id, `structure-id`),
            by = c("Tracer" = "id")) %>%
  rowwise() %>%
  mutate(
    major_division = find_major_division(`structure-id`,
                                         major_division_dict,
                                         aba_region_labels)
  ) %>%
  ungroup()

###repeat the same process for the included experiments
knox_joined_full <- knox_experiments_included %>%
  left_join(tracer_sites %>% select(id, `structure-id`),
            by = c("id" = "id")) %>%
  rowwise() %>%
  mutate(
    major_division = find_major_division(`structure-id`,
                                         major_division_dict,
                                         aba_region_labels)
  ) %>%
  ungroup()

##sort by major division
# Get the desired order from the column names of your major_division_dict
div_order <- colnames(major_division_dict)

# Sort knox_joined_full by major_division in that order
knox_joined_full_sorted <- knox_joined_full %>%
  mutate(major_division = factor(major_division, levels = div_order)) %>%
  arrange(major_division)

knox_joined_full <- knox_joined_full_sorted


# Keep this as the working df (you had "tracer_removal_df <- tracer_joined")
tracer_removal_df <- tracer_joined

# --- 2) Build color map for major_division ---
# If you also have fmd/fmd_col, you can include them in the union as you showed.
# If not, we just use the levels present in tracer_joined:
# --- Fixed legend/order from major_division_dict ---
# Preserve the column order exactly as written in your dict
div_lvls <- colnames(major_division_dict)

# Map colors to levels in dict order
div_cols <- setNames(iwant_hex[seq_along(div_lvls)], div_lvls)

# Make major_division a factor with dict order (controls legend order everywhere)
tracer_removal_df <- tracer_removal_df %>%
  mutate(major_division = factor(major_division, levels = div_lvls))

annotation_colors <- list(major_division = div_cols)


###########Sankey plot################
######################################
# --- 0. Identify QC / metric columns explicitly (fixes the error) ---
# Identify columns
auto_cols <- names(tracer_removal_df)[grepl("^Auto", names(tracer_removal_df))]
manual_cols <- names(tracer_removal_df)[grepl("^Manual", names(tracer_removal_df))]
metric_cols <- names(tracer_removal_df)[grepl("^(Auto|Manual)", names(tracer_removal_df))]
metric_cols
# sanity check: this should only include columns like "Auto OOB Proj.", "Manual Nonspecific Inj.", etc.

# --- 1. Compute first_removal robustly ---
first_removal_df <- tracer_removal_df %>%
  # keep other columns as-is, but compute first_removal based only on metric_cols
  rowwise() %>%
  mutate(
    # extract metrics for this row, coerce to numeric, handle NAs
    first_removal = {
      vals <- as.numeric(c_across(all_of(metric_cols)))
      # find first index where value == 1 (handle NAs)
      idx <- which(!is.na(vals) & vals == 1)
      if (length(idx) == 0) {
        "kept"
      } else {
        metric_cols[min(idx)]
      }
    }
  ) %>%
  ungroup()

# --- 2. Prepare aggregated alluvial dataframe ---
# Make sure major_division is character (or factor with desired levels)
first_removal_df <- first_removal_df %>%
  mutate(
    first_removal = factor(first_removal, levels = c(metric_cols, "kept")),
    major_division = as.character(major_division)  # convert factor -> character safely
  )

major_div_levels <- names(div_cols)    # from your earlier code
first_removal_df$major_division <- factor(first_removal_df$major_division,
                                          levels = major_div_levels)

# Aggregate counts
alluv_df_with_tracer <- first_removal_df %>%
  group_by(first_removal, major_division) %>%
  summarise(freq = n(), .groups = "drop") %>%
  filter(freq > 0)

# Aggregate counts
alluv_df <- first_removal_df %>%
  group_by(first_removal, major_division) %>%
  summarise(freq = n(), .groups = "drop") %>%
  filter(freq > 0)

# --- 3. Plot using ggalluvial (with n=... under left-hand major-division strata) ---

# Ensure major_division is a factor in the desired order
alluv_df <- alluv_df %>%
  mutate(major_division = factor(major_division, levels = names(div_cols)))

# --- Fix: reverse order so labels align correctly with strata positions ---
left_counts <- alluv_df %>%
  group_by(major_division) %>%
  summarise(n = sum(freq), .groups = "drop") %>%
  # reverse order so that first division (Isocortex) appears at top
  arrange(desc(factor(major_division, levels = names(div_cols)))) %>%
  mutate(
    bottom = lag(cumsum(n), default = 0),
    top    = bottom + n,
    center = (bottom + top) / 2,
    gap = max(n) * 0.1,
    label_y = center - gap
  )

left_counts <- left_counts[which(left_counts$n > 3),]

# --- Right-hand counts (QC failure mode totals) ---
right_counts <- alluv_df %>%
  group_by(first_removal) %>%
  summarise(n = sum(freq), .groups = "drop") %>%
  # reverse order so that first removal criterion appears at top
  arrange(desc(factor(first_removal, levels = c(metric_cols, "kept")))) %>%
  mutate(
    bottom = lag(cumsum(n), default = 0),
    top    = bottom + n,
    center = (bottom + top) / 2,
    gap = max(n) * 0.1,
    label_y = center - gap
  ) %>%
  filter(n > 3)   # same filtering rule as left_counts


stratum_labels <- ggplot_build(
  ggplot(alluv_df,
         aes(axis1 = major_division,
             axis2 = first_removal,
             y = freq)) +
    geom_stratum(width = 0.3)
)$data[[1]]


# Build the plot: left axis = Major division, right axis = first_removal (Removal criterion)
p <- ggplot(alluv_df,
            aes(axis1 = major_division,
                axis2 = first_removal,
                y = freq)) +

  geom_alluvium(aes(fill = major_division),
                width = 0.3,
                knot.pos = 0.5) +

  geom_stratum(width = 0.3,
               fill = "grey90",
               color = "black") +

  ## --- STRATUM LABELS (repelled, per ggalluvial vignette) ---
  geom_text_repel(
    data = stratum_labels,
    aes(x = x, y = y, label = stratum),
    inherit.aes = FALSE,
    size = 9,
    direction = "y",
    nudge_x = ifelse(stratum_labels$x == 1, -0.35, 0.35),
    segment.color = "grey40",
    box.padding = 0.3,
    min.segment.length = 0
  ) +

  ## --- LEFT n= labels ---
  geom_text_repel(
    data = left_counts,
    aes(x = 1, y = label_y, label = paste0("n=", n)),
    inherit.aes = FALSE,
    size = 7,
    nudge_x = -0.35,
    min.segment.length = 0,
    segment.color = NA


  ) +

  ## --- RIGHT n= labels ---
  geom_text_repel(
    data = right_counts,
    aes(x = 2, y = label_y, label = paste0("n=", n)),
    inherit.aes = FALSE,
    size = 7,
    nudge_x=0.35,
    min.segment.length = 0,
    segment.color = NA
  ) +

  scale_x_discrete(
    limits = c("Major division", "Removal criterion"),
    expand = c(.05, .05)
  ) +

  scale_fill_manual(values = div_cols, na.value = "grey50") +

  labs(
    y = "Number of tracers",
    x = NULL,
    title = "Sankey (alluvial) of tracers removed by first QC criterion, grouped by major division"
  ) +

  theme_void(base_size = 30) +
  theme(legend.position = "right")

p

ggsave("figures/fig_3a_removed_tracers_sankey.png", p, dpi = 300, width = 24, height = 14)



# =========================
# B) HEATMAP with Tracer labels colored by major_division
# =========================
# --- Fixes: define metric_levels and ensure plotting order matches tracer_levels_ordered ---

# ensure metric_levels exists (use metric_cols order)
metric_levels <- metric_cols

# compute first_removal per tracer if not already present (robust) 
if (!"first_removal" %in% names(tracer_removal_df))  { 
  tracer_removal_df <- tracer_removal_df %>% rowwise() %>%
    mutate( first_removal =
              { vals <- as.numeric(c_across(all_of(metric_cols))) 
              idx <- which(!is.na(vals) & vals == 1) 
              if (length(idx) == 0) "kept" else metric_cols[min(idx)] } ) %>% 
    ungroup() 
} 
# make ordering factors for sorting 
tracer_removal_df <- tracer_removal_df %>% mutate( major_division_f = 
                                                     factor(major_division, levels = major_div_levels), 
                                                   first_removal_f = factor(first_removal, levels = c(metric_cols, "kept")) ) 
# Build ordering: by major division (as in dict), then by first_removal (metric order), then by tracer id 
tracer_order_df <- tracer_removal_df %>% arrange(major_division_f, first_removal_f, Tracer) 

# final tracer level vector (in the order we want the rows) 
tracer_levels_ordered <- tracer_order_df$Tracer

# Build df_flagged and arrange rows by tracer_levels_ordered
df_flagged <- tracer_removal_df %>%
  mutate(
    only_manual = if_all(all_of(auto_cols),   ~ .x == 0) &
      if_any(all_of(manual_cols), ~ .x == 1)
  ) %>%
  arrange(match(Tracer, tracer_levels_ordered))

# For plotting with ggplot2 we must make the y variable a factor with levels in the desired order.
# Do this *only* for plotting (it doesn't harm other data manipulations, but it's required for ggplot).
df_flagged <- df_flagged %>%
  mutate(Tracer = factor(Tracer, levels = tracer_levels_ordered))

# pivot to long form
df_long <- df_flagged %>%
  pivot_longer(cols = -c(Tracer, only_manual, major_division, `structure-id`, first_removal,
                         major_division_f, first_removal_f),
               names_to = "Metric", values_to = "Value") %>%
  mutate(
    Metric = factor(Metric, levels = metric_levels),   # use the metric order you want
    type = dplyr::case_when(
      Value == 0 ~ "zero",
      Metric %in% auto_cols   & Value == 1 ~ "auto",
      Metric %in% manual_cols & Value == 1 &  only_manual ~ "manual_only",
      Metric %in% manual_cols & Value == 1 & !only_manual ~ "manual",
      TRUE ~ "zero"
    )
  )

# Build colored y-axis labels using tracer_order_df to keep consistent division colors
tracer_label_colors <- tracer_order_df %>%
  distinct(Tracer, major_division) %>%
  mutate(
    lbl = paste0("<span style='color:", div_cols[major_division], "'>", Tracer, "</span>")
  )

ylab_map <- setNames(tracer_label_colors$lbl, tracer_label_colors$Tracer)

# Now plotting should respect tracer_levels_ordered because Tracer is a factor with those levels
p <- ggplot(df_long, aes(x = Metric, y = Tracer, fill = type)) +
  geom_tile(color = "white") +
  scale_fill_manual(
    name = NULL,
    values = c(
      zero        = "white",
      auto        = "steelblue",
      manual      = "green",
      manual_only = "yellow"
    ),
    breaks = c("auto", "manual", "manual_only", "zero"),
    labels = c("Automated", "Manual", "Manual-only tracer", "0")
  ) +
  scale_y_discrete(labels = ylab_map, limits = rev(levels(df_long$Tracer))) +
  labs(x = NULL, y = "Tracer (colored by major division)",
       title = "Overall Tracers to Remove from Connectome") +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = ggtext::element_markdown(),  # enable colored labels
    panel.grid = element_blank()
  )

# Save if desired
ggsave("figures/supp_fig_1_heatmap.png", p, width = 24, height = 16, dpi = 300)

####################look at overall ABA regions with removed tracers###########################
df_with_structures <- tracer_removal_df 

# --- Count frequency of each structure-name ---
# structure_IDs has two columns: structure-id and n (the counts)
structure_IDs <- df_with_structures %>%
  count(`structure-id`) %>%
  arrange(desc(n))

# initialize with 0s
aba_label_file_missing_regions <- rep(0, length(aba_label_file))

# update with counts instead of 1
matched_idx <- match(aba_label_file, structure_IDs$`structure-id`)
aba_label_file_missing_regions[!is.na(matched_idx)] <- structure_IDs$n[matched_idx[!is.na(matched_idx)]]

#####repeat for overall tracer df
df_with_structures <- tracer_sites
structure_IDs <- df_with_structures %>%
  count(`structure-id`) %>%
  arrange(desc(n))

aba_label_file_knox_regions <- rep(0, length(aba_label_file))

# update with counts instead of 1
matched_idx <- match(aba_label_file, structure_IDs$`structure-id`)
aba_label_file_knox_regions[!is.na(matched_idx)] <- structure_IDs$n[matched_idx[!is.na(matched_idx)]]

####write out output volumes
vol_output_dir="../derivatives/excluded_tracer_aggregate_volumes/"
mincWriteVolume(aba_label_file_missing_regions, like=aba_label_filepath, paste0(vol_output_dir, "overall_excluded_tracer_regions.mnc"))
####save the RATIO of missing injections to total injections within each region
output_ratio_vol <- rep(0, length(aba_label_file))
output_ratio_vol[which(aba_label_file_knox_regions > 0)] <- 
  aba_label_file_missing_regions[which(aba_label_file_knox_regions > 0)]/aba_label_file_knox_regions[which(aba_label_file_knox_regions > 0)]
mincWriteVolume(output_ratio_vol, 
                like=aba_label_filepath, paste0(vol_output_dir, "overall_excluded_tracer_regions_ratio.mnc"))


####visualize using Yohan's visualization script
allen_template_path_50um <- paste0(allen_input_dir, "average_template_50.mnc")
allen_mask_path_50um <- paste0(allen_input_dir, "mask_50um.mnc")
allen_50um_template_df <- prepare_masked_anatomy(allen_template_path_50um, allen_mask_path_50um, "y", seq(-7.5, 5, 1))[[2]]
excluded_rgn_df <- prepare_masked_anatomy(paste0(vol_output_dir,"./overall_excluded_tracer_regions.mnc"), allen_mask_path_50um, "y", seq(-7.5, 5, 1))[[2]]

allen_50um_template_df <- allen_50um_template_df %>% filter(mask_value == 1)
excluded_rgn_df <- excluded_rgn_df %>% filter(mask_value == 1) %>% filter(intensity > 0)

# Combine "slice_world" with "y" for unique facet labels
excluded_rgn_df <- excluded_rgn_df %>% 
  mutate(facet_label = paste0("y = ", y))

p <- ggplot(data = allen_50um_template_df, mapping = aes(x = x, y = z)) +
  geom_raster(mapping = aes(fill = intensity), interpolate = TRUE) +
  scale_fill_gradient(
    low = 'black',
    high = 'white',
    oob = scales::squish,
    guide = 'none'
  ) +
  ggnewscale::new_scale_fill() + 
  geom_raster(data = excluded_rgn_df, mapping = aes(fill = factor(intensity))) +
  scale_fill_manual(
    values = c("1" = "yellow", "2" = "orange", "3" = "red"),
    breaks = c("1", "2", "3"),
    name = "Number of tracers removed",
    drop = FALSE
  ) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~ slice_world, ncol = 4, labeller = labeller(slice_world = function(x) paste0("Slice: ", x))) +  
  coord_fixed(ratio = 1) +
  theme_void(base_size=30) +
  theme(
    panel.spacing = unit(0, "npc"),                     # Keep spacing minimal
    strip.text = element_text(size = 16, face = "bold"), # Increase facet label text size
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5), # Make title bigger and center it
    plot.subtitle = element_text(size = 20, hjust = 0.5, face = "italic") # Make subtitle bigger and center it
  ) +
  labs(
    title = "Injection Regions with Excluded Tracers", # Add title
  )

ggsave("figures/supp_fig_2a_inj_rgn_counts.png", p, width = 24, height = 16, dpi = 300)  # adjust width/height as needed

#######REPEAT FOR RATIO######################
excluded_rgn_df <- prepare_masked_anatomy(paste0(vol_output_dir,"overall_excluded_tracer_regions_ratio.mnc"), allen_mask_path_50um, "y", seq(-7.5, 5, 1))[[2]]
excluded_rgn_df <- excluded_rgn_df %>% filter(mask_value == 1) %>% filter(intensity > 0)

# Combine "slice_world" with "y" for unique facet labels
excluded_rgn_df <- excluded_rgn_df %>% 
  mutate(facet_label = paste0("y = ", y))

p <- ggplot(data = allen_50um_template_df, mapping = aes(x = x, y = z)) +
  geom_raster(mapping = aes(fill = intensity), interpolate = TRUE) +
  scale_fill_gradient(
    low = 'black',
    high = 'white',
    oob = scales::squish,
    guide = 'none'
  ) +
  ggnewscale::new_scale_fill() + 
  geom_raster(data = excluded_rgn_df, mapping = aes(fill = intensity)) +
  scale_fill_gradient(
    name="Prop. exps.",
    low="yellow",
    high='red',
    guide="colourbar"
  ) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~ slice_world, ncol = 3, labeller = labeller(slice_world = function(x) paste0("Slice: ", x))) +  
  coord_fixed(ratio = 1) +
  theme_void(base_size=30) +
  theme(
    panel.spacing = unit(0, "npc"),                     # Keep spacing minimal
    strip.text = element_text(size = 16, face = "bold"), # Increase facet label text size
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5), # Make title bigger and center it
    plot.subtitle = element_text(size = 20, hjust = 0.5, face = "italic") # Make subtitle bigger and center it
  ) +
  labs(
    title = "Proportion of Removed Injections Per Region", # Add title
  )

ggsave("figures/figure_3b_inj_rgn_ratio.png", p, width = 24, height = 16, dpi = 300)  # adjust width/height as needed

###########Look at missing projection density (overall, binary > 0.1)##########
aba_missing_proj_vol <- rep(0, length(aba_label_file))
for(tracer in tracer_removal_df$Tracer)  {
  projection_density_bin0.1 <- mincArray(mincGetVolume(paste0(allen_tracer_dir, "/",tracer,"/projection_density_bin0.1.mnc")))
  aba_missing_proj_vol[which(projection_density_bin0.1 > 0)] <- aba_missing_proj_vol[which(projection_density_bin0.1 > 0)] +1
}

mincWriteVolume(aba_missing_proj_vol, paste0(vol_output_dir,"knox_lost_proj_density_bin0.1.mnc"), like=aba_label_filepath)
###count the overall number of experiments @ each voxel; divide the map above
aba_overall_proj_vol <- rep(0, length(aba_label_file))
for(tracer in knox_experiments_included$id)  {
  if(tracer == 310207648) { #missing data
    next
  }
  projection_density_bin0.1 <- mincArray(mincGetVolume(paste0(allen_tracer_dir, "/",tracer,"/projection_density_bin0.1.mnc")))
  aba_overall_proj_vol[which(projection_density_bin0.1 > 0)] <- aba_overall_proj_vol[which(projection_density_bin0.1 > 0)] +1
}

###calculate volume fraction of missing vs. non-missing experiments @ each voxel
aba_missing_proj_vol_ratio <- rep(0, length(aba_label_file))
aba_missing_proj_vol_ratio[which(aba_overall_proj_vol > 0)] <- aba_missing_proj_vol[which(aba_overall_proj_vol > 0)] / aba_overall_proj_vol[which(aba_overall_proj_vol > 0)]
mincWriteVolume(aba_missing_proj_vol_ratio, paste0(vol_output_dir,"knox_lost_proj_density_ratio_bin0.1.mnc"), like=aba_label_filepath)

###visualize missing proj. densities
excluded_rgn_df <- prepare_anatomy(paste0(vol_output_dir,"./knox_lost_proj_density_bin0.1.mnc"), "y", seq(-7.5, 5, 1))[[2]]
excluded_rgn_df <- excluded_rgn_df %>%filter(intensity > 0)

# Combine "slice_world" with "y" for unique facet labels
excluded_rgn_df <- excluded_rgn_df %>% 
  mutate(facet_label = paste0("y = ", y))

p <- ggplot(data = allen_50um_template_df, mapping = aes(x = x, y = z)) +
  geom_raster(mapping = aes(fill = intensity), interpolate = TRUE) +
  scale_fill_gradient(
    low = 'black',
    high = 'white',
    oob = scales::squish,
    guide = 'none'
  ) +
  ggnewscale::new_scale_fill() + 
  geom_raster(data = excluded_rgn_df, mapping = aes(fill = intensity)) +
  scale_fill_gradient(
    name="Num. tracers",
    low="yellow",
    high='red',
    guide="colourbar"
  ) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~ slice_world, ncol = 4, labeller = labeller(slice_world = function(x) paste0("Slice: ", x))) +  
  coord_fixed(ratio = 1) +
  theme_void(base_size=30) +
  theme(
    panel.spacing = unit(0, "npc"),                     # Keep spacing minimal
    strip.text = element_text(size = 16, face = "bold"), # Increase facet label text size
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5), # Make title bigger and center it
    plot.subtitle = element_text(size = 20, hjust = 0.5, face = "italic") # Make subtitle bigger and center it
  ) +
  labs(
    title = "Number of Experiments with Excluded Projection > 0.1", # Add title
  )
ggsave("figures/supp_fig_2b_proj_voxels.png", p, width=24, height=16, dpi=300)


###visualize lost projection density ratio @ each voxel
excluded_rgn_df <- prepare_anatomy(paste0(vol_output_dir,"./knox_lost_proj_density_ratio_bin0.1.mnc"), "y", seq(-7.5, 5, 1))[[2]]
excluded_rgn_df <- excluded_rgn_df %>%filter(intensity > 0)

# Combine "slice_world" with "y" for unique facet labels
excluded_rgn_df <- excluded_rgn_df %>% 
  mutate(facet_label = paste0("y = ", y))

p <- ggplot(data = allen_50um_template_df, mapping = aes(x = x, y = z)) +
  geom_raster(mapping = aes(fill = intensity), interpolate = TRUE) +
  scale_fill_gradient(
    low = 'black',
    high = 'white',
    oob = scales::squish,
    guide = 'none'
  ) +
  ggnewscale::new_scale_fill() + 
  geom_raster(data = excluded_rgn_df, mapping = aes(fill = intensity)) +
  scale_fill_gradient(
    name="Prop. exps.",
    low="yellow",
    high='red',
    guide="colourbar"
  ) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~ slice_world, ncol = 3, labeller = labeller(slice_world = function(x) paste0("Slice: ", x))) +  
  coord_fixed(ratio = 1) +
  theme_void(base_size=30) +
  theme(
    panel.spacing = unit(0, "npc"),                     # Keep spacing minimal
    strip.text = element_text(size = 16, face = "bold"), # Increase facet label text size
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5), # Make title bigger and center it
    plot.subtitle = element_text(size = 20, hjust = 0.5, face = "italic") # Make subtitle bigger and center it
  ) +
  labs(
    title = "Proportion of Experiments With Projection Voxels > 0.1 Excluded After QC", # Add title
  )
ggsave("figures/figure_3b_proj_voxel_ratio.png", p, width=24, height=16, dpi=300)

