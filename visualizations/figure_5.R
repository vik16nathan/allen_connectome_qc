setwd("..")
library("pacman")
library("RMINC")
pacman::p_load(tidyverse, pheatmap, dplyr, readr, readxl, patchwork, R.matlab, ggrepel, combinat)
source("./ggslicer/plotting_functions/plotting_functions.R")
source("./ggslicer/plotting_functions/plotting_functions_labels.R")

####Paths and parameters####
allen_input_dir <- "../preprocessed/allen_template_inputs/"

aba_label_filepath <- paste0(allen_input_dir, "AMBA_25um_resampled_50um_int.mnc")
aba_label_file <- round(mincGetVolume(aba_label_filepath))

knox_conn_old <- as.data.frame(read_csv("mouse_connectivity_models/paper/connectivity/voxel-standard-model/normalized_connection_strength_original.csv"))
knox_conn_region_numbers <- knox_conn_old[2:nrow(knox_conn_old), 1]
num_rgns <- length(knox_conn_region_numbers)

###major division dictionary (depth 4/5 in Allen region tree)
major_division_dict <- data.frame(
  Isocortex = 315, OLF = 698, HPF = 1089, CTXsp = 703, STR = 477, PAL = 803,
  Thal = 549, Hypothal = 1097, Midbrain = 313, Pons = 771,
  Medulla = 354, CB = 512
)

find_major_division <- function(label, major_division_dict, aba_region_labels) {
  sample_string <- aba_region_labels[which(aba_region_labels[, "structure ID"] == label), "structure_id_path"]
  region_hierarchy_list <- as.numeric(unlist(strsplit(sample_string, "/")))
  colnames(major_division_dict)[which(major_division_dict %in% region_hierarchy_list)]
}

###organize region numbers --> names --> major subdivisions
#load the dictionary from label numbers <--> region acronyms
aba_region_filepath <- paste0(allen_input_dir, "allen_ccfv3_tree_wang_2020_s2.xlsx")
aba_region_labels <- as.data.frame(read_excel(aba_region_filepath))
colnames(aba_region_labels) <- aba_region_labels[1, ]
aba_region_labels <- aba_region_labels[2:nrow(aba_region_labels), ]

#############FIND MAJOR DIVISIONS#############
knox_row_major_divisions <- c()
for (label in knox_conn_region_numbers) {
  knox_row_major_divisions <- c(knox_row_major_divisions,
                                find_major_division(label, major_division_dict, aba_region_labels))
}

#############HELPER FUNCTION - CONVERT KNOX --> LABEL FILE REGIONS############
## create a data frame converting Allen connectome regions --> label file regions
## one mesoscale region can correspond to multiple label file regions (children)
metadata_comm_rich_club_old <- as.data.frame(read_csv("../derivatives/community_louvain/metadata_community_rich_club_old_strength_pct80.csv"))

mesoscale_rgn_to_label_file_rgn_dict <- data.frame(matrix(ncol = 2))
colnames(mesoscale_rgn_to_label_file_rgn_dict) <- c("mesoscale region", "label file region")

for (rgn in unique(metadata_comm_rich_club_old$names)) {
  if (rgn %in% aba_label_file) {
    mesoscale_rgn_to_label_file_rgn_dict <- rbind(mesoscale_rgn_to_label_file_rgn_dict, c(rgn, rgn))
  } else {
    # find the region's children
    children <- as.numeric(aba_region_labels[which(aba_region_labels$parent_id == rgn), "structure ID"])
    mesoscale_rgn_to_label_file_rgn_dict <- rbind(mesoscale_rgn_to_label_file_rgn_dict,
                                                   c(rgn, paste(children, collapse = ",")))
  }
}

mesoscale_rgn_to_label_file_rgn_dict <- mesoscale_rgn_to_label_file_rgn_dict[2:nrow(mesoscale_rgn_to_label_file_rgn_dict), ]
write_csv(mesoscale_rgn_to_label_file_rgn_dict,
          "../derivatives/community_louvain/mesoscale_rgn_to_label_file_rgn_dict.csv")

## find the size of each region
mesoscale_rgn_to_size_dict <- data.frame(matrix(nrow = 0, ncol = 2))
for (i in seq_len(nrow(mesoscale_rgn_to_label_file_rgn_dict))) {
  rgn <- mesoscale_rgn_to_label_file_rgn_dict[i, "mesoscale region"]
  subregion_list <- as.numeric(unlist(strsplit(mesoscale_rgn_to_label_file_rgn_dict[i, "label file region"], ",")))
  total_num_voxels <- 0
  for (sr in subregion_list) {
    total_num_voxels <- total_num_voxels + length(which(aba_label_file == sr))
  }
  mesoscale_rgn_to_size_dict <- rbind(mesoscale_rgn_to_size_dict,
                                      c(as.numeric(rgn), as.numeric(total_num_voxels)))
}
colnames(mesoscale_rgn_to_size_dict) <- c("mesoscale region", "num voxels")

#####################FIG 5A: RICH CLUBS#######################################
############LOAD MAT FILES FOR OLD/NEW RICH CLUBS##############################

## calculate separately for connection strength vs. connection density
rich_club_dir <- "../derivatives/rich_club/"

for (suffix in c("density_pct80", "strength_pct80")) {

  num_rgns <- 291

  ### Rich club percentages for each node
  pct_rich_club_old <- readMat(paste0(rich_club_dir, "knox_conn_old_", suffix, "_richOrNot_pctg.mat"))
  pct_rich_club_new <- readMat(paste0(rich_club_dir, "knox_conn_new_", suffix, "_richOrNot_pctg.mat"))

  ### topological rich clubs using "average" significant rich club degree
  topo_rich_club_old <- readMat(paste0(rich_club_dir, "knox_conn_old_", suffix, "_richOrNot_topo.mat"))
  topo_rich_club_new <- readMat(paste0(rich_club_dir, "knox_conn_new_", suffix, "_richOrNot_topo.mat"))

  ### ASSERT symmetry between ipsi/contra
  pct_rich_club_old$richOrNot.mult.final.ids.pctg[1:num_rgns] == pct_rich_club_old$richOrNot.mult.final.ids.pctg[(num_rgns + 1):(2 * num_rgns)]
  pct_rich_club_new$richOrNot.mult.final.ids.pctg[1:num_rgns] == pct_rich_club_new$richOrNot.mult.final.ids.pctg[(num_rgns + 1):(2 * num_rgns)]
  topo_rich_club_old$richOrNot[1:num_rgns] == topo_rich_club_old$richOrNot[(num_rgns + 1):(2 * num_rgns)]
  topo_rich_club_new$richOrNot[1:num_rgns] == topo_rich_club_new$richOrNot[(num_rgns + 1):(2 * num_rgns)]

  rich_club_df <- tibble(region = knox_conn_region_numbers) %>%
    mutate(
      pct_rich_club_old = as.numeric(pct_rich_club_old$richOrNot.mult.final.ids.pctg[1:num_rgns]),
      pct_rich_club_new = as.numeric(pct_rich_club_new$richOrNot.mult.final.ids.pctg[1:num_rgns]),
      rich_or_not_old   = as.numeric(topo_rich_club_old$richOrNot[1:num_rgns]),
      rich_or_not_new   = as.numeric(topo_rich_club_new$richOrNot[1:num_rgns])
    )

  ### diff: whether a node is in the old vs. new rich club
  rich_club_df <- rich_club_df %>%
    mutate(diff = as.integer(rich_or_not_new - rich_or_not_old))

  # make sure structure ID is numeric to match
  aba_region_labels <- aba_region_labels %>%
    mutate(`structure ID` = as.numeric(`structure ID`))

  # join full structure name into rich_club_df
  rich_club_df <- rich_club_df %>%
    left_join(aba_region_labels %>%
                select(`structure ID`, `full structure name`),
              by = c("region" = "structure ID")) %>%
    mutate(label_full = paste(`full structure name`))

  ### highlight nodes with highest residuals relative to y = x
  rich_club_df <- rich_club_df %>%
    mutate(resid = abs(pct_rich_club_new - pct_rich_club_old))

  N <- 10
  top_resid <- rich_club_df %>%
    slice_max(order_by = resid, n = N)

  rich_club_df <- rich_club_df %>%
    mutate(is_highlight = label_full %in% top_resid$label_full)

  # add major_division to each row
  rich_club_df <- rich_club_df %>%
    mutate(
      region = as.numeric(region),
      major_division = vapply(region, function(lbl)
        find_major_division(lbl, major_division_dict, aba_region_labels), character(1))
    )

  # point colors with correct precedence
  rich_club_df <- rich_club_df %>%
    mutate(
      point_color = case_when(
        diff == -1     ~ "blue",
        diff ==  1     ~ "red",
        is_highlight   ~ as.character(major_division),
        TRUE           ~ "other"
      )
    )

  top_resid <- inner_join(top_resid, rich_club_df)

  # iWantHue palette
  iwant_hex <- c("#a83537", "#4fc79c", "#63348a", "#73c161", "#6280d6",
                 "#d3a046", "#ca78cd", "#4c792a", "#bb467a", "#aab248",
                 "#a96126", "#ff846b")

  md    <- knox_row_major_divisions
  fmd   <- if (is.factor(md)) factor(md, levels = levels(md)) else factor(md, levels = unique(md))
  div_lvls <- unique(as.character(fmd))
  div_cols <- setNames(iwant_hex[seq_along(div_lvls)], div_lvls)

  annotation_colors <- list(major_division = div_cols)
  pal      <- c("other" = "grey30", div_cols)
  pal_base <- c(pal, "red" = "red", "blue" = "blue")

  # join region sizes into plot df
  df_plot <- rich_club_df %>%
    left_join(mesoscale_rgn_to_size_dict %>%
                rename(region = `mesoscale region`), by = "region") %>%
    mutate(
      `num voxels` = ifelse(is.na(`num voxels`), 1, `num voxels`),
      plot_size    = `num voxels`
    )

  # base scatter plot
  p_base <- ggplot(df_plot, aes(x = pct_rich_club_old, y = pct_rich_club_new)) +
    geom_point(aes(
      color = point_color,
      size  = plot_size,
      alpha = ifelse(point_color == "other", 0.2, 1)
    )) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40") +
    scale_size_continuous(range = c(1, 6), guide = "none") +
    scale_alpha_identity() +
    scale_color_manual(values = pal_base, guide = "none") +
    scale_x_continuous(limits = c(-0.05, 1.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.05, 1.5), expand = c(0, 0)) +
    labs(x = "Old rich club percentage", y = "New rich club percentage") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none")

  df_plot_subset  <- df_plot[which(df_plot$diff != 0 & !df_plot$region %in% top_resid$region), ]
  combined_filt_df <- merge(x = df_plot_subset, y = top_resid, all = TRUE)

  p_full_boxed <- p_base +
    geom_text_repel(
      data               = combined_filt_df,
      aes(label = label_full, color = point_color),
      size               = 3.5,
      point.size         = 5,
      box.padding        = 1.2,
      point.padding      = 0.8,
      force              = 20,
      force_pull         = 0.1,
      min.segment.length = 0,
      segment.size       = 0.4,
      max.overlaps       = Inf,
      max.iter           = 20000,
      xlim               = c(-0.05, 1.45),
      ylim               = c(-0.05, 1.45),
      seed               = 42
    )

  ggsave(paste0("figures/figure_5a_rich_club_pct_", suffix, ".png"), p_full_boxed,
         width = 8, height = 6, units = "in")

  ###########5B: Visualize old vs. new rich club regions on brain slices#########

  aba_label_file_rich_club_regions <- rep(0, length(aba_label_file))

  shared_regions <- as.numeric(as.matrix(as.data.frame(
    rich_club_df[which(rich_club_df$rich_or_not_old == 1 & rich_club_df$rich_or_not_new == 1), "region"]
  )))
  old_regions <- as.numeric(as.matrix(as.data.frame(
    rich_club_df[which(rich_club_df$rich_or_not_old == 1 & rich_club_df$rich_or_not_new == 0), "region"]
  )))
  new_regions <- as.numeric(as.matrix(as.data.frame(
    rich_club_df[which(rich_club_df$rich_or_not_old == 0 & rich_club_df$rich_or_not_new == 1), "region"]
  )))

  aba_label_file_rich_club_regions[which(aba_label_file %in% shared_regions)] <- 1
  aba_label_file_rich_club_regions[which(aba_label_file %in% old_regions)]    <- 2
  aba_label_file_rich_club_regions[which(aba_label_file %in% new_regions)]    <- 3

  # fill in missing regions via children
  missing_rgn_indices      <- which(!shared_regions %in% aba_label_file)
  shared_regions_children  <- aba_region_labels[
    which(aba_region_labels$parent_id %in% shared_regions[missing_rgn_indices]), "structure ID"
  ]
  aba_label_file_rich_club_regions[which(aba_label_file %in% shared_regions_children)] <- 1
  which(!old_regions %in% aba_label_file)  ## empty
  which(!new_regions %in% aba_label_file)  ## empty

  rich_club_diff_shared_path <- paste0(rich_club_dir, "rich_club_strength_regions_shared_old_new_", suffix, ".mnc")
  mincWriteVolume(aba_label_file_rich_club_regions, like = aba_label_filepath, rich_club_diff_shared_path)

  allen_template_path_50um <- paste0(allen_input_dir, "average_template_50.mnc")
  allen_mask_path_50um     <- paste0(allen_input_dir, "mask_50um.mnc")

  allen_50um_template_df <- prepare_masked_anatomy(allen_template_path_50um, allen_mask_path_50um, "y", seq(-7.5, 5, 1))[[2]]
  rich_club_vis_df        <- prepare_masked_anatomy(rich_club_diff_shared_path, allen_mask_path_50um, "y", seq(-7.5, 5, 1))[[2]]

  allen_50um_template_df <- allen_50um_template_df %>% filter(mask_value == 1)
  rich_club_vis_df        <- rich_club_vis_df %>% filter(mask_value == 1) %>% filter(intensity > 0)
  rich_club_vis_df$intensity <- as.integer(rich_club_vis_df$intensity)

  p <- ggplot(data = allen_50um_template_df, mapping = aes(x = x, y = z)) +
    geom_raster(mapping = aes(fill = intensity), interpolate = TRUE) +
    scale_fill_gradient(low = "black", high = "white", oob = scales::squish, guide = "none") +
    ggnewscale::new_scale_fill() +
    geom_raster(data = rich_club_vis_df, aes(fill = factor(intensity, levels = c(1, 2, 3)))) +
    scale_fill_manual(
      values = c("1" = "darkgreen", "2" = "blue", "3" = "red"),
      breaks = c("1", "2", "3"),
      labels = c("shared", "lost", "gained"),
      name   = "Rich club status",
      drop   = FALSE
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_wrap(~ slice_world, ncol = 4,
               labeller = labeller(slice_world = function(x) paste0("Slice: ", x))) +
    coord_fixed(ratio = 1) +
    labs(title = "Rich Club Regions") +
    theme_void(base_size = 30) +
    theme(
      panel.spacing  = unit(0, "npc"),
      strip.text     = element_text(size = 16, face = "bold"),
      plot.title     = element_text(size = 24, face = "bold", hjust = 0.5),
      plot.subtitle  = element_text(size = 20, hjust = 0.5, face = "italic"),
      panel.grid     = element_blank()
    )

  ggsave(paste0("figures/figure_5b_rich_club_regions_", suffix, ".png"), p,
         width = 16, height = 16, dpi = 300)
}

########################5C: Louvain Community Analysis#########################
## load old vs. new community assignments

for (suffix in c("density_pct80", "strength_pct80")) {

  print(paste("Suffix:", suffix))

  metadata_comm_rich_club_old <- as.data.frame(read_csv(paste0("../derivatives/community_louvain/metadata_community_rich_club_old_", suffix, ".csv")))
  metadata_comm_rich_club_new <- as.data.frame(read_csv(paste0("../derivatives/community_louvain/metadata_community_rich_club_new_", suffix, ".csv")))

  ### write out .mnc files with community assignments in old/new connectomes
  aba_label_file_communities_old <- rep(0, length(aba_label_file))
  aba_label_file_communities_new <- rep(0, length(aba_label_file))

  for (comm in unique(metadata_comm_rich_club_old$comm_assign)) {
    metadata_comm_rgns <- as.numeric(metadata_comm_rich_club_old[
      which(metadata_comm_rich_club_old$comm_assign == comm), "names"
    ])
    label_file_comm_rgns <- paste0(
      mesoscale_rgn_to_label_file_rgn_dict[
        which(mesoscale_rgn_to_label_file_rgn_dict$`mesoscale region` %in% metadata_comm_rgns),
        "label file region"
      ],
      collapse = ","
    )
    label_file_comm_rgns <- as.numeric(unlist(strsplit(label_file_comm_rgns, split = ",")))
    print(sum(aba_label_file_communities_old[which(aba_label_file %in% label_file_comm_rgns)]))
    aba_label_file_communities_old[which(aba_label_file %in% label_file_comm_rgns)] <- comm
  }

  # confirm communities 1 and 5 are ipsi/contra duplicates
  metadata_comm_rgns_1 <- as.numeric(metadata_comm_rich_club_old[which(metadata_comm_rich_club_old$comm_assign == 1), "names"])
  metadata_comm_rgns_5 <- as.numeric(metadata_comm_rich_club_old[which(metadata_comm_rich_club_old$comm_assign == 5), "names"])
  metadata_comm_rgns_5[which(!metadata_comm_rgns_5 %in% metadata_comm_rgns_1)]
  metadata_comm_rgns_1[which(!metadata_comm_rgns_1 %in% metadata_comm_rgns_5)]

  for (comm in unique(metadata_comm_rich_club_new$comm_assign)) {
    metadata_comm_rgns <- as.numeric(metadata_comm_rich_club_new[
      which(metadata_comm_rich_club_new$comm_assign == comm), "names"
    ])
    label_file_comm_rgns <- paste0(
      mesoscale_rgn_to_label_file_rgn_dict[
        which(mesoscale_rgn_to_label_file_rgn_dict$`mesoscale region` %in% metadata_comm_rgns),
        "label file region"
      ],
      collapse = ","
    )
    label_file_comm_rgns <- as.numeric(unlist(strsplit(label_file_comm_rgns, split = ",")))
    aba_label_file_communities_new[which(aba_label_file %in% label_file_comm_rgns)] <- comm
  }

  mincWriteVolume(aba_label_file_communities_old,
                  paste0("../derivatives/community_louvain/aba_label_file_communities_old_", suffix, ".mnc"),
                  like = aba_label_filepath)

  new_comm_path <- paste0("../derivatives/community_louvain/aba_label_file_communities_new_", suffix, ".mnc")
  mincWriteVolume(aba_label_file_communities_new, new_comm_path, like = aba_label_filepath)

  ### Merge communities 1 and 5 (ipsi/contra are interchangeable)
  metadata_comm_rich_club_old[which(metadata_comm_rich_club_old$comm_assign == 1), "comm_assign"] <- 5
  metadata_comm_rich_club_new[which(metadata_comm_rich_club_new$comm_assign == 1), "comm_assign"] <- 5

  ### Bug fix: region 186 is duplicated across communities 2 and 5 in old clustering; assign to 5
  metadata_comm_rich_club_old[which(metadata_comm_rich_club_old$names == 186), "comm_assign"] <- 5



  # ── AUTO-REMAP: find the bijective label permutation that best aligns new → old ──
  shared_names <- intersect(metadata_comm_rich_club_old$names,
                            metadata_comm_rich_club_new$names)
  old_matched  <- metadata_comm_rich_club_old$comm_assign[
    match(shared_names, metadata_comm_rich_club_old$names)
  ]
  new_matched  <- metadata_comm_rich_club_new$comm_assign[
    match(shared_names, metadata_comm_rich_club_new$names)
  ]

  comm_labels <- sort(unique(metadata_comm_rich_club_old$comm_assign))  # e.g. 2, 3, 4, 5

  # co-occurrence matrix: rows = new label, cols = old label
  co_mat <- table(
    factor(new_matched, levels = comm_labels),
    factor(old_matched, levels = comm_labels)
  )

  # brute-force optimal bijective mapping over all permutations
  best_score   <- -1
  best_mapping <- comm_labels  # identity fallback

  for (perm in combinat::permn(comm_labels)) {
    perm <- unlist(perm)
    score <- sum(diag(co_mat[as.character(perm), as.character(comm_labels)]))
    if (score > best_score) {
      best_score   <- score
      best_mapping <- perm
    }
  }

  # mapping_vec: named vector where name = new label, value = remapped label
  mapping_vec <- setNames(comm_labels, as.character(best_mapping))

  cat(sprintf("[%s] Best remap (new → old): %s → %s  (%d/%d matched)\n",
              suffix,
              paste(names(mapping_vec), collapse = ","),
              paste(mapping_vec,        collapse = ","),
              best_score, length(shared_names)))

  # apply remapping to new community assignments
  metadata_comm_rich_club_new$comm_assign <-
    mapping_vec[as.character(metadata_comm_rich_club_new$comm_assign)]
  # ── END AUTO-REMAP ──────────────────────────────────────────────────────────

  ### identify regions whose community assignment changed between old and new
  diff_comm_rgns_old <- unique(metadata_comm_rich_club_old[
    which(metadata_comm_rich_club_old$comm_assign != metadata_comm_rich_club_new$comm_assign),
    c("names", "comm_assign")
  ])
  diff_comm_rgns_new <- unique(metadata_comm_rich_club_new[
    which(metadata_comm_rich_club_old$comm_assign != metadata_comm_rich_club_new$comm_assign),
    c("names", "comm_assign")
  ])

  ### If a region appears twice (ipsi and contra rows had different comm_assign values),
  ### keep only the first occurrence (ipsilateral) in each data frame.
  diff_comm_rgns_old <- diff_comm_rgns_old[!duplicated(diff_comm_rgns_old$names), ]
  diff_comm_rgns_new <- diff_comm_rgns_new[!duplicated(diff_comm_rgns_new$names), ]

  diff_rgn_names <- aba_region_labels$`full structure name`[
    match(diff_comm_rgns_old$names, aba_region_labels$`structure ID`)
  ]

  diff_comm_rgns <- cbind(diff_rgn_names, diff_comm_rgns_old, diff_comm_rgns_new$comm_assign)
  colnames(diff_comm_rgns) <- c("full structure name", "structure ID", "old community", "new community")
  print(diff_comm_rgns)

  ### write .mnc file flagging regions with changed community assignments
  diff_rgn_file <- rep(0, length(aba_label_file))
  label_file_diff_rgns <- paste0(
    mesoscale_rgn_to_label_file_rgn_dict[
      which(mesoscale_rgn_to_label_file_rgn_dict$`mesoscale region` %in% diff_comm_rgns$`structure ID`),
      "label file region"
    ],
    collapse = ","
  )
  label_file_diff_rgns <- as.numeric(unlist(strsplit(label_file_diff_rgns, split = ",")))
  diff_rgn_file[which(aba_label_file %in% label_file_diff_rgns)] <- 1

  diff_rgn_output_path <- paste0("../derivatives/community_louvain/rgns_with_diff_communities_rebuilt_", suffix, ".mnc")
  mincWriteVolume(diff_rgn_file, diff_rgn_output_path, like = aba_label_filepath)

  ### visualize communities + differences on brain slices
  allen_template_path_50um <- paste0(allen_input_dir, "average_template_50.mnc")
  allen_mask_path_50um     <- paste0(allen_input_dir, "mask_50um.mnc")

  allen_50um_template_df <- prepare_masked_anatomy(allen_template_path_50um, allen_mask_path_50um, "y", seq(-7, 5, 1))[[2]]
  community_vis_df        <- prepare_masked_anatomy(new_comm_path, allen_mask_path_50um, "y", seq(-7, 5, 1))[[2]]
  community_diff_df       <- prepare_masked_anatomy(diff_rgn_output_path, allen_mask_path_50um, "y", seq(-7, 5, 1))[[2]]

  allen_50um_template_df <- allen_50um_template_df %>% filter(mask_value == 1)
  community_vis_df        <- community_vis_df %>% filter(mask_value == 1) %>% filter(intensity > 0)
  community_diff_df       <- community_diff_df %>% filter(mask_value == 1) %>% filter(intensity > 0)

  p <- ggplot(data = allen_50um_template_df, mapping = aes(x = x, y = z)) +
    geom_raster(mapping = aes(fill = intensity), interpolate = TRUE) +
    scale_fill_gradient(low = "black", high = "white", oob = scales::squish, guide = "none") +
    ggnewscale::new_scale_fill() +
    geom_raster(
      data = community_vis_df,
      aes(fill = factor(intensity, levels = c(2, 3, 4, 5)))
    ) +
    scale_fill_manual(
      values = if (suffix == "strength_pct80") {
        c("2" = "#fe345b", "3" = "#8b3355", "4" = "#3ee094", "5" = "#baad00")
      } else {
        c("2" = "#fe345b", "3" = "#3ee094", "4" = "#8b3355", "5" = "#baad00")
      },
      breaks = c("2", "3", "4", "5"),
      labels = c("Brainstem/Hip", "Midbrain", "CB/Hindbrain", "CTX/Olf"),
      name   = "Community status",
      drop   = FALSE
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_wrap(~ slice_world, ncol = 4,
               labeller = labeller(slice_world = function(x) paste0("Slice: ", x))) +
    coord_fixed(ratio = 1) +
    labs(title = "Louvain Communities") +
    ggnewscale::new_scale_fill() +
    geom_tile(
      data = community_diff_df,
      aes(fill = factor(intensity, levels = c(1)))
    ) +
    scale_fill_manual(
      values = c("1" = "lavender"),
      breaks = c("1"),
      labels = c("diff rgn"),
      name   = "Community diff.",
      drop   = FALSE
    ) +
    theme_void(base_size = 30) +
    theme(
      panel.spacing  = unit(0, "npc"),
      strip.text     = element_text(size = 16, face = "bold"),
      plot.title     = element_text(size = 24, face = "bold", hjust = 0.5),
      plot.subtitle  = element_text(size = 20, hjust = 0.5, face = "italic"),
      panel.grid     = element_blank()
    )

  ggsave(paste0("figures/figure_5b_communities_", suffix, ".png"), p,
         width = 16, height = 16, dpi = 300)
}
