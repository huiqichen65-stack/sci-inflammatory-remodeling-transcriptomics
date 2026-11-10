# ============================================================
# 82_Figure1_GLIA_resident_microglia_landscape_v29_CFleft_Elabel_fix.R
# Figure 1 | Resident microglial purification and injury-state landscape
# v29 final spacing refinement: move the top legends of panels C/F further left
# and nudge only the 'Resident-derived microglia' label in Panel E further
# toward the lower-right for cleaner separation from the panel letter.
# VISUAL-ONLY REVISION: no QC, embedding, clustering, marker discovery, state scoring,
# composition/Ro-e calculation, biological labels, colors, or frozen coordinates are changed.
#
# Project: adult SCI resident microglia program
# Target journal: GLIA
#
# DESIGN PRINCIPLES
# - Read frozen RDS/CSV only.
# - Never rerun QC, Harmony, UMAP, clustering, marker discovery, or state scoring.
# - A/B use the SAME frozen broad-atlas UMAP coordinates and SAME point order.
# - E uses Step2C final lineage adjudication.
# - G/H use the Step3B cleaned resident-microglia object, with the exact
#   frozen Step2C resident lineage classes carried into Step3.
#   (Canonical / activated / IFN-responsive / cycling; no new DAM label is
#   invented here. DAM-like transition is reserved for Figure 2.)
# - Biological sample is the composition display unit; no cell-level P values.
# - UMAP point layers are rasterized, text/axes remain vector.
# - Arial font, white background, explicit colors, no default ggplot palette.
#
# OUTPUT
# <PROJECT_ROOT>/09_final_figures/GLIA/Figure1/
#   Figure1A.* ... Figure1H.*
#   Figure_1_complete.*
#   00_input_audit.csv
#   01_schema_audit.csv
#   02_plot_data_used.csv
#   03_panel_layout_audit.csv
#   source_file_md5.csv
#   sessionInfo.txt
# ============================================================

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
set.seed(20260811)

# ------------------------------------------------------------
# 0. USER SETTINGS
# ------------------------------------------------------------
project_root <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
out_dir <- file.path(
  project_root,
  "09_final_figures", "GLIA", "Figure1", "v29_CFleft_Elabel_fix"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Build all panels by default. For development, e.g. c("A","B","G","H")
panels_to_build <- LETTERS[1:8]

# Final manuscript canvas harmonized to Figure 4 v8
complete_width_mm  <- 180
complete_height_mm <- 230

# ------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------
required_pkgs <- c(
  "Seurat",
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "patchwork",
  "scales",
  "ggrastr",
  "ggrepel",
  "svglite",
  "ragg"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them first, then rerun."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(scales)
  library(ggrepel)
})

# ------------------------------------------------------------
# 2. FROZEN INPUTS
# ------------------------------------------------------------
input_files <- c(
  broad_rds = file.path(
    project_root,
    "04_results/step1F_GSE162610_broad_annotation",
    "GSE162610_broad_annotation_audit.rds"
  ),
  broad_marker_availability = file.path(
    project_root,
    "04_results/step1F_GSE162610_broad_annotation/tables",
    "00_marker_gene_availability.csv"
  ),
  broad_type_table = file.path(
    project_root,
    "04_results/step1F_GSE162610_broad_annotation/tables",
    "04_cluster_broad_type_suggestions.csv"
  ),
  myeloid_rds = file.path(
    project_root,
    "04_results/step2C_GSE162610_final_microglia_extraction",
    "GSE162610_step2_full_myeloid_with_final_lineage.rds"
  ),
  lineage_map = file.path(
    project_root,
    "04_results/step2C_GSE162610_final_microglia_extraction/tables",
    "01_final_step2_lineage_map.csv"
  ),
  microglia_rds = file.path(
    project_root,
    "04_results/step3B_GSE162610_state_scoring",
    "GSE162610_step3B_clean_state_scored.rds"
  ),
  cleanup_map = file.path(
    project_root,
    "04_results/step3B_GSE162610_state_scoring/tables",
    "01_step3B_residual_cleanup_map.csv"
  )
)

input_audit <- data.frame(
  key = names(input_files),
  path = unname(input_files),
  exists = file.exists(input_files),
  size_MB = ifelse(
    file.exists(input_files),
    round(file.info(input_files)$size / 1024^2, 3),
    NA_real_
  ),
  stringsAsFactors = FALSE
)

write.csv(
  input_audit,
  file.path(out_dir, "00_input_audit.csv"),
  row.names = FALSE
)

if (!all(input_audit$exists)) {
  print(input_audit)
  stop(
    "At least one required Figure 1 input is missing. ",
    "See 00_input_audit.csv."
  )
}

# MD5 audit. This reads the frozen files but does not modify them.
message("Computing source-file MD5...")
source_md5 <- data.frame(
  key = names(input_files),
  path = unname(input_files),
  md5 = unname(tools::md5sum(input_files)),
  stringsAsFactors = FALSE
)

write.csv(
  source_md5,
  file.path(out_dir, "source_file_md5.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 3. GLOBAL STYLE + FROZEN COLORS
# ------------------------------------------------------------

# Frozen timepoint colors from the plotting plan
time_cols <- c(
  "Uninjured" = "#4C78A8",  # A: Microglia/Myeloid blue
  "1 dpi"     = "#EDC948",  # A: Neuron yellow
  "3 dpi"     = "#F28E2B",  # A: T-cell orange
  "7 dpi"     = "#E45756"   # A: Neutrophil red
)

time_levels <- names(time_cols)

# Explicit broad-lineage palette
broad_cols <- c(
  "Microglia/Myeloid" = "#4C78A8",
  "Neutrophil"        = "#E45756",
  "Astrocyte"         = "#E7298A",
  "Oligodendrocyte"   = "#59A14F",
  "OPC"               = "#76B7B2",
  "Endothelial"       = "#1B9E77",
  "Fibroblast"        = "#7570B3",
  "Pericyte/VSMC"     = "#66A61E",
  "Neuron"            = "#EDC948",
  "T cell"            = "#F28E2B",
  "Ependymal"         = "#B279A2",
  "B cell"            = "#8CD17D",
  "Mast"              = "#FF9DA7",
  "Other"             = "#9D9D9D"
)

# Step2C final lineage display colors
lineage_cols <- c(
  "Resident-derived microglia"   = "#4C78A8",  # A blue
  "Infiltrating mono/macrophage" = "#E45756",  # A red
  "BAM-like macrophage"          = "#7570B3"   # A purple
)

# Figure-1 resident state classes: all are direct display collapses of
# frozen Step2C final lineage calls, not newly inferred state labels.
resident_state_cols <- c(
  "Canonical resident"      = "#4C78A8",  # A blue
  "Activated resident"      = "#E45756",  # A red
  "IFN-responsive resident" = "#1B9E77",  # A teal/green
  "Cycling resident"        = "#B279A2"   # A purple
)

dot_low  <- "#2166AC"
dot_mid  <- "#F7F7F7"
dot_high <- "#B2182B"

base_family <- "Arial"

theme_pub <- theme_classic(
  base_family = base_family,
  base_size = 8.0
) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(size = 8.0, color = "black"),
    axis.text = element_text(size = 7.2, color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 6.4),
    legend.key.height = grid::unit(3.2, "mm"),
    legend.key.width = grid::unit(3.2, "mm"),
    plot.margin = margin(4, 4, 5, 5, unit = "mm"),
    plot.tag = element_text(
      family = base_family,
      face = "bold",
      size = 16,
      color = "black"
    ),
    plot.tag.position = c(-0.04, 1.04)
  )

theme_umap <- theme_pub +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )

# ------------------------------------------------------------
# 4. HELPERS
# ------------------------------------------------------------

check_meta <- function(obj, required, object_name) {
  miss <- setdiff(required, colnames(obj[[]]))
  if (length(miss) > 0) {
    stop(
      object_name, " is missing metadata: ",
      paste(miss, collapse = ", ")
    )
  }
}

check_reduction <- function(obj, reduction, object_name) {
  if (!reduction %in% Reductions(obj)) {
    stop(
      object_name, " is missing frozen reduction '", reduction,
      "'. Available reductions: ",
      paste(Reductions(obj), collapse = ", ")
    )
  }
}

schema_row <- function(obj, object_name) {
  data.frame(
    object = object_name,
    class = paste(class(obj), collapse = ";"),
    n_cells = ncol(obj),
    n_features = nrow(obj),
    assays = paste(Assays(obj), collapse = ";"),
    reductions = paste(Reductions(obj), collapse = ";"),
    metadata_columns = paste(colnames(obj[[]]), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

extract_umap_df <- function(obj, reduction, extra_cols) {
  emb <- Embeddings(obj, reduction = reduction)
  if (ncol(emb) < 2) stop("Reduction ", reduction, " has <2 dimensions.")

  df <- data.frame(
    cell_id = rownames(emb),
    UMAP_1 = emb[, 1],
    UMAP_2 = emb[, 2],
    stringsAsFactors = FALSE
  )

  meta <- obj[[]]
  meta$cell_id <- rownames(meta)

  extra_cols <- unique(extra_cols)
  extra_cols <- extra_cols[extra_cols %in% colnames(meta)]

  if (length(extra_cols) > 0) {
    df <- left_join(
      df,
      meta[, c("cell_id", extra_cols), drop = FALSE],
      by = "cell_id"
    )
  }

  df
}

# Same point order is reused for A and B.
freeze_point_order <- function(df, seed) {
  set.seed(seed)
  df[sample(seq_len(nrow(df))), , drop = FALSE]
}

collapse_broad_type <- function(x) {
  x <- trimws(as.character(x))

  out <- dplyr::case_when(
    x %in% c("Microglia", "Myeloid_general", "Macrophage_Monocyte") ~ "Microglia/Myeloid",
    x == "Neutrophil"      ~ "Neutrophil",
    x == "Astrocyte"       ~ "Astrocyte",
    x == "Oligodendrocyte" ~ "Oligodendrocyte",
    x == "OPC"             ~ "OPC",
    x == "Endothelial"     ~ "Endothelial",
    x == "Fibroblast"      ~ "Fibroblast",
    x == "Pericyte_VSMC"   ~ "Pericyte/VSMC",
    x == "Neuron"          ~ "Neuron",
    x == "T_cell"          ~ "T cell",
    x == "B_cell"          ~ "B cell",
    x == "Mast"            ~ "Mast",
    x == "Ependymal"       ~ "Ependymal",
    TRUE                   ~ NA_character_
  )

  out
}
collapse_final_lineage <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x %in% c(
      "canonical_resident_microglia",
      "activated_resident_microglia",
      "cycling_S_microglia",
      "cycling_G2M_microglia",
      "IFN_responsive_resident_microglia"
    ) ~ "Resident-derived microglia",
    x == "BAM_like_macrophage" ~ "BAM-like macrophage",
    TRUE ~ "Infiltrating mono/macrophage"
  )
}

collapse_resident_state <- function(step2_cluster) {
  x <- as.character(step2_cluster)
  dplyr::case_when(
    x == "0"              ~ "Canonical resident",
    x %in% c("3", "4")    ~ "Activated resident",
    x == "9"              ~ "IFN-responsive resident",
    x %in% c("5", "7")    ~ "Cycling resident",
    TRUE                  ~ NA_character_
  )
}

# Explicit vector UMAP axes in lower-left. These remain visible in SVG/PDF.
add_umap_axes <- function(
  p,
  df,
  x_fraction = 0.132,
  y_fraction = 0.132
) {
  xr <- range(df$UMAP_1, na.rm = TRUE)
  yr <- range(df$UMAP_2, na.rm = TRUE)
  dx <- diff(xr)
  dy <- diff(yr)

  # Larger dedicated lower-left pocket plus longer arrows.
  x0 <- xr[1] - 0.135 * dx
  y0 <- yr[1] - 0.155 * dy
  x1 <- x0 + x_fraction * dx
  y1 <- y0 + y_fraction * dy

  x_limits <- c(xr[1] - 0.215 * dx, xr[2] + 0.020 * dx)
  y_limits <- c(yr[1] - 0.255 * dy, yr[2] + 0.045 * dy)

  p +
    scale_x_continuous(
      limits = x_limits,
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      limits = y_limits,
      expand = expansion(mult = 0)
    ) +
    coord_equal(clip = "off") +
    annotate(
      "segment",
      x = x0, xend = x1,
      y = y0, yend = y0,
      linewidth = 0.62,
      color = "black",
      lineend = "round",
      arrow = grid::arrow(
        type = "closed",
        length = grid::unit(1.6, "mm")
      )
    ) +
    annotate(
      "segment",
      x = x0, xend = x0,
      y = y0, yend = y1,
      linewidth = 0.62,
      color = "black",
      lineend = "round",
      arrow = grid::arrow(
        type = "closed",
        length = grid::unit(1.6, "mm")
      )
    ) +
    annotate(
      "text",
      x = (x0 + x1) / 2,
      y = y0 - 0.072 * dy,
      label = "UMAP_1",
      family = base_family,
      size = 2.05,
      hjust = 0.5,
      vjust = 1
    ) +
    annotate(
      "text",
      x = x0 - 0.078 * dx,
      y = (y0 + y1) / 2,
      label = "UMAP_2",
      family = base_family,
      size = 2.05,
      angle = 90,
      hjust = 0.5,
      vjust = 0.5
    )
}

label_centers <- function(df, group_col, subgroup_col = NULL) {
  tmp <- df %>%
    filter(!is.na(.data[[group_col]]))

  # When one displayed class contains several frozen clusters, place its label
  # inside the largest contributing frozen cluster rather than in an empty
  # geometric midpoint between islands.
  if (!is.null(subgroup_col) && subgroup_col %in% colnames(tmp)) {
    largest_subgroup <- tmp %>%
      filter(!is.na(.data[[subgroup_col]])) %>%
      count(
        .data[[group_col]],
        .data[[subgroup_col]],
        name = "n_cells"
      ) %>%
      group_by(.data[[group_col]]) %>%
      slice_max(
        order_by = n_cells,
        n = 1,
        with_ties = FALSE
      ) %>%
      ungroup() %>%
      select(all_of(c(group_col, subgroup_col)))

    tmp <- inner_join(
      tmp,
      largest_subgroup,
      by = c(group_col, subgroup_col)
    )
  }

  tmp %>%
    group_by(.data[[group_col]]) %>%
    summarise(
      label = first(as.character(.data[[group_col]])),
      UMAP_1 = median(UMAP_1, na.rm = TRUE),
      UMAP_2 = median(UMAP_2, na.rm = TRUE),
      .groups = "drop"
    )
}

make_umap <- function(
  df,
  color_col,
  color_values,
  point_size = 0.22,
  label_groups = TRUE,
  label_subgroup_col = NULL,
  label_map = NULL,
  label_size = 2.05,
  label_force = 0.55,
  label_nudge_x_frac = 0,
  label_nudge_y_frac = 0,
  label_extra_nudge_df = NULL,
  show_legend = TRUE,
  total_n_label = FALSE,
  seed_tag = "UMAP"
) {
  # Preserve only colors that are actually present.
  present <- unique(as.character(df[[color_col]]))
  present <- present[!is.na(present)]
  missing_colors <- setdiff(present, names(color_values))
  if (length(missing_colors) > 0) {
    stop(
      seed_tag, ": explicit colors are missing for: ",
      paste(missing_colors, collapse = ", ")
    )
  }

  p <- ggplot(
    df,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      color = .data[[color_col]]
    )
  ) +
    ggrastr::geom_point_rast(
      size = point_size,
      alpha = 0.86,
      raster.dpi = 400,
      stroke = 0
    ) +
    scale_color_manual(
      values = color_values,
      drop = FALSE
    ) +
    theme_umap +
    guides(
      color = guide_legend(
        override.aes = list(size = 2.2, alpha = 1),
        title.position = "top"
      )
    ) +
    theme(
      legend.position = if (show_legend) "right" else "none",
      legend.key.height = grid::unit(2.1, "mm"),
      legend.key.width = grid::unit(2.1, "mm"),
      legend.text = element_text(size = 6.2, colour = "black"),
      legend.spacing.y = grid::unit(0.25, "mm")
    )

  if (label_groups) {
    centers <- label_centers(
      df,
      group_col = color_col,
      subgroup_col = label_subgroup_col
    )

    if (!is.null(label_map)) {
      hit <- centers$label %in% names(label_map)
      centers$label[hit] <- unname(label_map[centers$label[hit]])
    }

    # Optional deterministic nudge for panel-specific label fine-tuning.
    xr_local <- range(df$UMAP_1, na.rm = TRUE)
    yr_local <- range(df$UMAP_2, na.rm = TRUE)

    if (!isTRUE(all.equal(label_nudge_x_frac, 0)) || !isTRUE(all.equal(label_nudge_y_frac, 0))) {
      centers$UMAP_1 <- centers$UMAP_1 + label_nudge_x_frac * diff(xr_local)
      centers$UMAP_2 <- centers$UMAP_2 + label_nudge_y_frac * diff(yr_local)
    }

    # Optional label-specific extra nudges after the global panel nudge.
    if (!is.null(label_extra_nudge_df) && nrow(label_extra_nudge_df) > 0) {
      required_nudge_cols <- c("label", "dx_frac", "dy_frac")
      if (!all(required_nudge_cols %in% colnames(label_extra_nudge_df))) {
        stop("label_extra_nudge_df must contain columns: label, dx_frac, dy_frac")
      }
      for (i in seq_len(nrow(label_extra_nudge_df))) {
        hit_i <- centers$label == as.character(label_extra_nudge_df$label[i])
        if (any(hit_i)) {
          centers$UMAP_1[hit_i] <- centers$UMAP_1[hit_i] +
            as.numeric(label_extra_nudge_df$dx_frac[i]) * diff(xr_local)
          centers$UMAP_2[hit_i] <- centers$UMAP_2[hit_i] +
            as.numeric(label_extra_nudge_df$dy_frac[i]) * diff(yr_local)
        }
      }
    }

    p <- p +
      ggrepel::geom_text_repel(
        data = centers,
        aes(x = UMAP_1, y = UMAP_2, label = label),
        inherit.aes = FALSE,
        family = base_family,
        size = label_size,
        color = "black",
        fontface = "plain",
        seed = 20260811,
        max.overlaps = Inf,
        box.padding = 0.20,
        point.padding = 0.05,
        min.segment.length = Inf,
        segment.color = NA,
        force = label_force,
        direction = "both"
      )
  }

  if (total_n_label) {
    xr <- range(df$UMAP_1, na.rm = TRUE)
    yr <- range(df$UMAP_2, na.rm = TRUE)

    p <- p +
      annotate(
        "text",
        x = xr[2],
        y = yr[2] + 0.065 * diff(yr),
        label = paste0("n = ", format(nrow(df), big.mark = ",")),
        family = base_family,
        size = 2.4,
        hjust = 1,
        vjust = 0.5
      )
  }

  add_umap_axes(p, df)
}

make_dotplot <- function(
  obj,
  group_col,
  marker_df,
  tag_letter,
  dot_max = 3.3,
  show_row_group_labels = TRUE,
  gene_text_size = 4.6,
  x_text_size = 4.9
) {
  check_meta(obj, group_col, paste0("DotPlot-", tag_letter))

  marker_df <- marker_df %>%
    filter(gene %in% rownames(obj)) %>%
    distinct(marker_group, gene)

  if (nrow(marker_df) == 0) {
    stop("Panel ", tag_letter, ": no requested marker genes found.")
  }

  # Freeze marker and marker-block order exactly as supplied.
  marker_group_levels <- unique(as.character(marker_df$marker_group))
  marker_df$marker_group <- factor(
    marker_df$marker_group,
    levels = marker_group_levels
  )
  marker_df$gene <- factor(
    marker_df$gene,
    levels = rev(unique(as.character(marker_df$gene)))
  )

  dp <- Seurat::DotPlot(
    object = obj,
    features = as.character(unique(marker_df$gene)),
    group.by = group_col,
    assay = "RNA",
    scale = TRUE
  )$data

  dp <- dp %>%
    mutate(
      gene = as.character(features.plot)
    ) %>%
    left_join(
      marker_df %>%
        mutate(gene = as.character(gene)),
      by = "gene"
    )

  dp$gene <- factor(
    dp$gene,
    levels = rev(as.character(unique(marker_df$gene)))
  )
  dp$marker_group <- factor(
    as.character(dp$marker_group),
    levels = marker_group_levels
  )

  # Preserve the frozen group factor order from Seurat metadata.
  group_vec <- obj[[]][[group_col]]
  if (is.factor(group_vec)) {
    group_levels <- levels(group_vec)
  } else {
    group_levels <- unique(as.character(group_vec))
  }
  group_levels <- group_levels[
    group_levels %in% unique(as.character(dp$id))
  ]
  dp$id <- factor(as.character(dp$id), levels = group_levels)

  ggplot(
    dp,
    aes(
      x = id,
      y = gene,
      size = pct.exp,
      color = avg.exp.scaled
    )
  ) +
    geom_point(alpha = 0.95) +
    scale_size_continuous(
      range = c(0.2, dot_max),
      limits = c(0, 100),
      breaks = c(25, 50, 75, 100),
      name = "Expressing\ncells (%)"
    ) +
    scale_color_gradient2(
      low = dot_low,
      mid = dot_mid,
      high = dot_high,
      midpoint = 0,
      name = "Mean\nexpression"
    ) +
    facet_grid(
      marker_group ~ .,
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ) +
    labs(x = NULL, y = NULL) +
    theme_pub +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = x_text_size
      ),
      axis.text.y = element_text(size = gene_text_size),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = if (show_row_group_labels) element_text(
        angle = 0,
        face = "bold",
        size = 5.4,
        hjust = 1
      ) else element_blank(),
      panel.spacing.y = grid::unit(0.42, "mm"),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.justification = c(0, 1),
      legend.title = element_text(
        size = 6.2,
        colour = "black",
        face = "plain"
      ),
      legend.text = element_text(
        size = 5.8,
        colour = "black"
      ),
      legend.key.height = grid::unit(2.3, "mm"),
      legend.key.width = grid::unit(4.8, "mm"),
      legend.spacing.x = grid::unit(0.8, "mm"),
      plot.tag.position = c(-0.10, 1.03),
      plot.margin = margin(4, 3, 4, 10, unit = "mm")
    ) +
    guides(
      size = guide_legend(
        order = 1,
        title.position = "top",
        nrow = 1,
        byrow = TRUE,
        override.aes = list(color = "black", alpha = 1)
      ),
      color = guide_colorbar(
        order = 2,
        title.position = "top",
        direction = "horizontal",
        barwidth = grid::unit(18, "mm"),
        barheight = grid::unit(2.5, "mm"),
        ticks = FALSE
      )
    ) +
    labs(tag = tag_letter)
}


make_dotplot_reference <- function(
  obj,
  group_col,
  marker_df,
  tag_letter,
  group_labels = NULL,
  dot_max = 3.2,
  x_text_size = 6.2,
  y_text_size = 6.5
) {
  check_meta(obj, group_col, paste0("ReferenceDotPlot-", tag_letter))

  marker_df <- marker_df %>%
    filter(gene %in% rownames(obj)) %>%
    distinct(marker_group, gene)

  if (nrow(marker_df) == 0) {
    stop("Panel ", tag_letter, ": no requested marker genes found.")
  }

  # Preserve the exact supplied marker order from left to right.
  gene_order <- unique(as.character(marker_df$gene))

  # Preserve the frozen group order from object metadata.
  group_vec <- obj[[]][[group_col]]
  if (is.factor(group_vec)) {
    group_levels <- levels(group_vec)
  } else {
    group_levels <- unique(as.character(group_vec))
  }

  dp <- Seurat::DotPlot(
    object = obj,
    features = gene_order,
    group.by = group_col,
    assay = "RNA",
    scale = TRUE
  )$data %>%
    mutate(
      gene = as.character(features.plot),
      group = as.character(id)
    )

  group_levels <- group_levels[
    group_levels %in% unique(dp$group)
  ]

  # Reference-style orientation:
  # first biological group is at the bottom, later groups above it.
  dp$gene <- factor(
    dp$gene,
    levels = gene_order
  )
  dp$group <- factor(
    dp$group,
    levels = group_levels
  )

  # Default y-axis labels use the frozen display names.
  y_labels <- setNames(group_levels, group_levels)
  if (!is.null(group_labels)) {
    common_names <- intersect(names(group_labels), names(y_labels))
    y_labels[common_names] <- group_labels[common_names]
  }

  ggplot(
    dp,
    aes(
      x = gene,
      y = group,
      size = pct.exp,
      color = avg.exp.scaled
    )
  ) +
    geom_point(
      alpha = 0.98,
      stroke = 0.12
    ) +
    scale_size_continuous(
      range = c(0.35, dot_max),
      limits = c(0, 100),
      breaks = c(25, 50, 75),
      name = "Percent Expressed"
    ) +
    scale_color_gradientn(
      colours = c(
        "#4575B4",  # blue: low
        "#FFFFBF",  # yellow: intermediate
        "#D73027"   # red: high
      ),
      values = scales::rescale(c(-1.5, 0, 1.5)),
      limits = c(-1.5, 1.5),
      breaks = c(-1, 0, 1),
      oob = scales::squish,
      name = "Average Expression"
    ) +
    scale_x_discrete(
      limits = gene_order,
      drop = FALSE,
      expand = expansion(add = c(0.72, 0.72))
    ) +
    scale_y_discrete(
      limits = group_levels,
      labels = y_labels,
      drop = FALSE,
      expand = expansion(add = c(0.62, 0.62))
    ) +
    labs(
      x = NULL,
      y = NULL,
      tag = tag_letter
    ) +
    theme_bw(
      base_family = base_family,
      base_size = 7.2
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        color = "#6A6A6A",
        linewidth = 0.42,
        fill = NA
      ),
      axis.line = element_blank(),
      axis.ticks = element_line(
        color = "#5F5F5F",
        linewidth = 0.28
      ),
      axis.text.x = element_text(
        family = base_family,
        size = x_text_size,
        color = "black",
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = element_text(
        family = base_family,
        size = y_text_size,
        color = "black",
        hjust = 1
      ),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.box.just = "left",
      legend.justification = c(0, 0.5),
      legend.margin = margin(0, 0, 1.0, 0, unit = "mm"),
      legend.box.margin = margin(0, 0, 0, 0, unit = "mm"),
      legend.title = element_text(
        family = base_family,
        size = 5.7,
        face = "plain",
        color = "black"
      ),
      legend.text = element_text(
        family = base_family,
        size = 5.5,
        color = "black"
      ),
      legend.key.height = grid::unit(2.5, "mm"),
      legend.key.width = grid::unit(3.0, "mm"),
      legend.spacing.x = grid::unit(0.80, "mm"),
      plot.tag = element_text(
        family = base_family,
        face = "bold",
        size = 16,
        color = "black"
      ),
      plot.tag.position = c(-0.04, 1.04),
      plot.margin = margin(4, 4, 5, 5, unit = "mm")
    ) +
    guides(
      size = guide_legend(
        order = 1,
        title.position = "left",
        nrow = 1,
        byrow = TRUE,
        keywidth = grid::unit(2.7, "mm"),
        keyheight = grid::unit(2.2, "mm"),
        override.aes = list(
          color = "black",
          alpha = 1
        )
      ),
      color = guide_colorbar(
        order = 2,
        title.position = "left",
        direction = "horizontal",
        barwidth = grid::unit(9.5, "mm"),
        barheight = grid::unit(2.2, "mm"),
        ticks = FALSE
      )
    )
}

add_tag <- function(p, letter) {
  p + labs(tag = letter) + theme_pub
}

save_plot_triplet <- function(
  p,
  stem,
  width_mm,
  height_mm,
  dpi = 600
) {
  pdf_file  <- file.path(out_dir, paste0(stem, ".pdf"))
  svg_file  <- file.path(out_dir, paste0(stem, ".svg"))
  tiff_file <- file.path(out_dir, paste0(stem, ".tiff"))

  ggsave(
    pdf_file,
    p,
    device = cairo_pdf,
    width = width_mm,
    height = height_mm,
    units = "mm",
    bg = "white"
  )

  ggsave(
    svg_file,
    p,
    device = svglite::svglite,
    width = width_mm,
    height = height_mm,
    units = "mm",
    bg = "white"
  )

  ggsave(
    tiff_file,
    p,
    device = ragg::agg_tiff,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = dpi,
    compression = "lzw",
    background = "white"
  )
}

# ------------------------------------------------------------
# 5. AUDIT CONTAINERS
# ------------------------------------------------------------
schema_audit <- list()
plot_data_audit <- list()
panel_plots <- list()

record_plot_data <- function(panel, source, n_rows, fields, note = "") {
  data.frame(
    panel = panel,
    source = source,
    n_rows = n_rows,
    fields = paste(fields, collapse = ";"),
    note = note,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# 6. PANELS A-D | BROAD ATLAS
# ------------------------------------------------------------
if (any(c("A", "B", "C", "D") %in% panels_to_build)) {

  message("Reading frozen broad atlas...")
  broad <- readRDS(input_files[["broad_rds"]])

  if (!inherits(broad, "Seurat")) {
    stop("Broad atlas RDS is not a Seurat object.")
  }

  DefaultAssay(broad) <- "RNA"

  check_meta(
    broad,
    c(
      "cluster_provisional",
      "sample_id",
      "timepoint",
      "broad_type_suggestion"
    ),
    "Broad atlas"
  )
  check_reduction(broad, "umap.harmony", "Broad atlas")

  schema_audit[["broad"]] <- schema_row(broad, "broad_atlas")

  # IMPORTANT v02:
  # Use the frozen Step1F cluster-to-broad-type table directly instead of
  # relying on a free-text metadata collapse. This fixes the v01 A/C/D bug
  # in which unmatched labels were displayed as a single "Other" class.
  broad_type_tbl <- read.csv(
    input_files[["broad_type_table"]],
    check.names = FALSE
  )

  required_broad_cols <- c("cluster", "broad_type_suggestion")
  if (!all(required_broad_cols %in% colnames(broad_type_tbl))) {
    stop(
      "Frozen broad-type table is missing: ",
      paste(
        setdiff(required_broad_cols, colnames(broad_type_tbl)),
        collapse = ", "
      )
    )
  }

  # IMPORTANT Step1F provenance repair:
  # The frozen table stores broad-cluster IDs as g0...g18, whereas the
  # Seurat metadata stores the same frozen clusters as 0...18.
  # Normalize ONLY this textual prefix before joining. This does not
  # recluster, relabel, or infer any biological annotation.
  normalize_step1F_cluster_id <- function(x) {
    sub(
      pattern = "^[gG]",
      replacement = "",
      x = trimws(as.character(x))
    )
  }

  broad_type_tbl$cluster_key_normalized <- normalize_step1F_cluster_id(
    broad_type_tbl$cluster
  )

  broad_cluster_to_type <- setNames(
    as.character(broad_type_tbl$broad_type_suggestion),
    broad_type_tbl$cluster_key_normalized
  )

  broad_cluster_key <- normalize_step1F_cluster_id(
    broad$cluster_provisional
  )

  broad_type_raw_frozen <- unname(
    broad_cluster_to_type[broad_cluster_key]
  )

  # Write an explicit mapping audit before any plot is drawn.
  broad_mapping_audit <- data.frame(
    cluster_provisional = as.character(broad$cluster_provisional),
    cluster_key_normalized = broad_cluster_key,
    frozen_broad_type = broad_type_raw_frozen,
    stringsAsFactors = FALSE
  ) %>%
    distinct() %>%
    arrange(as.numeric(cluster_key_normalized))

  write.csv(
    broad_mapping_audit,
    file.path(out_dir, "00B_broad_cluster_mapping_audit.csv"),
    row.names = FALSE
  )

  if (anyNA(broad_type_raw_frozen)) {
    bad_clusters <- sort(
      unique(
        as.character(
          broad$cluster_provisional[
            is.na(broad_type_raw_frozen)
          ]
        )
      )
    )
    stop(
      "Some normalized broad clusters were not found in the frozen Step1F map: ",
      paste(bad_clusters, collapse = ", "),
      ". See 00B_broad_cluster_mapping_audit.csv. No broad label was guessed."
    )
  }

  # All 19 frozen Step1F broad clusters must map exactly once.
  expected_clusters <- sort(unique(normalize_step1F_cluster_id(broad_type_tbl$cluster)))
  observed_clusters <- sort(unique(broad_cluster_key))

  if (!identical(expected_clusters, observed_clusters)) {
    stop(
      "Frozen Step1F cluster sets differ after normalization. ",
      "Expected: ", paste(expected_clusters, collapse = ", "),
      " | observed: ", paste(observed_clusters, collapse = ", "),
      ". No plot was generated."
    )
  }

  message(
    "Step1F broad-cluster mapping passed after ID normalization (g0..g18 -> 0..18)."
  )

  broad$major_cell_type_display <- collapse_broad_type(
    broad_type_raw_frozen
  )

  if (anyNA(broad$major_cell_type_display)) {
    bad_types <- sort(
      unique(
        broad_type_raw_frozen[
          is.na(broad$major_cell_type_display)
        ]
      )
    )
    stop(
      "Unrecognized frozen broad-type suggestion(s): ",
      paste(bad_types, collapse = ", "),
      ". Update the explicit display map rather than silently using Other."
    )
  }

  present_broad_types <- names(broad_cols)[
    names(broad_cols) %in% unique(broad$major_cell_type_display)
  ]

  broad$major_cell_type_display <- factor(
    broad$major_cell_type_display,
    levels = present_broad_types
  )

  broad_cols_present <- broad_cols[present_broad_types]

  broad$timepoint <- factor(
    as.character(broad$timepoint),
    levels = time_levels
  )

  if (any(is.na(broad$timepoint))) {
    bad <- unique(as.character(broad[[]]$timepoint[is.na(broad$timepoint)]))
    stop(
      "Unexpected broad-atlas timepoint value(s): ",
      paste(bad, collapse = ", ")
    )
  }

  broad_umap <- extract_umap_df(
    broad,
    reduction = "umap.harmony",
    extra_cols = c(
      "sample_id",
      "timepoint",
      "cluster_provisional",
      "major_cell_type_display"
    )
  )

  broad_umap <- freeze_point_order(
    broad_umap,
    seed = 20260811
  )

  # A
  if ("A" %in% panels_to_build) {
    pA <- make_umap(
      broad_umap,
      color_col = "major_cell_type_display",
      color_values = broad_cols_present,
      point_size = 0.18,
      label_groups = TRUE,
      label_subgroup_col = "cluster_provisional",
      label_map = c(
        "Microglia/Myeloid" = "Microglia/\nMyeloid",
        "Pericyte/VSMC" = "Pericyte/\nVSMC"
      ),
      label_size = 2.18,
      label_force = 0.72,
      show_legend = TRUE,
      total_n_label = TRUE,
      seed_tag = "Panel A"
    ) +
      labs(
        color = "Major cell type",
        tag = "A"
      ) +
      theme_pub +
      theme(
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "right"
      )

    panel_plots[["A"]] <- pA

    plot_data_audit[["A"]] <- record_plot_data(
      "A",
      input_files[["broad_rds"]],
      nrow(broad_umap),
      c("umap.harmony", "broad_type_suggestion"),
      "Frozen Step1F cluster-to-broad-type table mapped to display labels; no UMAP recalculation."
    )
  }

  # B
  if ("B" %in% panels_to_build) {
    pB <- make_umap(
      broad_umap,
      color_col = "timepoint",
      color_values = time_cols,
      point_size = 0.18,
      label_groups = FALSE,
      show_legend = TRUE,
      total_n_label = FALSE,
      seed_tag = "Panel B"
    ) +
      labs(
        color = NULL,
        tag = "B"
      ) +
      theme_pub +
      theme(
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "right"
      )

    panel_plots[["B"]] <- pB

    plot_data_audit[["B"]] <- record_plot_data(
      "B",
      input_files[["broad_rds"]],
      nrow(broad_umap),
      c("umap.harmony", "timepoint"),
      "Exactly the same embedding and point order as panel A."
    )
  }

  # C | compact frozen canonical marker panel
  if ("C" %in% panels_to_build) {

    marker_C <- tribble(
      ~marker_group,          ~gene,
      "Microglia/Myeloid",    "P2ry12",
      "Microglia/Myeloid",    "Tmem119",
      "Microglia/Myeloid",    "Lyz2",
      "Neutrophil",           "Ly6g",
      "Neutrophil",           "Retnlg",
      "Astrocyte",            "Aqp4",
      "Astrocyte",            "Gfap",
      "Oligodendrocyte",      "Mbp",
      "Oligodendrocyte",      "Plp1",
      "OPC",                  "Pdgfra",
      "OPC",                  "Cspg4",
      "Endothelial",          "Pecam1",
      "Endothelial",          "Cdh5",
      "Fibroblast",           "Col1a1",
      "Fibroblast",           "Dcn",
      "Pericyte/VSMC",        "Pdgfrb",
      "Pericyte/VSMC",        "Rgs5",
      "Neuron",               "Snap25",
      "Neuron",               "Rbfox3",
      "T cell",               "Cd3d",
      "T cell",               "Cd3e",
      "Ependymal",            "Foxj1",
      "Ependymal",            "Pifo"
    )

    pC <- make_dotplot_reference(
      broad,
      group_col = "major_cell_type_display",
      marker_df = marker_C,
      tag_letter = "C",
      group_labels = c(
        "Microglia/Myeloid" = "Microglia/Myeloid",
        "Neutrophil" = "Neutrophil",
        "Astrocyte" = "Astrocyte",
        "Oligodendrocyte" = "Oligodendrocyte",
        "OPC" = "OPC",
        "Endothelial" = "Endothelial",
        "Fibroblast" = "Fibroblast",
        "Pericyte/VSMC" = "Pericyte/VSMC",
        "Neuron" = "Neuron",
        "T cell" = "T cell",
        "Ependymal" = "Ependymal"
      ),
      dot_max = 3.05,
      x_text_size = 5.2,
      y_text_size = 6.4
    ) +
      theme(
        legend.box.margin = margin(0, 0, 0, -8.5, unit = "mm")
      )

    panel_plots[["C"]] <- pC

    plot_data_audit[["C"]] <- record_plot_data(
      "C",
      input_files[["broad_rds"]],
      length(unique(marker_C$gene)),
      c("major_cell_type_display", "RNA normalized data"),
      "Canonical genes are a fixed subset of the frozen Step1F marker panel; groups use the frozen Step1F cluster-to-broad-type map; no marker discovery."
    )
  }

  # D | sample-level major-cell-type composition
  if ("D" %in% panels_to_build) {

    # Remove display categories absent from the broad atlas.
    present_types <- present_broad_types

    comp_D <- broad[[]] %>%
      transmute(
        sample_id = as.character(sample_id),
        timepoint = factor(as.character(timepoint), levels = time_levels),
        major_cell_type_display = as.character(major_cell_type_display)
      ) %>%
      filter(
        !is.na(sample_id),
        !is.na(timepoint),
        !is.na(major_cell_type_display)
      ) %>%
      count(
        sample_id,
        timepoint,
        major_cell_type_display,
        name = "n_cells"
      ) %>%
      group_by(sample_id, timepoint) %>%
      tidyr::complete(
        major_cell_type_display = present_types,
        fill = list(n_cells = 0)
      ) %>%
      mutate(
        sample_total = sum(n_cells),
        proportion = ifelse(
          sample_total > 0,
          n_cells / sample_total,
          NA_real_
        )
      ) %>%
      ungroup() %>%
      filter(!is.na(proportion)) %>%
      mutate(
        major_cell_type_display = factor(
          major_cell_type_display,
          levels = present_types
        )
      )

    # Sample-level Ro/e calculation.
    # Expected proportion is the overall abundance of each major cell type
    # across the complete frozen broad atlas. For each biological sample:
    # Ro/e = observed within-sample proportion / expected overall proportion.
    # Ro/e > 1 = enrichment; Ro/e < 1 = depletion.

    expected_D <- broad[[]] %>%
      transmute(
        major_cell_type_display = as.character(major_cell_type_display)
      ) %>%
      filter(!is.na(major_cell_type_display)) %>%
      count(
        major_cell_type_display,
        name = "overall_cells"
      ) %>%
      mutate(
        expected_proportion = overall_cells / sum(overall_cells)
      )

    roe_D <- comp_D %>%
      left_join(
        expected_D,
        by = "major_cell_type_display"
      ) %>%
      mutate(
        roe = ifelse(
          expected_proportion > 0,
          proportion / expected_proportion,
          NA_real_
        )
      )

    summary_D <- roe_D %>%
      group_by(timepoint, major_cell_type_display) %>%
      summarise(
        median_roe = median(roe, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        timepoint = factor(
          as.character(timepoint),
          levels = rev(time_levels)
        ),
        major_cell_type_display = factor(
          as.character(major_cell_type_display),
          levels = present_types
        )
      )

    finite_roe_D <- summary_D$median_roe[
      is.finite(summary_D$median_roe)
    ]

    if (length(finite_roe_D) == 0) {
      stop("Panel D: no finite Ro/e values were available.")
    }

    roe_upper_D <- max(
      2.5,
      as.numeric(
        quantile(
          finite_roe_D,
          probs = 0.97,
          na.rm = TRUE
        )
      )
    )

    # Panel D redesign:
    # - Each pale point = one biological sample
    # - Large filled point = timepoint median
    # - Short vertical bar = interquartile range (IQR)
    # - NO lines connect different cell types because cell types are categorical,
    #   not an ordered biological trajectory.
    # - sqrt y-scale compresses the very high Microglia/Myeloid proportion while
    #   expanding the low-abundance cell types; tick labels remain true percentages.

    # Panel D | Reference-style Ro/e enrichment profile
    # The user-requested visual form follows the supplied reference:
    #   - x axis: major cell types
    #   - y axis: median sample-level Ro/e
    #   - one colored profile per SCI timepoint
    #   - circular points connected by dashed lines
    #   - horizontal reference line at Ro/e = 1
    #
    # Ro/e itself remains sample-level:
    # observed within-sample proportion / atlas-wide expected proportion.

    summary_D <- summary_D %>%
      mutate(
        timepoint = factor(
          as.character(timepoint),
          levels = time_levels
        ),
        major_cell_type_display = factor(
          as.character(major_cell_type_display),
          levels = present_types
        )
      )

    # Dynamic y range with room above the largest median Ro/e.
    finite_roe_D <- summary_D$median_roe[
      is.finite(summary_D$median_roe)
    ]

    y_max_D <- max(
      1.25,
      as.numeric(max(finite_roe_D, na.rm = TRUE)) * 1.10
    )

    # Use a modest lower bound of zero; reference line at 1 remains explicit.
    pD <- ggplot(
      summary_D,
      aes(
        x = major_cell_type_display,
        y = median_roe,
        color = timepoint,
        group = timepoint
      )
    ) +
      geom_hline(
        yintercept = 1,
        color = "#8A8A8A",
        linewidth = 0.34,
        linetype = "solid"
      ) +
      geom_line(
        linewidth = 0.62,
        linetype = "dashed",
        alpha = 0.95
      ) +
      geom_point(
        size = 2.15,
        shape = 16,
        alpha = 1
      ) +
      scale_color_manual(
        values = time_cols,
        breaks = time_levels,
        drop = FALSE,
        name = "Time after SCI"
      ) +
      scale_x_discrete(
        limits = present_types,
        drop = FALSE,
        expand = expansion(add = c(0.25, 0.25))
      ) +
      scale_y_continuous(
        limits = c(0, y_max_D),
        breaks = scales::pretty_breaks(n = 5),
        expand = expansion(mult = c(0, 0.03))
      ) +
      labs(
        x = NULL,
        y = "Ro/e",
        tag = "D"
      ) +
      theme_classic(
        base_family = base_family,
        base_size = 8.0
      ) +
      theme(
        plot.background = element_rect(
          fill = "white",
          color = NA
        ),
        panel.background = element_rect(
          fill = "white",
          color = NA
        ),
        axis.line = element_line(
          color = "#5F5F5F",
          linewidth = 0.38
        ),
        axis.ticks = element_line(
          color = "#5F5F5F",
          linewidth = 0.30
        ),
        axis.text.x = element_text(
          family = base_family,
          size = 6.2,
          color = "black",
          angle = 45,
          hjust = 1,
          vjust = 1
        ),
        axis.text.y = element_text(
          family = base_family,
          size = 6.5,
          color = "black"
        ),
        axis.title.y = element_text(
          family = base_family,
          size = 7.2,
          color = "black",
          margin = margin(r = 3)
        ),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.justification = c(0, 1),
        legend.margin = margin(0, 0, 0, -1.5, unit = "mm"),
        legend.title = element_text(
          family = base_family,
          size = 6.0,
          face = "plain"
        ),
        legend.text = element_text(
          family = base_family,
          size = 5.8
        ),
        legend.key.width = grid::unit(4.5, "mm"),
        legend.key.height = grid::unit(2.3, "mm"),
        legend.spacing.x = grid::unit(0.30, "mm"),
        plot.tag = element_text(
          family = base_family,
          face = "bold",
          size = 16,
          color = "black"
        ),
        plot.tag.position = c(-0.04, 1.04),
        plot.margin = margin(4, 4, 5, 5, unit = "mm")
      ) +
      guides(
        color = guide_legend(
          title.position = "left",
          nrow = 1,
          byrow = TRUE,
          override.aes = list(
            linewidth = 0.62,
            linetype = "dashed",
            size = 2.15
          )
        )
      )

    panel_plots[["D"]] <- pD

    plot_data_audit[["D"]] <- record_plot_data(
      "D",
      input_files[["broad_rds"]],
      nrow(comp_D),
      c("sample_id", "timepoint", "major_cell_type_display", "proportion"),
      "Reference-style Ro/e profile: each biological sample Ro/e equals observed within-sample proportion divided by atlas-wide expected proportion; panel shows timepoint median Ro/e as colored points connected by dashed lines; horizontal line marks Ro/e = 1; no cell-count P values."
    )
  }

  rm(broad)
  gc()
}

# ------------------------------------------------------------
# 7. PANELS E-F | MYELOID LINEAGE ADJUDICATION
# ------------------------------------------------------------
if (any(c("E", "F") %in% panels_to_build)) {

  message("Reading frozen Step2C full myeloid object...")
  myeloid <- readRDS(input_files[["myeloid_rds"]])

  if (!inherits(myeloid, "Seurat")) {
    stop("Step2C full myeloid RDS is not a Seurat object.")
  }

  DefaultAssay(myeloid) <- "RNA"

  check_meta(
    myeloid,
    c(
      "myeloid_cluster_r04",
      "sample_id",
      "timepoint",
      "final_step2_lineage"
    ),
    "Step2C full myeloid"
  )
  check_reduction(
    myeloid,
    "umap.myeloid",
    "Step2C full myeloid"
  )

  schema_audit[["myeloid"]] <- schema_row(
    myeloid,
    "step2C_full_myeloid"
  )

  myeloid$lineage_display <- collapse_final_lineage(
    myeloid$final_step2_lineage
  )

  myeloid$lineage_display <- factor(
    myeloid$lineage_display,
    levels = names(lineage_cols)
  )

  myeloid_umap <- extract_umap_df(
    myeloid,
    reduction = "umap.myeloid",
    extra_cols = c(
      "sample_id",
      "timepoint",
      "myeloid_cluster_r04",
      "final_step2_lineage",
      "lineage_display"
    )
  )

  myeloid_umap <- freeze_point_order(
    myeloid_umap,
    seed = 20260812
  )

  # E
  if ("E" %in% panels_to_build) {
    pE <- make_umap(
      myeloid_umap,
      color_col = "lineage_display",
      color_values = lineage_cols,
      point_size = 0.23,
      label_groups = TRUE,
      label_subgroup_col = "myeloid_cluster_r04",
      label_map = c(
        "Resident-derived microglia" = "Resident-derived\nmicroglia",
        "Infiltrating mono/macrophage" = "Infiltrating\nmono/macrophage",
        "BAM-like macrophage" = "BAM-like\nmacrophage"
      ),
      label_size = 2.18,
      label_force = 0.70,
      label_nudge_x_frac = 0.025,
      label_nudge_y_frac = -0.020,
      label_extra_nudge_df = data.frame(
        label = "Resident-derived\nmicroglia",
        dx_frac = 0.020,
        dy_frac = -0.018,
        stringsAsFactors = FALSE
      ),
      show_legend = TRUE,
      total_n_label = FALSE,
      seed_tag = "Panel E"
    ) +
      labs(
        color = "Final lineage",
        tag = "E"
      ) +
      theme_pub +
      theme(
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "right"
      )

    panel_plots[["E"]] <- pE

    plot_data_audit[["E"]] <- record_plot_data(
      "E",
      input_files[["myeloid_rds"]],
      nrow(myeloid_umap),
      c("umap.myeloid", "final_step2_lineage"),
      "Step2C final lineage; excluded nonresident populations are not reclassified."
    )
  }

  # F
  if ("F" %in% panels_to_build) {
    marker_F <- tribble(
      ~marker_group,     ~gene,
      "Resident",        "P2ry12",
      "Resident",        "Tmem119",
      "Resident",        "Fcrls",
      "Resident",        "Hexb",
      "Infiltrating",    "Ccr2",
      "Infiltrating",    "Ly6c2",
      "Infiltrating",    "Arg1",
      "Infiltrating",    "F13a1"
    )

    pF <- make_dotplot_reference(
      myeloid,
      group_col = "lineage_display",
      marker_df = marker_F,
      tag_letter = "F",
      group_labels = c(
        "Resident-derived microglia" = "Resident-derived\nmicroglia",
        "Infiltrating mono/macrophage" = "Infiltrating\nmono/macrophage",
        "BAM-like macrophage" = "BAM-like\nmacrophage"
      ),
      dot_max = 4.10,
      x_text_size = 6.0,
      y_text_size = 6.6
    ) +
      theme(
        legend.box.margin = margin(0, 0, 0, -8.5, unit = "mm")
      )

    panel_plots[["F"]] <- pF

    plot_data_audit[["F"]] <- record_plot_data(
      "F",
      input_files[["myeloid_rds"]],
      length(unique(marker_F$gene)),
      c("lineage_display", "RNA normalized data"),
      "Frozen resident-vs-infiltrating marker panel from the Figure1 plan; no marker discovery."
    )
  }

  rm(myeloid)
  gc()
}

# ------------------------------------------------------------
# 8. PANELS G-H | CLEAN RESIDENT MICROGLIA
# ------------------------------------------------------------
if (any(c("G", "H") %in% panels_to_build)) {

  message("Reading frozen Step3B cleaned resident-microglia object...")
  micro <- readRDS(input_files[["microglia_rds"]])

  if (!inherits(micro, "Seurat")) {
    stop("Step3B clean state-scored RDS is not a Seurat object.")
  }

  DefaultAssay(micro) <- "RNA"

  check_meta(
    micro,
    c(
      "step3A_cluster_r04",
      "sample_id",
      "timepoint",
      "step2_final_cluster"
    ),
    "Step3B clean microglia"
  )
  check_reduction(
    micro,
    "umap.microglia.step3A",
    "Step3B clean microglia"
  )

  schema_audit[["microglia"]] <- schema_row(
    micro,
    "step3B_clean_microglia"
  )

  # Exact Step2C frozen lineage classes, carried into Step3.
  # The two cell-cycle phases are combined only for display.
  micro$resident_state_display <- collapse_resident_state(
    micro$step2_final_cluster
  )

  if (any(is.na(micro$resident_state_display))) {
    bad_clusters <- sort(
      unique(
        as.character(
          micro$step2_final_cluster[
            is.na(micro$resident_state_display)
          ]
        )
      )
    )
    stop(
      "Unexpected step2_final_cluster value(s) in Step3B clean object: ",
      paste(bad_clusters, collapse = ", "),
      ". No state label was guessed."
    )
  }

  micro$resident_state_display <- factor(
    micro$resident_state_display,
    levels = names(resident_state_cols)
  )

  micro$timepoint <- factor(
    as.character(micro$timepoint),
    levels = time_levels
  )

  micro_umap <- extract_umap_df(
    micro,
    reduction = "umap.microglia.step3A",
    extra_cols = c(
      "sample_id",
      "timepoint",
      "step2_final_cluster",
      "resident_state_display"
    )
  )

  micro_umap <- freeze_point_order(
    micro_umap,
    seed = 20260813
  )

  # G
  if ("G" %in% panels_to_build) {
    pG <- make_umap(
      micro_umap,
      color_col = "resident_state_display",
      color_values = resident_state_cols,
      point_size = 0.24,
      label_groups = TRUE,
      label_subgroup_col = "step2_final_cluster",
      label_map = c(
        "Canonical resident" = "Canonical\nresident",
        "Activated resident" = "Activated\nresident",
        "IFN-responsive resident" = "IFN-responsive\nresident",
        "Cycling resident" = "Cycling\nresident"
      ),
      label_size = 2.22,
      label_force = 0.68,
      show_legend = TRUE,
      total_n_label = FALSE,
      seed_tag = "Panel G"
    ) +
      labs(
        color = "Resident state class",
        tag = "G"
      ) +
      theme_pub +
      theme(
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "right"
      )

    panel_plots[["G"]] <- pG

    plot_data_audit[["G"]] <- record_plot_data(
      "G",
      input_files[["microglia_rds"]],
      nrow(micro_umap),
      c(
        "umap.microglia.step3A",
        "step2_final_cluster",
        "resident_state_display"
      ),
      paste(
        "Step3B residual-cleaned resident microglia.",
        "Display classes are exact Step2C frozen lineage calls;",
        "no DAM/acute-inflammatory state is newly assigned in Figure1."
      )
    )
  }

  # H
  if ("H" %in% panels_to_build) {

    comp_H <- micro[[]] %>%
      transmute(
        sample_id = as.character(sample_id),
        timepoint = factor(as.character(timepoint), levels = time_levels),
        resident_state_display = as.character(resident_state_display)
      ) %>%
      filter(
        !is.na(sample_id),
        !is.na(timepoint),
        !is.na(resident_state_display)
      ) %>%
      count(
        sample_id,
        timepoint,
        resident_state_display,
        name = "n_cells"
      ) %>%
      group_by(sample_id, timepoint) %>%
      tidyr::complete(
        resident_state_display = names(resident_state_cols),
        fill = list(n_cells = 0)
      ) %>%
      mutate(
        sample_total = sum(n_cells),
        proportion = ifelse(
          sample_total > 0,
          n_cells / sample_total,
          NA_real_
        )
      ) %>%
      ungroup() %>%
      filter(!is.na(proportion))

    # Atlas-wide expected state proportion across the frozen clean-microglia set.
    expected_H <- micro[[]] %>%
      transmute(
        resident_state_display = as.character(resident_state_display)
      ) %>%
      filter(!is.na(resident_state_display)) %>%
      count(
        resident_state_display,
        name = "overall_cells"
      ) %>%
      mutate(
        expected_proportion = overall_cells / sum(overall_cells)
      )

    roe_H <- comp_H %>%
      left_join(
        expected_H,
        by = "resident_state_display"
      ) %>%
      mutate(
        roe = ifelse(
          expected_proportion > 0,
          proportion / expected_proportion,
          NA_real_
        ),
        resident_state_display = factor(
          resident_state_display,
          levels = names(resident_state_cols)
        )
      )

    summary_H <- roe_H %>%
      group_by(timepoint, resident_state_display) %>%
      summarise(
        median_roe = median(roe, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        timepoint = factor(
          as.character(timepoint),
          levels = time_levels
        ),
        resident_state_display = factor(
          as.character(resident_state_display),
          levels = names(resident_state_cols)
        )
      )

    finite_roe_H <- summary_H$median_roe[
      is.finite(summary_H$median_roe)
    ]

    if (length(finite_roe_H) == 0) {
      stop("Panel H: no finite Ro/e values were available.")
    }

    y_max_H <- max(
      1.25,
      as.numeric(max(finite_roe_H, na.rm = TRUE)) * 1.10
    )

    # D-style Ro/e profile: timepoint-colored circular points + dashed lines.
    pH <- ggplot(
      summary_H,
      aes(
        x = resident_state_display,
        y = median_roe,
        color = timepoint,
        group = timepoint
      )
    ) +
      geom_hline(
        yintercept = 1,
        color = "#8A8A8A",
        linewidth = 0.34,
        linetype = "solid"
      ) +
      geom_line(
        linewidth = 0.62,
        linetype = "dashed",
        alpha = 0.95
      ) +
      geom_point(
        size = 2.15,
        shape = 16,
        alpha = 1
      ) +
      scale_color_manual(
        values = time_cols,
        breaks = time_levels,
        drop = FALSE,
        name = "Time after SCI"
      ) +
      scale_x_discrete(
        labels = c(
          "Canonical resident" = "Canonical\nresident",
          "Activated resident" = "Activated\nresident",
          "IFN-responsive resident" = "IFN-responsive\nresident",
          "Cycling resident" = "Cycling\nresident"
        ),
        drop = FALSE,
        expand = expansion(add = c(0.25, 0.25))
      ) +
      scale_y_continuous(
        limits = c(0, y_max_H),
        breaks = scales::pretty_breaks(n = 5),
        expand = expansion(mult = c(0, 0.03))
      ) +
      labs(
        x = NULL,
        y = "Ro/e",
        tag = "H"
      ) +
      theme_classic(
        base_family = base_family,
        base_size = 8.0
      ) +
      theme(
        plot.background = element_rect(
          fill = "white",
          color = NA
        ),
        panel.background = element_rect(
          fill = "white",
          color = NA
        ),
        axis.line = element_line(
          color = "#5F5F5F",
          linewidth = 0.38
        ),
        axis.ticks = element_line(
          color = "#5F5F5F",
          linewidth = 0.30
        ),
        axis.text.x = element_text(
          family = base_family,
          size = 6.3,
          color = "black",
          angle = 25,
          hjust = 1,
          vjust = 1,
          lineheight = 0.92
        ),
        axis.text.y = element_text(
          family = base_family,
          size = 6.5,
          color = "black"
        ),
        axis.title.y = element_text(
          family = base_family,
          size = 7.2,
          color = "black",
          margin = margin(r = 3)
        ),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.justification = c(0, 1),
        legend.margin = margin(0, 0, 0, -1.5, unit = "mm"),
        legend.title = element_text(
          family = base_family,
          size = 6.0,
          face = "plain"
        ),
        legend.text = element_text(
          family = base_family,
          size = 5.8
        ),
        legend.key.width = grid::unit(4.5, "mm"),
        legend.key.height = grid::unit(2.3, "mm"),
        legend.spacing.x = grid::unit(0.30, "mm"),
        plot.tag = element_text(
          family = base_family,
          face = "bold",
          size = 16,
          color = "black"
        ),
        plot.tag.position = c(-0.04, 1.04),
        plot.margin = margin(4, 4, 5, 5, unit = "mm")
      ) +
      guides(
        color = guide_legend(
          title.position = "left",
          nrow = 1,
          byrow = TRUE,
          override.aes = list(
            linewidth = 0.62,
            linetype = "dashed",
            size = 2.15
          )
        )
      )

    panel_plots[["H"]] <- pH

    plot_data_audit[["H"]] <- record_plot_data(
      "H",
      input_files[["microglia_rds"]],
      nrow(roe_H),
      c(
        "sample_id",
        "timepoint",
        "resident_state_display",
        "proportion",
        "expected_proportion",
        "roe"
      ),
      paste(
        "D-style sample-level Ro/e profile.",
        "Each sample Ro/e = observed within-sample state proportion / overall expected state proportion;",
        "colored points and dashed lines show timepoint medians;",
        "horizontal line marks Ro/e = 1; no cell-count P values."
      )
    )
  }

  rm(micro)
  gc()
}

# ------------------------------------------------------------
# 9. WRITE AUDITS
# ------------------------------------------------------------
schema_df <- bind_rows(schema_audit)
write.csv(
  schema_df,
  file.path(out_dir, "01_schema_audit.csv"),
  row.names = FALSE
)

plot_data_df <- bind_rows(plot_data_audit)
write.csv(
  plot_data_df,
  file.path(out_dir, "02_plot_data_used.csv"),
  row.names = FALSE
)

panel_layout <- data.frame(
  panel = LETTERS[1:8],
  row = c(1,1,2,2,3,3,4,4),
  column = c(1,2,1,2,1,2,1,2),
  individual_width_mm = c(92,92,92,92,92,92,92,92),
  individual_height_mm = c(62,62,70,70,66,66,62,62),
  complete_canvas = paste0(complete_width_mm, " mm x ", complete_height_mm, " mm"),
  stringsAsFactors = FALSE
)

write.csv(
  panel_layout,
  file.path(out_dir, "03_panel_layout_audit.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. SAVE INDIVIDUAL PANELS
# ------------------------------------------------------------
panel_dims <- list(
  A = c(92, 62),
  B = c(92, 62),
  C = c(98, 70),
  D = c(92, 70),
  E = c(92, 66),
  F = c(98, 66),
  G = c(92, 62),
  H = c(92, 66)
)

for (letter in intersect(names(panel_plots), panels_to_build)) {
  dims <- panel_dims[[letter]]

  message("Saving Figure1", letter, " ...")
  save_plot_triplet(
    panel_plots[[letter]],
    stem = paste0("Figure1", letter),
    width_mm = dims[1],
    height_mm = dims[2],
    dpi = 600
  )
}

# ------------------------------------------------------------
# 11. COMPLETE FIGURE
# ------------------------------------------------------------
if (all(LETTERS[1:8] %in% names(panel_plots))) {

  # v24:
  # Return to the exact v18 4 x 2 layout.
  #
  # The ONLY change is the composite panel-letter system:
  # - A, C, E, G use the same fixed left-column x/y anchor.
  # - B, D, F, H use the same fixed right-column x/y anchor.
  #
  # We remove each ggplot's internal tag only in the COMPLETE figure,
  # then overlay the panel letter in the same fixed position inside each
  # equal-width panel container. This avoids differences caused by UMAP,
  # DotPlot, and Ro/e plot margins/legends.
  #
  # Standalone Figure1A.pdf ... Figure1H.pdf are NOT changed by this block.

  make_aligned_complete_panel <- function(letter) {

    # Remove the original internal tag for composite assembly only.
    p_no_tag <- panel_plots[[letter]] +
      labs(tag = NULL) +
      theme(
        plot.tag = element_blank()
      )

    # Identical letter coordinate for every panel container.
    tag_overlay <- ggplot() +
      xlim(0, 1) +
      ylim(0, 1) +
      annotate(
        "text",
        x = 0.012,
        y = 0.992,
        label = letter,
        hjust = 0,
        vjust = 1,
        family = base_family,
        fontface = "bold",
        size = 16 / .pt,
        color = "black"
      ) +
      theme_void()

    patchwork::wrap_elements(
      full = p_no_tag,
      clip = FALSE
    ) +
      patchwork::inset_element(
        tag_overlay,
        left = 0,
        bottom = 0,
        right = 1,
        top = 1,
        align_to = "full",
        clip = FALSE
      )
  }

  fixed_panels <- lapply(
    LETTERS[1:8],
    make_aligned_complete_panel
  )
  names(fixed_panels) <- LETTERS[1:8]

  # Exact v18 4 x 2 geometry.
  p_complete <-
    (fixed_panels$A | fixed_panels$B) /
    (fixed_panels$C | fixed_panels$D) /
    (fixed_panels$E | fixed_panels$F) /
    (fixed_panels$G | fixed_panels$H) +
    plot_layout(
      widths = c(1, 1),
      heights = c(1.00, 1.15, 1.08, 1.00)
    ) +
    plot_annotation(
      theme = theme(
        plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(3, 3, 3, 3, unit = "mm")
      )
    )

  message("Saving complete Figure 1 v29 spacing refinement (180 x 230 mm; Figure-4-matched panel letters) ...")
  save_plot_triplet(
    p_complete,
    stem = "Figure_1_complete_v29_CFleft_Elabel_fix",
    width_mm = complete_width_mm,
    height_mm = complete_height_mm,
    dpi = 600
  )
}

# v29 visual-standard audit
write.csv(
  data.frame(
    setting = c(
      "complete_canvas",
      "composite_panel_letter",
      "font_family",
      "C_x_text_pt",
      "F_x_text_pt",
      "CF_legend_left_shift_mm",
      "E_global_label_nudge_x_frac",
      "E_global_label_nudge_y_frac",
      "E_resident_derived_extra_dx_frac",
      "E_resident_derived_extra_dy_frac"
    ),
    value = c(
      paste0(complete_width_mm, " x ", complete_height_mm, " mm"),
      "16 pt Arial Bold",
      base_family,
      "5.2",
      "6.0",
      "-8.5",
      "0.025",
      "-0.020",
      "0.020",
      "-0.018"
    ),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "04_v29_GLIA_visual_standard_audit.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 12. SESSION INFO + FINAL REPORT
# ------------------------------------------------------------
capture.output(
  sessionInfo(),
  file = file.path(out_dir, "sessionInfo.txt")
)

message("")
message("============================================================")
message("Figure 1 generation finished.")
message("Output directory:")
message(out_dir)
message("")
message("Key review files:")
message("1) Figure_1_complete_v29_CFleft_Elabel_fix.pdf")
message("2) Figure1A.pdf")
message("3) Figure1B.pdf")
message("4) Figure1G.pdf")
message("5) Figure1H.pdf")
message("6) 01_schema_audit.csv")
message("7) 02_plot_data_used.csv")
message("============================================================")
