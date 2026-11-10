# ============================================================
# v4 GLIA readability revision: larger axis/legend text only; communication data and visual encoding unchanged.
# Figure 4C v3 | Reference-style bubble plot
# Signaling from Target program
#
# GOAL
#   Redraw Figure 4C so it visually follows the reference panel:
#     - single source: Target program
#     - x-axis: source -> recipient combinations
#     - y-axis: curated ligand-receptor pairs
#     - color: communication probability
#     - size : p-value class
#
# IMPORTANT
#   - Uses the ORIGINAL frozen Step8B nboot100 communication table
#   - DOES NOT rerun CellChat analysis
#   - Outputs one-page PDF + PNG
#
# INPUT
#   04_results/
#   step8B_CONFIRM_GSE162610_7dpi_targeted_fullcell_nboot100/
#   tables/03_fullcell_targeted_all_communications.tsv
# ============================================================

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------
pkgs <- c("ggplot2", "dplyr", "readr", "stringr", "grid", "scales")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(miss) > 0) {
  stop(
    "Missing package(s): ", paste(miss, collapse = ", "),
    "\nInstall with:\ninstall.packages(c(",
    paste(sprintf('\"%s\"', miss), collapse = ", "), "))"
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(grid)
  library(scales)
})

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------
ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
COMM_FILE_PRIMARY <- file.path(
  ROOT,
  "04_results",
  "step8B_CONFIRM_GSE162610_7dpi_targeted_fullcell_nboot100",
  "tables",
  "03_fullcell_targeted_all_communications.tsv"
)

find_exact_file <- function(root, filename) {
  hits <- list.files(
    root,
    pattern = paste0("^", gsub("\\.", "\\\\.", filename), "$"),
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

if (file.exists(COMM_FILE_PRIMARY)) {
  COMM_FILE <- COMM_FILE_PRIMARY
} else {
  COMM_FILE <- find_exact_file(ROOT, "03_fullcell_targeted_all_communications.tsv")
}

if (is.na(COMM_FILE) || !file.exists(COMM_FILE)) {
  stop(
    "Could not find the frozen communication table.\n\nExpected:\n",
    COMM_FILE_PRIMARY
  )
}

OUTDIR <- file.path(
  ROOT,
  "09_final_figures",
  "GLIA",
  "Figure4",
  "single_panels",
  "Figure4C"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

PDF_FILE <- file.path(
  OUTDIR,
  "Figure4C_signaling_from_target_program_reference_style_v4_GLIA_readability.pdf"
)

PNG_FILE <- file.path(
  OUTDIR,
  "Figure4C_signaling_from_target_program_reference_style_v4_GLIA_readability.png"
)

if (file.exists(PDF_FILE)) unlink(PDF_FILE)
if (file.exists(PNG_FILE)) unlink(PNG_FILE)

# ------------------------------------------------------------
# 2. Read the original frozen table
# ------------------------------------------------------------
comm <- readr::read_tsv(COMM_FILE, show_col_types = FALSE)

required_cols <- c(
  "timepoint", "source", "target", "ligand", "receptor",
  "interaction_name", "pathway_name", "probability", "p_value"
)

missing_cols <- setdiff(required_cols, colnames(comm))
if (length(missing_cols) > 0) {
  stop(
    "Missing required column(s): ",
    paste(missing_cols, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 3. Fixed cell order (same 10 compartments as Figure 4A/B)
# ------------------------------------------------------------
recipient_raw <- c(
  "Microglia_Arm_2",
  "Other_myeloid",
  "Microglia_other",
  "Neutrophil",
  "T_cell",
  "Astrocyte",
  "Endothelial",
  "Fibroblast",
  "Pericyte_VSMC",
  "OPC"
)

recipient_display <- c(
  "Microglia_Arm_2" = "Target program",
  "Other_myeloid"   = "Other myeloid",
  "Microglia_other" = "Microglia other",
  "Neutrophil"      = "Neutrophil",
  "T_cell"          = "T cell",
  "Astrocyte"       = "Astrocyte",
  "Endothelial"     = "Endothelial",
  "Fibroblast"      = "Fibroblast",
  "Pericyte_VSMC"   = "Pericyte/VSMC",
  "OPC"             = "OPC"
)

# ------------------------------------------------------------
# 4. Keep ONLY 7 dpi outgoing communication from Target program
# ------------------------------------------------------------
plot_df <- comm %>%
  filter(
    timepoint == "7dpi",
    source == "Microglia_Arm_2",
    target %in% recipient_raw,
    is.finite(probability),
    is.finite(p_value),
    p_value <= 0.05
  )

if (nrow(plot_df) == 0) {
  stop("No significant outgoing communication remained after filtering.")
}

# ------------------------------------------------------------
# 5. Curated LR order
#
# To make the panel closer to the reference figure, use a fixed
# biologically meaningful pathway order and then LR order within it.
# ------------------------------------------------------------
interaction_order <- c(
  "TNF_TNFRSF1A",
  "TNF_TNFRSF1B",
  "ICAM1_ITGAM_ITGB2",
  "ICAM1_ITGAM_ITGB2L",
  "CCL3_CCR1",
  "CCL2_CCR2",
  "OSM_LIFR_IL6ST",
  "OSM_OSMR_IL6ST",
  "SPP1_CD44",
  "SPP1_ITGAV_ITGB1",
  "SPP1_ITGAV_ITGB5",
  "SPP1_ITGA5_ITGB1",
  "TGFB1_TGFBR1_TGFBR2",
  "ENTPD1_ADORA1",
  "ENTPD1_ADORA2B",
  "PPIA_BSG"
)

# Keep only interactions that are actually present
interaction_order_present <- interaction_order[
  interaction_order %in% unique(plot_df$interaction_name)
]

# Any unexpected extra interactions are appended at the end
extra_interactions <- setdiff(
  unique(plot_df$interaction_name),
  interaction_order_present
)

interaction_order_final <- c(
  interaction_order_present,
  sort(extra_interactions)
)

# ------------------------------------------------------------
# 6. Labels
# ------------------------------------------------------------
pretty_lr_label <- function(ligand, receptor) {
  rec <- receptor %>%
    str_replace_all("_", "+")
  lig <- str_to_title(ligand)
  paste0(lig, " \u2013 ", rec)
}

label_map <- plot_df %>%
  distinct(interaction_name, ligand, receptor) %>%
  mutate(
    lr_label = pretty_lr_label(ligand, receptor)
  )

plot_df <- plot_df %>%
  left_join(
    label_map %>% select(interaction_name, lr_label),
    by = "interaction_name"
  ) %>%
  mutate(
    recipient_label = unname(recipient_display[target]),
    x_label = paste0("Target program \u2192 ", recipient_label),
    p_class = case_when(
      p_value < 0.01 ~ "p < 0.01",
      p_value <= 0.05 ~ "0.01 < p \u2264 0.05",
      TRUE ~ NA_character_
    )
  )

# x-axis order fixed
x_levels <- paste0(
  "Target program \u2192 ",
  unname(recipient_display[recipient_raw])
)

plot_df$x_label <- factor(plot_df$x_label, levels = x_levels)

# y-axis order fixed
y_levels <- label_map %>%
  mutate(interaction_name = as.character(interaction_name)) %>%
  filter(interaction_name %in% interaction_order_final) %>%
  arrange(match(interaction_name, interaction_order_final)) %>%
  pull(lr_label)

plot_df$lr_label <- factor(plot_df$lr_label, levels = rev(unique(y_levels)))

plot_df$p_class <- factor(
  plot_df$p_class,
  levels = c("0.01 < p \u2264 0.05", "p < 0.01")
)

write.csv(
  plot_df,
  file.path(OUTDIR, "Figure4C_v4_plot_data.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Font
# ------------------------------------------------------------
FONT <- "Arial"

if (.Platform$OS.type == "windows") {
  try(
    grDevices::windowsFonts(
      Arial = grDevices::windowsFont("Arial")
    ),
    silent = TRUE
  )
}

# ------------------------------------------------------------
# 8. Probability scale
#
# To avoid almost all dots collapsing into the blue end, use
# quantile-spaced colour anchors while still labeling min/max.
# This keeps the plot more visually similar to the reference.
# ------------------------------------------------------------
prob_min <- min(plot_df$probability, na.rm = TRUE)
prob_q25 <- as.numeric(quantile(plot_df$probability, 0.25, na.rm = TRUE))
prob_q50 <- as.numeric(quantile(plot_df$probability, 0.50, na.rm = TRUE))
prob_q75 <- as.numeric(quantile(plot_df$probability, 0.75, na.rm = TRUE))
prob_max <- max(plot_df$probability, na.rm = TRUE)

colour_values <- rescale(
  c(prob_min, prob_q25, prob_q50, prob_q75, prob_max),
  to = c(0, 1)
)

# ------------------------------------------------------------
# 9. Plot
# ------------------------------------------------------------
p <- ggplot(
  plot_df,
  aes(
    x = x_label,
    y = lr_label,
    colour = probability,
    size = p_class
  )
) +
  geom_point(
    alpha = 0.98,
    na.rm = TRUE
  ) +

  scale_colour_gradientn(
    colours = c(
      "#313695",  # deep blue
      "#2C7BB6",
      "#66C2A5",
      "#FEE08B",
      "#D73027"   # red
    ),
    values = colour_values,
    limits = c(prob_min, prob_max),
    breaks = c(prob_min, prob_max),
    labels = c("min", "max"),
    name = "Commun. Prob."
  ) +

  scale_size_manual(
    values = c(
      "0.01 < p \u2264 0.05" = 2.8,
      "p < 0.01" = 4.5
    ),
    name = "p-value",
    drop = FALSE
  ) +

  labs(
    x = NULL,
    y = NULL,
    title = "Signaling from Target program"
  ) +

  theme_bw(base_family = FONT) +

  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 15,
      face = "plain",
      margin = margin(b = 5)
    ),

    panel.grid.major = element_line(
      colour = "#D3D3D3",
      linewidth = 0.38,
      linetype = "dotted"
    ),

    panel.grid.minor = element_blank(),

    panel.border = element_rect(
      colour = "#555555",
      linewidth = 0.6
    ),

    axis.text.x = element_text(
      angle = 48,
      hjust = 1,
      vjust = 1,
      size = 9.5,
      colour = "#17A2B8"
    ),

    axis.text.y = element_text(
      size = 9.5,
      colour = "black"
    ),

    axis.ticks = element_line(
      colour = "#666666",
      linewidth = 0.35
    ),

    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.18, "cm"),

    legend.title = element_text(
      size = 9.2,
      face = "plain"
    ),

    legend.text = element_text(
      size = 8.4
    ),

    plot.margin = margin(
      t = 18, r = 20, b = 24, l = 18
    )
  ) +

  guides(
    colour = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(2.7, "cm"),
      barwidth = unit(0.34, "cm"),
      order = 1
    ),
    size = guide_legend(
      title.position = "top",
      order = 2,
      override.aes = list(
        colour = "black",
        alpha = 1
      )
    )
  )

# ------------------------------------------------------------
# 10. Save helper
# ------------------------------------------------------------
draw_with_panel_letter <- function(plot_obj) {
  grid.newpage()
  grid.draw(ggplotGrob(plot_obj))

  grid.text(
    "C",
    x = unit(0.014, "npc"),
    y = unit(0.985, "npc"),
    just = c("left", "top"),
    gp = gpar(
      fontsize = 27,
      fontface = "bold",
      fontfamily = FONT
    )
  )
}

# ------------------------------------------------------------
# 11. Save PDF and PNG
# ------------------------------------------------------------
if (capabilities("cairo")) {
  grDevices::cairo_pdf(
    filename = PDF_FILE,
    width = 8.9,
    height = 7.5,
    family = FONT,
    onefile = FALSE,
    bg = "white"
  )
} else {
  grDevices::pdf(
    file = PDF_FILE,
    width = 8.9,
    height = 7.5,
    family = "Helvetica",
    onefile = FALSE,
    useDingbats = FALSE,
    bg = "white"
  )
}

draw_with_panel_letter(p)
dev.off()

grDevices::png(
  filename = PNG_FILE,
  width = 2850,
  height = 2400,
  res = 320,
  bg = "white"
)

draw_with_panel_letter(p)
dev.off()

# ------------------------------------------------------------
# 12. Final audit
# ------------------------------------------------------------
write.csv(
  data.frame(
    input_table = COMM_FILE,
    n_plot_rows = nrow(plot_df),
    n_unique_lr = length(unique(plot_df$interaction_name)),
    n_unique_targets = length(unique(plot_df$target)),
    pdf_exists = file.exists(PDF_FILE),
    png_exists = file.exists(PNG_FILE),
    stringsAsFactors = FALSE
  ),
  file.path(OUTDIR, "Figure4C_v3_output_audit.csv"),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("Figure 4C v3 completed successfully.\n")
cat("Reference-style bubble plot rebuilt from the original Step8B table.\n")
cat("Source fixed at Target program (Microglia_Arm_2), 7 dpi.\n")
cat("No CellChat analysis was rerun.\n\n")
cat("Input table:\n  ", COMM_FILE, "\n", sep = "")
cat("PDF:\n  ", PDF_FILE, "\n", sep = "")
cat("PNG:\n  ", PNG_FILE, "\n", sep = "")
cat("============================================================\n")
