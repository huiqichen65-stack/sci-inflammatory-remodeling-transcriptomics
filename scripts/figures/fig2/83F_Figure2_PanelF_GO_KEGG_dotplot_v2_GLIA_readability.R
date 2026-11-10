# ==============================================================================
# 83F_Figure2_PanelF_GO_KEGG_dotplot_v2_GLIA_readability.R
# v2 visual-only revision: enlarge pathway/axis/legend text for final 180-mm master.
# Figure 2F | Target-trajectory GO BP + KEGG enrichment dot plot
#
# Goal
# ----
# Draw a publication-style two-panel enrichment dot plot similar to the
# user-provided reference:
#   - upper panel: GO BP
#   - lower panel: KEGG
#   - x-axis: Gene ratio
#   - point size: Gene count
#   - point color: -log10(FDR)
#
# Frozen-data contract
# --------------------
# - Read frozen Step6CR GO/KEGG tables only.
# - Do NOT rerun enrichment analysis.
# - Select representative terms from the already-frozen significant results.
# ==============================================================================

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 0. USER SETTINGS
# ------------------------------------------------------------------------------
project_root <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
out_dir <- file.path(
  project_root,
  "09_final_figures", "GLIA", "Figure2_panelF_v2_GLIA_readability"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

panel_width_mm  <- 148
panel_height_mm <- 78

top_n_go   <- 6
top_n_kegg <- 4

# ------------------------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------------------------
required_pkgs <- c(
  "ggplot2",
  "dplyr",
  "stringr",
  "scales",
  "svglite",
  "ragg"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them first, then rerun."
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(scales)
})

# ------------------------------------------------------------------------------
# 2. INPUT
# ------------------------------------------------------------------------------
input_files <- c(
  step6CR_GO = file.path(
    project_root,
    "04_results/step6CR_GSE162610_branch_DEG_GO_KEGG",
    "tables", "01_all_branch_GO_results.csv"
  ),
  step6CR_KEGG = file.path(
    project_root,
    "04_results/step6CR_GSE162610_branch_DEG_GO_KEGG",
    "tables", "02_all_branch_KEGG_results.csv"
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
  stop("Required frozen input missing. See 00_input_audit.csv.")
}

message("Reading frozen Step6CR GO / KEGG tables...")
go_tbl <- read.csv(input_files[["step6CR_GO"]], check.names = FALSE)
kegg_tbl <- read.csv(input_files[["step6CR_KEGG"]], check.names = FALSE)

# ------------------------------------------------------------------------------
# 3. HELPERS
# ------------------------------------------------------------------------------
first_existing_col <- function(df, candidates, required = TRUE, context = "table") {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) {
    if (required) {
      stop(
        context, " is missing expected column(s): ",
        paste(candidates, collapse = ", "),
        ". Observed columns: ",
        paste(colnames(df), collapse = ", ")
      )
    } else {
      return(NULL)
    }
  }
  hit[1]
}

safe_numeric <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  if (all(is.na(out))) return(rep(NA_real_, length(x)))
  out
}

parse_ratio <- function(x) {
  x <- as.character(x)
  out <- rep(NA_real_, length(x))
  slash_idx <- grepl("/", x, fixed = TRUE)

  # numeric already?
  numeric_try <- suppressWarnings(as.numeric(x))
  out[!is.na(numeric_try)] <- numeric_try[!is.na(numeric_try)]

  # ratio like 5/12
  if (any(slash_idx)) {
    parts <- strsplit(x[slash_idx], "/", fixed = TRUE)
    vals <- vapply(parts, function(p) {
      if (length(p) != 2) return(NA_real_)
      a <- suppressWarnings(as.numeric(p[1]))
      b <- suppressWarnings(as.numeric(p[2]))
      if (is.na(a) || is.na(b) || b == 0) return(NA_real_)
      a / b
    }, numeric(1))
    out[slash_idx] <- vals
  }

  out
}

choose_target_branch <- function(df) {
  if (!"branch" %in% colnames(df)) return(df)
  branches <- unique(as.character(df$branch))
  target_like <- branches[grepl("Arm_?2|Target", branches, ignore.case = TRUE)]
  if (length(target_like) == 0) return(df)
  df %>% filter(branch %in% target_like)
}

wrap_term <- function(x, width = 29) {
  stringr::str_wrap(as.character(x), width = width)
}

pick_terms <- function(df, top_n, source_name) {
  if (nrow(df) == 0) return(df)

  df2 <- df %>%
    mutate(
      wrapped_term = wrap_term(Description, width = 29),
      description_length = nchar(Description)
    ) %>%
    arrange(
      FDR_value,
      desc(GeneRatio_value),
      desc(Count_value),
      description_length,
      Description
    )

  # Keep only one entry per exact description
  df2 <- df2 %>% distinct(Description, .keep_all = TRUE)

  # Take top N
  df2 <- df2 %>% slice_head(n = top_n)

  # Order terms from top to bottom visually like the reference
  df2 <- df2 %>%
    arrange(GeneRatio_value, FDR_value) %>%
    mutate(term_factor = factor(wrapped_term, levels = wrapped_term))

  df2$source <- source_name
  df2
}

# ------------------------------------------------------------------------------
# 4. STANDARDIZE GO / KEGG TABLES
# ------------------------------------------------------------------------------
standardize_enrich_tbl <- function(df, source_name) {

  desc_col <- first_existing_col(
    df, c("Description"),
    context = paste(source_name, "table")
  )

  fdr_col <- first_existing_col(
    df, c("best_FDR", "p.adjust", "qvalue"),
    required = FALSE,
    context = paste(source_name, "table")
  )

  count_col <- first_existing_col(
    df, c("Count", "maximum_count"),
    required = FALSE,
    context = paste(source_name, "table")
  )

  ratio_col <- first_existing_col(
    df, c("GeneRatio", "gene_ratio", "RichFactor", "FoldEnrichment"),
    required = FALSE,
    context = paste(source_name, "table")
  )

  sig_col <- first_existing_col(
    df, c("significant"),
    required = FALSE,
    context = paste(source_name, "table")
  )

  out <- df
  out$Description <- as.character(out[[desc_col]])

  if (!is.null(fdr_col)) {
    out$FDR_value <- safe_numeric(out[[fdr_col]])
  } else {
    out$FDR_value <- NA_real_
  }

  if (!is.null(count_col)) {
    out$Count_value <- safe_numeric(out[[count_col]])
  } else {
    out$Count_value <- NA_real_
  }

  if (!is.null(ratio_col)) {
    out$GeneRatio_value <- parse_ratio(out[[ratio_col]])
  } else {
    out$GeneRatio_value <- NA_real_
  }

  if (!is.null(sig_col)) {
    out$significant_flag <- as.logical(out[[sig_col]])
  } else {
    out$significant_flag <- TRUE
  }

  out$source_name <- source_name
  out
}

go_std   <- standardize_enrich_tbl(go_tbl, "GO BP")
kegg_std <- standardize_enrich_tbl(kegg_tbl, "KEGG")

# Restrict to target / Arm_2 branch if the column exists
go_std   <- choose_target_branch(go_std)
kegg_std <- choose_target_branch(kegg_std)

# Keep significant rows only, and require usable numbers
go_sig <- go_std %>%
  filter(significant_flag | is.na(significant_flag)) %>%
  filter(!is.na(FDR_value), !is.na(Count_value), !is.na(GeneRatio_value))

kegg_sig <- kegg_std %>%
  filter(significant_flag | is.na(significant_flag)) %>%
  filter(!is.na(FDR_value), !is.na(Count_value), !is.na(GeneRatio_value))

if (nrow(go_sig) == 0) {
  stop("No usable significant GO BP rows found in frozen Step6CR table.")
}
if (nrow(kegg_sig) == 0) {
  stop("No usable significant KEGG rows found in frozen Step6CR table.")
}

go_plot   <- pick_terms(go_sig, top_n_go, "GO BP")
kegg_plot <- pick_terms(kegg_sig, top_n_kegg, "KEGG")

plot_df <- bind_rows(go_plot, kegg_plot) %>%
  mutate(
    neglog10FDR = -log10(FDR_value),
    source = factor(source, levels = c("GO BP", "KEGG"))
  )

if (nrow(plot_df) == 0) {
  stop("No terms selected for Figure 2F.")
}

write.csv(
  plot_df,
  file.path(out_dir, "01_selected_terms.csv"),
  row.names = FALSE
)

# Shared x-axis limits
x_min <- min(plot_df$GeneRatio_value, na.rm = TRUE)
x_max <- max(plot_df$GeneRatio_value, na.rm = TRUE)
x_pad <- max(0.02, 0.06 * (x_max - x_min))
x_limits <- c(max(0, x_min - x_pad), min(1.05, x_max + x_pad))

# Nice size scaling
count_breaks <- sort(unique(plot_df$Count_value))
if (length(count_breaks) > 4) {
  count_breaks <- pretty(range(plot_df$Count_value, na.rm = TRUE), n = 4)
  count_breaks <- count_breaks[count_breaks >= min(plot_df$Count_value) &
                                 count_breaks <= max(plot_df$Count_value)]
}
count_breaks <- unique(round(count_breaks))

# ------------------------------------------------------------------------------
# 5. DISPLAY STYLE
# ------------------------------------------------------------------------------
BASE_FAMILY <- "Arial"

theme_pub <- theme_bw(base_family = BASE_FAMILY, base_size = 9.2) +
  theme(
    panel.grid.major = element_line(color = "#D9D9D9", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "#333333", linewidth = 0.65),
    axis.title = element_text(size = 8.6, color = "black"),
    axis.text = element_text(size = 7.8, color = "black"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.title = element_text(size = 8.0, color = "black"),
    legend.text = element_text(size = 7.4, color = "black"),
    plot.tag = element_text(face = "bold", size = 16, color = "black"),
    plot.tag.position = c(0.01, 0.99),
    plot.margin = margin(7, 5, 4, 7, "mm")
  )

pF <- ggplot(
  plot_df,
  aes(
    x = GeneRatio_value,
    y = term_factor,
    size = Count_value,
    colour = neglog10FDR
  )
) +
  geom_point(alpha = 0.95) +
  facet_grid(
    source ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_size_continuous(
    name = "Gene count",
    range = c(1.8, 6.6),
    breaks = count_breaks
  ) +
  scale_colour_gradient(
    name = expression(-log[10](FDR)),
    low  = "#4A76B8",
    high = "#DE352D"
  ) +
  scale_x_continuous(
    name = "Gene ratio",
    limits = x_limits,
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    y = NULL,
    tag = "F"
  ) +
  theme_pub +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(
      angle = 90,
      face = "bold",
      size = 9.3,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 7.4,
      lineheight = 0.94,
      margin = margin(r = 1.5, unit = "mm")
    ),
    legend.position = "right",
    legend.box = "vertical",
    legend.box.just = "top",
    legend.box.margin = margin(t = 4, r = 1, b = 0, l = 2, unit = "mm"),
    legend.margin = margin(t = 2, r = 0, b = 0, l = 0, unit = "mm"),
    legend.spacing.y = unit(3, "mm")
  ) +
  guides(
    size = guide_legend(
      order = 1,
      title.position = "top",
      title.hjust = 0.5,
      keyheight = unit(4.5, "mm"),
      keywidth = unit(4.5, "mm"),
      override.aes = list(alpha = 1)
    ),
    colour = guide_colorbar(
      order = 2,
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(22, "mm"),
      barwidth = unit(4.2, "mm")
    )
  )

# We want visible strip labels like the reference.
pF <- pF +
  theme(
    strip.background = element_blank(),
    strip.text.y.left = element_text(
      angle = 90,
      face = "bold",
      size = 8.5,
      color = "black"
    )
  )

# ------------------------------------------------------------------------------
# 6. SAVE
# ------------------------------------------------------------------------------
stem <- "Figure2F_GO_KEGG_dotplot_v2_GLIA_readability"

message(
  "Rendering with font family: ", BASE_FAMILY,
  "; Cairo support: ", capabilities("cairo")
)

if (!capabilities("cairo")) {
  stop(
    "This R installation does not report Cairo support. ",
    "Figure 2F requires cairo_pdf() for reliable Arial rendering."
  )
}

# PDF
grDevices::cairo_pdf(
  filename = file.path(out_dir, paste0(stem, ".pdf")),
  width = panel_width_mm / 25.4,
  height = panel_height_mm / 25.4,
  family = BASE_FAMILY,
  onefile = FALSE,
  bg = "white"
)
print(pF)
dev.off()

# SVG
svglite::svglite(
  file = file.path(out_dir, paste0(stem, ".svg")),
  width = panel_width_mm / 25.4,
  height = panel_height_mm / 25.4,
  bg = "white"
)
print(pF)
dev.off()

# TIFF
ragg::agg_tiff(
  filename = file.path(out_dir, paste0(stem, ".tiff")),
  width = panel_width_mm,
  height = panel_height_mm,
  units = "mm",
  res = 600,
  compression = "lzw",
  background = "white"
)
print(pF)
dev.off()

# PNG
ragg::agg_png(
  filename = file.path(out_dir, paste0(stem, ".png")),
  width = panel_width_mm,
  height = panel_height_mm,
  units = "mm",
  res = 300,
  background = "white"
)
print(pF)
dev.off()

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "sessionInfo.txt")
)

message("")
message("============================================================")
message("Figure 2F GO/KEGG dot plot v2 GLIA readability COMPLETE")
message("Primary review file:")
message(file.path(out_dir, paste0(stem, ".pdf")))
message("============================================================")
