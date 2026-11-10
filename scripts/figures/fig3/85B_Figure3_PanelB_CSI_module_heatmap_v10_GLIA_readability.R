# ==============================================================================
# 85B_Figure3_PanelB_CSI_module_heatmap_v10_GLIA_readability.R
# v10 visual-only revision for final 180-mm Figure 3:
# enlarge right annotation/table and CSI key text without changing frozen CSI/modules.
#
# Figure 3B | CSI regulon-module landscape (compact dendrogram, purple palette)
#
# WHY THIS VERSION
# ----------------
# v2 was scientifically correct but visually too different from the reference:
# it lacked the dendrogram/annotation-strip language and the large right-hand
# table dominated the panel.
#
# v3 therefore follows the reference Fig. 5E visual grammar much more closely:
#   1) square CSI heatmap
#   2) frozen hierarchical dendrogram
#   3) thin module annotation strips
#   4) red boxes around frozen module blocks
#   5) module labels around the diagonal
#   6) compact three-column annotation at right
#   7) horizontal CSI color key immediately below the matrix
#
# IMPORTANT SCIENTIFIC LIMIT
# --------------------------
# The reference paper contains many CSI modules, whereas our Step7E analysis
# deliberately froze k = 3 modules. Therefore the NUMBER/SIZE of blocks must
# remain different from the reference. This script imitates the visual language,
# not the reference paper's biological module count.
#
# NO REANALYSIS
# -------------
# Uses only frozen Step7E:
#   CSI
#   CSI_hclust
#   CSI_module_assignment
#   CSI_module_summary
#
# Historical labels are changed only for display:
#   Arm_2 -> Target program
#   Arm_1 -> Alternative branch
#
# OUTPUT
# ------
# <PROJECT_ROOT>/09_final_figures/GLIA/Figure3B_v9/
# ==============================================================================

rm(list = ls())
gc()
options(stringsAsFactors = FALSE, scipen = 999)

# ------------------------------------------------------------------------------
# 0. USER SETTINGS
# ------------------------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
OUT_DIR <- file.path(
  PROJECT_ROOT,
  "09_final_figures", "GLIA", "Figure3B_v10_GLIA_readability"
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FONT_FAMILY <- "Arial"

# Panel dimensions
FIG_WIDTH_MM  <- 176
FIG_HEIGHT_MM <- 88

# Colors
COL_TEXT <- "#222222"
COL_BOX  <- "#D94B35"

MODULE_COLORS <- c(
  "M1" = "#46A7C6",
  "M2" = "#E39A28",
  "M3" = "#35A56C"
)

STATE_DISPLAY <- c(
  "Shared_prebranch" = "Shared prebranch",
  "Arm_2" = "Target program",
  "Arm_1" = "Alternative branch"
)

STATE_COLORS <- c(
  "Shared prebranch" = "#7A6AA6",
  "Target program" = "#B64B45",
  "Alternative branch" = "#7F8C8D"
)

# CSI palette: low -> high; deliberately close to the reference
# CSI palette sampled directly from the user-provided reference Fig. 5E
# horizontal color key (low -> high). These anchors preserve the reference's
# cream -> gray-lavender -> deep blue-purple appearance.
CSI_COLORS <- grDevices::colorRampPalette(
  c(
    "#F8F4DC",
    "#E3DED0",
    "#C8C3C4",
    "#AEA9B6",
    "#928CA7",
    "#787399",
    "#5B5A8C",
    "#3E437E",
    "#203173"
  )
)(256)

# ------------------------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------------------------
required_pkgs <- c("data.table", "svglite", "ragg")

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0L) {
  stop(
    "Missing required R package(s): ",
    paste(missing_pkgs, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
})

# ------------------------------------------------------------------------------
# 2. FROZEN INPUTS
# ------------------------------------------------------------------------------
STEP7E_RDS <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step7E_GSE162610_RSS_CSI_regulon_specificity_modules",
  "RDS",
  "GSE162610_step7E_RSS_CSI_results.rds"
)

STEP7E_MODULE_SUMMARY <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step7E_GSE162610_RSS_CSI_regulon_specificity_modules",
  "tables",
  "11_CSI_module_summary.tsv"
)

input_audit <- data.frame(
  file = c("Step7E RDS", "Step7E module summary"),
  path = c(STEP7E_RDS, STEP7E_MODULE_SUMMARY),
  exists = c(
    file.exists(STEP7E_RDS),
    file.exists(STEP7E_MODULE_SUMMARY)
  ),
  size_bytes = c(
    if (file.exists(STEP7E_RDS)) file.info(STEP7E_RDS)$size else NA_real_,
    if (file.exists(STEP7E_MODULE_SUMMARY)) file.info(STEP7E_MODULE_SUMMARY)$size else NA_real_
  ),
  stringsAsFactors = FALSE
)

write.csv(
  input_audit,
  file.path(OUT_DIR, "00_input_audit.csv"),
  row.names = FALSE
)

if (!all(input_audit$exists)) {
  stop("Figure 3B input missing. Inspect 00_input_audit.csv.", call. = FALSE)
}

step7e <- readRDS(STEP7E_RDS)
module_summary <- fread(STEP7E_MODULE_SUMMARY, data.table = FALSE)

# ------------------------------------------------------------------------------
# 3. STRICT FROZEN SCHEMA
# ------------------------------------------------------------------------------
required_elements <- c(
  "CSI",
  "CSI_hclust",
  "CSI_module_assignment"
)

schema_audit <- data.frame(
  source = "GSE162610_step7E_RSS_CSI_results.rds",
  item = required_elements,
  present = required_elements %in% names(step7e),
  stringsAsFactors = FALSE
)

write.csv(
  schema_audit,
  file.path(OUT_DIR, "01_schema_audit.csv"),
  row.names = FALSE
)

if (!all(schema_audit$present)) {
  stop(
    "Frozen Step7E RDS schema does not match expected save block. ",
    "Inspect 01_schema_audit.csv.",
    call. = FALSE
  )
}

CSI <- as.matrix(step7e$CSI)
HC  <- step7e$CSI_hclust
ASSIGN <- as.data.frame(
  step7e$CSI_module_assignment,
  stringsAsFactors = FALSE
)

if (!inherits(HC, "hclust")) {
  stop("CSI_hclust is not an hclust object.", call. = FALSE)
}

# Detect exact assignment columns robustly
reg_col <- intersect(
  c("regulon", "Regulon", "tf_regulon"),
  colnames(ASSIGN)
)[1]

mod_col <- intersect(
  c("CSI_module", "module", "Module"),
  colnames(ASSIGN)
)[1]

if (is.na(reg_col) || is.na(mod_col)) {
  stop(
    "Cannot identify regulon/module columns in CSI_module_assignment.",
    call. = FALSE
  )
}

ASSIGN$regulon_use <- as.character(ASSIGN[[reg_col]])
ASSIGN$module_use  <- as.character(ASSIGN[[mod_col]])

# ------------------------------------------------------------------------------
# 4. EXACT FROZEN DENDROGRAM ORDER
# ------------------------------------------------------------------------------
hc_labels <- HC$labels
if (is.null(hc_labels)) {
  hc_labels <- rownames(CSI)
}

reg_order <- hc_labels[HC$order]
reg_order <- reg_order[
  reg_order %in% rownames(CSI) &
    reg_order %in% ASSIGN$regulon_use
]

if (length(reg_order) < 2L) {
  stop("Cannot recover frozen CSI_hclust regulon order.", call. = FALSE)
}

CSI_ord <- CSI[
  reg_order,
  reg_order,
  drop = FALSE
]

module_vec <- ASSIGN$module_use[
  match(reg_order, ASSIGN$regulon_use)
]

if (anyNA(module_vec)) {
  stop("Frozen module assignment incomplete after ordering.", call. = FALSE)
}

N <- length(reg_order)

order_audit <- data.frame(
  order_index = seq_len(N),
  regulon = reg_order,
  CSI_module = module_vec,
  stringsAsFactors = FALSE
)

write.csv(
  order_audit,
  file.path(OUT_DIR, "02_Figure3B_regulon_order.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 5. CONTIGUOUS MODULE BLOCKS
# ------------------------------------------------------------------------------
runs <- rle(module_vec)

run_end <- cumsum(runs$lengths)
run_start <- c(1L, head(run_end, -1L) + 1L)

boundaries <- data.frame(
  CSI_module = runs$values,
  start = run_start,
  end = run_end,
  size = runs$lengths,
  stringsAsFactors = FALSE
)

# Each frozen module must be one contiguous dendrogram block.
if (anyDuplicated(boundaries$CSI_module)) {
  stop(
    "A frozen CSI module is split into multiple dendrogram blocks. ",
    "Do not redraw until this is audited.",
    call. = FALSE
  )
}

boundaries$mid <- (boundaries$start + boundaries$end) / 2

write.csv(
  boundaries,
  file.path(OUT_DIR, "03_Figure3B_module_boundaries.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. MODULE SUMMARY
# ------------------------------------------------------------------------------
required_summary <- c(
  "CSI_module",
  "module_size",
  "preferred_state_by_module_RSS",
  "preferred_module_RSS",
  "key_TFs"
)

missing_summary <- setdiff(required_summary, colnames(module_summary))
if (length(missing_summary) > 0L) {
  stop(
    "Missing module summary field(s): ",
    paste(missing_summary, collapse = ", "),
    call. = FALSE
  )
}

summary_used <- module_summary[
  module_summary$CSI_module %in% unique(module_vec),
  required_summary,
  drop = FALSE
]

# order M1 M2 M3 for annotation
summary_used <- summary_used[
  order(match(summary_used$CSI_module, c("M1", "M2", "M3"))),
  ,
  drop = FALSE
]

summary_used$state_display <- unname(
  STATE_DISPLAY[
    as.character(summary_used$preferred_state_by_module_RSS)
  ]
)

short_tfs <- function(x, n = 4L) {
  vapply(as.character(x), function(z) {
    p <- trimws(strsplit(z, "|", fixed = TRUE)[[1]])
    p <- p[nzchar(p)]
    paste(head(p, n), collapse = ", ")
  }, FUN.VALUE = character(1))
}

summary_used$key_TFs_short <- short_tfs(summary_used$key_TFs, 4L)

write.csv(
  summary_used,
  file.path(OUT_DIR, "04_Figure3B_module_summary_used.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. DRAWING HELPERS
# ------------------------------------------------------------------------------
module_col_vec <- unname(MODULE_COLORS[module_vec])
if (anyNA(module_col_vec)) {
  stop("Missing module color for one or more frozen modules.", call. = FALSE)
}

# Base graphics image() expects z[x, y]; construct a matrix so displayed rows
# preserve the same top-to-bottom frozen order.
z_heat <- t(CSI_ord[N:1, , drop = FALSE])

# Color mapping helper
map_to_col_index <- function(x, ncol = length(CSI_COLORS)) {
  x <- pmin(pmax(x, 0), 1)
  pmax(1L, pmin(ncol, floor(x * (ncol - 1)) + 1L))
}

# ------------------------------------------------------------------------------
# 8. MAIN DRAW FUNCTION
# ------------------------------------------------------------------------------
draw_figure3B <- function() {

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  # Layout closely follows reference:
  #
  #             top dendrogram
  #     left    square heatmap       compact annotation
  #     dendro
  #             CSI color key
  #
  # Reserve a real OUTER top margin for the panel label.
  # This is safer than putting B inside a tiny layout cell and prevents clipping.
  par(oma = c(0, 0, 1.8, 0))

  layout(
    matrix(
      c(
        0,  1,  0,
        2,  3,  4,
        0,  5,  0
      ),
      nrow = 3,
      byrow = TRUE
    ),
    # compact left dendrogram | heatmap | right annotation
    widths = c(0.20, 4.80, 4.00),
    # make the CSI key row thinner so it sits immediately under the heatmap
    heights = c(0.20, 5.05, 0.34)
  )

  # ---------------------------------------------------------------------------
  # 8.1 TOP DENDROGRAM
  # ---------------------------------------------------------------------------
  par(
    mar = c(0, 0.20, 0.02, 0.02),
    family = FONT_FAMILY,
    xaxs = "i",
    yaxs = "i",
    lwd = 0.28
  )

  plot(
    as.dendrogram(HC),
    leaflab = "none",
    axes = FALSE,
    xaxs = "i",
    yaxs = "i",
    horiz = FALSE
  )

  # ---------------------------------------------------------------------------
  # 8.2 LEFT DENDROGRAM
  # ---------------------------------------------------------------------------
  par(
    mar = c(0.20, 0.02, 0.02, 0),
    family = FONT_FAMILY,
    xaxs = "i",
    yaxs = "i",
    lwd = 0.28
  )

  plot(
    as.dendrogram(HC),
    leaflab = "none",
    axes = FALSE,
    xaxs = "i",
    yaxs = "i",
    horiz = TRUE
  )

  # ---------------------------------------------------------------------------
  # 8.3 HEATMAP
  # ---------------------------------------------------------------------------
  par(
    mar = c(0.06, 0.55, 0.05, 0.05),
    family = FONT_FAMILY,
    xaxs = "i",
    yaxs = "i"
  )

  image(
    x = seq_len(N),
    y = seq_len(N),
    z = z_heat,
    col = CSI_COLORS,
    zlim = c(0, 1),
    axes = FALSE,
    xlab = "",
    ylab = "",
    useRaster = TRUE
  )

  # Top module annotation strip
  strip_h <- max(0.9, N * 0.008)

  for (i in seq_len(N)) {
    rect(
      xleft = i - 0.5,
      ybottom = N + 0.5,
      xright = i + 0.5,
      ytop = N + 0.5 + strip_h,
      col = module_col_vec[i],
      border = NA,
      xpd = NA
    )
  }

  # Left module annotation strip
  strip_w <- max(0.9, N * 0.008)

  for (i in seq_len(N)) {
    # top matrix row i is display y = N-i+1
    yy <- N - i + 1
    rect(
      xleft = 0.5 - strip_w,
      ybottom = yy - 0.5,
      xright = 0.5,
      ytop = yy + 0.5,
      col = module_col_vec[i],
      border = NA,
      xpd = NA
    )
  }

  # Thin module boxes along diagonal
  for (j in seq_len(nrow(boundaries))) {

    b <- boundaries[j, ]

    # x coordinates follow original order
    x1 <- b$start - 0.5
    x2 <- b$end + 0.5

    # convert original top-down row positions into image y coordinates
    y1 <- N - b$end + 0.5
    y2 <- N - b$start + 1.5

    rect(
      x1, y1, x2, y2,
      border = COL_BOX,
      lwd = 0.75
    )

    # Module label near upper-right edge of each block
    text(
      x = min(x2 - 1.2, b$end),
      y = y2 - 1.6,
      labels = b$CSI_module,
      col = MODULE_COLORS[[b$CSI_module]],
      cex = 0.82,
      font = 2,
      adj = c(1, 1),
      xpd = NA
    )
  }

  box(col = "#B7B7B7", lwd = 0.45)

  # ---------------------------------------------------------------------------
  # 8.4 RIGHT COMPACT ANNOTATION
  # ---------------------------------------------------------------------------
  par(
    mar = c(0.55, 0.42, 0.30, 0.15),
    family = FONT_FAMILY,
    xaxs = "i",
    yaxs = "i"
  )

  plot.new()
  plot.window(
    xlim = c(0, 1),
    ylim = c(0, 1)
  )

  # Reference-like column headers and gray rules
  headers <- c("Module", "Key TFs", "Preferred state")
  hx <- c(0.02, 0.22, 0.69)

  for (i in seq_along(headers)) {
    text(
      hx[i], 0.91,
      headers[i],
      adj = c(0, 0.5),
      col = "#999999",
      cex = 1.10,
      font = 2
    )
  }

  segments(0.02, 0.865, 0.18, 0.865, col = "#9D9D9D", lwd = 1.0)
  segments(0.23, 0.865, 0.65, 0.865, col = "#9D9D9D", lwd = 1.0)
  segments(0.70, 0.865, 0.99, 0.865, col = "#9D9D9D", lwd = 1.0)

  y_rows <- seq(0.72, 0.28, length.out = nrow(summary_used))

  for (i in seq_len(nrow(summary_used))) {

    mm <- summary_used$CSI_module[i]
    st <- summary_used$state_display[i]

    text(
      0.02, y_rows[i],
      mm,
      adj = c(0, 0.5),
      col = MODULE_COLORS[[mm]],
      cex = 1.10,
      font = 2
    )

    text(
      0.23, y_rows[i],
      summary_used$key_TFs_short[i],
      adj = c(0, 0.5),
      col = "#5A3F86",
      cex = 1.04
    )

    text(
      0.70, y_rows[i],
      st,
      adj = c(0, 0.5),
      col = STATE_COLORS[[st]],
      cex = 1.04
    )
  }

  # ---------------------------------------------------------------------------
  # 8.5 HORIZONTAL CSI COLOR KEY
  # ---------------------------------------------------------------------------
  par(
    mar = c(0.02, 0.55, 0.02, 0.05),
    family = FONT_FAMILY,
    xaxs = "i",
    yaxs = "i"
  )

  plot.new()
  plot.window(
    xlim = c(0, 1),
    ylim = c(0, 1)
  )

  vals <- seq(0, 1, length.out = length(CSI_COLORS))

  for (i in seq_len(length(CSI_COLORS) - 1L)) {
    rect(
      vals[i], 0.66,
      vals[i + 1], 0.92,
      col = CSI_COLORS[i],
      border = NA
    )
  }

  text(
    0.00, 0.38,
    "low",
    adj = c(0, 0.5),
    cex = 0.82,
    col = COL_TEXT
  )

  text(
    1.00, 0.38,
    "high",
    adj = c(1, 0.5),
    cex = 0.82,
    col = COL_TEXT
  )

  text(
    0.50, 0.38,
    "Connection specificity index (CSI)",
    adj = c(0.5, 0.5),
    cex = 0.88,
    font = 2,
    col = COL_TEXT
  )



  # ---------------------------------------------------------------------------
  # 8.6 PANEL LABEL "B" -- TRUE OUTER MARGIN
  # ---------------------------------------------------------------------------
  # Draw B in the reserved outer top margin after every layout panel is drawn.
  # A generous top outer margin plus line=0.15 keeps the full glyph inside
  # PDF/SVG/TIFF/PNG device bounds.
  mtext(
    "B",
    side = 3,
    outer = TRUE,
    at = 0.012,
    adj = 0,
    line = 0.15,
    cex = 1.25,
    font = 2,
    family = FONT_FAMILY,
    col = COL_TEXT
  )


}

# ------------------------------------------------------------------------------
# 9. EXPORT DEVICES
# ------------------------------------------------------------------------------
WIDTH_IN  <- FIG_WIDTH_MM / 25.4
HEIGHT_IN <- FIG_HEIGHT_MM / 25.4

PDF_FILE <- file.path(
  OUT_DIR,
  "Figure3B_CSI_module_heatmap_v10_GLIA_readability.pdf"
)

SVG_FILE <- file.path(
  OUT_DIR,
  "Figure3B_CSI_module_heatmap_v10_GLIA_readability.svg"
)

TIFF_FILE <- file.path(
  OUT_DIR,
  "Figure3B_CSI_module_heatmap_v10_GLIA_readability.tiff"
)

PNG_FILE <- file.path(
  OUT_DIR,
  "Figure3B_CSI_module_heatmap_v10_GLIA_readability.png"
)

# PDF
grDevices::cairo_pdf(
  filename = PDF_FILE,
  width = WIDTH_IN,
  height = HEIGHT_IN,
  family = FONT_FAMILY,
  onefile = TRUE
)
draw_figure3B()
grDevices::dev.off()

# SVG
svglite::svglite(
  file = SVG_FILE,
  width = WIDTH_IN,
  height = HEIGHT_IN,
  system_fonts = list(sans = FONT_FAMILY)
)
draw_figure3B()
grDevices::dev.off()

# TIFF 600 dpi
ragg::agg_tiff(
  filename = TIFF_FILE,
  width = WIDTH_IN,
  height = HEIGHT_IN,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)
draw_figure3B()
grDevices::dev.off()

# PNG 300 dpi
ragg::agg_png(
  filename = PNG_FILE,
  width = WIDTH_IN,
  height = HEIGHT_IN,
  units = "in",
  res = 300,
  background = "white"
)
draw_figure3B()
grDevices::dev.off()

# ------------------------------------------------------------------------------
# 10. MD5 / SESSION
# ------------------------------------------------------------------------------
source_md5 <- data.frame(
  file = c(STEP7E_RDS, STEP7E_MODULE_SUMMARY),
  md5 = unname(
    tools::md5sum(
      c(STEP7E_RDS, STEP7E_MODULE_SUMMARY)
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  source_md5,
  file.path(OUT_DIR, "source_file_md5.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_DIR, "sessionInfo.txt")
)

# ------------------------------------------------------------------------------
# 11. CONSOLE SUMMARY
# ------------------------------------------------------------------------------
cat("\n============================================================\n")
cat("Figure 3B v10 GLIA readability completed.\n")
cat("Output directory:\n", OUT_DIR, "\n", sep = "")
cat("\nVisual design:\n")
cat("- square CSI matrix with reference-sampled CSI palette\n")
cat("- top + left frozen dendrograms reduced to reference-like depth (~4% of heatmap)\n")
cat("- panel label B is drawn in a true outer top margin, fully visible and non-overlapping\n")
cat("- module annotation strips top + left\n")
cat("- red diagonal module boxes\n")
cat("- compact reference-like annotation at right\n")
cat("- CSI color key moved tightly against the heatmap\n")
cat("- note: two modules can legitimately share the same preferred state\n")
cat("\nScientific constraint:\n")
cat("- M1 and M3 both display as Alternative branch because the frozen Step7E module summary assigns both modules that same preferred state; this is biological/frozen output, not a plotting bug.\n")
cat("- Our frozen Step7E has k=3 modules; reference paper has many modules.\n")
cat("- Therefore block number/size should NOT be forced to match the reference.\n")
cat("- No re-clustering performed.\n")
cat("============================================================\n")
