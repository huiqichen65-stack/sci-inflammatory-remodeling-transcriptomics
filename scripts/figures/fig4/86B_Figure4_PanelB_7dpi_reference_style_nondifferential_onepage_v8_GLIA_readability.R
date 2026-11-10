# ============================================================
# v8 GLIA readability revision: enlarged labels/axes/legend only; no CellChat re-analysis.
# Figure 4B v7 | MANUAL GRID LAYOUT — guaranteed no legend overlap
# 7 dpi, NON-differential CellChat network
#
# WHY THIS VERSION:
# Previous ComplexHeatmap layout kept allowing the legend to intrude
# into the right heatmap and created excessive left whitespace.
# This version does NOT let ComplexHeatmap control the page geometry.
# Every component is drawn at fixed page coordinates with base grid.
#
# Left  : Number of interactions
# Right : Interaction strength
# Heatmap colors = within-7dpi relative values
#                  z-score(log1p(raw network value))
# Top bars        = incoming totals (column sums)
# Right bars      = outgoing totals (row sums)
#
# One PDF page only.
# No CellChat re-analysis.
# ============================================================

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# 0. Package
# ------------------------------------------------------------
if (!requireNamespace("CellChat", quietly = TRUE)) {
  stop("CellChat is required to read the frozen CellChat RDS object.")
}

suppressPackageStartupMessages({
  library(CellChat)
  library(grid)
})

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------
ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
STEP8B_PREVIEW <- file.path(
  ROOT,
  "04_results",
  "step8B_PREVIEW_GSE162610_CellChat_nboot20"
)

RDS_PRIMARY <- file.path(
  STEP8B_PREVIEW,
  "RDS",
  "timepoint_checkpoints",
  "GSE162610_step8B_PREVIEW_CellChat_7dpi_checkpoint.rds"
)

RDS_FALLBACK1 <- file.path(
  STEP8B_PREVIEW,
  "RDS",
  "GSE162610_step8B_PREVIEW_CellChat_results.rds"
)

RDS_FALLBACK2 <- file.path(
  STEP8B_PREVIEW,
  "RDS",
  "GSE162610_step8B_merged_CellChat_1_3_7dpi.rds"
)

OUTDIR <- file.path(
  ROOT,
  "09_final_figures",
  "GLIA",
  "Figure4",
  "single_panels",
  "Figure4B"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

PDF_FILE <- file.path(
  OUTDIR,
  "Figure4B_7dpi_number_strength_reference_style_onepage_v8_GLIA_readability.pdf"
)

if (file.exists(PDF_FILE)) unlink(PDF_FILE)

# ------------------------------------------------------------
# 2. Read frozen CellChat object
# ------------------------------------------------------------
find_CellChat_object <- function(x) {
  if (methods::is(x, "CellChat")) return(x)

  if (is.list(x)) {
    for (i in seq_along(x)) {
      ans <- try(find_CellChat_object(x[[i]]), silent = TRUE)
      if (
        !inherits(ans, "try-error") &&
        !is.null(ans) &&
        methods::is(ans, "CellChat")
      ) {
        return(ans)
      }
    }
  }

  NULL
}

load_CellChat_from_candidates <- function(paths) {
  existing <- paths[file.exists(paths)]

  if (length(existing) == 0) {
    stop(
      "No frozen Step8B CellChat RDS was found.\nChecked:\n",
      paste(paths, collapse = "\n")
    )
  }

  for (f in existing) {
    message("Trying frozen CellChat RDS: ", f)
    obj <- readRDS(f)
    cc <- find_CellChat_object(obj)

    if (!is.null(cc)) {
      message("CellChat object found in: ", f)
      return(cc)
    }
  }

  stop("RDS files exist, but no CellChat object was found.")
}

cellchat7 <- load_CellChat_from_candidates(
  c(RDS_PRIMARY, RDS_FALLBACK1, RDS_FALLBACK2)
)

if (is.null(cellchat7@net$count) || is.null(cellchat7@net$weight)) {
  stop("Frozen CellChat object lacks @net$count and/or @net$weight.")
}

count_raw  <- as.matrix(cellchat7@net$count)
weight_raw <- as.matrix(cellchat7@net$weight)

# ------------------------------------------------------------
# 3. Same 10 groups used in Figure 4A/B
# ------------------------------------------------------------
raw_groups <- c(
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

display_map <- c(
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

missing_groups <- setdiff(
  raw_groups,
  Reduce(
    intersect,
    list(
      rownames(count_raw),
      colnames(count_raw),
      rownames(weight_raw),
      colnames(weight_raw)
    )
  )
)

if (length(missing_groups) > 0) {
  stop(
    "Missing required group(s): ",
    paste(missing_groups, collapse = ", ")
  )
}

count_mat <- count_raw[raw_groups, raw_groups, drop = FALSE]
weight_mat <- weight_raw[raw_groups, raw_groups, drop = FALSE]

display_groups <- unname(display_map[raw_groups])

rownames(count_mat) <- colnames(count_mat) <- display_groups
rownames(weight_mat) <- colnames(weight_mat) <- display_groups

# ------------------------------------------------------------
# 4. Relative values used ONLY for heatmap color
# ------------------------------------------------------------
Z_CAP <- 1.8

relative_matrix <- function(mat, cap = Z_CAP) {
  x <- log1p(as.numeric(mat))
  mu <- mean(x, na.rm = TRUE)
  sig <- stats::sd(x, na.rm = TRUE)

  if (!is.finite(sig) || sig == 0) {
    z <- rep(0, length(x))
  } else {
    z <- (x - mu) / sig
  }

  z <- pmax(-cap, pmin(cap, z))

  matrix(
    z,
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dimnames(mat)
  )
}

count_rel  <- relative_matrix(count_mat)
weight_rel <- relative_matrix(weight_mat)

# ------------------------------------------------------------
# 5. Marginal totals
# ------------------------------------------------------------
count_in  <- colSums(count_mat, na.rm = TRUE)
count_out <- rowSums(count_mat, na.rm = TRUE)

weight_in  <- colSums(weight_mat, na.rm = TRUE)
weight_out <- rowSums(weight_mat, na.rm = TRUE)

# ------------------------------------------------------------
# 6. Colors
# ------------------------------------------------------------
group_cols <- c(
  "Target program"  = "#B13A36",
  "Other myeloid"   = "#D95F02",
  "Microglia other" = "#7663B0",
  "Neutrophil"      = "#D9A300",
  "T cell"          = "#149A99",
  "Astrocyte"       = "#D81B80",
  "Endothelial"     = "#159570",
  "Fibroblast"      = "#6B61B8",
  "Pericyte/VSMC"   = "#5E9E1B",
  "OPC"             = "#3C6EA8"
)

cell_cols <- unname(group_cols[display_groups])

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
# 7. Red-white-blue interpolation
# ------------------------------------------------------------
hex_to_rgb <- function(hex) {
  grDevices::col2rgb(hex) / 255
}

interp_col <- function(v, vmin = -Z_CAP, vmid = 0, vmax = Z_CAP) {
  blue  <- hex_to_rgb("#2166AC")
  white <- hex_to_rgb("#F7F7F7")
  red   <- hex_to_rgb("#B2182B")

  v <- max(vmin, min(vmax, v))

  if (v <= vmid) {
    t <- (v - vmin) / (vmid - vmin)
    rgbv <- blue * (1 - t) + white * t
  } else {
    t <- (v - vmid) / (vmax - vmid)
    rgbv <- white * (1 - t) + red * t
  }

  grDevices::rgb(rgbv[1], rgbv[2], rgbv[3])
}

# ------------------------------------------------------------
# 8. Nice marginal-axis limits
# ------------------------------------------------------------
nice_max <- function(x) {
  m <- max(x, na.rm = TRUE)
  if (!is.finite(m) || m <= 0) return(1)
  pretty(c(0, m), n = 3)[length(pretty(c(0, m), n = 3))]
}

# ------------------------------------------------------------
# 9. Draw one reference-style panel at FIXED coordinates
# ------------------------------------------------------------
draw_panel <- function(
  rel_mat,
  raw_mat,
  incoming,
  outgoing,
  x0,
  y0,
  hm_w,
  hm_h,
  title,
  show_row_labels = FALSE,
  show_y_title = FALSE
) {

  nr <- nrow(rel_mat)
  nc <- ncol(rel_mat)

  cw <- hm_w / nc
  ch <- hm_h / nr

  top_gap <- 0.055
  top_h   <- 0.28

  right_gap <- 0.065
  right_w   <- 0.30

  # -------------------------
  # Heatmap cells
  # -------------------------
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {

      cx <- x0 + (j - 0.5) * cw
      cy <- y0 + hm_h - (i - 0.5) * ch

      grid.rect(
        x = unit(cx, "in"),
        y = unit(cy, "in"),
        width = unit(cw, "in"),
        height = unit(ch, "in"),
        gp = gpar(
          fill = interp_col(rel_mat[i, j]),
          col = "white",
          lwd = 0.55
        )
      )
    }
  }

  # -------------------------
  # Column labels
  # -------------------------
  for (j in seq_len(nc)) {
    grid.text(
      colnames(rel_mat)[j],
      x = unit(x0 + (j - 0.5) * cw, "in"),
      y = unit(y0 - 0.055, "in"),
      rot = 90,
      just = c("right", "center"),
      gp = gpar(
        fontsize = 8.6,
        fontfamily = FONT
      )
    )
  }

  # -------------------------
  # Row labels only on left panel
  # -------------------------
  if (show_row_labels) {
    for (i in seq_len(nr)) {
      grid.text(
        rownames(rel_mat)[i],
        x = unit(x0 - 0.085, "in"),
        y = unit(y0 + hm_h - (i - 0.5) * ch, "in"),
        just = c("right", "center"),
        gp = gpar(
          fontsize = 10.0,
          fontfamily = FONT
        )
      )
    }
  }

  # -------------------------
  # Top incoming bars
  # -------------------------
  in_max <- nice_max(incoming)

  # baseline
  grid.lines(
    x = unit(c(x0, x0 + hm_w), "in"),
    y = unit(c(y0 + hm_h + top_gap, y0 + hm_h + top_gap), "in"),
    gp = gpar(col = "#555555", lwd = 0.45)
  )

  for (j in seq_len(nc)) {

    bh <- if (in_max > 0) top_h * incoming[j] / in_max else 0

    grid.rect(
      x = unit(x0 + (j - 0.5) * cw, "in"),
      y = unit(y0 + hm_h + top_gap + bh / 2, "in"),
      width = unit(cw * 0.68, "in"),
      height = unit(max(bh, 0.002), "in"),
      gp = gpar(
        fill = cell_cols[j],
        col = NA
      )
    )
  }

  # top marginal y axis
  ax_x <- x0 - 0.04
  grid.lines(
    x = unit(c(ax_x, ax_x), "in"),
    y = unit(
      c(y0 + hm_h + top_gap, y0 + hm_h + top_gap + top_h),
      "in"
    ),
    gp = gpar(col = "#555555", lwd = 0.45)
  )

  ticks <- c(0, in_max / 2, in_max)
  for (tv in ticks) {
    yy <- y0 + hm_h + top_gap + top_h * tv / in_max

    grid.lines(
      x = unit(c(ax_x - 0.025, ax_x), "in"),
      y = unit(c(yy, yy), "in"),
      gp = gpar(col = "#555555", lwd = 0.4)
    )

    grid.text(
      format(round(tv, 1), trim = TRUE),
      x = unit(ax_x - 0.04, "in"),
      y = unit(yy, "in"),
      just = c("right", "center"),
      gp = gpar(
        fontsize = 5.4,
        fontfamily = FONT,
        col = "#444444"
      )
    )
  }

  # -------------------------
  # Right outgoing bars
  # -------------------------
  out_max <- nice_max(outgoing)

  bar_x0 <- x0 + hm_w + right_gap

  for (i in seq_len(nr)) {

    bw <- if (out_max > 0) right_w * outgoing[i] / out_max else 0
    cy <- y0 + hm_h - (i - 0.5) * ch

    grid.rect(
      x = unit(bar_x0 + bw / 2, "in"),
      y = unit(cy, "in"),
      width = unit(max(bw, 0.002), "in"),
      height = unit(ch * 0.58, "in"),
      gp = gpar(
        fill = cell_cols[i],
        col = NA
      )
    )
  }

  # right marginal x-axis
  axis_y <- y0 - 0.03

  grid.lines(
    x = unit(c(bar_x0, bar_x0 + right_w), "in"),
    y = unit(c(axis_y, axis_y), "in"),
    gp = gpar(col = "#555555", lwd = 0.4)
  )

  rticks <- c(0, out_max / 2, out_max)

  for (tv in rticks) {
    xx <- bar_x0 + right_w * tv / out_max

    grid.lines(
      x = unit(c(xx, xx), "in"),
      y = unit(c(axis_y, axis_y - 0.025), "in"),
      gp = gpar(col = "#555555", lwd = 0.4)
    )

    grid.text(
      format(round(tv, 1), trim = TRUE),
      x = unit(xx, "in"),
      y = unit(axis_y - 0.045, "in"),
      rot = 90,
      just = c("right", "center"),
      gp = gpar(
        fontsize = 5.1,
        fontfamily = FONT,
        col = "#444444"
      )
    )
  }

  # -------------------------
  # Panel title
  # -------------------------
  grid.text(
    title,
    x = unit(x0 + hm_w / 2, "in"),
    y = unit(y0 + hm_h + top_gap + top_h + 0.16, "in"),
    gp = gpar(
      fontsize = 11.2,
      fontface = "bold",
      fontfamily = FONT
    )
  )

  # -------------------------
  # Y title (left panel only)
  # -------------------------
  if (show_y_title) {
    grid.text(
      "Sources (Sender)",
      x = unit(x0 - 1.22, "in"),
      y = unit(y0 + hm_h / 2, "in"),
      rot = 90,
      gp = gpar(
        fontsize = 9.0,
        fontfamily = FONT
      )
    )
  }
}

# ------------------------------------------------------------
# 10. Draw legend at FIXED coordinates
# ------------------------------------------------------------
draw_relative_legend <- function(x0, y0, w = 0.14, h = 1.10) {

  nseg <- 120
  vals <- seq(-Z_CAP, Z_CAP, length.out = nseg)

  for (k in seq_len(nseg)) {
    yy <- y0 + (k - 1) * h / nseg

    grid.rect(
      x = unit(x0 + w / 2, "in"),
      y = unit(yy + h / (2 * nseg), "in"),
      width = unit(w, "in"),
      height = unit(h / nseg + 0.002, "in"),
      gp = gpar(
        fill = interp_col(vals[k]),
        col = NA
      )
    )
  }

  grid.rect(
    x = unit(x0 + w / 2, "in"),
    y = unit(y0 + h / 2, "in"),
    width = unit(w, "in"),
    height = unit(h, "in"),
    gp = gpar(
      fill = NA,
      col = "#777777",
      lwd = 0.45
    )
  )

  at_vals <- c(-1.8, -0.9, 0, 0.9, 1.8)

  for (v in at_vals) {
    yy <- y0 + h * (v + Z_CAP) / (2 * Z_CAP)

    grid.lines(
      x = unit(c(x0 + w, x0 + w + 0.035), "in"),
      y = unit(c(yy, yy), "in"),
      gp = gpar(col = "#555555", lwd = 0.4)
    )

    grid.text(
      format(v, trim = TRUE),
      x = unit(x0 + w + 0.055, "in"),
      y = unit(yy, "in"),
      just = c("left", "center"),
      gp = gpar(
        fontsize = 6.0,
        fontfamily = FONT
      )
    )
  }

  grid.text(
    "Relative values",
    x = unit(x0 - 0.12, "in"),
    y = unit(y0 + h / 2, "in"),
    rot = 90,
    gp = gpar(
      fontsize = 7.0,
      fontfamily = FONT
    )
  )
}

# ------------------------------------------------------------
# 11. Exact page geometry
#
# KEY FIXES:
#   - left heatmap starts at x = 1.72 in, not in the middle of page
#   - right heatmap starts at x = 4.90 in
#   - legend starts at x = 8.15 in
#   => legend can NEVER touch the right heatmap
# ------------------------------------------------------------
PAGE_W <- 9.00
PAGE_H <- 4.20

HM_W <- 2.28
HM_H <- 2.28

LEFT_X  <- 1.72
RIGHT_X <- 4.90
HM_Y    <- 1.02

LEGEND_X <- 8.18
LEGEND_Y <- 1.58

# ------------------------------------------------------------
# 12. One-page PDF
# ------------------------------------------------------------
if (capabilities("cairo")) {
  grDevices::cairo_pdf(
    filename = PDF_FILE,
    width = PAGE_W,
    height = PAGE_H,
    family = FONT,
    onefile = FALSE,
    bg = "white"
  )
} else {
  grDevices::pdf(
    file = PDF_FILE,
    width = PAGE_W,
    height = PAGE_H,
    family = "Helvetica",
    onefile = FALSE,
    useDingbats = FALSE,
    bg = "white"
  )
}

grid.newpage()

# Panel B — close to the figure, not separated by a huge blank area
grid.text(
  "B",
  x = unit(0.10, "in"),
  y = unit(PAGE_H - 0.10, "in"),
  just = c("left", "top"),
  gp = gpar(
    fontsize = 20,
    fontface = "bold",
    fontfamily = FONT
  )
)

# Left panel
draw_panel(
  rel_mat = count_rel,
  raw_mat = count_mat,
  incoming = count_in,
  outgoing = count_out,
  x0 = LEFT_X,
  y0 = HM_Y,
  hm_w = HM_W,
  hm_h = HM_H,
  title = "Number of interactions",
  show_row_labels = TRUE,
  show_y_title = TRUE
)

# Right panel
draw_panel(
  rel_mat = weight_rel,
  raw_mat = weight_mat,
  incoming = weight_in,
  outgoing = weight_out,
  x0 = RIGHT_X,
  y0 = HM_Y,
  hm_w = HM_W,
  hm_h = HM_H,
  title = "Interaction strength",
  show_row_labels = FALSE,
  show_y_title = FALSE
)

# Legend — completely separate fixed region
draw_relative_legend(
  x0 = LEGEND_X,
  y0 = LEGEND_Y,
  w = 0.14,
  h = 1.10
)

dev.off()

# ------------------------------------------------------------
# 13. Check
# ------------------------------------------------------------
if (
  !file.exists(PDF_FILE) ||
  is.na(file.info(PDF_FILE)$size) ||
  file.info(PDF_FILE)$size < 1000
) {
  stop("PDF was not written correctly.")
}

cat("\n============================================================\n")
cat("Figure 4B v8 GLIA readability revision completed.\n")
cat("Manual fixed-coordinate grid layout.\n")
cat("Legend is in a physically separate right-side region.\n")
cat("Left blank space has been removed substantially.\n")
cat("One PDF page only.\n")
cat("No CellChat analysis was rerun.\n\n")
cat("PDF:\n")
cat("  ", PDF_FILE, "\n", sep = "")
cat("============================================================\n")
