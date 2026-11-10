# ==============================================================================
# 86F_Figure4_PanelF_ONE_FILE_v8_TNF_NFKB_BH.R
#
# Figure 4F | GSE234774 biological-library-level exploratory association
#
# ONE FILE = calculation + correlation + plotting + export
#
# LEFT:
#   Target-program AUCell score in resident microglia
#        vs
#   Endothelial TNF/NF-kB response score
#
# RIGHT:
#   Reference-style triangular correlation matrix:
#   fixed grey circles + red/blue pie sectors proportional to |Spearman R|
#   Significance is determined using BH-adjusted P values across 36 unique pairs.
#
# IMPORTANT
# ---------
# - No CellChat rerun.
# - No reclustering.
# - Target-program library scores are read from the FROZEN Step11B library table:
#     record_type == "LIBRARY_SCORE"
#     method      == "AUCell_library_median"
#     program     == "Arm_2_specific"
#   The historical Arm_2 name is only used to read the frozen result; the plot
#   displays "Target program".
# - Endothelial response programs are scored with the same deterministic
#   gene-wise z-score logic used in Step8C.
# - Step10D broad lineage uses major_lineage_candidate. This is appropriate here
#   only for broad compartments such as Endothelial/Fibroblast/OPC etc.;
#   no fine new biological labels are inferred.
# - Statistical unit = biological library.
# - The original endothelial-activation association remains in the matrix.
# - The 36 unique off-diagonal correlations form one BH correction family.
# - This is exploratory association analysis; it does not establish causality.
# ==============================================================================

# Do not clear the user workspace when sourced interactively.
gc()
options(stringsAsFactors = FALSE, scipen = 999)

# ------------------------------------------------------------------------------
# 0. PROJECT PATHS
# ------------------------------------------------------------------------------
ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
STEP10D_DIR <- file.path(
  ROOT,
  "04_results",
  "step10D_GSE234774_marker_first_major_lineage_annotation"
)

STEP10D_RDS <- file.path(
  STEP10D_DIR,
  "RDS",
  "GSE234774_step10D_major_lineage_annotation_evidence_BPCells.rds"
)

STEP10D_REVIEW_TABLE <- file.path(
  STEP10D_DIR,
  "tables",
  "05_cluster_annotation_review_table.csv"
)

STEP11B_DIR <- file.path(
  ROOT,
  "04_results",
  "step11B_GSE234774_refined_main_Arm_specific_shared_AUCell_ssGSEA"
)

STEP11B_LIBRARY_TABLE <- file.path(
  STEP11B_DIR,
  "tables",
  "05_refined_main_library_and_timepoint_scores.csv"
)

STEP8C_GENESET_TSV <- file.path(
  ROOT,
  "04_results",
  "step8C_GSE162610_recipient_program_integration",
  "tables",
  "02_recipient_program_gene_coverage.tsv"
)

OUTDIR <- file.path(
  ROOT,
  "09_final_figures",
  "GLIA",
  "Figure4",
  "single_panels",
  "Figure4F"
)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

PDF_FILE <- file.path(
  OUTDIR,
  "Figure4F_GSE234774_library_level_target_TNF_NFKB_correlation_v8_BH.pdf"
)
PNG_FILE <- file.path(
  OUTDIR,
  "Figure4F_GSE234774_library_level_target_TNF_NFKB_correlation_v8_BH.png"
)
SUMMARY_CSV <- file.path(
  OUTDIR,
  "Figure4F_v8_library_level_summary.csv"
)
CORR_CSV <- file.path(
  OUTDIR,
  "Figure4F_v8_library_level_correlation_matrix_BH.csv"
)
AUDIT_CSV <- file.path(
  OUTDIR,
  "Figure4F_v8_audit.csv"
)

FONT_FAMILY <- "Arial"

if (.Platform$OS.type == "windows") {
  try(
    grDevices::windowsFonts(
      Arial = grDevices::windowsFont("Arial")
    ),
    silent = TRUE
  )
}

# ------------------------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------------------------
required_pkgs <- c(
  "Seurat",
  "SeuratObject",
  "ggplot2",
  "dplyr",
  "tibble",
  "tidyr",
  "stringr",
  "Matrix",
  "patchwork",
  "data.table"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall with:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
    "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(Matrix)
  library(patchwork)
  library(data.table)
})

# ------------------------------------------------------------------------------
# 2. HELPERS
# ------------------------------------------------------------------------------
find_exact_file <- function(root, filename) {
  hits <- list.files(
    root,
    pattern = paste0("^", gsub("\\.", "\\\\.", filename), "$"),
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(hits) == 0L) return(NA_character_)
  hits[1]
}

resolve_file <- function(primary, filename) {
  if (file.exists(primary)) return(primary)
  hit <- find_exact_file(file.path(ROOT, "04_results"), filename)
  if (!is.na(hit) && file.exists(hit)) return(hit)
  stop("Could not find required file:\n", primary, call. = FALSE)
}

standardize_timepoint <- function(x) {
  z <- tolower(trimws(as.character(x)))
  z <- gsub("_", " ", z)
  z <- gsub("-", " ", z)
  z <- gsub("\\s+", " ", z)

  out <- dplyr::case_when(
    z %in% c("uninjured", "naive", "control", "sham") ~ "Uninjured",
    z %in% c("1 dpi", "1dpi", "d1", "day 1") ~ "1 dpi",
    z %in% c("4 dpi", "4dpi", "d4", "day 4") ~ "4 dpi",
    z %in% c("7 dpi", "7dpi", "d7", "day 7") ~ "7 dpi",
    z %in% c("14 dpi", "14dpi", "d14", "day 14") ~ "14 dpi",
    z %in% c("1 month", "1month", "1 mo", "1m", "30 dpi", "30dpi") ~ "1 month",
    z %in% c("2 months", "2 month", "2months", "2month",
             "2 mo", "2m", "60 dpi", "60dpi") ~ "2 months",
    TRUE ~ as.character(x)
  )

  factor(
    out,
    levels = c(
      "Uninjured", "1 dpi", "4 dpi", "7 dpi",
      "14 dpi", "1 month", "2 months"
    )
  )
}

match_genes <- function(genes, available) {
  mp <- setNames(available, tolower(available))
  out <- unname(mp[tolower(unique(genes))])
  unique(out[!is.na(out)])
}

get_normalized_subset <- function(obj, genes, cells) {
  DefaultAssay(obj) <- "RNA"

  available_genes <- rownames(obj[["RNA"]])
  available_cells <- colnames(obj)

  genes <- match_genes(genes, available_genes)
  cells <- intersect(cells, available_cells)

  if (length(genes) == 0L) {
    stop("No requested program genes were found in the RNA assay.", call. = FALSE)
  }
  if (length(cells) == 0L) {
    stop("No requested cells were found in the atlas object.", call. = FALSE)
  }

  x <- try(
    SeuratObject::LayerData(
      object = obj,
      assay = "RNA",
      layer = "data",
      features = genes,
      cells = cells
    ),
    silent = TRUE
  )

  if (inherits(x, "try-error")) {
    x <- try(
      Seurat::GetAssayData(
        object = obj,
        assay = "RNA",
        layer = "data"
      )[genes, cells, drop = FALSE],
      silent = TRUE
    )
  }

  if (inherits(x, "try-error")) {
    x <- Seurat::GetAssayData(
      object = obj,
      assay = "RNA",
      slot = "data"
    )[genes, cells, drop = FALSE]
  }

  x
}

# Same Step8C logic: deterministic gene-wise z scoring.
program_score <- function(expr_sub) {
  z <- t(scale(t(as.matrix(expr_sub))))
  z[!is.finite(z)] <- 0
  colMeans(z, na.rm = TRUE)
}

library_mean_score <- function(meta, sample_col, cells, score, out_name) {
  x <- tibble(
    cell = cells,
    score = as.numeric(score)
  )

  meta %>%
    rownames_to_column("cell") %>%
    transmute(
      cell = cell,
      biological_library = as.character(.data[[sample_col]])
    ) %>%
    inner_join(x, by = "cell") %>%
    filter(!is.na(biological_library), is.finite(score)) %>%
    group_by(biological_library) %>%
    summarise(
      !!out_name := mean(score, na.rm = TRUE),
      .groups = "drop"
    )
}

library_prop <- function(meta, sample_col, type_col, regex, out_name) {
  meta %>%
    rownames_to_column("cell") %>%
    transmute(
      biological_library = as.character(.data[[sample_col]]),
      broad_type = as.character(.data[[type_col]])
    ) %>%
    filter(!is.na(biological_library), !is.na(broad_type)) %>%
    mutate(hit = grepl(regex, tolower(broad_type))) %>%
    group_by(biological_library) %>%
    summarise(
      !!out_name := mean(hit),
      .groups = "drop"
    )
}

spearman_one <- function(x, y) {
  ok <- complete.cases(x, y) & is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]

  if (length(x) < 4L || sd(x) == 0 || sd(y) == 0) {
    return(list(rho = NA_real_, p = NA_real_, n = length(x)))
  }

  ct <- suppressWarnings(
    stats::cor.test(x, y, method = "spearman", exact = FALSE)
  )

  list(
    rho = unname(ct$estimate),
    p = ct$p.value,
    n = length(x)
  )
}

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return(format(p, scientific = TRUE, digits = 2))
  sprintf("%.4f", p)
}

# ------------------------------------------------------------------------------
# 3. RESOLVE INPUT FILES
# ------------------------------------------------------------------------------
STEP10D_RDS <- resolve_file(
  STEP10D_RDS,
  "GSE234774_step10D_major_lineage_annotation_evidence_BPCells.rds"
)

STEP10D_REVIEW_TABLE <- resolve_file(
  STEP10D_REVIEW_TABLE,
  "05_cluster_annotation_review_table.csv"
)

STEP11B_LIBRARY_TABLE <- resolve_file(
  STEP11B_LIBRARY_TABLE,
  "05_refined_main_library_and_timepoint_scores.csv"
)

STEP8C_GENESET_TSV <- resolve_file(
  STEP8C_GENESET_TSV,
  "02_recipient_program_gene_coverage.tsv"
)

cat("\nResolved inputs:\n")
cat("Step10D atlas       : ", STEP10D_RDS, "\n", sep = "")
cat("Step10D review table: ", STEP10D_REVIEW_TABLE, "\n", sep = "")
cat("Step11B library data: ", STEP11B_LIBRARY_TABLE, "\n", sep = "")
cat("Step8C gene sets    : ", STEP8C_GENESET_TSV, "\n\n", sep = "")

# ------------------------------------------------------------------------------
# 4. LOAD STEP10D ATLAS AND DEFINE BROAD LINEAGE COLUMN
# ------------------------------------------------------------------------------
atlas <- readRDS(STEP10D_RDS)

if (!inherits(atlas, "Seurat")) {
  stop("Step10D RDS is not a Seurat object.", call. = FALSE)
}

meta_atlas <- atlas[[]]

# Use biological_library explicitly whenever available.
sample_candidates <- c(
  "biological_library",
  "library_clean",
  "library",
  "orig.ident"
)

sample_col <- sample_candidates[sample_candidates %in% colnames(meta_atlas)][1]

if (is.na(sample_col) || length(sample_col) == 0L) {
  stop(
    "Could not find a biological-library field in Step10D metadata.",
    call. = FALSE
  )
}

time_candidates <- c(
  "timepoint",
  "timepoint_standard",
  "injury_time"
)
time_col <- time_candidates[time_candidates %in% colnames(meta_atlas)][1]

if (is.na(time_col) || length(time_col) == 0L) {
  stop("Could not find Step10D timepoint field.", call. = FALSE)
}

# Step10D saved this exact review-only broad lineage field.
# If it is unexpectedly absent, reconstruct it from frozen cluster + review table.
if ("major_lineage_candidate" %in% colnames(meta_atlas)) {

  broad_type_col <- "major_lineage_candidate"

} else {

  if (!"hcon_res_0_2" %in% colnames(meta_atlas)) {
    stop(
      "Neither major_lineage_candidate nor hcon_res_0_2 is available in Step10D metadata.",
      call. = FALSE
    )
  }

  review_tbl <- read.csv(
    STEP10D_REVIEW_TABLE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!all(c("cluster", "top_panel") %in% colnames(review_tbl))) {
    stop(
      "Step10D review table must contain cluster and top_panel.",
      call. = FALSE
    )
  }

  candidate_map <- setNames(
    as.character(review_tbl$top_panel),
    as.character(review_tbl$cluster)
  )

  meta_atlas$major_lineage_candidate_rebuilt <- unname(
    candidate_map[
      as.character(meta_atlas$hcon_res_0_2)
    ]
  )

  if (anyNA(meta_atlas$major_lineage_candidate_rebuilt)) {
    stop(
      "Failed to reconstruct major-lineage candidate labels for every nucleus.",
      call. = FALSE
    )
  }

  broad_type_col <- "major_lineage_candidate_rebuilt"
}

cat("Using biological-library field: ", sample_col, "\n", sep = "")
cat("Using timepoint field         : ", time_col, "\n", sep = "")
cat("Using broad-lineage field     : ", broad_type_col, "\n\n", sep = "")

# ------------------------------------------------------------------------------
# 5. READ FROZEN STEP11B TARGET-PROGRAM LIBRARY SCORES
# ------------------------------------------------------------------------------
step11b <- data.table::fread(
  STEP11B_LIBRARY_TABLE,
  data.table = FALSE,
  showProgress = FALSE
)

required_step11b <- c(
  "record_type",
  "biological_library",
  "timepoint",
  "method",
  "program",
  "score"
)

miss11 <- setdiff(required_step11b, colnames(step11b))
if (length(miss11) > 0L) {
  stop(
    "Step11B library table missing columns: ",
    paste(miss11, collapse = ", "),
    call. = FALSE
  )
}

target_lib <- step11b %>%
  filter(
    record_type == "LIBRARY_SCORE",
    method == "AUCell_library_median",
    program == "Arm_2_specific"
  ) %>%
  transmute(
    biological_library = as.character(biological_library),
    timepoint = standardize_timepoint(timepoint),
    target_program = as.numeric(score),
    n_resident = if ("n_cells" %in% colnames(step11b)) {
      as.numeric(n_cells)
    } else {
      NA_real_
    }
  ) %>%
  filter(
    !is.na(biological_library),
    is.finite(target_program)
  ) %>%
  distinct(
    biological_library,
    timepoint,
    .keep_all = TRUE
  )

if (nrow(target_lib) < 4L) {
  stop(
    "Too few frozen Step11B target-program biological libraries were recovered.",
    call. = FALSE
  )
}

cat("Frozen Step11B target-program libraries: ", nrow(target_lib), "\n\n", sep = "")

# ------------------------------------------------------------------------------
# 6. READ FROZEN STEP8C RECIPIENT PROGRAMS
# ------------------------------------------------------------------------------
gs_tbl <- data.table::fread(
  STEP8C_GENESET_TSV,
  data.table = FALSE,
  showProgress = FALSE
)

if (!all(c("response_program", "present_genes") %in% colnames(gs_tbl))) {
  stop(
    "Step8C gene-set table must contain response_program and present_genes.",
    call. = FALSE
  )
}

extract_genes <- function(program_name) {
  x <- gs_tbl %>%
    filter(response_program == program_name)

  if (nrow(x) != 1L) {
    stop(
      "Could not uniquely identify Step8C program: ",
      program_name,
      call. = FALSE
    )
  }

  genes <- unlist(
    strsplit(
      as.character(x$present_genes[1]),
      "\\s*\\|\\s*"
    )
  )

  genes <- trimws(genes)
  unique(genes[nchar(genes) > 0])
}

genes_endo_activation <- extract_genes(
  "ENDOTHELIAL_ACTIVATION"
)

genes_endo_tnf <- extract_genes(
  "TNF_NFKB_RESPONSE"
)

# ------------------------------------------------------------------------------
# 7. ENDOTHELIAL CELLS + PROGRAM SCORING
# ------------------------------------------------------------------------------
endo_cells <- rownames(meta_atlas)[
  grepl(
    "^endothelial$|endothelial",
    tolower(
      as.character(
        meta_atlas[[broad_type_col]]
      )
    )
  )
]

if (length(endo_cells) == 0L) {
  stop(
    "No Endothelial cells found from Step10D broad-lineage labels.",
    call. = FALSE
  )
}

cat("Endothelial nuclei: ", length(endo_cells), "\n", sep = "")

expr_endo_activation <- get_normalized_subset(
  atlas,
  genes_endo_activation,
  endo_cells
)

expr_endo_tnf <- get_normalized_subset(
  atlas,
  genes_endo_tnf,
  endo_cells
)

endo_activation_scores <- program_score(
  expr_endo_activation
)

endo_tnf_scores <- program_score(
  expr_endo_tnf
)

endo_activation_lib <- library_mean_score(
  meta = meta_atlas,
  sample_col = sample_col,
  cells = colnames(expr_endo_activation),
  score = endo_activation_scores,
  out_name = "endothelial_activation"
)

endo_tnf_lib <- library_mean_score(
  meta = meta_atlas,
  sample_col = sample_col,
  cells = colnames(expr_endo_tnf),
  score = endo_tnf_scores,
  out_name = "endothelial_tnf_nfkB"
)

# ------------------------------------------------------------------------------
# 8. BROAD-COMPARTMENT LIBRARY PROPORTIONS
# ------------------------------------------------------------------------------
prop_endothelial <- library_prop(
  meta_atlas, sample_col, broad_type_col,
  "endothelial",
  "prop_endothelial"
)

prop_fibroblast <- library_prop(
  meta_atlas, sample_col, broad_type_col,
  "fibroblast",
  "prop_fibroblast"
)

prop_pericyte <- library_prop(
  meta_atlas, sample_col, broad_type_col,
  "pericyte",
  "prop_pericyte"
)

prop_astrocyte <- library_prop(
  meta_atlas, sample_col, broad_type_col,
  "astrocyte",
  "prop_astrocyte"
)

prop_opc <- library_prop(
  meta_atlas, sample_col, broad_type_col,
  "opc",
  "prop_opc"
)

prop_microglia <- library_prop(
  meta_atlas, sample_col, broad_type_col,
  "resident_microglia",
  "prop_microglia"
)

# ------------------------------------------------------------------------------
# 9. MERGE LIBRARY-LEVEL DATA
# ------------------------------------------------------------------------------
library_df <- target_lib %>%
  left_join(
    endo_activation_lib,
    by = "biological_library"
  ) %>%
  left_join(
    endo_tnf_lib,
    by = "biological_library"
  ) %>%
  left_join(
    prop_endothelial,
    by = "biological_library"
  ) %>%
  left_join(
    prop_fibroblast,
    by = "biological_library"
  ) %>%
  left_join(
    prop_pericyte,
    by = "biological_library"
  ) %>%
  left_join(
    prop_astrocyte,
    by = "biological_library"
  ) %>%
  left_join(
    prop_opc,
    by = "biological_library"
  ) %>%
  left_join(
    prop_microglia,
    by = "biological_library"
  ) %>%
  arrange(
    timepoint,
    biological_library
  )

write.csv(
  library_df,
  SUMMARY_CSV,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 10. LEFT PANEL STATISTICS
# ------------------------------------------------------------------------------
scatter_df <- library_df %>%
  filter(
    is.finite(target_program),
    is.finite(endothelial_tnf_nfkB)
  )

left_stat <- spearman_one(
  scatter_df$target_program,
  scatter_df$endothelial_tnf_nfkB
)

# ------------------------------------------------------------------------------
# 11. RIGHT PANEL CORRELATION MATRIX
#     Reference-style fixed grey circles + red/blue pie sectors
# ------------------------------------------------------------------------------
var_map <- c(
  "Microglia" = "prop_microglia",
  "OPC" = "prop_opc",
  "Astrocyte" = "prop_astrocyte",
  "Pericyte" = "prop_pericyte",
  "Fibroblast" = "prop_fibroblast",
  "Endothelial" = "prop_endothelial",
  "Endothelial TNF/NFκB" = "endothelial_tnf_nfkB",
  "Endothelial activation" = "endothelial_activation",
  "Target program" = "target_program"
)

valid <- vapply(
  unname(var_map),
  function(v) {
    if (!v %in% colnames(library_df)) return(FALSE)
    x <- library_df[[v]]
    x <- x[is.finite(x)]
    length(x) >= 4L && sd(x) > 0
  },
  logical(1)
)

var_map <- var_map[valid]

if (length(var_map) < 4L) {
  stop(
    "Too few valid library-level variables for Figure4F matrix.",
    call. = FALSE
  )
}

# Keep Target program at the top row of the reference-style triangle.
if ("Target program" %in% names(var_map)) {
  var_map <- c(
    var_map[names(var_map) != "Target program"],
    var_map["Target program"]
  )
}

matrix_levels <- names(var_map)
n_vars <- length(matrix_levels)

corr_list <- list()
k <- 1L

# Reference-style triangle:
# row i contains correlations with columns 1:i, including diagonal self-correlation.
for (i in seq_along(var_map)) {
  for (j in seq_len(i)) {

    row_name <- names(var_map)[i]
    col_name <- names(var_map)[j]

    if (i == j) {

      st <- list(
        rho = 1,
        p = 0,
        n = sum(is.finite(library_df[[unname(var_map[i])]]))
      )

    } else {

      st <- spearman_one(
        library_df[[unname(var_map[j])]],
        library_df[[unname(var_map[i])]]
      )
    }

    corr_list[[k]] <- tibble(
      row_name = row_name,
      col_name = col_name,
      row_index = i,
      col_index = j,
      rho = st$rho,
      p = st$p,
      n = st$n
    )

    k <- k + 1L
  }
}

corr_df <- bind_rows(corr_list) %>%
  mutate(
    abs_rho = pmin(abs(rho), 1),
    p_adj = NA_real_
  )

# Correct the 36 unique off-diagonal pairwise tests as one family.
# Diagonal self-correlations are descriptive and are excluded from correction.
idx_pairwise <- which(
  corr_df$row_name != corr_df$col_name & is.finite(corr_df$p)
)
if (length(idx_pairwise) != 36L) {
  stop(
    "Expected 36 unique off-diagonal correlations for nine variables; found ",
    length(idx_pairwise),
    ". Check the variables included in Figure 4F.",
    call. = FALSE
  )
}
corr_df$p_adj[idx_pairwise] <- stats::p.adjust(
  corr_df$p[idx_pairwise],
  method = "BH"
)

corr_df <- corr_df %>%
  mutate(
    status = case_when(
      row_name == col_name ~ "Pos",
      is.na(rho) | is.na(p_adj) ~ "Unsig",
      p_adj >= 0.05 ~ "Unsig",
      rho > 0 ~ "Pos",
      rho < 0 ~ "Neg",
      TRUE ~ "Unsig"
    )
  )

# The left-panel BH-adjusted P value is retrieved from the same 36-test family.
left_pair <- corr_df %>%
  filter(
    (row_name == "Target program" & col_name == "Endothelial TNF/NFκB") |
      (row_name == "Endothelial TNF/NFκB" & col_name == "Target program")
  )
if (nrow(left_pair) != 1L) {
  stop("Could not uniquely recover Target program–Endothelial TNF/NFκB pair.", call. = FALSE)
}
left_p_adj <- left_pair$p_adj[1]
left_label <- paste0(
  "\u03C1 = ", sprintf("%.3f", left_stat$rho),
  "\nBH-adjusted P = ", format_p(left_p_adj)
)

write.csv(
  corr_df,
  CORR_CSV,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 12. LEFT SCATTER
# ------------------------------------------------------------------------------
time_cols <- c(
  "Uninjured" = "#6D6D6D",
  "1 dpi" = "#3A923A",
  "4 dpi" = "#E17C05",
  "7 dpi" = "#7B62A3",
  "14 dpi" = "#D64A8A",
  "1 month" = "#3A8FB7",
  "2 months" = "#B68D19"
)

x_rng <- range(
  scatter_df$target_program,
  na.rm = TRUE
)

y_rng <- range(
  scatter_df$endothelial_tnf_nfkB,
  na.rm = TRUE
)

p_left <- ggplot(
  scatter_df,
  aes(
    x = target_program,
    y = endothelial_tnf_nfkB,
    colour = timepoint
  )
) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    colour = "black",
    fill = "grey75",
    alpha = 0.35,
    linewidth = 0.75,
    formula = y ~ x
  ) +
  geom_point(
    size = 3.2,
    alpha = 0.95
  ) +
  scale_colour_manual(
    values = time_cols,
    drop = FALSE
  ) +
  scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0.03, 0.18))
  ) +
  annotate(
    "text",
    x = x_rng[1] + 0.03 * diff(x_rng),
    y = y_rng[2] - 0.02 * diff(y_rng),
    label = left_label,
    hjust = 0,
    vjust = 1,
    family = FONT_FAMILY,
    fontface = "italic",
    size = 4.0
  ) +
  labs(
    title = "GSE234774",
    x = "Target-program AUCell score",
    y = "Endothelial TNF/NF-κB response score"
  ) +
  theme_classic(
    base_family = FONT_FAMILY
  ) +
  theme(
    plot.title = element_text(
      size = 11.5,
      hjust = 0.5,
      face = "plain"
    ),
    axis.title = element_text(
      size = 9.4,
      colour = "black"
    ),
    axis.text = element_text(
      size = 8.3,
      colour = "black"
    ),
    axis.line = element_line(
      linewidth = 0.45,
      colour = "black"
    ),
    legend.position = "none",
    plot.margin = margin(10, 3, 2, 2)
  )

# ------------------------------------------------------------------------------
# 13. RIGHT REFERENCE-STYLE PIE-SECTOR CORRELATION MATRIX
# ------------------------------------------------------------------------------
RADIUS <- 0.40
N_CIRCLE_SEGMENTS <- 100L

make_circle_polygon <- function(cx, cy, radius, group_id) {
  theta <- seq(0, 2 * pi, length.out = N_CIRCLE_SEGMENTS + 1L)

  tibble(
    x = cx + radius * cos(theta),
    y = cy + radius * sin(theta),
    group = group_id
  )
}

make_sector_polygon <- function(cx, cy, radius, fraction, group_id) {

  fraction <- max(0, min(1, fraction))

  if (fraction <= 0) {
    return(NULL)
  }

  if (fraction >= 0.999) {
    theta <- seq(
      0,
      2 * pi,
      length.out = N_CIRCLE_SEGMENTS + 1L
    )

    return(
      tibble(
        x = cx + radius * cos(theta),
        y = cy + radius * sin(theta),
        group = group_id
      )
    )
  }

  # Start at 12 o'clock and fill clockwise, as in the reference figure.
  theta_start <- pi / 2
  theta_end <- theta_start - 2 * pi * fraction

  n_arc <- max(
    8L,
    ceiling(N_CIRCLE_SEGMENTS * fraction)
  )

  theta <- seq(
    theta_start,
    theta_end,
    length.out = n_arc
  )

  tibble(
    x = c(
      cx,
      cx + radius * cos(theta),
      cx
    ),
    y = c(
      cy,
      cy + radius * sin(theta),
      cy
    ),
    group = group_id
  )
}

# Fixed grey background circles for every triangle position.
circle_polys <- lapply(
  seq_len(nrow(corr_df)),
  function(i) {
    make_circle_polygon(
      cx = corr_df$col_index[i],
      cy = corr_df$row_index[i],
      radius = RADIUS,
      group_id = paste0("circle_", i)
    )
  }
)

circle_df <- bind_rows(circle_polys)

# Colored sectors only for significant correlations and diagonal self-correlations.
sector_rows <- corr_df %>%
  filter(
    status %in% c("Pos", "Neg"),
    is.finite(abs_rho),
    abs_rho > 0
  )

sector_polys <- lapply(
  seq_len(nrow(sector_rows)),
  function(i) {

    poly <- make_sector_polygon(
      cx = sector_rows$col_index[i],
      cy = sector_rows$row_index[i],
      radius = RADIUS,
      fraction = sector_rows$abs_rho[i],
      group_id = paste0("sector_", i)
    )

    if (is.null(poly)) return(NULL)

    poly$status <- sector_rows$status[i]
    poly
  }
)

sector_df <- bind_rows(sector_polys)

legend_dummy <- tibble(
  x = rep(1, 3),
  y = rep(1, 3),
  status = factor(
    c("Pos", "Neg", "Unsig"),
    levels = c("Pos", "Neg", "Unsig")
  )
)

p_right <- ggplot() +
  geom_polygon(
    data = circle_df,
    aes(
      x = x,
      y = y,
      group = group
    ),
    fill = "#D9D9D9",
    colour = NA
  ) +
  geom_polygon(
    data = sector_df,
    aes(
      x = x,
      y = y,
      group = group,
      fill = status
    ),
    colour = NA
  ) +
  # Invisible dummy layer guarantees the 3-item legend, including Unsig.
  geom_point(
    data = legend_dummy,
    aes(
      x = x,
      y = y,
      fill = status
    ),
    shape = 22,
    size = 4,
    alpha = 0
  ) +
  scale_fill_manual(
    values = c(
      "Pos" = "#E31A1C",
      "Neg" = "#24579B",
      "Unsig" = "#D9D9D9"
    ),
    breaks = c(
      "Pos",
      "Neg",
      "Unsig"
    ),
    limits = c(
      "Pos",
      "Neg",
      "Unsig"
    ),
    drop = FALSE,
    name = "correlation"
  ) +
  scale_x_continuous(
    breaks = seq_len(n_vars),
    labels = matrix_levels,
    position = "top",
    limits = c(0.45, n_vars + 1.75),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq_len(n_vars),
    labels = matrix_levels,
    limits = c(0.45, n_vars + 0.55),
    expand = c(0, 0)
  ) +
  coord_fixed(
    clip = "off"
  ) +
  labs(
    title = "GSE234774",
    x = NULL,
    y = NULL
  ) +
  theme_void(
    base_family = FONT_FAMILY
  ) +
  theme(
    plot.title = element_text(
      size = 11.2,
      hjust = 0.5,
      face = "plain",
      margin = margin(b = 1)
    ),
    axis.text.x = element_text(
      size = 8.1,
      angle = 45,
      hjust = 0,
      vjust = 0,
      colour = "black",
      margin = margin(b = 0)
    ),
    axis.text.y = element_text(
      size = 8.3,
      colour = "black",
      margin = margin(r = 4)
    ),
    legend.title = element_text(
      size = 8.6
    ),
    legend.text = element_text(
      size = 8.0
    ),
    legend.key.size = unit(0.42, "cm"),
    plot.margin = margin(1, 1, 1, 1)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        alpha = 1,
        shape = 22,
        size = 5,
        colour = NA
      )
    )
  )

# Put the legend immediately to the right of the matrix, inside the right-panel
# viewport, so it cannot create a large outer blank region.
if (utils::packageVersion("ggplot2") >= "3.5.0") {
  p_right <- p_right +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.895, 0.49),
      legend.justification = c(0, 0.5),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
} else {
  p_right <- p_right +
    theme(
      legend.position = c(0.895, 0.49),
      legend.justification = c(0, 0.5),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
}

# ------------------------------------------------------------------------------
# 14. FINAL LAYOUT
#     Manual viewports avoid patchwork height-alignment whitespace.
# ------------------------------------------------------------------------------
draw_figure4F <- function() {

  grid::grid.newpage()

  # TRUE visible panel letter: fixed inside the page, never off-canvas.
  grid::grid.text(
    "F",
    x = grid::unit(0.012, "npc"),
    y = grid::unit(0.972, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(
      fontfamily = FONT_FAMILY,
      fontface = "bold",
      fontsize = 28,
      col = "black"
    )
  )

  # LEFT: independent viewport, so the right-panel rotated labels no longer
  # force extra whitespace above the scatter/title.
  print(
    p_left,
    vp = grid::viewport(
      x = 0.190,
      y = 0.475,
      width = 0.335,
      height = 0.850,
      just = c("center", "center")
    ),
    newpage = FALSE
  )

  # RIGHT: starts almost immediately after the left panel.
  # Legend is positioned inside the deliberately narrow right-side strip.
  print(
    p_right,
    vp = grid::viewport(
      x = 0.665,
      y = 0.490,
      width = 0.595,
      height = 0.910,
      just = c("center", "center")
    ),
    newpage = FALSE
  )
}

# ------------------------------------------------------------------------------
# 15. EXPORT: ONE PAGE
# ------------------------------------------------------------------------------
if (file.exists(PDF_FILE)) unlink(PDF_FILE)
if (file.exists(PNG_FILE)) unlink(PNG_FILE)

if (capabilities("cairo")) {
  grDevices::cairo_pdf(
    filename = PDF_FILE,
    width = 11.4,
    height = 4.55,
    family = FONT_FAMILY,
    onefile = FALSE,
    bg = "white"
  )
} else {
  grDevices::pdf(
    file = PDF_FILE,
    width = 11.4,
    height = 4.55,
    onefile = FALSE,
    useDingbats = FALSE,
    bg = "white"
  )
}

draw_figure4F()
grDevices::dev.off()

grDevices::png(
  filename = PNG_FILE,
  width = 3420,
  height = 1365,
  res = 300,
  bg = "white"
)
draw_figure4F()
grDevices::dev.off()

# ------------------------------------------------------------------------------
# 16. AUDIT
# ------------------------------------------------------------------------------
audit <- tibble(
  item = c(
    "Step10D atlas",
    "Step10D review table",
    "Step11B frozen library table",
    "Step8C gene-set table",
    "Step10D library field",
    "Step10D broad-lineage field",
    "Step11B target filter",
    "Endothelial nuclei",
    "Target libraries",
    "Libraries in scatter",
    "Left outcome",
    "Left Spearman R",
    "Left nominal p",
    "Left BH-adjusted p",
    "BH pairwise test family size",
    "PDF",
    "PNG"
  ),
  value = c(
    STEP10D_RDS,
    STEP10D_REVIEW_TABLE,
    STEP11B_LIBRARY_TABLE,
    STEP8C_GENESET_TSV,
    sample_col,
    broad_type_col,
    "LIBRARY_SCORE | AUCell_library_median | Arm_2_specific",
    as.character(length(endo_cells)),
    as.character(nrow(target_lib)),
    as.character(nrow(scatter_df)),
    "Target program vs endothelial TNF/NF-kB response",
    as.character(signif(left_stat$rho, 5)),
    as.character(signif(left_stat$p, 5)),
    as.character(signif(left_p_adj, 5)),
    as.character(length(idx_pairwise)),
    PDF_FILE,
    PNG_FILE
  )
)

write.csv(
  audit,
  AUDIT_CSV,
  row.names = FALSE
)

cat("\n============================================================\n")
cat("Figure 4F v8 completed.\n")
cat("Left panel: Target program vs endothelial TNF/NF-kB response.\n")
cat("Correlation matrix significance: BH-adjusted across 36 unique pairwise tests.\n\n")
cat("Target-program source:\n")
cat("  Frozen Step11B biological-library AUCell median\n")
cat("Recipient-program scoring:\n")
cat("  Step8C deterministic gene-wise z-score logic\n\n")
cat("Left panel:\n")
cat(
  "  Spearman rho = ",
  signif(left_stat$rho, 4),
  "; nominal p = ",
  signif(left_stat$p, 4),
  "; BH-adjusted p = ",
  signif(left_p_adj, 4),
  "; n = ",
  left_stat$n,
  "\n",
  sep = ""
)
cat("\nPDF:\n  ", PDF_FILE, "\n", sep = "")
cat("PNG:\n  ", PNG_FILE, "\n", sep = "")
cat("Summary CSV:\n  ", SUMMARY_CSV, "\n", sep = "")
cat("Correlation CSV:\n  ", CORR_CSV, "\n", sep = "")
cat("============================================================\n")
