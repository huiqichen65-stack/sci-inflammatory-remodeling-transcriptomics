# ==============================================================================
# 24_step5C_GSE162610_Y29_branch_dependent_genes_v2.R
# ==============================================================================
# STEP 5C v2
# Y_29 branch-dependent genes + branch smoothing
#
# IMPORTANT CORRECTIONS FROM v1
#
# 1) Formal compare_models() inference now uses NEGATIVE BINOMIAL.
#    Quasi-Poisson has no true likelihood and is inappropriate for formal
#    likelihood-ratio comparison with compare_models().
#
# 2) Formal branch-dependence test is BEAM-like:
#
#      FULL:
#        expression ~ sample
#                   + branch
#                   + ns(pseudotime, df=3)
#                   + branch:ns(pseudotime, df=3)
#
#      NULL:
#        expression ~ sample
#                   + ns(pseudotime, df=3)
#
#    This asks whether allowing branch-specific expression trajectories gives
#    significantly better fit than one common trajectory.
#
# 3) A stricter SECONDARY shape-only interaction test is also calculated:
#
#      FULL versus ADDITIVE:
#        expression ~ sample
#                   + branch
#                   + ns(pseudotime, df=3)
#
#    Thus:
#      branch_dependence_q = primary Step5C branch-gene evidence
#      shape_interaction_q = secondary evidence for different curve shapes
#
# 4) Windows is explicitly run with cores=1 for Monocle model fitting.
#
# 5) If zero genes pass a threshold, the script DOES NOT abort before saving
#    diagnostics. It produces explicit null-result figures/tables instead.
#
# FIXED BIOLOGY
#   - Branch topology: corrected Step5B-R2 Y_29 Arm_1 / Arm_2
#   - Dynamics pseudotime: Step5B-R balanced 225-cell-root pseudotime
#   - Failed Y_19 single-node pseudotime is NOT used.
#   - repair/scar candidate scores are NOT used for gene discovery.
#
# This is still STEP 5. GO/KEGG/GSEA and formal pseudobulk DEGs are Step 6.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(20260731)

# ==============================================================================
# 0. PATHS / PARAMETERS
# ==============================================================================

PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
STEP5BR_DIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step5BR_GSE162610_injury_response_trajectory_refinement"
)

STEP5BR2_DIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step5BR2_GSE162610_single_root_branch_audit"
)

INPUT_CDS <- file.path(
  STEP5BR_DIR,
  "GSE162610_step5BR_injury_monocle3_cds.rds"
)

INPUT_SEURAT <- file.path(
  STEP5BR2_DIR,
  "GSE162610_step5BR2_single_root_seurat.rds"
)

OUTDIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step5C_GSE162610_Y29_branch_dependent_genes_v2"
)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

CACHE_GRAPH <- file.path(
  OUTDIR,
  "cache_graph_test.rds"
)

OUT_FIG1 <- file.path(
  OUTDIR,
  "Fig5C_v2_01_branch_gene_discovery.pdf"
)

OUT_FIG2 <- file.path(
  OUTDIR,
  "Fig5C_v2_02_branch_dynamic_heatmap.pdf"
)

OUT_FIG3 <- file.path(
  OUTDIR,
  "Fig5C_v2_03_top_branch_gene_smoothing.pdf"
)

OUT_FIG4 <- file.path(
  OUTDIR,
  "Fig5C_v2_04_branch_module_replicate_audit.pdf"
)

OUT_XLSX <- file.path(
  OUTDIR,
  "Table5C_v2_Y29_branch_dependent_genes.xlsx"
)

OUT_MODULE_RDS <- file.path(
  OUTDIR,
  "GSE162610_step5C_v2_branch_gene_modules.rds"
)

OUT_SEURAT <- file.path(
  OUTDIR,
  "GSE162610_step5C_v2_branch_scored_seurat.rds"
)

OUT_REP <- file.path(
  OUTDIR,
  "Step5C_v2_report.txt"
)

OUT_SESS <- file.path(
  OUTDIR,
  "Step5C_v2_sessionInfo.txt"
)

TIME_COL <- "timepoint"
SAMPLE_COL <- "sample_id"
PSEUDOTIME_COL <- "step5BR_pseudotime"
BRANCH_ARM_COL <- "step5BR2_branch_arm"
PRIMARY_BRANCHPOINT_COL <- "step5BR2_primary_branchpoint"

EXPECTED_BRANCHPOINT <- "Y_29"

MIN_DETECTION_FRACTION <- 0.05
GRAPH_Q_MAX <- 0.05
GRAPH_MORAN_MIN <- 0.05

BRANCH_DEP_Q_MAX <- 0.05
SHAPE_INTERACTION_Q_MAX <- 0.05

MAX_GENES_FOR_NB_MODEL <- 5000

N_PSEUDOTIME_BINS <- 20
TOP_GENES_PER_ARM_HEATMAP <- 35
TOP_GENES_PER_ARM_SMOOTH <- 6

MIN_ABS_7DPI_MEAN_DIFF <- 0.10

# pbmclapply parallelism used internally by Monocle is not supported on Windows.
CORES <- if (.Platform$OS.type == "windows") {
  1L
} else {
  max(
    1L,
    min(
      4L,
      parallel::detectCores(logical = TRUE) - 1L
    )
  )
}

# ==============================================================================
# 1. PACKAGES
# ==============================================================================

required_pkgs <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "monocle3",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "S4Vectors",
  "splines",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "patchwork",
  "openxlsx"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "缺少以下R包：",
      paste(missing_pkgs, collapse = ", ")
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(monocle3)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(splines)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# ==============================================================================
# 2. HELPERS
# ==============================================================================

msg <- function(...) {
  cat(
    sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
    ...,
    "\n",
    sep = ""
  )
}

capture_table <- function(x) {
  paste(
    capture.output(
      print(
        tibble::as_tibble(x),
        n = Inf,
        width = Inf
      )
    ),
    collapse = "\n"
  )
}

get_normalized_matrix <- function(obj) {

  assay <- DefaultAssay(obj)

  x <- tryCatch(
    Seurat::GetAssayData(
      obj,
      assay = assay,
      layer = "data"
    ),
    error = function(e) NULL
  )

  if (!is.null(x)) {
    return(x)
  }

  Seurat::GetAssayData(
    obj,
    assay = assay,
    slot = "data"
  )
}

safe_add_sheet <- function(
    wb,
    sheet_name,
    df) {

  openxlsx::addWorksheet(
    wb,
    sheetName = sheet_name
  )

  ddf <- tryCatch(
    as.data.frame(df),
    error = function(e) NULL
  )

  if (
    is.null(ddf) ||
    nrow(ddf) == 0 ||
    ncol(ddf) == 0
  ) {

    openxlsx::writeData(
      wb,
      sheet = sheet_name,
      x = data.frame(
        note =
          "No results generated under the predefined Step5C-v2 criteria.",
        stringsAsFactors = FALSE
      )
    )

  } else {

    openxlsx::writeDataTable(
      wb,
      sheet = sheet_name,
      x = ddf,
      tableStyle = "TableStyleMedium2"
    )
  }

  openxlsx::freezePane(
    wb,
    sheet = sheet_name,
    firstRow = TRUE
  )
}

standardize_compare <- function(
    x,
    prefix) {

  x <- tibble::as_tibble(
    x,
    rownames = "gene_id_from_rownames"
  )

  if (!"gene_id" %in% colnames(x)) {

    if ("id" %in% colnames(x)) {
      x$gene_id <- as.character(x$id)
    } else {
      x$gene_id <- as.character(
        x$gene_id_from_rownames
      )
    }
  }

  p_candidates <- c(
    "p_value",
    "pval",
    "p.value",
    "lrtest_p_value"
  )

  p_col <- p_candidates[
    p_candidates %in%
      colnames(x)
  ][1]

  if (
    length(p_col) == 0 ||
    is.na(p_col)
  ) {
    stop(
      paste0(
        "compare_models output lacks a recognized p-value column. Columns: ",
        paste(colnames(x), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  p <- as.numeric(
    x[[p_col]]
  )

  q_candidates <- c(
    "q_value",
    "qval",
    "q.value"
  )

  q_col <- q_candidates[
    q_candidates %in%
      colnames(x)
  ][1]

  if (
    length(q_col) > 0 &&
    !is.na(q_col)
  ) {
    q <- as.numeric(
      x[[q_col]]
    )
  } else {
    q <- stats::p.adjust(
      p,
      method = "BH"
    )
  }

  out <- tibble::tibble(
    gene_id =
      as.character(
        x$gene_id
      )
  )

  out[[paste0(prefix, "_p_value")]] <- p
  out[[paste0(prefix, "_q_value")]] <- q

  out
}

row_zscore <- function(mat) {

  mu <- rowMeans(
    mat,
    na.rm = TRUE
  )

  ss <- apply(
    mat,
    1,
    stats::sd,
    na.rm = TRUE
  )

  ss[
    !is.finite(ss) |
      ss == 0
  ] <- 1

  z <- sweep(
    mat,
    1,
    mu,
    "-"
  )

  sweep(
    z,
    1,
    ss,
    "/"
  )
}

null_plot <- function(label) {

  ggplot2::ggplot(
    tibble::tibble(
      x = 1,
      y = 1
    ),
    ggplot2::aes(x, y)
  ) +
    ggplot2::annotate(
      "text",
      x = 1,
      y = 1,
      label = label,
      size = 5
    ) +
    ggplot2::theme_void()
}

# ==============================================================================
# 3. READ INPUTS
# ==============================================================================

msg("====================================================")
msg("GSE162610 Step 5C v2 | corrected likelihood branch-gene test")
msg("====================================================")

if (!file.exists(INPUT_CDS)) {
  stop(
    paste0(
      "找不到Step5B-R CDS：",
      INPUT_CDS
    ),
    call. = FALSE
  )
}

if (!file.exists(INPUT_SEURAT)) {
  stop(
    paste0(
      "找不到Step5B-R2 Seurat object：",
      INPUT_SEURAT
    ),
    call. = FALSE
  )
}

cds <- readRDS(
  INPUT_CDS
)

obj <- readRDS(
  INPUT_SEURAT
)

DefaultAssay(obj) <- "RNA"

required_meta <- c(
  TIME_COL,
  SAMPLE_COL,
  PSEUDOTIME_COL,
  BRANCH_ARM_COL,
  PRIMARY_BRANCHPOINT_COL
)

missing_meta <- required_meta[
  !required_meta %in%
    colnames(obj[[]])
]

if (length(missing_meta) > 0) {
  stop(
    paste0(
      "缺少Step5C-v2所需metadata：",
      paste(
        missing_meta,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

# ==============================================================================
# 4. FREEZE Y_29 / ARM_1 / ARM_2
# ==============================================================================

bp_values <- unique(
  stats::na.omit(
    as.character(
      obj[[PRIMARY_BRANCHPOINT_COL, drop = TRUE]]
    )
  )
)

if (length(bp_values) != 1) {
  stop(
    "Step5B-R2 primary branchpoint is not uniquely defined.",
    call. = FALSE
  )
}

PRIMARY_BRANCHPOINT <- bp_values[1]

msg(
  "Frozen corrected branchpoint: ",
  PRIMARY_BRANCHPOINT
)

if (
  PRIMARY_BRANCHPOINT !=
    EXPECTED_BRANCHPOINT
) {
  warning(
    paste0(
      "Stored branchpoint is ",
      PRIMARY_BRANCHPOINT,
      " rather than expected Y_29."
    )
  )
}

obj_cells <- colnames(
  obj
)

cds_cells <- colnames(
  cds
)

if (!all(cds_cells %in% obj_cells)) {
  stop(
    "Some CDS cells are absent from Step5B-R2 Seurat object.",
    call. = FALSE
  )
}

branch_arm_all <- as.character(
  obj[[BRANCH_ARM_COL, drop = TRUE]]
)
names(branch_arm_all) <- obj_cells

pt_all <- as.numeric(
  obj[[PSEUDOTIME_COL, drop = TRUE]]
)
names(pt_all) <- obj_cells

time_all <- as.character(
  obj[[TIME_COL, drop = TRUE]]
)
names(time_all) <- obj_cells

sample_all <- as.character(
  obj[[SAMPLE_COL, drop = TRUE]]
)
names(sample_all) <- obj_cells

arm_counts_7dpi <- tibble::tibble(
  cell = obj_cells,
  branch_arm =
    branch_arm_all[
      obj_cells
    ],
  timepoint =
    time_all[
      obj_cells
    ]
) |>
  dplyr::filter(
    timepoint == "7 dpi",
    !is.na(branch_arm),
    grepl(
      "^Arm_",
      branch_arm
    )
  ) |>
  dplyr::count(
    branch_arm,
    name = "n_7dpi"
  ) |>
  dplyr::arrange(
    dplyr::desc(n_7dpi)
  )

if (nrow(arm_counts_7dpi) < 2) {
  stop(
    "Fewer than two corrected downstream arms contain 7 dpi cells.",
    call. = FALSE
  )
}

MAJOR_ARMS <- arm_counts_7dpi$branch_arm[
  1:2
]

arm1 <- MAJOR_ARMS[1]
arm2 <- MAJOR_ARMS[2]

msg(
  "Major arms: ",
  arm1,
  " and ",
  arm2
)

branch_cells <- obj_cells[
  branch_arm_all[
    obj_cells
  ] %in%
    MAJOR_ARMS &
    time_all[
      obj_cells
    ] %in%
    c(
      "3 dpi",
      "7 dpi"
    ) &
    is.finite(
      pt_all[
        obj_cells
      ]
    )
]

branch_meta <- tibble::tibble(
  cell = branch_cells,
  branch_arm =
    branch_arm_all[
      branch_cells
    ],
  timepoint =
    time_all[
      branch_cells
    ],
  sample =
    sample_all[
      branch_cells
    ],
  pseudotime_step5BR =
    pt_all[
      branch_cells
    ]
)

branch_cell_summary <- branch_meta |>
  dplyr::count(
    branch_arm,
    timepoint,
    sample,
    name = "n_cells"
  )

msg(
  "Branch-model cells: ",
  length(
    branch_cells
  )
)

print(
  branch_meta |>
    dplyr::count(
      branch_arm,
      timepoint,
      name = "n_cells"
    )
)

# ==============================================================================
# 5. GENE PREFILTER + graph_test
# ==============================================================================

counts <- SummarizedExperiment::assay(
  cds,
  "counts"
)

counts_branch <- counts[
  ,
  branch_cells,
  drop = FALSE
]

detected_n <- Matrix::rowSums(
  counts_branch > 0
)

gene_prefilter_tbl <- tibble::tibble(
  gene_id =
    rownames(
      counts_branch
    ),
  detected_cells =
    as.integer(
      detected_n
    ),
  detection_fraction =
    as.numeric(
      detected_n /
        ncol(
          counts_branch
        )
    )
)

prefilter_genes <- gene_prefilter_tbl$gene_id[
  gene_prefilter_tbl$detection_fraction >=
    MIN_DETECTION_FRACTION
]

msg(
  "Genes detected in >=5% of branch cells: ",
  length(
    prefilter_genes
  )
)

# Cache graph_test so reruns do not repeat the ~9 minute calculation.
if (file.exists(CACHE_GRAPH)) {

  msg(
    "Loading cached graph_test results..."
  )

  graph_test_tbl <- readRDS(
    CACHE_GRAPH
  )

} else {

  msg(
    "Running graph_test on ",
    length(
      prefilter_genes
    ),
    " genes; cores = ",
    CORES,
    "..."
  )

  cds_graph <- cds[
    prefilter_genes,
    ,
    drop = FALSE
  ]

  graph_test_res <- monocle3::graph_test(
    cds_graph,
    neighbor_graph =
      "principal_graph",
    cores = CORES
  )

  graph_test_tbl <- tibble::as_tibble(
    graph_test_res,
    rownames = "gene_id"
  ) |>
    dplyr::left_join(
      gene_prefilter_tbl,
      by = "gene_id"
    ) |>
    dplyr::mutate(
      graph_dynamic =
        !is.na(q_value) &
        q_value <=
          GRAPH_Q_MAX &
        morans_I >=
          GRAPH_MORAN_MIN
    ) |>
    dplyr::arrange(
      q_value,
      dplyr::desc(
        morans_I
      )
    )

  saveRDS(
    graph_test_tbl,
    CACHE_GRAPH
  )
}

graph_dynamic_genes <- graph_test_tbl$gene_id[
  graph_test_tbl$graph_dynamic
]

msg(
  "Graph-dynamic genes: ",
  length(
    graph_dynamic_genes
  )
)

fit_gene_tbl <- graph_test_tbl |>
  dplyr::filter(
    graph_dynamic
  ) |>
  dplyr::arrange(
    q_value,
    dplyr::desc(
      morans_I
    )
  )

if (
  nrow(fit_gene_tbl) >
    MAX_GENES_FOR_NB_MODEL
) {

  fit_gene_tbl <- fit_gene_tbl |>
    dplyr::slice_head(
      n =
        MAX_GENES_FOR_NB_MODEL
    )
}

fit_genes <- fit_gene_tbl$gene_id

if (length(fit_genes) < 10) {
  stop(
    "Too few graph-dynamic genes for branch modeling.",
    call. = FALSE
  )
}

msg(
  "Genes entering NEGATIVE-BINOMIAL branch model: ",
  length(
    fit_genes
  )
)

# ==============================================================================
# 6. PREPARE CDS FOR LIKELIHOOD-BASED BRANCH TEST
# ==============================================================================

cds_branch <- cds[
  fit_genes,
  branch_cells,
  drop = FALSE
]

meta_ordered <- branch_meta[
  match(
    colnames(cds_branch),
    branch_meta$cell
  ),
  ,
  drop = FALSE
]

if (
  !identical(
    as.character(
      meta_ordered$cell
    ),
    colnames(
      cds_branch
    )
  )
) {
  stop(
    "Branch metadata order mismatch.",
    call. = FALSE
  )
}

pt_min <- min(
  meta_ordered$pseudotime_step5BR
)

pt_max <- max(
  meta_ordered$pseudotime_step5BR
)

meta_ordered$pseudotime_scaled <-
  (
    meta_ordered$pseudotime_step5BR -
      pt_min
  ) /
  (
    pt_max -
      pt_min
  )

SummarizedExperiment::colData(
  cds_branch
)$step5C_sample <-
  factor(
    meta_ordered$sample
  )

SummarizedExperiment::colData(
  cds_branch
)$step5C_branch_arm <-
  factor(
    meta_ordered$branch_arm,
    levels =
      MAJOR_ARMS
  )

SummarizedExperiment::colData(
  cds_branch
)$step5C_pseudotime_scaled <-
  meta_ordered$pseudotime_scaled

# ==============================================================================
# 7. CORRECTED NEGATIVE-BINOMIAL MODELS
# ==============================================================================

# BEAM-like formal branch-dependence test:
# full = branch-specific trajectory
# null = one common trajectory
FULL_FORMULA <- paste0(
  "~ step5C_sample + ",
  "step5C_branch_arm + ",
  "splines::ns(step5C_pseudotime_scaled, df = 3) + ",
  "step5C_branch_arm:splines::ns(step5C_pseudotime_scaled, df = 3)"
)

NULL_FORMULA <- paste0(
  "~ step5C_sample + ",
  "splines::ns(step5C_pseudotime_scaled, df = 3)"
)

# Secondary test:
# Does adding the interaction improve a model that already knows branch identity?
ADDITIVE_FORMULA <- paste0(
  "~ step5C_sample + ",
  "step5C_branch_arm + ",
  "splines::ns(step5C_pseudotime_scaled, df = 3)"
)

msg(
  "Formal family: NEGATIVE BINOMIAL"
)

msg(
  "Full model: ",
  FULL_FORMULA
)

msg(
  "Null model: ",
  NULL_FORMULA
)

msg(
  "Additive model: ",
  ADDITIVE_FORMULA
)

msg(
  "Fitting FULL negative-binomial model. Windows cores = ",
  CORES,
  "..."
)

full_fit <- monocle3::fit_models(
  cds_branch,
  model_formula_str =
    FULL_FORMULA,
  expression_family =
    "negbinomial",
  clean_model = TRUE,
  cores = CORES
)

msg(
  "Fitting NULL negative-binomial model..."
)

null_fit <- monocle3::fit_models(
  cds_branch,
  model_formula_str =
    NULL_FORMULA,
  expression_family =
    "negbinomial",
  clean_model = TRUE,
  cores = CORES
)

msg(
  "Fitting ADDITIVE negative-binomial model..."
)

additive_fit <- monocle3::fit_models(
  cds_branch,
  model_formula_str =
    ADDITIVE_FORMULA,
  expression_family =
    "negbinomial",
  clean_model = TRUE,
  cores = CORES
)

fit_status_tbl <- tibble::tibble(
  gene_id = fit_genes,
  full_status =
    as.character(
      full_fit$status[
        match(
          fit_genes,
          full_fit$gene_id
        )
      ]
    ),
  null_status =
    as.character(
      null_fit$status[
        match(
          fit_genes,
          null_fit$gene_id
        )
      ]
    ),
  additive_status =
    as.character(
      additive_fit$status[
        match(
          fit_genes,
          additive_fit$gene_id
        )
      ]
    )
)

msg(
  "FULL model OK: ",
  sum(
    fit_status_tbl$full_status ==
      "OK",
    na.rm = TRUE
  ),
  "/",
  nrow(
    fit_status_tbl
  )
)

msg(
  "NULL model OK: ",
  sum(
    fit_status_tbl$null_status ==
      "OK",
    na.rm = TRUE
  ),
  "/",
  nrow(
    fit_status_tbl
  )
)

msg(
  "ADDITIVE model OK: ",
  sum(
    fit_status_tbl$additive_status ==
      "OK",
    na.rm = TRUE
  ),
  "/",
  nrow(
    fit_status_tbl
  )
)

# ==============================================================================
# 8. LIKELIHOOD-RATIO TESTS
# ==============================================================================

msg(
  "Formal BEAM-like comparison: FULL vs NULL..."
)

branch_dep_raw <- monocle3::compare_models(
  full_fit,
  null_fit
)

branch_dep_tbl <- standardize_compare(
  branch_dep_raw,
  "branch_dependence"
)

msg(
  "Secondary interaction-only comparison: FULL vs ADDITIVE..."
)

shape_raw <- monocle3::compare_models(
  full_fit,
  additive_fit
)

shape_tbl <- standardize_compare(
  shape_raw,
  "shape_interaction"
)

# Diagnostics that matter if no gene survives BH correction.
branch_dep_valid <- branch_dep_tbl |>
  dplyr::filter(
    is.finite(
      branch_dependence_p_value
    )
  )

shape_valid <- shape_tbl |>
  dplyr::filter(
    is.finite(
      shape_interaction_p_value
    )
  )

msg(
  "Valid branch-dependence LRT p-values: ",
  nrow(
    branch_dep_valid
  )
)

msg(
  "Minimum branch-dependence p = ",
  signif(
    min(
      branch_dep_valid$branch_dependence_p_value,
      na.rm = TRUE
    ),
    4
  )
)

msg(
  "Minimum branch-dependence q = ",
  signif(
    min(
      branch_dep_valid$branch_dependence_q_value,
      na.rm = TRUE
    ),
    4
  )
)

msg(
  "Raw branch-dependence p<0.05 genes = ",
  sum(
    branch_dep_valid$branch_dependence_p_value <
      0.05,
    na.rm = TRUE
  )
)

# ==============================================================================
# 9. 7 dpi EFFECTS + REPLICATE DIRECTION
# ==============================================================================

msg(
  "Calculating descriptive 7 dpi Arm_1 vs Arm_2 effects..."
)

norm_mat <- get_normalized_matrix(
  obj
)

available_fit_genes <- intersect(
  fit_genes,
  rownames(
    norm_mat
  )
)

seven_cells <- branch_cells[
  time_all[
    branch_cells
  ] ==
    "7 dpi"
]

seven_arm1 <- seven_cells[
  branch_arm_all[
    seven_cells
  ] ==
    arm1
]

seven_arm2 <- seven_cells[
  branch_arm_all[
    seven_cells
  ] ==
    arm2
]

mean_arm1 <- Matrix::rowMeans(
  norm_mat[
    available_fit_genes,
    seven_arm1,
    drop = FALSE
  ]
)

mean_arm2 <- Matrix::rowMeans(
  norm_mat[
    available_fit_genes,
    seven_arm2,
    drop = FALSE
  ]
)

terminal_effect_tbl <- tibble::tibble(
  gene_id =
    available_fit_genes,
  mean_7dpi_arm1 =
    as.numeric(
      mean_arm1
    ),
  mean_7dpi_arm2 =
    as.numeric(
      mean_arm2
    ),
  mean_diff_arm1_minus_arm2 =
    as.numeric(
      mean_arm1 -
        mean_arm2
    )
)

seven_samples <- sort(
  unique(
    sample_all[
      seven_cells
    ]
  )
)

replicate_effect_tbl <- dplyr::bind_rows(
  lapply(
    seven_samples,
    function(smp) {

      c1 <- seven_arm1[
        sample_all[
          seven_arm1
        ] ==
          smp
      ]

      c2 <- seven_arm2[
        sample_all[
          seven_arm2
        ] ==
          smp
      ]

      if (
        length(c1) == 0 ||
        length(c2) == 0
      ) {
        return(NULL)
      }

      d1 <- Matrix::rowMeans(
        norm_mat[
          available_fit_genes,
          c1,
          drop = FALSE
        ]
      )

      d2 <- Matrix::rowMeans(
        norm_mat[
          available_fit_genes,
          c2,
          drop = FALSE
        ]
      )

      tibble::tibble(
        gene_id =
          available_fit_genes,
        sample = smp,
        n_arm1 =
          length(c1),
        n_arm2 =
          length(c2),
        mean_diff_arm1_minus_arm2 =
          as.numeric(
            d1 -
              d2
          )
      )
    }
  )
)

replicate_direction_tbl <- replicate_effect_tbl |>
  dplyr::group_by(
    gene_id
  ) |>
  dplyr::summarise(
    n_7dpi_samples =
      dplyr::n(),
    n_positive =
      sum(
        mean_diff_arm1_minus_arm2 > 0
      ),
    n_negative =
      sum(
        mean_diff_arm1_minus_arm2 < 0
      ),
    direction_consistent =
      (
        n_positive ==
          n_7dpi_samples
      ) |
      (
        n_negative ==
          n_7dpi_samples
      ),
    .groups = "drop"
  )

# ==============================================================================
# 10. MASTER TABLE
# ==============================================================================

branch_gene_tbl <- fit_gene_tbl |>
  dplyr::select(
    gene_id,
    graph_morans_I =
      morans_I,
    graph_p_value =
      p_value,
    graph_q_value =
      q_value,
    detection_fraction
  ) |>
  dplyr::left_join(
    fit_status_tbl,
    by = "gene_id"
  ) |>
  dplyr::left_join(
    branch_dep_tbl,
    by = "gene_id"
  ) |>
  dplyr::left_join(
    shape_tbl,
    by = "gene_id"
  ) |>
  dplyr::left_join(
    terminal_effect_tbl,
    by = "gene_id"
  ) |>
  dplyr::left_join(
    replicate_direction_tbl,
    by = "gene_id"
  ) |>
  dplyr::mutate(
    branch_dependent =
      is.finite(
        branch_dependence_q_value
      ) &
      branch_dependence_q_value <=
        BRANCH_DEP_Q_MAX,

    shape_interaction_significant =
      is.finite(
        shape_interaction_q_value
      ) &
      shape_interaction_q_value <=
        SHAPE_INTERACTION_Q_MAX,

    terminal_branch_bias =
      dplyr::case_when(
        mean_diff_arm1_minus_arm2 > 0 ~
          arm1,
        mean_diff_arm1_minus_arm2 < 0 ~
          arm2,
        TRUE ~
          "Neutral"
      ),

    abs_7dpi_mean_diff =
      abs(
        mean_diff_arm1_minus_arm2
      ),

    robust_display_gene =
      branch_dependent &
      direction_consistent &
      abs_7dpi_mean_diff >=
        MIN_ABS_7DPI_MEAN_DIFF
  ) |>
  dplyr::arrange(
    branch_dependence_q_value,
    dplyr::desc(
      abs_7dpi_mean_diff
    )
  )

branch_sig_tbl <- branch_gene_tbl |>
  dplyr::filter(
    branch_dependent
  )

shape_sig_tbl <- branch_gene_tbl |>
  dplyr::filter(
    shape_interaction_significant
  )

robust_tbl <- branch_gene_tbl |>
  dplyr::filter(
    robust_display_gene
  )

arm1_genes <- branch_sig_tbl$gene_id[
  branch_sig_tbl$terminal_branch_bias ==
    arm1
]

arm2_genes <- branch_sig_tbl$gene_id[
  branch_sig_tbl$terminal_branch_bias ==
    arm2
]

msg(
  "Formal BEAM-like branch-dependent genes q<=0.05: ",
  nrow(
    branch_sig_tbl
  )
)

msg(
  "Secondary shape-interaction genes q<=0.05: ",
  nrow(
    shape_sig_tbl
  )
)

msg(
  arm1,
  "-biased formal genes: ",
  length(
    arm1_genes
  )
)

msg(
  arm2,
  "-biased formal genes: ",
  length(
    arm2_genes
  )
)

msg(
  "Robust display genes: ",
  nrow(
    robust_tbl
  )
)

# ==============================================================================
# 11. CHOOSE GENES FOR DYNAMIC VISUALIZATION
# ==============================================================================

choose_top <- function(
    tbl,
    arm_name,
    n) {

  x <- tbl |>
    dplyr::filter(
      branch_dependent,
      terminal_branch_bias ==
        arm_name
    ) |>
    dplyr::arrange(
      dplyr::desc(
        direction_consistent
      ),
      branch_dependence_q_value,
      dplyr::desc(
        abs_7dpi_mean_diff
      )
    )

  head(
    x$gene_id,
    n
  )
}

top_arm1_heat <- choose_top(
  branch_gene_tbl,
  arm1,
  TOP_GENES_PER_ARM_HEATMAP
)

top_arm2_heat <- choose_top(
  branch_gene_tbl,
  arm2,
  TOP_GENES_PER_ARM_HEATMAP
)

heatmap_genes <- unique(
  c(
    top_arm1_heat,
    top_arm2_heat
  )
)

top_arm1_smooth <- head(
  top_arm1_heat,
  TOP_GENES_PER_ARM_SMOOTH
)

top_arm2_smooth <- head(
  top_arm2_heat,
  TOP_GENES_PER_ARM_SMOOTH
)

smooth_genes <- unique(
  c(
    top_arm1_smooth,
    top_arm2_smooth
  )
)

# ==============================================================================
# 12. BINNED BRANCH SMOOTHING
# ==============================================================================

smooth_tbl <- tibble::tibble()
module_bin_summary <- tibble::tibble()
module_7dpi_sample_summary <- tibble::tibble()

if (length(heatmap_genes) >= 2) {

  msg(
    "Building branch smoothing across ",
    N_PSEUDOTIME_BINS,
    " pseudotime bins..."
  )

  branch_model_meta <- branch_meta |>
    dplyr::mutate(
      pseudotime_scaled =
        (
          pseudotime_step5BR -
            min(
              pseudotime_step5BR
            )
        ) /
        (
          max(
            pseudotime_step5BR
          ) -
            min(
              pseudotime_step5BR
            )
        ),

      pt_bin =
        cut(
          pseudotime_scaled,
          breaks =
            seq(
              0,
              1,
              length.out =
                N_PSEUDOTIME_BINS + 1
            ),
          include.lowest = TRUE,
          labels = FALSE
        )
    )

  bin_meta <- branch_model_meta |>
    dplyr::group_by(
      branch_arm,
      pt_bin
    ) |>
    dplyr::summarise(
      n_cells =
        dplyr::n(),
      pt_scaled_median =
        median(
          pseudotime_scaled
        ),
      pt_raw_median =
        median(
          pseudotime_step5BR
        ),
      .groups = "drop"
    )

  smooth_expr_genes <- intersect(
    heatmap_genes,
    rownames(
      norm_mat
    )
  )

  smooth_list <- list()

  for (aa in MAJOR_ARMS) {

    bins_aa <- sort(
      unique(
        branch_model_meta$pt_bin[
          branch_model_meta$branch_arm ==
            aa
        ]
      )
    )

    for (bb in bins_aa) {

      cc <- branch_model_meta$cell[
        branch_model_meta$branch_arm ==
          aa &
        branch_model_meta$pt_bin ==
          bb
      ]

      if (length(cc) == 0) {
        next
      }

      vv <- Matrix::rowMeans(
        norm_mat[
          smooth_expr_genes,
          cc,
          drop = FALSE
        ]
      )

      smooth_list[[
        paste0(
          aa,
          "__",
          bb
        )
      ]] <- tibble::tibble(
        gene_id =
          smooth_expr_genes,
        branch_arm = aa,
        pt_bin =
          as.integer(bb),
        mean_expression =
          as.numeric(vv),
        n_cells =
          length(cc)
      )
    }
  }

  smooth_tbl <- dplyr::bind_rows(
    smooth_list
  ) |>
    dplyr::left_join(
      bin_meta,
      by = c(
        "branch_arm",
        "pt_bin"
      )
    )
}

# ==============================================================================
# 13. MODULE SCORES
# ==============================================================================

score_module <- function(
    genes,
    cells) {

  genes <- intersect(
    genes,
    rownames(
      norm_mat
    )
  )

  if (length(genes) == 0) {
    return(
      rep(
        NA_real_,
        length(cells)
      )
    )
  }

  as.numeric(
    Matrix::colMeans(
      norm_mat[
        genes,
        cells,
        drop = FALSE
      ]
    )
  )
}

obj$step5C_v2_Arm1_branch_module_score <-
  NA_real_

obj$step5C_v2_Arm2_branch_module_score <-
  NA_real_

idx_score <- match(
  branch_cells,
  colnames(obj)
)

obj$step5C_v2_Arm1_branch_module_score[
  idx_score
] <- score_module(
  arm1_genes,
  branch_cells
)

obj$step5C_v2_Arm2_branch_module_score[
  idx_score
] <- score_module(
  arm2_genes,
  branch_cells
)

if (
  length(arm1_genes) > 0 ||
  length(arm2_genes) > 0
) {

  module_score_tbl <- tibble::tibble(
    cell =
      branch_cells,
    branch_arm =
      branch_arm_all[
        branch_cells
      ],
    timepoint =
      time_all[
        branch_cells
      ],
    sample =
      sample_all[
        branch_cells
      ],
    pseudotime =
      pt_all[
        branch_cells
      ],
    Arm1_module_score =
      obj$step5C_v2_Arm1_branch_module_score[
        idx_score
      ],
    Arm2_module_score =
      obj$step5C_v2_Arm2_branch_module_score[
        idx_score
      ]
  ) |>
    dplyr::mutate(
      pseudotime_scaled =
        (
          pseudotime -
            min(pseudotime)
        ) /
        (
          max(pseudotime) -
            min(pseudotime)
        ),
      pt_bin =
        cut(
          pseudotime_scaled,
          breaks =
            seq(
              0,
              1,
              length.out =
                N_PSEUDOTIME_BINS + 1
            ),
          include.lowest = TRUE,
          labels = FALSE
        )
    )

  module_bin_summary <- module_score_tbl |>
    dplyr::group_by(
      branch_arm,
      pt_bin
    ) |>
    dplyr::summarise(
      n_cells =
        dplyr::n(),
      pseudotime_scaled =
        median(
          pseudotime_scaled
        ),
      Arm1_module_score =
        mean(
          Arm1_module_score,
          na.rm = TRUE
        ),
      Arm2_module_score =
        mean(
          Arm2_module_score,
          na.rm = TRUE
        ),
      .groups = "drop"
    ) |>
    tidyr::pivot_longer(
      cols = c(
        Arm1_module_score,
        Arm2_module_score
      ),
      names_to = "module",
      values_to =
        "module_score"
    )

  module_7dpi_sample_summary <- module_score_tbl |>
    dplyr::filter(
      timepoint == "7 dpi"
    ) |>
    dplyr::group_by(
      sample,
      branch_arm
    ) |>
    dplyr::summarise(
      n_cells =
        dplyr::n(),
      Arm1_module_median =
        median(
          Arm1_module_score,
          na.rm = TRUE
        ),
      Arm2_module_median =
        median(
          Arm2_module_score,
          na.rm = TRUE
        ),
      .groups = "drop"
    )
}

# ==============================================================================
# 14. FIGURE 1 -- DISCOVERY / LRT DIAGNOSTICS
# ==============================================================================

msg(
  "Creating Fig5C_v2_01_branch_gene_discovery.pdf..."
)

graph_plot_tbl <- graph_test_tbl |>
  dplyr::mutate(
    neglog10_q =
      -log10(
        pmax(
          q_value,
          1e-300
        )
      )
  )

p1a <- ggplot2::ggplot(
  graph_plot_tbl,
  ggplot2::aes(
    x = morans_I,
    y = neglog10_q
  )
) +
  ggplot2::geom_point(
    size = 0.55,
    alpha = 0.35
  ) +
  ggplot2::geom_vline(
    xintercept =
      GRAPH_MORAN_MIN,
    linetype = 2,
    linewidth = 0.4
  ) +
  ggplot2::geom_hline(
    yintercept =
      -log10(
        GRAPH_Q_MAX
      ),
    linetype = 2,
    linewidth = 0.4
  ) +
  ggplot2::labs(
    title =
      "Principal-graph dynamic genes",
    x = "Moran's I",
    y = "-log10(q)"
  ) +
  ggplot2::theme_classic(
    base_size = 11
  )

lrt_plot_tbl <- branch_gene_tbl |>
  dplyr::filter(
    is.finite(
      branch_dependence_q_value
    ),
    is.finite(
      mean_diff_arm1_minus_arm2
    )
  ) |>
  dplyr::mutate(
    neglog10_branch_q =
      -log10(
        pmax(
          branch_dependence_q_value,
          1e-300
        )
      )
  )

p1b <- ggplot2::ggplot(
  lrt_plot_tbl,
  ggplot2::aes(
    x =
      mean_diff_arm1_minus_arm2,
    y =
      neglog10_branch_q
  )
) +
  ggplot2::geom_point(
    size = 0.65,
    alpha = 0.40
  ) +
  ggplot2::geom_hline(
    yintercept =
      -log10(
        BRANCH_DEP_Q_MAX
      ),
    linetype = 2,
    linewidth = 0.4
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = 0.35
  ) +
  ggplot2::labs(
    title =
      "BEAM-like branch dependence",
    subtitle =
      "Negative-binomial likelihood-ratio test: FULL versus common-trajectory NULL",
    x =
      paste0(
        "7 dpi mean expression: ",
        arm1,
        " - ",
        arm2
      ),
    y =
      "-log10(branch-dependence q)"
  ) +
  ggplot2::theme_classic(
    base_size = 11
  )

shape_plot_tbl <- branch_gene_tbl |>
  dplyr::filter(
    is.finite(
      shape_interaction_q_value
    )
  ) |>
  dplyr::mutate(
    neglog10_shape_q =
      -log10(
        pmax(
          shape_interaction_q_value,
          1e-300
        )
      )
  )

p1c <- ggplot2::ggplot(
  shape_plot_tbl,
  ggplot2::aes(
    x =
      branch_dependence_q_value,
    y =
      shape_interaction_q_value
  )
) +
  ggplot2::geom_point(
    size = 0.65,
    alpha = 0.35
  ) +
  ggplot2::geom_hline(
    yintercept =
      SHAPE_INTERACTION_Q_MAX,
    linetype = 2,
    linewidth = 0.4
  ) +
  ggplot2::geom_vline(
    xintercept =
      BRANCH_DEP_Q_MAX,
    linetype = 2,
    linewidth = 0.4
  ) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    title =
      "Primary branch dependence vs shape-only interaction",
    x =
      "BEAM-like branch-dependence q",
    y =
      "Interaction-only q"
  ) +
  ggplot2::theme_classic(
    base_size = 10
  )

gene_count_tbl <- tibble::tibble(
  category = c(
    "Graph-dynamic",
    "BEAM-like branch-dependent",
    "Shape-interaction",
    paste0(
      arm1,
      "-biased"
    ),
    paste0(
      arm2,
      "-biased"
    ),
    "Robust display"
  ),
  n_genes = c(
    length(
      graph_dynamic_genes
    ),
    nrow(
      branch_sig_tbl
    ),
    nrow(
      shape_sig_tbl
    ),
    length(
      arm1_genes
    ),
    length(
      arm2_genes
    ),
    nrow(
      robust_tbl
    )
  )
)

p1d <- ggplot2::ggplot(
  gene_count_tbl,
  ggplot2::aes(
    x = category,
    y = n_genes
  )
) +
  ggplot2::geom_col(
    width = 0.70
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = n_genes
    ),
    vjust = -0.4,
    size = 3.4
  ) +
  ggplot2::labs(
    title =
      "Corrected Step5C evidence cascade",
    x = NULL,
    y = "Genes"
  ) +
  ggplot2::theme_classic(
    base_size = 9
  ) +
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 25,
        hjust = 1
      )
  )

grDevices::pdf(
  OUT_FIG1,
  width = 14,
  height = 10,
  useDingbats = FALSE
)

print(
  (p1a | p1b) /
    (p1c | p1d) +
    patchwork::plot_annotation(
      title =
        paste0(
          "Step 5C v2 | ",
          PRIMARY_BRANCHPOINT,
          " corrected likelihood-based branch-gene discovery"
        )
    )
)

grDevices::dev.off()

# ==============================================================================
# 15. FIGURE 2 -- DYNAMIC HEATMAP
# ==============================================================================

msg(
  "Creating Fig5C_v2_02_branch_dynamic_heatmap.pdf..."
)

if (
  length(heatmap_genes) >= 2 &&
  nrow(smooth_tbl) > 0
) {

  heat_wide <- smooth_tbl |>
    dplyr::filter(
      gene_id %in%
        heatmap_genes
    ) |>
    dplyr::mutate(
      arm_bin =
        paste0(
          branch_arm,
          "__",
          sprintf(
            "%02d",
            pt_bin
          )
        )
    ) |>
    dplyr::select(
      gene_id,
      arm_bin,
      mean_expression
    ) |>
    tidyr::pivot_wider(
      names_from =
        arm_bin,
      values_from =
        mean_expression
    )

  heat_mat <- as.matrix(
    heat_wide[
      ,
      setdiff(
        colnames(
          heat_wide
        ),
        "gene_id"
      ),
      drop = FALSE
    ]
  )

  rownames(
    heat_mat
  ) <- heat_wide$gene_id

  heat_mat[
    is.na(
      heat_mat
    )
  ] <- 0

  heat_z <- row_zscore(
    heat_mat
  )

  heat_long <- as.data.frame(
    heat_z
  ) |>
    tibble::rownames_to_column(
      "gene_id"
    ) |>
    tidyr::pivot_longer(
      cols = -gene_id,
      names_to =
        "arm_bin",
      values_to =
        "z_expression"
    ) |>
    dplyr::mutate(
      branch_arm =
        ifelse(
          grepl(
            paste0(
              "^",
              arm1,
              "__"
            ),
            arm_bin
          ),
          arm1,
          arm2
        ),

      pt_bin =
        as.integer(
          sub(
            "^.*__([0-9]+)$",
            "\\1",
            arm_bin
          )
        ),

      gene_module =
        ifelse(
          gene_id %in%
            top_arm1_heat,
          paste0(
            arm1,
            "-biased"
          ),
          paste0(
            arm2,
            "-biased"
          )
        )
    )

  gene_order <- unique(
    c(
      rev(
        top_arm1_heat
      ),
      rev(
        top_arm2_heat
      )
    )
  )

  heat_long$gene_id <- factor(
    heat_long$gene_id,
    levels =
      gene_order
  )

  heat_long$branch_arm <- factor(
    heat_long$branch_arm,
    levels =
      MAJOR_ARMS
  )

  p2 <- ggplot2::ggplot(
    heat_long,
    ggplot2::aes(
      x = pt_bin,
      y = gene_id,
      fill = z_expression
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::facet_grid(
      gene_module ~ branch_arm,
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::scale_fill_gradient2(
      midpoint = 0,
      name =
        "Row z-score"
    ) +
    ggplot2::labs(
      title =
        paste0(
          PRIMARY_BRANCHPOINT,
          " branch-dependent gene dynamics"
        ),
      subtitle =
        "Genes selected by negative-binomial BEAM-like branch-dependence test",
      x =
        "Post-branch pseudotime bin",
      y =
        "Gene"
    ) +
    ggplot2::theme_classic(
      base_size = 8
    ) +
    ggplot2::theme(
      axis.text.y =
        ggplot2::element_text(
          size = 5.5
        )
    )

} else {

  p2 <- null_plot(
    paste0(
      "No genes passed the formal branch-dependence criterion (q <= ",
      BRANCH_DEP_Q_MAX,
      ").\nThis is retained as a valid null result."
    )
  )
}

grDevices::pdf(
  OUT_FIG2,
  width = 13,
  height = 14,
  useDingbats = FALSE
)

print(p2)

grDevices::dev.off()

# ==============================================================================
# 16. FIGURE 3 -- TOP GENE SMOOTHING
# ==============================================================================

msg(
  "Creating Fig5C_v2_03_top_branch_gene_smoothing.pdf..."
)

if (
  length(smooth_genes) >= 2 &&
  nrow(smooth_tbl) > 0
) {

  top_gene_plot_tbl <- smooth_tbl |>
    dplyr::filter(
      gene_id %in%
        smooth_genes
    ) |>
    dplyr::mutate(
      gene_module =
        ifelse(
          gene_id %in%
            top_arm1_smooth,
          paste0(
            arm1,
            "-biased"
          ),
          paste0(
            arm2,
            "-biased"
          )
        )
    )

  p3 <- ggplot2::ggplot(
    top_gene_plot_tbl,
    ggplot2::aes(
      x =
        pt_scaled_median,
      y =
        mean_expression,
      group =
        branch_arm,
      linetype =
        branch_arm
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.75
    ) +
    ggplot2::geom_point(
      size = 1.4
    ) +
    ggplot2::facet_wrap(
      gene_module ~ gene_id,
      scales = "free_y",
      ncol = 3
    ) +
    ggplot2::labs(
      title =
        "Top branch-dependent genes along retained Step5B-R pseudotime",
      subtitle =
        "Neutral Arm_1/Arm_2 labels; repair/scar scores were not used",
      x =
        "Scaled pseudotime",
      y =
        "Mean normalized expression",
      linetype =
        "Graph arm"
    ) +
    ggplot2::theme_classic(
      base_size = 9
    )

} else {

  p3 <- null_plot(
    "No formal branch-dependent genes available for smoothing."
  )
}

grDevices::pdf(
  OUT_FIG3,
  width = 13,
  height = 11,
  useDingbats = FALSE
)

print(p3)

grDevices::dev.off()

# ==============================================================================
# 17. FIGURE 4 -- MODULE / REPLICATE AUDIT
# ==============================================================================

msg(
  "Creating Fig5C_v2_04_branch_module_replicate_audit.pdf..."
)

if (nrow(module_bin_summary) > 0) {

  p4a <- ggplot2::ggplot(
    module_bin_summary,
    ggplot2::aes(
      x =
        pseudotime_scaled,
      y =
        module_score,
      group =
        branch_arm,
      linetype =
        branch_arm
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.9
    ) +
    ggplot2::geom_point(
      size = 1.6
    ) +
    ggplot2::facet_wrap(
      ~module,
      nrow = 1,
      scales = "free_y"
    ) +
    ggplot2::labs(
      title =
        "Branch-gene module dynamics",
      x =
        "Scaled retained pseudotime",
      y =
        "Mean module score",
      linetype =
        "Graph arm"
    ) +
    ggplot2::theme_classic(
      base_size = 10
    )

} else {

  p4a <- null_plot(
    "No formal branch module could be constructed."
  )
}

if (
  nrow(
    module_7dpi_sample_summary
  ) > 0
) {

  rep_long <- module_7dpi_sample_summary |>
    tidyr::pivot_longer(
      cols = c(
        Arm1_module_median,
        Arm2_module_median
      ),
      names_to = "module",
      values_to =
        "sample_median_score"
    )

  p4b <- ggplot2::ggplot(
    rep_long,
    ggplot2::aes(
      x = branch_arm,
      y =
        sample_median_score,
      shape = sample
    )
  ) +
    ggplot2::geom_point(
      size = 3
    ) +
    ggplot2::facet_wrap(
      ~module,
      nrow = 1,
      scales = "free_y"
    ) +
    ggplot2::labs(
      title =
        "7 dpi module scores by biological sample",
      x =
        "Corrected graph arm",
      y =
        "Sample-level median module score",
      shape =
        "7 dpi sample"
    ) +
    ggplot2::theme_classic(
      base_size = 10
    )

} else {

  p4b <- null_plot(
    "No 7 dpi branch-module replicate summary."
  )
}

top_effect_plot <- branch_gene_tbl |>
  dplyr::filter(
    branch_dependent
  ) |>
  dplyr::arrange(
    branch_dependence_q_value,
    dplyr::desc(
      abs_7dpi_mean_diff
    )
  ) |>
  dplyr::group_by(
    terminal_branch_bias
  ) |>
  dplyr::slice_head(
    n = 15
  ) |>
  dplyr::ungroup()

if (nrow(top_effect_plot) > 0) {

  p4c <- ggplot2::ggplot(
    top_effect_plot,
    ggplot2::aes(
      x =
        mean_diff_arm1_minus_arm2,
      y =
        reorder(
          gene_id,
          mean_diff_arm1_minus_arm2
        )
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_vline(
      xintercept = 0,
      linewidth = 0.35
    ) +
    ggplot2::labs(
      title =
        "Largest 7 dpi differences among formal branch genes",
      x =
        paste0(
          arm1,
          " - ",
          arm2,
          " mean expression"
        ),
      y = "Gene"
    ) +
    ggplot2::theme_classic(
      base_size = 9
    )

  top_gene_ids <- unique(
    as.character(
      top_effect_plot$gene_id
    )
  )

  rep_gene_top <- replicate_effect_tbl |>
    dplyr::filter(
      gene_id %in%
        top_gene_ids
    )

  p4d <- ggplot2::ggplot(
    rep_gene_top,
    ggplot2::aes(
      x = gene_id,
      y =
        mean_diff_arm1_minus_arm2,
      shape = sample
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.35
    ) +
    ggplot2::geom_point(
      size = 2.3
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title =
        "Top branch genes in each 7 dpi sample",
      x = "Gene",
      y =
        paste0(
          arm1,
          " - ",
          arm2
        ),
      shape =
        "7 dpi sample"
    ) +
    ggplot2::theme_classic(
      base_size = 8
    )

} else {

  p4c <- null_plot(
    "No formal branch-dependent genes."
  )

  p4d <- null_plot(
    "No replicate gene audit because no formal branch-dependent genes passed."
  )
}

grDevices::pdf(
  OUT_FIG4,
  width = 14,
  height = 10,
  useDingbats = FALSE
)

print(
  (p4a | p4b) /
    (p4c | p4d) +
    patchwork::plot_annotation(
      title =
        "Step 5C v2 | Branch-module dynamics and replicate audit"
    )
)

grDevices::dev.off()

# ==============================================================================
# 18. SAVE GENE SETS / OBJECT
# ==============================================================================

module_object <- list(
  primary_branchpoint =
    PRIMARY_BRANCHPOINT,
  major_arms =
    MAJOR_ARMS,
  retained_pseudotime =
    PSEUDOTIME_COL,
  formal_family =
    "negative binomial",
  full_formula =
    FULL_FORMULA,
  null_formula =
    NULL_FORMULA,
  additive_formula =
    ADDITIVE_FORMULA,
  graph_dynamic_genes =
    graph_dynamic_genes,
  branch_dependent_genes =
    branch_sig_tbl$gene_id,
  shape_interaction_genes =
    shape_sig_tbl$gene_id,
  Arm1_biased_genes =
    arm1_genes,
  Arm2_biased_genes =
    arm2_genes,
  robust_display_genes =
    robust_tbl$gene_id,
  branch_gene_table =
    branch_gene_tbl
)

saveRDS(
  module_object,
  OUT_MODULE_RDS,
  compress = TRUE
)

saveRDS(
  obj,
  OUT_SEURAT,
  compress = TRUE
)

writeLines(
  arm1_genes,
  file.path(
    OUTDIR,
    "Step5C_v2_Arm1_branch_dependent_genes.txt"
  )
)

writeLines(
  arm2_genes,
  file.path(
    OUTDIR,
    "Step5C_v2_Arm2_branch_dependent_genes.txt"
  )
)

# ==============================================================================
# 19. EXCEL
# ==============================================================================

msg(
  "Writing Table5C_v2_Y29_branch_dependent_genes.xlsx..."
)

parameter_tbl <- tibble::tibble(
  parameter = c(
    "Primary branchpoint",
    "Arm 1",
    "Arm 2",
    "Pseudotime used",
    "Branch-model cells",
    "Genes detected >=5%",
    "Graph-dynamic genes",
    "Genes entering NB models",
    "Model family",
    "Valid branch-dependence LRT p-values",
    "Minimum branch-dependence p",
    "Minimum branch-dependence q",
    "Raw branch-dependence p<0.05 genes",
    "Formal branch-dependent q<=0.05 genes",
    "Shape-interaction q<=0.05 genes",
    paste0(
      arm1,
      "-biased genes"
    ),
    paste0(
      arm2,
      "-biased genes"
    ),
    "Robust display genes",
    "R cores"
  ),
  value = c(
    PRIMARY_BRANCHPOINT,
    arm1,
    arm2,
    PSEUDOTIME_COL,
    length(
      branch_cells
    ),
    length(
      prefilter_genes
    ),
    length(
      graph_dynamic_genes
    ),
    length(
      fit_genes
    ),
    "negbinomial",
    nrow(
      branch_dep_valid
    ),
    min(
      branch_dep_valid$branch_dependence_p_value,
      na.rm = TRUE
    ),
    min(
      branch_dep_valid$branch_dependence_q_value,
      na.rm = TRUE
    ),
    sum(
      branch_dep_valid$branch_dependence_p_value <
        0.05,
      na.rm = TRUE
    ),
    nrow(
      branch_sig_tbl
    ),
    nrow(
      shape_sig_tbl
    ),
    length(
      arm1_genes
    ),
    length(
      arm2_genes
    ),
    nrow(
      robust_tbl
    ),
    CORES
  )
)

wb <- openxlsx::createWorkbook()

safe_add_sheet(
  wb,
  "Parameters",
  parameter_tbl
)

safe_add_sheet(
  wb,
  "Branch_cell_summary",
  branch_cell_summary
)

safe_add_sheet(
  wb,
  "Gene_prefilter",
  gene_prefilter_tbl
)

safe_add_sheet(
  wb,
  "Graph_test",
  graph_test_tbl
)

safe_add_sheet(
  wb,
  "Model_fit_status",
  fit_status_tbl
)

safe_add_sheet(
  wb,
  "Branch_dependence_LRT",
  branch_dep_tbl
)

safe_add_sheet(
  wb,
  "Shape_interaction_LRT",
  shape_tbl
)

safe_add_sheet(
  wb,
  "Branch_gene_master",
  branch_gene_tbl
)

safe_add_sheet(
  wb,
  "Formal_branch_genes",
  branch_sig_tbl
)

safe_add_sheet(
  wb,
  "Shape_interaction_genes",
  shape_sig_tbl
)

safe_add_sheet(
  wb,
  "Robust_display_genes",
  robust_tbl
)

safe_add_sheet(
  wb,
  "7dpi_replicate_effects",
  replicate_effect_tbl
)

safe_add_sheet(
  wb,
  "7dpi_direction_summary",
  replicate_direction_tbl
)

safe_add_sheet(
  wb,
  "Smoothing_binned",
  smooth_tbl
)

safe_add_sheet(
  wb,
  "Module_7dpi_samples",
  module_7dpi_sample_summary
)

openxlsx::saveWorkbook(
  wb,
  OUT_XLSX,
  overwrite = TRUE
)

# ==============================================================================
# 20. REPORT
# ==============================================================================

msg(
  "Writing Step5C_v2_report.txt..."
)

top_arm1_report <- branch_gene_tbl |>
  dplyr::filter(
    branch_dependent,
    terminal_branch_bias ==
      arm1
  ) |>
  dplyr::select(
    gene_id,
    graph_morans_I,
    graph_q_value,
    branch_dependence_q_value,
    shape_interaction_q_value,
    mean_diff_arm1_minus_arm2,
    direction_consistent
  ) |>
  dplyr::slice_head(
    n = 20
  )

top_arm2_report <- branch_gene_tbl |>
  dplyr::filter(
    branch_dependent,
    terminal_branch_bias ==
      arm2
  ) |>
  dplyr::select(
    gene_id,
    graph_morans_I,
    graph_q_value,
    branch_dependence_q_value,
    shape_interaction_q_value,
    mean_diff_arm1_minus_arm2,
    direction_consistent
  ) |>
  dplyr::slice_head(
    n = 20
  )

report_lines <- c(
  "======================================================================",
  "GSE162610 | STEP 5C v2 REPORT",
  "Corrected likelihood-based Y_29 branch-dependent gene analysis",
  "======================================================================",
  "",
  "WHY v2:",
  "v1 used quasipoisson fit_models followed by compare_models likelihood-ratio",
  "testing. Quasi-Poisson has no true likelihood, so v1's zero significant",
  "interaction genes is not treated as a biological negative result.",
  "",
  "v2 uses negative-binomial models for all likelihood-ratio comparisons.",
  "",
  "--------------------------- FIXED TOPOLOGY -----------------------------",
  paste0(
    "Corrected branchpoint: ",
    PRIMARY_BRANCHPOINT
  ),
  paste0(
    "Major arms: ",
    arm1,
    " vs ",
    arm2
  ),
  paste0(
    "Pseudotime used: ",
    PSEUDOTIME_COL
  ),
  "The failed Y_19 single-node pseudotime is not used.",
  "",
  "-------------------------- MODEL DESIGN -------------------------------",
  "PRIMARY BEAM-like branch-dependence test:",
  paste0(
    "  FULL: ",
    FULL_FORMULA
  ),
  paste0(
    "  NULL: ",
    NULL_FORMULA
  ),
  "",
  "SECONDARY shape-only interaction test:",
  paste0(
    "  FULL: ",
    FULL_FORMULA
  ),
  paste0(
    "  ADDITIVE: ",
    ADDITIVE_FORMULA
  ),
  "",
  "Expression family: negative binomial",
  paste0(
    "Windows model cores: ",
    CORES
  ),
  "",
  "---------------------------- CELLS ------------------------------------",
  paste0(
    "Branch-model cells: ",
    length(
      branch_cells
    )
  ),
  capture_table(
    branch_meta |>
      dplyr::count(
        branch_arm,
        timepoint,
        sample,
        name = "n_cells"
      )
  ),
  "",
  "-------------------------- graph_test ---------------------------------",
  paste0(
    "Genes detected >=5%: ",
    length(
      prefilter_genes
    )
  ),
  paste0(
    "Graph-dynamic genes: ",
    length(
      graph_dynamic_genes
    )
  ),
  paste0(
    "Genes entering negative-binomial models: ",
    length(
      fit_genes
    )
  ),
  "",
  "---------------------- LIKELIHOOD DIAGNOSTICS -------------------------",
  paste0(
    "Valid formal branch-dependence p-values: ",
    nrow(
      branch_dep_valid
    )
  ),
  paste0(
    "Minimum formal branch-dependence p-value: ",
    signif(
      min(
        branch_dep_valid$branch_dependence_p_value,
        na.rm = TRUE
      ),
      6
    )
  ),
  paste0(
    "Minimum formal branch-dependence q-value: ",
    signif(
      min(
        branch_dep_valid$branch_dependence_q_value,
        na.rm = TRUE
      ),
      6
    )
  ),
  paste0(
    "Raw p<0.05 branch-dependence genes: ",
    sum(
      branch_dep_valid$branch_dependence_p_value <
        0.05,
      na.rm = TRUE
    )
  ),
  "",
  "----------------------- FORMAL STEP5C RESULT --------------------------",
  paste0(
    "BEAM-like branch-dependent genes q<=",
    BRANCH_DEP_Q_MAX,
    ": ",
    nrow(
      branch_sig_tbl
    )
  ),
  paste0(
    "Secondary shape-interaction genes q<=",
    SHAPE_INTERACTION_Q_MAX,
    ": ",
    nrow(
      shape_sig_tbl
    )
  ),
  paste0(
    arm1,
    "-biased formal genes: ",
    length(
      arm1_genes
    )
  ),
  paste0(
    arm2,
    "-biased formal genes: ",
    length(
      arm2_genes
    )
  ),
  paste0(
    "Replicate-consistent + effect-filtered display genes: ",
    nrow(
      robust_tbl
    )
  ),
  "",
  paste0(
    "Top ",
    arm1,
    "-biased genes:"
  ),
  capture_table(
    top_arm1_report
  ),
  "",
  paste0(
    "Top ",
    arm2,
    "-biased genes:"
  ),
  capture_table(
    top_arm2_report
  ),
  "",
  "-------------------------- INTERPRETATION ------------------------------",
  "Formal Step5C branch-dependent genes are defined by the BEAM-like",
  "negative-binomial likelihood-ratio test (FULL vs NULL) after independent",
  "graph_test filtering.",
  "",
  "The stricter FULL vs ADDITIVE comparison asks specifically whether the",
  "Arm1/Arm2 curves differ in SHAPE beyond a constant branch offset.",
  "",
  "Therefore it is possible and biologically coherent to have many formal",
  "branch-dependent genes but fewer shape-interaction genes.",
  "",
  "If even the corrected negative-binomial BEAM-like test yields no significant",
  "genes, that null result should be taken seriously and Step5C should not be",
  "forced by relaxing q-value thresholds.",
  "",
  "Step5C remains part of Step 5. Step5D Monocle2 DDRTree sensitivity follows",
  "only after evaluating these corrected results.",
  "",
  "Formal outputs:",
  paste0("  ", OUT_FIG1),
  paste0("  ", OUT_FIG2),
  paste0("  ", OUT_FIG3),
  paste0("  ", OUT_FIG4),
  paste0("  ", OUT_XLSX),
  paste0("  ", OUT_MODULE_RDS),
  paste0("  ", OUT_SEURAT),
  paste0("  ", OUT_REP),
  paste0("  ", OUT_SESS)
)

writeLines(
  report_lines,
  OUT_REP,
  useBytes = TRUE
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  OUT_SESS,
  useBytes = TRUE
)

# ==============================================================================
# 21. FINISH
# ==============================================================================

msg("====================================================")
msg("STEP 5C v2 COMPLETE")
msg("====================================================")

msg(
  "Graph-dynamic genes: ",
  length(
    graph_dynamic_genes
  )
)

msg(
  "Formal branch-dependent genes: ",
  nrow(
    branch_sig_tbl
  )
)

msg(
  "Shape-interaction genes: ",
  nrow(
    shape_sig_tbl
  )
)

msg(
  arm1,
  "-biased genes: ",
  length(
    arm1_genes
  )
)

msg(
  arm2,
  "-biased genes: ",
  length(
    arm2_genes
  )
)

msg(
  "Robust display genes: ",
  nrow(
    robust_tbl
  )
)

msg("")
msg("请上传以下6个结果判断Step5C：")
msg("  1) Fig5C_v2_01_branch_gene_discovery.pdf")
msg("  2) Fig5C_v2_02_branch_dynamic_heatmap.pdf")
msg("  3) Fig5C_v2_03_top_branch_gene_smoothing.pdf")
msg("  4) Fig5C_v2_04_branch_module_replicate_audit.pdf")
msg("  5) Table5C_v2_Y29_branch_dependent_genes.xlsx")
msg("  6) Step5C_v2_report.txt")
