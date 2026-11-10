# ==============================================================================
# 49_step8D_GSE162610_mechanism_closure_QC_v2_METADATA_FIX_ONE_FILE.R
# ==============================================================================
# GSE162610 adult SCI activated/DAM-like microglia project
#
# STEP 8D | MECHANISM CLOSURE AND QC
#
# PURPOSE
#   A) Apply the pre-specified >=5% recipient receptor-positive threshold to the
#      completed Step8C recipient-program analysis and generate strict mechanism
#      tiers without overwriting Step8C.
#
#   B) Audit the proposed upstream regulon-to-ligand links:
#         Nfkb1 / Nfe2l2
#             ->
#         Tnf / Icam1 / Ccl2 / Ccl3
#
#      The audit combines:
#        1. Step7C motif-supported direct TF-target membership and target weight.
#        2. Step7D continuous regulon AUC.
#        3. Full-atlas normalized ligand expression.
#        4. Two separate 7 dpi biological samples.
#        5. Arm_1 versus Arm_2 state concordance.
#        6. Within-sample single-cell Spearman association.
#        7. A branch/sample-adjusted residual association.
#
# IMPORTANT INTERPRETATION BOUNDARY
#   - A motif-supported target plus AUC-expression association is an integrated
#     regulatory association, not perturbational proof of direct transcription.
#   - Cell-level correlations and P values are descriptive; the two 7 dpi
#     samples are used for direction replication, not for formal n=2 inference.
#   - Failure to identify a direct motif edge does not invalidate the TF as a
#     state regulator; it downgrades the TF-ligand connection to indirect/state
#     association.
#   - Arm_1 and Arm_2 remain neutral trajectory labels.
#
# RSTUDIO
# source("<DOWNLOADS>/49_step8D_GSE162610_mechanism_closure_QC_v2_METADATA_FIX_ONE_FILE.R")
# ==============================================================================

rm(list = ls())
gc()

options(
  stringsAsFactors = FALSE,
  width = 300,
  warn = 1,
  timeout = 86400,
  future.globals.maxSize = 32 * 1024^3
)

set.seed(20260803L)

cat(
  "R parser check PASS: the complete Step8D file was parsed by source().\n"
)

# ==============================================================================
# 0. FROZEN SETTINGS
# ==============================================================================

PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
RECEPTOR_POSITIVE_PERCENT_THRESHOLD <- 5
MIN_RECEPTOR_POSITIVE_CELLS <- 10L
MIN_RECEPTOR_NEGATIVE_CELLS <- 10L
MIN_TARGET_CELLS_PER_SAMPLE <- 20L

FOCUS_TFS <- c(
  "Nfkb1",
  "Nfe2l2"
)

FOCUS_LIGANDS <- c(
  "Tnf",
  "Icam1",
  "Ccl2",
  "Ccl3"
)

FORMAL_7DPI_STATES <- c(
  "Arm_1",
  "Arm_2"
)

EXPECTED_AUC_REGULONS <- 292L
EXPECTED_FROZEN_CELLS <- 4203L
EXPECTED_7DPI_ARM1_CELLS <- 1399L
EXPECTED_7DPI_ARM2_CELLS <- 1142L
EXPECTED_7DPI_SAMPLES <- 2L

MIN_ABS_SAMPLE_RHO_FOR_SUPPORT <- 0.05
MIN_ABS_ADJUSTED_RHO_FOR_SUPPORT <- 0.05

# ==============================================================================
# 1. INPUT PATHS
# ==============================================================================

ATLAS_RDS_FILE <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step1F_GSE162610_broad_annotation",
  "GSE162610_broad_annotation_audit.rds"
)

STEP7A_METADATA_FILE <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step7A_GSE162610_pySCENIC_input_environment_audit",
  "input",
  "GSE162610_step7A_frozen_4203_cell_metadata.csv"
)

STEP7C_DIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step7C_GSE162610_cisTarget_motif_pruning_regulon_inference"
)

STEP7C_TARGET_FILE <- file.path(
  STEP7C_DIR,
  "tables",
  "05_regulon_targets_long.tsv"
)

STEP7C_DECISION_FILE <- file.path(
  STEP7C_DIR,
  "tables",
  "08_final_decision.csv"
)

STEP7D_DIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step7D_GSE162610_AUCell_regulon_activity_scoring"
)

STEP7D_AUC_FILE <- file.path(
  STEP7D_DIR,
  "RDS",
  "GSE162610_step7D_AUC_regulons_by_cells.rds"
)

STEP7D_DECISION_FILE <- file.path(
  STEP7D_DIR,
  "tables",
  "08_final_decision.csv"
)

STEP7G_DIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step7G_GSE162610_regulon_target_program_closure"
)

STEP7G_CLOSURE_FILE <- file.path(
  STEP7G_DIR,
  "tables",
  "03_timing_candidate_target_closure.tsv"
)

STEP7G_DECISION_FILE <- file.path(
  STEP7G_DIR,
  "tables",
  "08_final_decision.csv"
)

STEP8A_DIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step8A_GSE162610_CellChat_input_preflight"
)

STEP8A_MANIFEST_FILE <- file.path(
  STEP8A_DIR,
  "RDS",
  "GSE162610_step8A_CellChat_input_manifest.rds"
)

STEP8A_CANDIDATE_LIGAND_FILE <- file.path(
  STEP8A_DIR,
  "tables",
  "06_candidate_regulon_ligand_catalog.tsv"
)

STEP8A_DECISION_FILE <- file.path(
  STEP8A_DIR,
  "tables",
  "08_final_decision.csv"
)

STEP8C_DIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step8C_GSE162610_recipient_program_integration"
)

STEP8C_RECEPTOR_AUDIT_FILE <- file.path(
  STEP8C_DIR,
  "tables",
  "05_7dpi_axis_program_receptor_audit_long.tsv"
)

STEP8C_AXIS_PROGRAM_FILE <- file.path(
  STEP8C_DIR,
  "tables",
  "06_axis_program_integrated_evidence.tsv"
)

STEP8C_AXIS_FILE <- file.path(
  STEP8C_DIR,
  "tables",
  "07_prioritized_mechanism_axes.tsv"
)

STEP8C_DECISION_FILE <- file.path(
  STEP8C_DIR,
  "tables",
  "09_final_decision.csv"
)

# ==============================================================================
# 2. OUTPUT PATHS
# ==============================================================================

OUTDIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step8D_GSE162610_mechanism_closure_QC"
)

TABDIR <- file.path(
  OUTDIR,
  "tables"
)

FIGDIR <- file.path(
  OUTDIR,
  "figures"
)

RDSDIR <- file.path(
  OUTDIR,
  "RDS"
)

REPORTDIR <- file.path(
  OUTDIR,
  "reports"
)

for (directory in c(
  OUTDIR,
  TABDIR,
  FIGDIR,
  RDSDIR,
  REPORTDIR
)) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

INPUT_AUDIT_FILE <- file.path(
  TABDIR,
  "00_input_file_audit.csv"
)

STRICT_RECEPTOR_AUDIT_FILE <- file.path(
  TABDIR,
  "01_receptor_5pct_strict_sample_audit.tsv"
)

STRICT_AXIS_PROGRAM_FILE <- file.path(
  TABDIR,
  "02_axis_program_strict_evidence.tsv"
)

STRICT_MECHANISM_AXES_FILE <- file.path(
  TABDIR,
  "03_strict_prioritized_mechanism_axes.tsv"
)

TF_LIGAND_MOTIF_AUDIT_FILE <- file.path(
  TABDIR,
  "04_focus_TF_ligand_motif_target_audit.tsv"
)

TF_LIGAND_CORRELATION_LONG_FILE <- file.path(
  TABDIR,
  "05_focus_TF_ligand_two_sample_correlation_long.tsv"
)

TF_LIGAND_CLOSURE_SUMMARY_FILE <- file.path(
  TABDIR,
  "06_focus_TF_ligand_closure_summary.tsv"
)

FINAL_MECHANISM_MODEL_FILE <- file.path(
  TABDIR,
  "07_strict_final_mechanism_model.tsv"
)

FINAL_DECISION_FILE <- file.path(
  TABDIR,
  "08_final_decision.csv"
)

RESULT_RDS_FILE <- file.path(
  RDSDIR,
  "GSE162610_step8D_mechanism_closure_QC_results.rds"
)

REPORT_FILE <- file.path(
  REPORTDIR,
  "Step8D_mechanism_closure_QC_report.txt"
)

ERROR_FILE <- file.path(
  REPORTDIR,
  "Step8D_mechanism_closure_QC_ERROR.txt"
)

FIGURE_STRICT_TIERS_FILE <- file.path(
  FIGDIR,
  "Fig8D_01_strict_mechanism_tiers.pdf"
)

FIGURE_TF_LIGAND_FILE <- file.path(
  FIGDIR,
  "Fig8D_02_TF_ligand_association_summary.pdf"
)

# ==============================================================================
# 3. PACKAGES
# ==============================================================================

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "data.table",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {
  stop(
    paste0(
      "Missing required package(s): ",
      paste(
        missing_packages,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

# ==============================================================================
# 4. HELPERS
# ==============================================================================

timestamp <- function() {
  format(
    Sys.time(),
    "%Y-%m-%d %H:%M:%S"
  )
}

msg <- function(...) {
  cat(
    sprintf(
      "[%s] %s\n",
      timestamp(),
      paste0(
        ...,
        collapse = ""
      )
    )
  )
  flush.console()
}

safe_path <- function(path) {
  normalizePath(
    path,
    winslash = "/",
    mustWork = FALSE
  )
}

write_csv_utf8 <- function(data, path) {
  data.table::fwrite(
    data,
    file = path,
    sep = ",",
    quote = TRUE,
    na = "NA"
  )
}

write_tsv <- function(data, path) {
  data.table::fwrite(
    data,
    file = path,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
}

abort_step <- function(text) {
  writeLines(
    c(
      "======================================================================",
      "GSE162610 | STEP 8D ERROR",
      "======================================================================",
      paste0(
        "Time: ",
        timestamp()
      ),
      "",
      text
    ),
    con = ERROR_FILE,
    useBytes = TRUE
  )

  stop(
    text,
    call. = FALSE
  )
}

as_logical_safe <- function(x) {
  if (is.logical(x)) {
    return(x)
  }

  normalized <- toupper(
    trimws(
      as.character(
        x
      )
    )
  )

  output <- rep(
    NA,
    length(
      normalized
    )
  )

  output[
    normalized %in%
      c(
        "TRUE",
        "T",
        "1",
        "PASS",
        "YES"
      )
  ] <- TRUE

  output[
    normalized %in%
      c(
        "FALSE",
        "F",
        "0",
        "FAIL",
        "NO"
      )
  ] <- FALSE

  output
}

decision_pass <- function(
  decision_table,
  domain_name
) {
  if (
    is.null(
      decision_table
    ) ||
    !"domain" %in%
      colnames(
        decision_table
      )
  ) {
    return(FALSE)
  }

  index <- which(
    as.character(
      decision_table$domain
    ) ==
      domain_name
  )

  if (length(index) == 0L) {
    return(FALSE)
  }

  if (
    "pass" %in%
      colnames(
        decision_table
      )
  ) {
    value <- as_logical_safe(
      decision_table$pass[
        index[[1L]]
      ]
    )

    return(
      isTRUE(
        value[[1L]]
      )
    )
  }

  if (
    "status" %in%
      colnames(
        decision_table
      )
  ) {
    return(
      identical(
        toupper(
          as.character(
            decision_table$status[
              index[[1L]]
            ]
          )
        ),
        "PASS"
      )
    )
  }

  FALSE
}

add_decision_row <- function(
  decision_list,
  domain,
  pass,
  detail = "",
  status_if_fail = "FAIL"
) {
  decision_list[[length( decision_list ) + 1L]] <- data.frame(
    domain = domain,
    pass = isTRUE(
      pass
    ),
    status = if (
      isTRUE(
        pass
      )
    ) {
      "PASS"
    } else {
      status_if_fail
    },
    detail = as.character(
      detail
    ),
    stringsAsFactors = FALSE
  )

  decision_list
}

collapse_unique <- function(x) {
  values <- unique(
    as.character(
      x
    )
  )

  values <- values[
    !is.na(
      values
    ) &
      nzchar(
        values
      )
  ]

  if (length(values) == 0L) {
    return(
      NA_character_
    )
  }

  paste(
    sort(
      values
    ),
    collapse = " | "
  )
}

safe_mean_finite <- function(x) {
  values <- x[
    is.finite(
      x
    )
  ]

  if (length(values) == 0L) {
    return(
      NA_real_
    )
  }

  mean(
    values
  )
}

safe_min_finite <- function(x) {
  values <- x[
    is.finite(
      x
    )
  ]

  if (length(values) == 0L) {
    return(
      NA_real_
    )
  }

  min(
    values
  )
}

safe_max_finite <- function(x) {
  values <- x[
    is.finite(
      x
    )
  ]

  if (length(values) == 0L) {
    return(
      NA_real_
    )
  }

  max(
    values
  )
}

safe_spearman <- function(
  x,
  y
) {
  valid <- is.finite(
    x
  ) &
    is.finite(
      y
    )

  x <- x[
    valid
  ]

  y <- y[
    valid
  ]

  if (
    length(
      x
    ) <
      10L ||
    length(
      unique(
        x
      )
    ) <
      2L ||
    length(
      unique(
        y
      )
    ) <
      2L
  ) {
    return(
      c(
        rho = NA_real_,
        p_value = NA_real_,
        n = length(
          x
        )
      )
    )
  }

  result <- tryCatch(
    suppressWarnings(
      stats::cor.test(
        x,
        y,
        method = "spearman",
        exact = FALSE
      )
    ),
    error = function(e) {
      NULL
    }
  )

  if (is.null(result)) {
    return(
      c(
        rho = NA_real_,
        p_value = NA_real_,
        n = length(
          x
        )
      )
    )
  }

  c(
    rho = unname(
      result$estimate
    ),
    p_value = result$p.value,
    n = length(
      x
    )
  )
}

safe_lm_residuals <- function(
  response,
  state,
  sample_id
) {
  input <- data.frame(
    response = as.numeric(
      response
    ),
    state = factor(
      state,
      levels = FORMAL_7DPI_STATES
    ),
    sample_id = factor(
      sample_id
    ),
    stringsAsFactors = FALSE
  )

  valid <- is.finite(
    input$response
  ) &
    !is.na(
      input$state
    ) &
    !is.na(
      input$sample_id
    )

  output <- rep(
    NA_real_,
    nrow(
      input
    )
  )

  if (
    sum(
      valid
    ) <
      20L
  ) {
    return(
      output
    )
  }

  fit <- tryCatch(
    stats::lm(
      response ~ state + sample_id,
      data = input[
        valid,
        ,
        drop = FALSE
      ]
    ),
    error = function(e) {
      NULL
    }
  )

  if (is.null(fit)) {
    return(
      output
    )
  }

  output[
    valid
  ] <- stats::residuals(
    fit
  )

  output
}

get_rna_data <- function(object) {
  output <- tryCatch(
    SeuratObject::LayerData(
      object = object,
      assay = "RNA",
      layer = "data"
    ),
    error = function(e) {
      NULL
    }
  )

  if (!is.null(output)) {
    return(
      output
    )
  }

  output <- tryCatch(
    Seurat::GetAssayData(
      object = object,
      assay = "RNA",
      slot = "data"
    ),
    error = function(e) {
      NULL
    }
  )

  if (is.null(output)) {
    stop(
      "Unable to retrieve normalized RNA data.",
      call. = FALSE
    )
  }

  output
}

standardize_state <- function(x) {
  value <- trimws(
    as.character(
      x
    )
  )

  normalized <- tolower(
    gsub(
      "[^a-z0-9]+",
      "",
      value
    )
  )

  output <- rep(
    NA_character_,
    length(
      normalized
    )
  )

  output[
    grepl(
      "arm1",
      normalized,
      fixed = TRUE
    ) |
    normalized %in%
      c(
        "arm01"
      )
  ] <- "Arm_1"

  output[
    grepl(
      "arm2",
      normalized,
      fixed = TRUE
    ) |
    normalized %in%
      c(
        "arm02"
      )
  ] <- "Arm_2"

  output[
    grepl(
      "sharedprebranch",
      normalized,
      fixed = TRUE
    ) |
    normalized %in%
      c(
        "shared",
        "prebranch"
      )
  ] <- "Shared_prebranch"

  output
}

standardize_timepoint <- function(
  x,
  sample_id = NULL
) {
  value <- trimws(
    as.character(
      x
    )
  )

  normalized <- tolower(
    gsub(
      "[^a-z0-9]+",
      "",
      value
    )
  )

  output <- rep(
    NA_character_,
    length(
      normalized
    )
  )

  output[
    normalized %in%
      c(
        "0",
        "0dpi",
        "day0",
        "d0",
        "uninjured",
        "naive",
        "control",
        "sham"
      ) |
    grepl(
      "uninjured",
      normalized,
      fixed = TRUE
    )
  ] <- "uninjured"

  output[
    normalized %in%
      c(
        "1",
        "1dpi",
        "day1",
        "d1",
        "24h",
        "24hr",
        "24hrs"
      ) |
    grepl(
      "(^|[^0-9])1dpi",
      normalized
    )
  ] <- "1dpi"

  output[
    normalized %in%
      c(
        "3",
        "3dpi",
        "day3",
        "d3",
        "72h",
        "72hr",
        "72hrs"
      ) |
    grepl(
      "(^|[^0-9])3dpi",
      normalized
    )
  ] <- "3dpi"

  output[
    normalized %in%
      c(
        "7",
        "7dpi",
        "day7",
        "d7",
        "1w",
        "1week"
      ) |
    grepl(
      "7dpi",
      normalized,
      fixed = TRUE
    ) |
    grepl(
      "day7",
      normalized,
      fixed = TRUE
    )
  ] <- "7dpi"

  if (!is.null(sample_id)) {
    sample_normalized <- tolower(
      gsub(
        "[^a-z0-9]+",
        "",
        as.character(
          sample_id
        )
      )
    )

    unresolved <- is.na(
      output
    )

    output[
      unresolved &
      grepl(
        "7dpi",
        sample_normalized,
        fixed = TRUE
      )
    ] <- "7dpi"

    output[
      unresolved &
      grepl(
        "3dpi",
        sample_normalized,
        fixed = TRUE
      )
    ] <- "3dpi"

    output[
      unresolved &
      grepl(
        "1dpi",
        sample_normalized,
        fixed = TRUE
      )
    ] <- "1dpi"

    output[
      unresolved &
      grepl(
        "uninjured",
        sample_normalized,
        fixed = TRUE
      )
    ] <- "uninjured"
  }

  output
}

state_from_cellchat_group <- function(x) {
  value <- as.character(
    x
  )

  output <- rep(
    NA_character_,
    length(
      value
    )
  )

  output[
    value ==
      "Microglia_Arm_1"
  ] <- "Arm_1"

  output[
    value ==
      "Microglia_Arm_2"
  ] <- "Arm_2"

  output[
    value ==
      "Microglia_Shared_prebranch"
  ] <- "Shared_prebranch"

  output
}

require_columns <- function(
  data,
  columns,
  data_name
) {
  missing_columns <- setdiff(
    columns,
    colnames(
      data
    )
  )

  if (length(missing_columns) > 0L) {
    stop(
      paste0(
        data_name,
        " is missing required columns: ",
        paste(
          missing_columns,
          collapse = ", "
        ),
        ". Current columns: ",
        paste(
          colnames(
            data
          ),
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }
}

resolve_auc_matrix <- function(object) {
  if (
    is.matrix(
      object
    ) ||
    inherits(
      object,
      "Matrix"
    )
  ) {
    return(
      as.matrix(
        object
      )
    )
  }

  if (is.data.frame(object)) {
    matrix_object <- as.matrix(
      object
    )

    storage.mode(
      matrix_object
    ) <- "numeric"

    return(
      matrix_object
    )
  }

  if (is.list(object)) {
    candidate_names <- c(
      "auc_regulons_by_cells",
      "AUC_regulons_by_cells",
      "auc_matrix",
      "AUC",
      "auc"
    )

    for (candidate_name in candidate_names) {
      if (
        candidate_name %in%
          names(
            object
          )
      ) {
        candidate <- object[[candidate_name]]

        if (
          is.matrix(
            candidate
          ) ||
          inherits(
            candidate,
            "Matrix"
          ) ||
          is.data.frame(
            candidate
          )
        ) {
          return(
            as.matrix(
              candidate
            )
          )
        }
      }
    }
  }

  stop(
    "Unable to resolve a regulons-by-cells AUC matrix from Step7D RDS.",
    call. = FALSE
  )
}

# ==============================================================================
# 5. MAIN
# ==============================================================================

main <- function() {
  start_time <- Sys.time()

  msg(
    "=============================================================="
  )
  msg(
    "GSE162610 Step 8D v2 | mechanism closure and QC"
  )
  msg(
    "=============================================================="
  )

  required_files <- c(
    ATLAS_RDS_FILE,
    STEP7A_METADATA_FILE,
    STEP7C_TARGET_FILE,
    STEP7C_DECISION_FILE,
    STEP7D_AUC_FILE,
    STEP7D_DECISION_FILE,
    STEP7G_CLOSURE_FILE,
    STEP7G_DECISION_FILE,
    STEP8A_MANIFEST_FILE,
    STEP8A_CANDIDATE_LIGAND_FILE,
    STEP8A_DECISION_FILE,
    STEP8C_RECEPTOR_AUDIT_FILE,
    STEP8C_AXIS_PROGRAM_FILE,
    STEP8C_AXIS_FILE,
    STEP8C_DECISION_FILE
  )

  input_audit <- data.frame(
    role = c(
      "Full normalized atlas",
      "Step7A frozen branch metadata",
      "Step7C motif-supported targets",
      "Step7C final decision",
      "Step7D continuous regulon AUC",
      "Step7D final decision",
      "Step7G candidate closure table",
      "Step7G final decision",
      "Step8A authoritative CellChat manifest",
      "Step8A candidate regulon-ligand catalog",
      "Step8A final decision",
      "Step8C sample receptor audit",
      "Step8C axis-program evidence",
      "Step8C prioritized mechanism axes",
      "Step8C final decision"
    ),
    path = safe_path(
      required_files
    ),
    exists = file.exists(
      required_files
    ),
    size_bytes = ifelse(
      file.exists(
        required_files
      ),
      file.info(
        required_files
      )$size,
      NA_real_
    ),
    stringsAsFactors = FALSE
  )

  write_csv_utf8(
    input_audit,
    INPUT_AUDIT_FILE
  )

  missing_files <- required_files[
    !file.exists(
      required_files
    )
  ]

  if (length(missing_files) > 0L) {
    abort_step(
      paste0(
        "Missing required file(s):\n",
        paste(
          safe_path(
            missing_files
          ),
          collapse = "\n"
        )
      )
    )
  }

  step7c_decision <- data.table::fread(
    STEP7C_DECISION_FILE,
    data.table = FALSE,
    showProgress = FALSE
  )

  step7d_decision <- data.table::fread(
    STEP7D_DECISION_FILE,
    data.table = FALSE,
    showProgress = FALSE
  )

  step7g_decision <- data.table::fread(
    STEP7G_DECISION_FILE,
    data.table = FALSE,
    showProgress = FALSE
  )

  step8a_decision <- data.table::fread(
    STEP8A_DECISION_FILE,
    data.table = FALSE,
    showProgress = FALSE
  )

  step8c_decision <- data.table::fread(
    STEP8C_DECISION_FILE,
    data.table = FALSE,
    showProgress = FALSE
  )

  if (
    !decision_pass(
      step7c_decision,
      "Overall Step 7C"
    )
  ) {
    abort_step(
      "Step7C Overall PASS was not found."
    )
  }

  if (
    !decision_pass(
      step7d_decision,
      "Overall Step 7D"
    )
  ) {
    abort_step(
      "Step7D Overall PASS was not found."
    )
  }

  if (
    !decision_pass(
      step7g_decision,
      "Overall Step 7G"
    )
  ) {
    abort_step(
      "Step7G Overall PASS was not found."
    )
  }

  if (
    !decision_pass(
      step8a_decision,
      "Overall Step 8A"
    )
  ) {
    abort_step(
      "Step8A Overall PASS was not found."
    )
  }

  if (
    !decision_pass(
      step8c_decision,
      "Overall Step 8C technical"
    )
  ) {
    abort_step(
      "Step8C technical PASS was not found."
    )
  }

  # ============================================================================
  # PART A. STRICT >=5% RECEPTOR-POSITIVE QC
  # ============================================================================

  msg(
    "Applying the strict >=5% receptor-positive threshold..."
  )

  receptor_audit <- data.table::fread(
    STEP8C_RECEPTOR_AUDIT_FILE,
    data.table = FALSE,
    showProgress = FALSE
  )

  axis_program_original <- data.table::fread(
    STEP8C_AXIS_PROGRAM_FILE,
    data.table = FALSE,
    showProgress = FALSE
  )

  axis_original <- data.table::fread(
    STEP8C_AXIS_FILE,
    data.table = FALSE,
    showProgress = FALSE
  )

  required_receptor_columns <- c(
    "sample_id",
    "target",
    "ligand",
    "receptor",
    "interaction_name",
    "pathway_name",
    "dominant_branch",
    "response_program",
    "program_class",
    "expected_direction",
    "expected_sign",
    "n_target_cells",
    "n_receptor_positive",
    "n_receptor_negative",
    "receptor_positive_percent",
    "effect_receptor_positive_minus_negative",
    "effect_in_expected_direction",
    "spearman_rho_receptor_vs_program",
    "evaluable",
    "expected_effect_match",
    "expected_correlation_match",
    "sample_support"
  )

  require_columns(
    receptor_audit,
    required_receptor_columns,
    "Step8C receptor audit"
  )

  receptor_audit$evaluable <- as_logical_safe(
    receptor_audit$evaluable
  )

  receptor_audit$expected_effect_match <- as_logical_safe(
    receptor_audit$expected_effect_match
  )

  receptor_audit$expected_correlation_match <- as_logical_safe(
    receptor_audit$expected_correlation_match
  )

  receptor_audit$sample_support <- as_logical_safe(
    receptor_audit$sample_support
  )

  receptor_audit$receptor_positive_5pct_pass <-
    is.finite(
      receptor_audit$receptor_positive_percent
    ) &
    receptor_audit$receptor_positive_percent >=
      RECEPTOR_POSITIVE_PERCENT_THRESHOLD

  receptor_audit$cell_count_strict_pass <-
    receptor_audit$n_target_cells >=
      MIN_TARGET_CELLS_PER_SAMPLE &
    receptor_audit$n_receptor_positive >=
      MIN_RECEPTOR_POSITIVE_CELLS &
    receptor_audit$n_receptor_negative >=
      MIN_RECEPTOR_NEGATIVE_CELLS

  receptor_audit$evaluable_strict <-
    receptor_audit$evaluable %in%
      TRUE &
    receptor_audit$receptor_positive_5pct_pass &
    receptor_audit$cell_count_strict_pass

  receptor_audit$sample_support_strict <-
    receptor_audit$sample_support %in%
      TRUE &
    receptor_audit$evaluable_strict

  receptor_audit$strict_qc_reason <- ifelse(
    !receptor_audit$receptor_positive_5pct_pass,
    "RECEPTOR_POSITIVE_PERCENT_BELOW_5",
    ifelse(
      !receptor_audit$cell_count_strict_pass,
      "INSUFFICIENT_POSITIVE_NEGATIVE_OR_TOTAL_CELLS",
      ifelse(
        receptor_audit$sample_support_strict,
        "STRICT_SAMPLE_SUPPORT",
        "STRICTLY_EVALUABLE_BUT_PROGRAM_NOT_SUPPORTED"
      )
    )
  )

  write_tsv(
    receptor_audit,
    STRICT_RECEPTOR_AUDIT_FILE
  )

  strict_summary <- data.table::as.data.table(
    receptor_audit
  )[
    ,
    .(
      n_samples = data.table::uniqueN(
        sample_id
      ),
      n_samples_evaluable_strict = sum(
        evaluable_strict,
        na.rm = TRUE
      ),
      n_samples_receptor_5pct_pass = sum(
        receptor_positive_5pct_pass,
        na.rm = TRUE
      ),
      n_samples_effect_direction_match = sum(
        expected_effect_match,
        na.rm = TRUE
      ),
      n_samples_correlation_direction_match = sum(
        expected_correlation_match,
        na.rm = TRUE
      ),
      n_samples_supported_strict = sum(
        sample_support_strict,
        na.rm = TRUE
      ),
      minimum_receptor_positive_percent_strict =
        safe_min_finite(
          receptor_positive_percent
        ),
      mean_effect_in_expected_direction_strict =
        safe_mean_finite(
          effect_in_expected_direction
        ),
      mean_spearman_rho_strict =
        safe_mean_finite(
          spearman_rho_receptor_vs_program
        ),
      strict_sample_details = paste(
        paste0(
          sample_id,
          ":pct=",
          round(
            receptor_positive_percent,
            2
          ),
          ",support=",
          sample_support_strict
        ),
        collapse = " | "
      )
    ),
    by = .(
      target,
      ligand,
      receptor,
      interaction_name,
      pathway_name,
      dominant_branch,
      response_program,
      program_class,
      expected_direction,
      expected_sign
    )
  ]

  strict_summary <- as.data.frame(
    strict_summary,
    stringsAsFactors = FALSE
  )

  strict_summary$recipient_program_status_strict <- ifelse(
    strict_summary$n_samples_evaluable_strict ==
      EXPECTED_7DPI_SAMPLES &
    strict_summary$n_samples_supported_strict ==
      EXPECTED_7DPI_SAMPLES,
    "SUPPORTED_IN_BOTH_7DPI_SAMPLES_STRICT",
    ifelse(
      strict_summary$n_samples_supported_strict >
        0L,
      "SUPPORTED_IN_ONE_7DPI_SAMPLE_STRICT",
      ifelse(
        strict_summary$n_samples_effect_direction_match >
          0L &
        strict_summary$n_samples_receptor_5pct_pass >
          0L,
        "DIRECTION_ONLY_STRICT",
        ifelse(
          strict_summary$n_samples_evaluable_strict ==
            0L,
          "NOT_EVALUABLE_AFTER_5PCT_QC",
          "NOT_SUPPORTED_STRICT"
        )
      )
    )
  )

  axis_program_keys <- c(
    "target",
    "ligand",
    "receptor",
    "interaction_name",
    "pathway_name",
    "dominant_branch",
    "response_program",
    "program_class",
    "expected_direction",
    "expected_sign"
  )

  strict_axis_program <- merge(
    axis_program_original,
    strict_summary,
    by = axis_program_keys,
    all.x = TRUE,
    sort = FALSE
  )

  strict_axis_program$program_evidence_tier_strict <- ifelse(
    strict_axis_program$recipient_program_status_strict ==
      "SUPPORTED_IN_BOTH_7DPI_SAMPLES_STRICT",
    "PROGRAM_TIER_1_TWO_SAMPLE_STRICT",
    ifelse(
      strict_axis_program$recipient_program_status_strict ==
        "SUPPORTED_IN_ONE_7DPI_SAMPLE_STRICT",
      "PROGRAM_TIER_2_ONE_SAMPLE_STRICT",
      ifelse(
        strict_axis_program$recipient_program_status_strict ==
          "DIRECTION_ONLY_STRICT",
        "PROGRAM_TIER_3_DIRECTION_ONLY_STRICT",
        "PROGRAM_TIER_4_NO_SUPPORT_STRICT"
      )
    )
  )

  strict_axis_program <- strict_axis_program[
    order(
      match(
        strict_axis_program$program_evidence_tier_strict,
        c(
          "PROGRAM_TIER_1_TWO_SAMPLE_STRICT",
          "PROGRAM_TIER_2_ONE_SAMPLE_STRICT",
          "PROGRAM_TIER_3_DIRECTION_ONLY_STRICT",
          "PROGRAM_TIER_4_NO_SUPPORT_STRICT"
        )
      ),
      -strict_axis_program$mean_effect_in_expected_direction_strict
    ),
    ,
    drop = FALSE
  ]

  write_tsv(
    strict_axis_program,
    STRICT_AXIS_PROGRAM_FILE
  )

  strict_axis_level <- data.table::as.data.table(
    strict_axis_program
  )[
    ,
    .(
      n_programs_tested = .N,
      n_programs_two_sample_strict = sum(
        recipient_program_status_strict ==
          "SUPPORTED_IN_BOTH_7DPI_SAMPLES_STRICT",
        na.rm = TRUE
      ),
      n_programs_one_sample_strict = sum(
        recipient_program_status_strict ==
          "SUPPORTED_IN_ONE_7DPI_SAMPLE_STRICT",
        na.rm = TRUE
      ),
      n_programs_direction_only_strict = sum(
        recipient_program_status_strict ==
          "DIRECTION_ONLY_STRICT",
        na.rm = TRUE
      ),
      n_programs_not_evaluable_after_5pct = sum(
        recipient_program_status_strict ==
          "NOT_EVALUABLE_AFTER_5PCT_QC",
        na.rm = TRUE
      ),
      minimum_receptor_positive_percent =
        safe_min_finite(
          minimum_receptor_positive_percent_strict
        ),
      top_effect_in_expected_direction_strict =
        safe_max_finite(
          mean_effect_in_expected_direction_strict
        ),
      supported_programs_two_sample_strict = collapse_unique(
        response_program[
          recipient_program_status_strict ==
            "SUPPORTED_IN_BOTH_7DPI_SAMPLES_STRICT"
        ]
      ),
      supported_programs_any_strict = collapse_unique(
        response_program[
          recipient_program_status_strict %in%
            c(
              "SUPPORTED_IN_BOTH_7DPI_SAMPLES_STRICT",
              "SUPPORTED_IN_ONE_7DPI_SAMPLE_STRICT"
            )
        ]
      ),
      Step7G_candidate_ligand = any(
        as_logical_safe(
          Step7G_candidate_ligand
        ),
        na.rm = TRUE
      ),
      Step7G_candidate_TFs = collapse_unique(
        Step7G_candidate_TFs
      ),
      mechanistic_priority_CellChat = collapse_unique(
        mechanistic_priority_CellChat
      ),
      original_program_tiers = collapse_unique(
        program_evidence_tier
      )
    ),
    by = .(
      target,
      ligand,
      receptor,
      interaction_name,
      pathway_name,
      dominant_branch
    )
  ]

  strict_axis_level <- as.data.frame(
    strict_axis_level,
    stringsAsFactors = FALSE
  )

  strict_axis_level$recipient_mechanism_status_strict <- ifelse(
    strict_axis_level$n_programs_two_sample_strict >
      0L,
    "RECIPIENT_PROGRAM_CONFIRMED_TWO_SAMPLE_STRICT",
    ifelse(
      strict_axis_level$n_programs_one_sample_strict >
        0L,
      "RECIPIENT_PROGRAM_PARTIAL_STRICT",
      ifelse(
        strict_axis_level$n_programs_direction_only_strict >
          0L,
        "RECIPIENT_PROGRAM_DIRECTION_ONLY_STRICT",
        "COMMUNICATION_ONLY_AFTER_STRICT_QC"
      )
    )
  )

  strict_axis_level$final_mechanism_tier_strict <- ifelse(
    strict_axis_level$mechanistic_priority_CellChat ==
      "TIER_1_CORE" &
    strict_axis_level$recipient_mechanism_status_strict ==
      "RECIPIENT_PROGRAM_CONFIRMED_TWO_SAMPLE_STRICT",
    "MECHANISM_TIER_1_STRICT",
    ifelse(
      strict_axis_level$recipient_mechanism_status_strict %in%
        c(
          "RECIPIENT_PROGRAM_CONFIRMED_TWO_SAMPLE_STRICT",
          "RECIPIENT_PROGRAM_PARTIAL_STRICT"
        ),
      "MECHANISM_TIER_2_STRICT",
      ifelse(
        strict_axis_level$recipient_mechanism_status_strict ==
          "RECIPIENT_PROGRAM_DIRECTION_ONLY_STRICT",
        "MECHANISM_TIER_3_STRICT",
        "COMMUNICATION_ONLY_STRICT"
      )
    )
  )

  original_tier_lookup <- unique(
    axis_original[
      ,
      c(
        "target",
        "ligand",
        "receptor",
        "interaction_name",
        "pathway_name",
        "dominant_branch",
        "final_mechanism_tier",
        "recipient_mechanism_status",
        "mechanism_chain"
      ),
      drop = FALSE
    ]
  )

  strict_axis_level <- merge(
    strict_axis_level,
    original_tier_lookup,
    by = c(
      "target",
      "ligand",
      "receptor",
      "interaction_name",
      "pathway_name",
      "dominant_branch"
    ),
    all.x = TRUE,
    sort = FALSE
  )

  strict_axis_level$tier_changed_by_5pct_QC <-
    paste0(
      strict_axis_level$final_mechanism_tier,
      "_STRICT"
    ) !=
      strict_axis_level$final_mechanism_tier_strict

  strict_axis_level$strict_mechanism_chain <- paste0(
    strict_axis_level$dominant_branch,
    " -> ",
    strict_axis_level$ligand,
    " -> ",
    strict_axis_level$receptor,
    " -> ",
    strict_axis_level$target,
    " -> ",
    ifelse(
      is.na(
        strict_axis_level$supported_programs_two_sample_strict
      ),
      strict_axis_level$supported_programs_any_strict,
      strict_axis_level$supported_programs_two_sample_strict
    )
  )

  strict_axis_level <- strict_axis_level[
    order(
      match(
        strict_axis_level$final_mechanism_tier_strict,
        c(
          "MECHANISM_TIER_1_STRICT",
          "MECHANISM_TIER_2_STRICT",
          "MECHANISM_TIER_3_STRICT",
          "COMMUNICATION_ONLY_STRICT"
        )
      ),
      -strict_axis_level$n_programs_two_sample_strict,
      -strict_axis_level$top_effect_in_expected_direction_strict
    ),
    ,
    drop = FALSE
  ]

  write_tsv(
    strict_axis_level,
    STRICT_MECHANISM_AXES_FILE
  )

  # ============================================================================
  # PART B. Nfkb1/Nfe2l2 -> ligand regulatory closure
  # ============================================================================

  msg(
    "Loading Step7 motif targets, AUC, metadata and normalized expression..."
  )

  regulon_targets <- data.table::fread(
    STEP7C_TARGET_FILE,
    data.table = FALSE,
    showProgress = FALSE
  )

  require_columns(
    regulon_targets,
    c(
      "regulon",
      "TF",
      "sign",
      "target",
      "weight"
    ),
    "Step7C regulon-target table"
  )

  regulon_targets$regulon <- as.character(
    regulon_targets$regulon
  )

  regulon_targets$TF <- as.character(
    regulon_targets$TF
  )

  regulon_targets$target <- as.character(
    regulon_targets$target
  )

  regulon_targets$weight <- suppressWarnings(
    as.numeric(
      regulon_targets$weight
    )
  )

  regulon_targets$weight_rank_within_regulon <- ave(
    -regulon_targets$weight,
    regulon_targets$regulon,
    FUN = function(values) {
      rank(
        values,
        ties.method = "min"
      )
    }
  )

  regulon_targets$weight_percentile_within_regulon <- ave(
    regulon_targets$weight,
    regulon_targets$regulon,
    FUN = function(values) {
      if (
        length(
          values
        ) ==
          1L
      ) {
        return(
          1
        )
      }

      rank(
        values,
        ties.method = "average"
      ) /
        length(
          values
        )
    }
  )

  focus_regulons <- unique(
    regulon_targets[
      regulon_targets$TF %in%
        FOCUS_TFS,
      c(
        "regulon",
        "TF",
        "sign"
      ),
      drop = FALSE
    ]
  )

  missing_focus_tfs <- setdiff(
    FOCUS_TFS,
    focus_regulons$TF
  )

  if (length(missing_focus_tfs) > 0L) {
    abort_step(
      paste0(
        "Focus TF(s) absent from Step7C regulons: ",
        paste(
          missing_focus_tfs,
          collapse = ", "
        )
      )
    )
  }

  auc_object <- readRDS(
    STEP7D_AUC_FILE
  )

  auc_matrix <- resolve_auc_matrix(
    auc_object
  )

  storage.mode(
    auc_matrix
  ) <- "numeric"

  if (
    nrow(
      auc_matrix
    ) !=
      EXPECTED_AUC_REGULONS ||
    ncol(
      auc_matrix
    ) !=
      EXPECTED_FROZEN_CELLS
  ) {
    abort_step(
      paste0(
        "Unexpected Step7D AUC dimensions: ",
        nrow(
          auc_matrix
        ),
        " x ",
        ncol(
          auc_matrix
        ),
        "; expected ",
        EXPECTED_AUC_REGULONS,
        " x ",
        EXPECTED_FROZEN_CELLS,
        "."
      )
    )
  }

  missing_focus_regulons <- setdiff(
    focus_regulons$regulon,
    rownames(
      auc_matrix
    )
  )

  if (length(missing_focus_regulons) > 0L) {
    abort_step(
      paste0(
        "Focus regulon(s) absent from AUC matrix: ",
        paste(
          missing_focus_regulons,
          collapse = ", "
        )
      )
    )
  }

  metadata <- data.table::fread(
    STEP7A_METADATA_FILE,
    data.table = FALSE,
    showProgress = FALSE
  )

  require_columns(
    metadata,
    c(
      "cell",
      "sample_id",
      "timepoint",
      "step5CR_state"
    ),
    "Step7A frozen metadata"
  )

  metadata$cell <- as.character(
    metadata$cell
  )

  metadata$sample_id <- as.character(
    metadata$sample_id
  )

  metadata$timepoint_raw <- as.character(
    metadata$timepoint
  )

  metadata$state_raw <- as.character(
    metadata$step5CR_state
  )

  metadata$timepoint_standard <- standardize_timepoint(
    metadata$timepoint_raw,
    sample_id = metadata$sample_id
  )

  metadata$step5CR_state <- standardize_state(
    metadata$state_raw
  )

  metadata_index <- match(
    colnames(
      auc_matrix
    ),
    metadata$cell
  )

  if (anyNA(metadata_index)) {
    abort_step(
      "Some AUC cells are absent from Step7A metadata."
    )
  }

  metadata <- metadata[
    metadata_index,
    ,
    drop = FALSE
  ]

  if (
    !identical(
      metadata$cell,
      colnames(
        auc_matrix
      )
    )
  ) {
    abort_step(
      "AUC matrix and metadata cell order do not match exactly."
    )
  }

  terminal_7dpi <- metadata[
    metadata$timepoint_standard ==
      "7dpi" &
    metadata$step5CR_state %in%
      FORMAL_7DPI_STATES,
    ,
    drop = FALSE
  ]

  observed_arm1 <- sum(
    terminal_7dpi$step5CR_state ==
      "Arm_1"
  )

  observed_arm2 <- sum(
    terminal_7dpi$step5CR_state ==
      "Arm_2"
  )

  observed_samples <- sort(
    unique(
      terminal_7dpi$sample_id
    )
  )

  metadata_source_used <- "Step7A_frozen_metadata_after_standardization"

  if (
    observed_arm1 !=
      EXPECTED_7DPI_ARM1_CELLS ||
    observed_arm2 !=
      EXPECTED_7DPI_ARM2_CELLS ||
    length(
      observed_samples
    ) !=
      EXPECTED_7DPI_SAMPLES
  ) {
    msg(
      "Step7A timepoint/state labels did not reproduce the frozen 7 dpi ",
      "branch counts. Falling back to the verified Step8A manifest."
    )

    step8a_manifest <- readRDS(
      STEP8A_MANIFEST_FILE
    )

    if (
      !is.list(
        step8a_manifest
      ) ||
      !"cell_metadata" %in%
        names(
          step8a_manifest
        )
    ) {
      abort_step(
        "Step8A manifest does not contain cell_metadata."
      )
    }

    manifest_metadata <- as.data.frame(
      step8a_manifest$cell_metadata,
      stringsAsFactors = FALSE
    )

    require_columns(
      manifest_metadata,
      c(
        "cell",
        "sample_id",
        "timepoint_standard",
        "cellchat_group"
      ),
      "Step8A manifest cell_metadata"
    )

    manifest_metadata$cell <- as.character(
      manifest_metadata$cell
    )

    manifest_metadata$sample_id <- as.character(
      manifest_metadata$sample_id
    )

    manifest_metadata$timepoint_standard <- standardize_timepoint(
      manifest_metadata$timepoint_standard,
      sample_id = manifest_metadata$sample_id
    )

    manifest_metadata$step5CR_state <- state_from_cellchat_group(
      manifest_metadata$cellchat_group
    )

    manifest_index <- match(
      colnames(
        auc_matrix
      ),
      manifest_metadata$cell
    )

    if (anyNA(manifest_index)) {
      abort_step(
        paste0(
          "Step8A manifest is missing ",
          sum(
            is.na(
              manifest_index
            )
          ),
          " AUC cells."
        )
      )
    }

    metadata <- manifest_metadata[
      manifest_index,
      c(
        "cell",
        "sample_id",
        "timepoint_standard",
        "step5CR_state"
      ),
      drop = FALSE
    ]

    if (
      !identical(
        metadata$cell,
        colnames(
          auc_matrix
        )
      )
    ) {
      abort_step(
        "AUC matrix and Step8A manifest cell order do not match exactly."
      )
    }

    terminal_7dpi <- metadata[
      metadata$timepoint_standard ==
        "7dpi" &
      metadata$step5CR_state %in%
        FORMAL_7DPI_STATES,
      ,
      drop = FALSE
    ]

    observed_arm1 <- sum(
      terminal_7dpi$step5CR_state ==
        "Arm_1"
    )

    observed_arm2 <- sum(
      terminal_7dpi$step5CR_state ==
        "Arm_2"
    )

    observed_samples <- sort(
      unique(
        terminal_7dpi$sample_id
      )
    )

    metadata_source_used <- "Step8A_verified_manifest"
  }

  if (
    observed_arm1 !=
      EXPECTED_7DPI_ARM1_CELLS ||
    observed_arm2 !=
      EXPECTED_7DPI_ARM2_CELLS ||
    length(
      observed_samples
    ) !=
      EXPECTED_7DPI_SAMPLES
  ) {
    raw_time_values <- collapse_unique(
      metadata$timepoint_standard
    )

    raw_state_values <- collapse_unique(
      metadata$step5CR_state
    )

    abort_step(
      paste0(
        "Unexpected 7 dpi terminal branch structure after metadata fallback: ",
        "Arm_1=",
        observed_arm1,
        "; Arm_2=",
        observed_arm2,
        "; samples=",
        length(
          observed_samples
        ),
        "; metadata_source=",
        metadata_source_used,
        "; standardized_timepoints=",
        raw_time_values,
        "; standardized_states=",
        raw_state_values,
        "."
      )
    )
  }

  msg(
    "Resolved 7 dpi branch metadata from ",
    metadata_source_used,
    ": Arm_1=",
    observed_arm1,
    "; Arm_2=",
    observed_arm2,
    "; samples=",
    paste(
      observed_samples,
      collapse = " | "
    ),
    "."
  )

  atlas <- readRDS(
    ATLAS_RDS_FILE
  )

  if (!inherits(atlas, "Seurat")) {
    abort_step(
      "The full atlas input is not a Seurat object."
    )
  }

  Seurat::DefaultAssay(
    atlas
  ) <- "RNA"

  expression_data <- get_rna_data(
    atlas
  )

  missing_terminal_cells <- setdiff(
    terminal_7dpi$cell,
    colnames(
      expression_data
    )
  )

  if (length(missing_terminal_cells) > 0L) {
    abort_step(
      paste0(
        "Normalized RNA data are missing ",
        length(
          missing_terminal_cells
        ),
        " 7 dpi branch cells."
      )
    )
  }

  missing_focus_ligands <- setdiff(
    FOCUS_LIGANDS,
    rownames(
      expression_data
    )
  )

  if (length(missing_focus_ligands) > 0L) {
    abort_step(
      paste0(
        "Focus ligand(s) absent from normalized RNA data: ",
        paste(
          missing_focus_ligands,
          collapse = ", "
        )
      )
    )
  }

  expression_focus <- as.matrix(
    expression_data[
      FOCUS_LIGANDS,
      terminal_7dpi$cell,
      drop = FALSE
    ]
  )

  rm(
    atlas,
    expression_data
  )
  gc()

  auc_focus <- auc_matrix[
    focus_regulons$regulon,
    terminal_7dpi$cell,
    drop = FALSE
  ]

  candidate_ligand_catalog <- data.table::fread(
    STEP8A_CANDIDATE_LIGAND_FILE,
    data.table = FALSE,
    showProgress = FALSE
  )

  step7g_closure <- data.table::fread(
    STEP7G_CLOSURE_FILE,
    data.table = FALSE,
    showProgress = FALSE
  )

  motif_rows <- list()
  motif_index <- 0L

  for (
    regulon_index in
    seq_len(
      nrow(
        focus_regulons
      )
    )
  ) {
    current_regulon <- focus_regulons$regulon[
      regulon_index
    ]

    current_tf <- focus_regulons$TF[
      regulon_index
    ]

    current_sign <- focus_regulons$sign[
      regulon_index
    ]

    current_targets <- regulon_targets[
      regulon_targets$regulon ==
        current_regulon,
      ,
      drop = FALSE
    ]

    for (ligand in FOCUS_LIGANDS) {
      direct_row <- current_targets[
        current_targets$target ==
          ligand,
        ,
        drop = FALSE
      ]

      in_step8a_catalog <- FALSE
      step8a_priority <- NA_character_
      step8a_expression_direction <- NA_character_

      if (
        all(
          c(
            "TF",
            "gene"
          ) %in%
            colnames(
              candidate_ligand_catalog
            )
        )
      ) {
        catalog_hit <- candidate_ligand_catalog[
          candidate_ligand_catalog$TF ==
            current_tf &
          candidate_ligand_catalog$gene ==
            ligand,
          ,
          drop = FALSE
        ]

        in_step8a_catalog <- nrow(
          catalog_hit
        ) >
          0L

        if (
          in_step8a_catalog &&
          "mechanistic_priority" %in%
            colnames(
              catalog_hit
            )
        ) {
          step8a_priority <- collapse_unique(
            catalog_hit$mechanistic_priority
          )
        }

        if (
          in_step8a_catalog &&
          "expression_direction_7dpi" %in%
            colnames(
              catalog_hit
            )
        ) {
          step8a_expression_direction <- collapse_unique(
            catalog_hit$expression_direction_7dpi
          )
        }
      }

      motif_index <- motif_index + 1L

      motif_rows[[motif_index]] <- data.frame(
        TF = current_tf,
        regulon = current_regulon,
        regulon_sign = current_sign,
        ligand = ligand,
        motif_supported_direct_target =
          nrow(
            direct_row
          ) >
            0L,
        target_weight = if (
          nrow(
            direct_row
          ) >
            0L
        ) {
          safe_max_finite(
            direct_row$weight
          )
        } else {
          NA_real_
        },
        target_weight_rank_within_regulon = if (
          nrow(
            direct_row
          ) >
            0L
        ) {
          safe_min_finite(
            direct_row$weight_rank_within_regulon
          )
        } else {
          NA_real_
        },
        target_weight_percentile_within_regulon = if (
          nrow(
            direct_row
          ) >
            0L
        ) {
          safe_max_finite(
            direct_row$weight_percentile_within_regulon
          )
        } else {
          NA_real_
        },
        present_in_Step8A_candidate_ligand_catalog =
          in_step8a_catalog,
        Step8A_mechanistic_priority =
          step8a_priority,
        Step8A_expression_direction_7dpi =
          step8a_expression_direction,
        stringsAsFactors = FALSE
      )
    }
  }

  motif_audit <- do.call(
    rbind,
    motif_rows
  )

  write_tsv(
    motif_audit,
    TF_LIGAND_MOTIF_AUDIT_FILE
  )

  correlation_rows <- list()
  correlation_index <- 0L

  for (
    motif_index in
    seq_len(
      nrow(
        motif_audit
      )
    )
  ) {
    current <- motif_audit[
      motif_index,
      ,
      drop = FALSE
    ]

    auc_values <- as.numeric(
      auc_focus[
        current$regulon,
        terminal_7dpi$cell,
        drop = TRUE
      ]
    )

    expression_values <- as.numeric(
      expression_focus[
        current$ligand,
        terminal_7dpi$cell,
        drop = TRUE
      ]
    )

    adjusted_auc <- safe_lm_residuals(
      response = auc_values,
      state = terminal_7dpi$step5CR_state,
      sample_id = terminal_7dpi$sample_id
    )

    adjusted_expression <- safe_lm_residuals(
      response = expression_values,
      state = terminal_7dpi$step5CR_state,
      sample_id = terminal_7dpi$sample_id
    )

    pooled_result <- safe_spearman(
      auc_values,
      expression_values
    )

    adjusted_result <- safe_spearman(
      adjusted_auc,
      adjusted_expression
    )

    correlation_index <- correlation_index + 1L

    correlation_rows[[correlation_index]] <- data.frame(
      TF = current$TF,
      regulon = current$regulon,
      ligand = current$ligand,
      analysis_level = "POOLED_7DPI_TERMINAL",
      sample_id = "ALL_7DPI_SAMPLES",
      n_cells = pooled_result[["n"]],
      spearman_rho = pooled_result[["rho"]],
      spearman_p_value = pooled_result[["p_value"]],
      mean_AUC_Arm_1 = safe_mean_finite(
        auc_values[
          terminal_7dpi$step5CR_state ==
            "Arm_1"
        ]
      ),
      mean_AUC_Arm_2 = safe_mean_finite(
        auc_values[
          terminal_7dpi$step5CR_state ==
            "Arm_2"
        ]
      ),
      delta_AUC_Arm2_minus_Arm1 =
        safe_mean_finite(
          auc_values[
            terminal_7dpi$step5CR_state ==
              "Arm_2"
          ]
        ) -
        safe_mean_finite(
          auc_values[
            terminal_7dpi$step5CR_state ==
              "Arm_1"
          ]
        ),
      mean_ligand_Arm_1 = safe_mean_finite(
        expression_values[
          terminal_7dpi$step5CR_state ==
            "Arm_1"
        ]
      ),
      mean_ligand_Arm_2 = safe_mean_finite(
        expression_values[
          terminal_7dpi$step5CR_state ==
            "Arm_2"
        ]
      ),
      delta_ligand_Arm2_minus_Arm1 =
        safe_mean_finite(
          expression_values[
            terminal_7dpi$step5CR_state ==
              "Arm_2"
          ]
        ) -
        safe_mean_finite(
          expression_values[
            terminal_7dpi$step5CR_state ==
              "Arm_1"
          ]
        ),
      TF_AUC_higher_in_Arm2 =
        safe_mean_finite(
          auc_values[
            terminal_7dpi$step5CR_state ==
              "Arm_2"
          ]
        ) >
        safe_mean_finite(
          auc_values[
            terminal_7dpi$step5CR_state ==
              "Arm_1"
          ]
        ),
      ligand_higher_in_Arm2 =
        safe_mean_finite(
          expression_values[
            terminal_7dpi$step5CR_state ==
              "Arm_2"
          ]
        ) >
        safe_mean_finite(
          expression_values[
            terminal_7dpi$step5CR_state ==
              "Arm_1"
          ]
        ),
      adjusted_spearman_rho = adjusted_result[["rho"]],
      adjusted_spearman_p_value = adjusted_result[["p_value"]],
      stringsAsFactors = FALSE
    )

    for (sample_id in observed_samples) {
      sample_index <- terminal_7dpi$sample_id ==
        sample_id

      sample_auc <- auc_values[
        sample_index
      ]

      sample_expression <- expression_values[
        sample_index
      ]

      sample_state <- terminal_7dpi$step5CR_state[
        sample_index
      ]

      sample_result <- safe_spearman(
        sample_auc,
        sample_expression
      )

      mean_auc_arm1 <- safe_mean_finite(
        sample_auc[
          sample_state ==
            "Arm_1"
        ]
      )

      mean_auc_arm2 <- safe_mean_finite(
        sample_auc[
          sample_state ==
            "Arm_2"
        ]
      )

      mean_ligand_arm1 <- safe_mean_finite(
        sample_expression[
          sample_state ==
            "Arm_1"
        ]
      )

      mean_ligand_arm2 <- safe_mean_finite(
        sample_expression[
          sample_state ==
            "Arm_2"
        ]
      )

      correlation_index <- correlation_index + 1L

      correlation_rows[[correlation_index]] <- data.frame(
        TF = current$TF,
        regulon = current$regulon,
        ligand = current$ligand,
        analysis_level = "WITHIN_7DPI_SAMPLE",
        sample_id = sample_id,
        n_cells = sample_result[["n"]],
        spearman_rho = sample_result[["rho"]],
        spearman_p_value = sample_result[["p_value"]],
        mean_AUC_Arm_1 = mean_auc_arm1,
        mean_AUC_Arm_2 = mean_auc_arm2,
        delta_AUC_Arm2_minus_Arm1 =
          mean_auc_arm2 -
          mean_auc_arm1,
        mean_ligand_Arm_1 = mean_ligand_arm1,
        mean_ligand_Arm_2 = mean_ligand_arm2,
        delta_ligand_Arm2_minus_Arm1 =
          mean_ligand_arm2 -
          mean_ligand_arm1,
        TF_AUC_higher_in_Arm2 =
          mean_auc_arm2 >
          mean_auc_arm1,
        ligand_higher_in_Arm2 =
          mean_ligand_arm2 >
          mean_ligand_arm1,
        adjusted_spearman_rho = NA_real_,
        adjusted_spearman_p_value = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }

  correlation_long <- do.call(
    rbind,
    correlation_rows
  )

  correlation_long$spearman_fdr <- stats::p.adjust(
    correlation_long$spearman_p_value,
    method = "BH"
  )

  pooled_rows <- correlation_long[
    correlation_long$analysis_level ==
      "POOLED_7DPI_TERMINAL",
    ,
    drop = FALSE
  ]

  sample_rows <- correlation_long[
    correlation_long$analysis_level ==
      "WITHIN_7DPI_SAMPLE",
    ,
    drop = FALSE
  ]

  sample_rows$sample_state_concordance <-
    sample_rows$TF_AUC_higher_in_Arm2 %in%
      TRUE &
    sample_rows$ligand_higher_in_Arm2 %in%
      TRUE

  sample_rows$sample_positive_correlation <-
    is.finite(
      sample_rows$spearman_rho
    ) &
    sample_rows$spearman_rho >=
      MIN_ABS_SAMPLE_RHO_FOR_SUPPORT

  correlation_long$sample_state_concordance <- NA
  correlation_long$sample_positive_correlation <- NA

  sample_match <- match(
    paste(
      sample_rows$TF,
      sample_rows$regulon,
      sample_rows$ligand,
      sample_rows$analysis_level,
      sample_rows$sample_id,
      sep = "||"
    ),
    paste(
      correlation_long$TF,
      correlation_long$regulon,
      correlation_long$ligand,
      correlation_long$analysis_level,
      correlation_long$sample_id,
      sep = "||"
    )
  )

  correlation_long$sample_state_concordance[
    sample_match
  ] <- sample_rows$sample_state_concordance

  correlation_long$sample_positive_correlation[
    sample_match
  ] <- sample_rows$sample_positive_correlation

  write_tsv(
    correlation_long,
    TF_LIGAND_CORRELATION_LONG_FILE
  )

  sample_summary <- data.table::as.data.table(
    sample_rows
  )[
    ,
    .(
      n_samples = data.table::uniqueN(
        sample_id
      ),
      n_samples_state_concordant = sum(
        sample_state_concordance,
        na.rm = TRUE
      ),
      n_samples_positive_correlation = sum(
        sample_positive_correlation,
        na.rm = TRUE
      ),
      mean_sample_rho = safe_mean_finite(
        spearman_rho
      ),
      minimum_sample_rho = safe_min_finite(
        spearman_rho
      ),
      sample_rhos = paste(
        paste0(
          sample_id,
          ":",
          round(
            spearman_rho,
            4
          )
        ),
        collapse = " | "
      ),
      sample_AUC_deltas = paste(
        paste0(
          sample_id,
          ":",
          round(
            delta_AUC_Arm2_minus_Arm1,
            5
          )
        ),
        collapse = " | "
      ),
      sample_ligand_deltas = paste(
        paste0(
          sample_id,
          ":",
          round(
            delta_ligand_Arm2_minus_Arm1,
            5
          )
        ),
        collapse = " | "
      )
    ),
    by = .(
      TF,
      regulon,
      ligand
    )
  ]

  sample_summary <- as.data.frame(
    sample_summary,
    stringsAsFactors = FALSE
  )

  pooled_keep <- pooled_rows[
    ,
    c(
      "TF",
      "regulon",
      "ligand",
      "spearman_rho",
      "spearman_p_value",
      "spearman_fdr",
      "adjusted_spearman_rho",
      "adjusted_spearman_p_value",
      "delta_AUC_Arm2_minus_Arm1",
      "delta_ligand_Arm2_minus_Arm1",
      "TF_AUC_higher_in_Arm2",
      "ligand_higher_in_Arm2"
    ),
    drop = FALSE
  ]

  colnames(
    pooled_keep
  )[
    colnames(
      pooled_keep
    ) ==
      "spearman_rho"
  ] <- "pooled_spearman_rho"

  colnames(
    pooled_keep
  )[
    colnames(
      pooled_keep
    ) ==
      "spearman_p_value"
  ] <- "pooled_spearman_p_value"

  colnames(
    pooled_keep
  )[
    colnames(
      pooled_keep
    ) ==
      "spearman_fdr"
  ] <- "pooled_spearman_fdr"

  colnames(
    pooled_keep
  )[
    colnames(
      pooled_keep
    ) ==
      "delta_AUC_Arm2_minus_Arm1"
  ] <- "pooled_delta_AUC_Arm2_minus_Arm1"

  colnames(
    pooled_keep
  )[
    colnames(
      pooled_keep
    ) ==
      "delta_ligand_Arm2_minus_Arm1"
  ] <- "pooled_delta_ligand_Arm2_minus_Arm1"

  colnames(
    pooled_keep
  )[
    colnames(
      pooled_keep
    ) ==
      "TF_AUC_higher_in_Arm2"
  ] <- "pooled_TF_AUC_higher_in_Arm2"

  colnames(
    pooled_keep
  )[
    colnames(
      pooled_keep
    ) ==
      "ligand_higher_in_Arm2"
  ] <- "pooled_ligand_higher_in_Arm2"

  closure_summary <- merge(
    motif_audit,
    sample_summary,
    by = c(
      "TF",
      "regulon",
      "ligand"
    ),
    all.x = TRUE,
    sort = FALSE
  )

  closure_summary <- merge(
    closure_summary,
    pooled_keep,
    by = c(
      "TF",
      "regulon",
      "ligand"
    ),
    all.x = TRUE,
    sort = FALSE
  )

  closure_summary$adjusted_positive_association <-
    is.finite(
      closure_summary$adjusted_spearman_rho
    ) &
    closure_summary$adjusted_spearman_rho >=
      MIN_ABS_ADJUSTED_RHO_FOR_SUPPORT

  closure_summary$two_sample_state_concordance <-
    closure_summary$n_samples_state_concordant ==
      EXPECTED_7DPI_SAMPLES

  closure_summary$two_sample_positive_correlation <-
    closure_summary$n_samples_positive_correlation ==
      EXPECTED_7DPI_SAMPLES

  closure_summary$regulatory_closure_status <- ifelse(
    closure_summary$motif_supported_direct_target &
      closure_summary$two_sample_state_concordance &
      closure_summary$two_sample_positive_correlation &
      closure_summary$adjusted_positive_association,
    "DIRECT_MOTIF_TARGET_STRONG_TWO_SAMPLE_CLOSURE",
    ifelse(
      closure_summary$motif_supported_direct_target &
        closure_summary$two_sample_state_concordance &
        (
          closure_summary$two_sample_positive_correlation |
          closure_summary$adjusted_positive_association
        ),
      "DIRECT_MOTIF_TARGET_PARTIAL_CORRELATION_CLOSURE",
      ifelse(
        closure_summary$motif_supported_direct_target &
          closure_summary$two_sample_state_concordance,
        "DIRECT_MOTIF_TARGET_STATE_CONCORDANCE_ONLY",
        ifelse(
          closure_summary$motif_supported_direct_target,
          "DIRECT_MOTIF_TARGET_WITHOUT_STABLE_STATE_CLOSURE",
          ifelse(
            closure_summary$two_sample_state_concordance &
              closure_summary$two_sample_positive_correlation,
            "INDIRECT_TWO_SAMPLE_STATE_AND_CORRELATION_ASSOCIATION",
            ifelse(
              closure_summary$two_sample_state_concordance,
              "INDIRECT_TWO_SAMPLE_STATE_ASSOCIATION",
              "NO_REGULON_LIGAND_CLOSURE"
            )
          )
        )
      )
    )
  )

  closure_summary$manuscript_interpretation <- ifelse(
    closure_summary$regulatory_closure_status ==
      "DIRECT_MOTIF_TARGET_STRONG_TWO_SAMPLE_CLOSURE",
    "Direct motif-supported candidate with replicated state and cell-level association",
    ifelse(
      grepl(
        "^DIRECT_MOTIF_TARGET",
        closure_summary$regulatory_closure_status
      ),
      "Direct motif-supported candidate, but association strength is incomplete",
      ifelse(
        grepl(
          "^INDIRECT",
          closure_summary$regulatory_closure_status
        ),
        "State-associated candidate without a motif-supported direct ligand edge",
        "No evidence for a direct TF-ligand closure in the current data"
      )
    )
  )

  closure_rank <- match(
    closure_summary$regulatory_closure_status,
    c(
      "DIRECT_MOTIF_TARGET_STRONG_TWO_SAMPLE_CLOSURE",
      "DIRECT_MOTIF_TARGET_PARTIAL_CORRELATION_CLOSURE",
      "DIRECT_MOTIF_TARGET_STATE_CONCORDANCE_ONLY",
      "DIRECT_MOTIF_TARGET_WITHOUT_STABLE_STATE_CLOSURE",
      "INDIRECT_TWO_SAMPLE_STATE_AND_CORRELATION_ASSOCIATION",
      "INDIRECT_TWO_SAMPLE_STATE_ASSOCIATION",
      "NO_REGULON_LIGAND_CLOSURE"
    )
  )

  closure_rank[
    is.na(
      closure_rank
    )
  ] <- 99L

  closure_summary <- closure_summary[
    order(
      closure_rank,
      closure_summary$TF,
      closure_summary$ligand
    ),
    ,
    drop = FALSE
  ]

  write_tsv(
    closure_summary,
    TF_LIGAND_CLOSURE_SUMMARY_FILE
  )

  # ============================================================================
  # PART C. FINAL STRICT MECHANISM MODEL
  # ============================================================================

  focus_closure_by_ligand <- data.table::as.data.table(
    closure_summary
  )[
    ,
    .(
      upstream_focus_TFs = collapse_unique(
        TF
      ),
      upstream_closure_statuses = collapse_unique(
        regulatory_closure_status
      ),
      upstream_direct_motif_TFs = collapse_unique(
        TF[
          motif_supported_direct_target
        ]
      ),
      upstream_strong_closure_TFs = collapse_unique(
        TF[
          regulatory_closure_status ==
            "DIRECT_MOTIF_TARGET_STRONG_TWO_SAMPLE_CLOSURE"
        ]
      ),
      upstream_best_closure_rank = min(
        match(
          regulatory_closure_status,
          c(
            "DIRECT_MOTIF_TARGET_STRONG_TWO_SAMPLE_CLOSURE",
            "DIRECT_MOTIF_TARGET_PARTIAL_CORRELATION_CLOSURE",
            "DIRECT_MOTIF_TARGET_STATE_CONCORDANCE_ONLY",
            "DIRECT_MOTIF_TARGET_WITHOUT_STABLE_STATE_CLOSURE",
            "INDIRECT_TWO_SAMPLE_STATE_AND_CORRELATION_ASSOCIATION",
            "INDIRECT_TWO_SAMPLE_STATE_ASSOCIATION",
            "NO_REGULON_LIGAND_CLOSURE"
          )
        ),
        na.rm = TRUE
      )
    ),
    by = ligand
  ]

  focus_closure_by_ligand <- as.data.frame(
    focus_closure_by_ligand,
    stringsAsFactors = FALSE
  )

  final_model <- merge(
    strict_axis_level,
    focus_closure_by_ligand,
    by = "ligand",
    all.x = TRUE,
    sort = FALSE
  )

  final_model$upstream_closure_scope <- ifelse(
    final_model$ligand %in%
      FOCUS_LIGANDS,
    "Nfkb1_Nfe2l2_FOCUS_AUDITED",
    "UPSTREAM_NOT_TESTED_IN_STEP8D"
  )

  final_model$final_model_status <- ifelse(
    final_model$final_mechanism_tier_strict ==
      "MECHANISM_TIER_1_STRICT" &
    final_model$upstream_best_closure_rank <=
      2,
    "FULL_ASSOCIATIVE_CHAIN_WITH_DIRECT_UPSTREAM_SUPPORT",
    ifelse(
      final_model$final_mechanism_tier_strict ==
        "MECHANISM_TIER_1_STRICT",
      "STRICT_LIGAND_RECEPTOR_RECIPIENT_CHAIN",
      ifelse(
        final_model$final_mechanism_tier_strict ==
          "MECHANISM_TIER_2_STRICT",
        "PARTIAL_RECIPIENT_PROGRAM_CHAIN",
        "COMMUNICATION_OR_DIRECTION_ONLY"
      )
    )
  )

  final_model$final_chain <- paste0(
    ifelse(
      is.na(
        final_model$upstream_direct_motif_TFs
      ),
      "[upstream unresolved]",
      final_model$upstream_direct_motif_TFs
    ),
    " -> ",
    final_model$dominant_branch,
    " ",
    final_model$ligand,
    " -> ",
    final_model$receptor,
    " -> ",
    final_model$target,
    " -> ",
    ifelse(
      is.na(
        final_model$supported_programs_two_sample_strict
      ),
      final_model$supported_programs_any_strict,
      final_model$supported_programs_two_sample_strict
    )
  )

  final_model <- final_model[
    order(
      match(
        final_model$final_model_status,
        c(
          "FULL_ASSOCIATIVE_CHAIN_WITH_DIRECT_UPSTREAM_SUPPORT",
          "STRICT_LIGAND_RECEPTOR_RECIPIENT_CHAIN",
          "PARTIAL_RECIPIENT_PROGRAM_CHAIN",
          "COMMUNICATION_OR_DIRECTION_ONLY"
        )
      ),
      match(
        final_model$final_mechanism_tier_strict,
        c(
          "MECHANISM_TIER_1_STRICT",
          "MECHANISM_TIER_2_STRICT",
          "MECHANISM_TIER_3_STRICT",
          "COMMUNICATION_ONLY_STRICT"
        )
      )
    ),
    ,
    drop = FALSE
  ]

  write_tsv(
    final_model,
    FINAL_MECHANISM_MODEL_FILE
  )

  # ============================================================================
  # FIGURES
  # ============================================================================

  strict_figure_data <- strict_axis_level

  strict_figure_data$axis_label <- paste0(
    strict_figure_data$interaction_name,
    " -> ",
    strict_figure_data$target
  )

  strict_figure_data$axis_label <- factor(
    strict_figure_data$axis_label,
    levels = rev(
      strict_figure_data$axis_label
    )
  )

  strict_plot <- ggplot2::ggplot(
    strict_figure_data,
    ggplot2::aes(
      x = n_programs_two_sample_strict,
      y = axis_label
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::labs(
      x = "Programs supported in both 7 dpi samples after >=5% receptor QC",
      y = NULL,
      title = "Strict recipient-program mechanism support"
    ) +
    ggplot2::theme_bw(
      base_size = 8
    )

  ggplot2::ggsave(
    filename = FIGURE_STRICT_TIERS_FILE,
    plot = strict_plot,
    width = 10,
    height = 8,
    units = "in"
  )

  closure_figure_data <- closure_summary

  closure_figure_data$pair_label <- paste0(
    closure_figure_data$TF,
    " -> ",
    closure_figure_data$ligand
  )

  closure_figure_data$pair_label <- factor(
    closure_figure_data$pair_label,
    levels = rev(
      closure_figure_data$pair_label
    )
  )

  closure_plot <- ggplot2::ggplot(
    closure_figure_data,
    ggplot2::aes(
      x = mean_sample_rho,
      y = pair_label
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = 2
    ) +
    ggplot2::labs(
      x = "Mean within-sample Spearman rho: regulon AUC vs ligand expression",
      y = NULL,
      title = "Nfkb1/Nfe2l2 to inflammatory-ligand closure audit"
    ) +
    ggplot2::theme_bw(
      base_size = 9
    )

  ggplot2::ggsave(
    filename = FIGURE_TF_LIGAND_FILE,
    plot = closure_plot,
    width = 8,
    height = 5,
    units = "in"
  )

  # ============================================================================
  # FINAL DECISION
  # ============================================================================

  decision_rows <- list()

  decision_rows <- add_decision_row(
    decision_rows,
    "Required input files",
    all(
      input_audit$exists
    ),
    detail = paste0(
      sum(
        input_audit$exists
      ),
      "/",
      nrow(
        input_audit
      ),
      " required files present"
    )
  )

  decision_rows <- add_decision_row(
    decision_rows,
    "Step7C/7D/7G/8A/8C PASS gates",
    TRUE,
    detail = "All upstream formal decision gates passed"
  )

  decision_rows <- add_decision_row(
    decision_rows,
    "Strict 5% receptor-positive QC",
    nrow(
      receptor_audit
    ) >
      0L &&
    all(
      c(
        "receptor_positive_5pct_pass",
        "sample_support_strict"
      ) %in%
        colnames(
          receptor_audit
        )
    ),
    detail = paste0(
      sum(
        receptor_audit$receptor_positive_5pct_pass,
        na.rm = TRUE
      ),
      "/",
      nrow(
        receptor_audit
      ),
      " sample-program rows passed >=5%"
    )
  )

  n_strict_tier1 <- sum(
    strict_axis_level$final_mechanism_tier_strict ==
      "MECHANISM_TIER_1_STRICT",
    na.rm = TRUE
  )

  n_strict_tier2 <- sum(
    strict_axis_level$final_mechanism_tier_strict ==
      "MECHANISM_TIER_2_STRICT",
    na.rm = TRUE
  )

  n_changed <- sum(
    strict_axis_level$tier_changed_by_5pct_QC,
    na.rm = TRUE
  )

  decision_rows <- add_decision_row(
    decision_rows,
    "Strict recipient mechanism tiers",
    n_strict_tier1 >
      0L,
    detail = paste0(
      n_strict_tier1,
      " strict Tier1; ",
      n_strict_tier2,
      " strict Tier2; ",
      n_changed,
      " axes changed tier after 5% QC"
    ),
    status_if_fail = "REVIEW_REQUIRED"
  )

  decision_rows <- add_decision_row(
    decision_rows,
    "AUC and metadata integrity",
    nrow(
      auc_matrix
    ) ==
      EXPECTED_AUC_REGULONS &&
    ncol(
      auc_matrix
    ) ==
      EXPECTED_FROZEN_CELLS &&
    observed_arm1 ==
      EXPECTED_7DPI_ARM1_CELLS &&
    observed_arm2 ==
      EXPECTED_7DPI_ARM2_CELLS,
    detail = paste0(
      "AUC=",
      nrow(
        auc_matrix
      ),
      "x",
      ncol(
        auc_matrix
      ),
      "; 7dpi Arm_1=",
      observed_arm1,
      "; Arm_2=",
      observed_arm2
    )
  )

  expected_focus_pairs <- length(
    FOCUS_TFS
  ) *
    length(
      FOCUS_LIGANDS
    )

  decision_rows <- add_decision_row(
    decision_rows,
    "Focus TF-ligand motif audit",
    nrow(
      motif_audit
    ) ==
      expected_focus_pairs,
    detail = paste0(
      nrow(
        motif_audit
      ),
      "/",
      expected_focus_pairs,
      " TF-ligand pairs audited; ",
      sum(
        motif_audit$motif_supported_direct_target
      ),
      " direct motif-supported edges"
    )
  )

  decision_rows <- add_decision_row(
    decision_rows,
    "Two-sample regulon-ligand audit",
    nrow(
      sample_rows
    ) ==
      expected_focus_pairs *
        EXPECTED_7DPI_SAMPLES,
    detail = paste0(
      nrow(
        sample_rows
      ),
      " sample-pair rows across ",
      length(
        observed_samples
      ),
      " 7 dpi samples"
    )
  )

  n_direct_strong <- sum(
    closure_summary$regulatory_closure_status ==
      "DIRECT_MOTIF_TARGET_STRONG_TWO_SAMPLE_CLOSURE",
    na.rm = TRUE
  )

  n_direct_any <- sum(
    grepl(
      "^DIRECT_MOTIF_TARGET",
      closure_summary$regulatory_closure_status
    ),
    na.rm = TRUE
  )

  biological_closure_status <- if (
    n_direct_strong >
      0L
  ) {
    "STRONG_DIRECT_TF_LIGAND_CLOSURE_IDENTIFIED"
  } else if (
    n_direct_any >
      0L
  ) {
    "DIRECT_TF_LIGAND_CANDIDATES_WITH_INCOMPLETE_CLOSURE"
  } else {
    "NO_DIRECT_NFKB1_NFE2L2_LIGAND_CLOSURE"
  }

  decision_rows <- add_decision_row(
    decision_rows,
    "Direct upstream regulatory closure",
    n_direct_any >
      0L,
    detail = paste0(
      biological_closure_status,
      "; strong=",
      n_direct_strong,
      "; any_direct=",
      n_direct_any
    ),
    status_if_fail = "REVIEW_REQUIRED"
  )

  technical_pass <- all(
    vapply(
      decision_rows[
        seq_len(
          min(
            7L,
            length(
              decision_rows
            )
          )
        )
      ],
      function(row) {
        isTRUE(
          row$pass[[1L]]
        )
      },
      FUN.VALUE = logical(1)
    )
  )

  decision_rows <- add_decision_row(
    decision_rows,
    "Overall Step 8D technical",
    technical_pass,
    detail = if (
      technical_pass
    ) {
      "PASS"
    } else {
      "At least one technical domain failed"
    }
  )

  decision_rows <- add_decision_row(
    decision_rows,
    "Overall Step 8D biological result",
    TRUE,
    detail = biological_closure_status
  )

  final_decision <- do.call(
    rbind,
    decision_rows
  )

  rownames(
    final_decision
  ) <- NULL

  write_csv_utf8(
    final_decision,
    FINAL_DECISION_FILE
  )

  result_object <- list(
    status = if (
      technical_pass
    ) {
      "PASS"
    } else {
      "FAIL"
    },
    biological_closure_status =
      biological_closure_status,
    script_version = "Step8D_v2_metadata_fix",
    metadata_source_used =
      metadata_source_used,
    generated = timestamp(),
    parameters = list(
      receptor_positive_percent_threshold =
        RECEPTOR_POSITIVE_PERCENT_THRESHOLD,
      focus_TFs =
        FOCUS_TFS,
      focus_ligands =
        FOCUS_LIGANDS,
      minimum_sample_rho =
        MIN_ABS_SAMPLE_RHO_FOR_SUPPORT,
      minimum_adjusted_rho =
        MIN_ABS_ADJUSTED_RHO_FOR_SUPPORT
    ),
    strict_receptor_audit =
      receptor_audit,
    strict_axis_program =
      strict_axis_program,
    strict_mechanism_axes =
      strict_axis_level,
    motif_audit =
      motif_audit,
    correlation_long =
      correlation_long,
    closure_summary =
      closure_summary,
    final_mechanism_model =
      final_model,
    final_decision =
      final_decision
  )

  saveRDS(
    result_object,
    RESULT_RDS_FILE,
    compress = FALSE
  )

  elapsed_minutes <- as.numeric(
    difftime(
      Sys.time(),
      start_time,
      units = "mins"
    )
  )

  strict_axis_text <- paste(
    capture.output(
      print(
        strict_axis_level[
          ,
          intersect(
            c(
              "target",
              "ligand",
              "receptor",
              "interaction_name",
              "pathway_name",
              "dominant_branch",
              "final_mechanism_tier",
              "final_mechanism_tier_strict",
              "tier_changed_by_5pct_QC",
              "minimum_receptor_positive_percent",
              "supported_programs_two_sample_strict",
              "strict_mechanism_chain"
            ),
            colnames(
              strict_axis_level
            )
          ),
          drop = FALSE
        ],
        row.names = FALSE
      )
    ),
    collapse = "\n"
  )

  closure_text <- paste(
    capture.output(
      print(
        closure_summary[
          ,
          intersect(
            c(
              "TF",
              "regulon",
              "ligand",
              "motif_supported_direct_target",
              "target_weight",
              "n_samples_state_concordant",
              "n_samples_positive_correlation",
              "mean_sample_rho",
              "adjusted_spearman_rho",
              "regulatory_closure_status",
              "manuscript_interpretation"
            ),
            colnames(
              closure_summary
            )
          ),
          drop = FALSE
        ],
        row.names = FALSE
      )
    ),
    collapse = "\n"
  )

  final_model_text <- paste(
    capture.output(
      print(
        final_model[
          ,
          intersect(
            c(
              "target",
              "ligand",
              "receptor",
              "dominant_branch",
              "final_mechanism_tier_strict",
              "upstream_direct_motif_TFs",
              "upstream_closure_statuses",
              "final_model_status",
              "final_chain"
            ),
            colnames(
              final_model
            )
          ),
          drop = FALSE
        ],
        row.names = FALSE
      )
    ),
    collapse = "\n"
  )

  report_lines <- c(
    "======================================================================",
    "GSE162610 | STEP 8D v2 MECHANISM CLOSURE AND QC",
    "======================================================================",
    "",
    paste0(
      "Generated: ",
      timestamp()
    ),
    paste0(
      "Technical status: ",
      if (
        technical_pass
      ) {
        "PASS"
      } else {
        "FAIL"
      }
    ),
    paste0(
      "Biological closure result: ",
      biological_closure_status
    ),
    paste0(
      "Elapsed time: ",
      round(
        elapsed_minutes,
        2
      ),
      " min"
    ),
    paste0(
      "7 dpi branch metadata source: ",
      metadata_source_used
    ),
    "",
    "PART A | STRICT >=5% RECIPIENT RECEPTOR QC:",
    paste0(
      "Strict MECHANISM_TIER_1 axes: ",
      n_strict_tier1
    ),
    paste0(
      "Strict MECHANISM_TIER_2 axes: ",
      n_strict_tier2
    ),
    paste0(
      "Axes whose tier changed after strict QC: ",
      n_changed
    ),
    "",
    strict_axis_text,
    "",
    "PART B | Nfkb1/Nfe2l2 -> Tnf/Icam1/Ccl2/Ccl3 CLOSURE:",
    paste0(
      "Motif-supported direct TF-ligand edges: ",
      sum(
        motif_audit$motif_supported_direct_target
      ),
      "/",
      nrow(
        motif_audit
      )
    ),
    paste0(
      "Strong direct two-sample closures: ",
      n_direct_strong
    ),
    "",
    closure_text,
    "",
    "PART C | STRICT FINAL MECHANISM MODEL:",
    final_model_text,
    "",
    "FINAL DECISION:",
    paste(
      capture.output(
        print(
          final_decision,
          row.names = FALSE
        )
      ),
      collapse = "\n"
    ),
    "",
    "INTERPRETATION BOUNDARY:",
    paste0(
      "The >=5% threshold is applied separately to every recipient sample. ",
      "Axes below this threshold are downgraded rather than deleted."
    ),
    paste0(
      "A direct motif-supported edge plus replicated AUC-expression association ",
      "is still an associative regulatory closure, not perturbational proof."
    ),
    paste0(
      "If Nfkb1 lacks a direct ligand edge, describe it as an Arm_2 state ",
      "regulator rather than as a proven direct driver of that ligand."
    ),
    paste0(
      "If Nfe2l2 is linked to an inflammatory ligand, do not automatically ",
      "interpret Nfe2l2 as pro-inflammatory; it may represent stress adaptation ",
      "coexisting with inflammatory output."
    ),
    "",
    "NEXT ACTION:",
    paste0(
      "Use the strict final mechanism model to freeze the GSE162610 discovery ",
      "signatures and axes. Then begin GSE142849 external long-term DAM/branch ",
      "validation without rerunning GSE162610 Steps 1-8D."
    ),
    "",
    "UPLOAD THESE FILES:",
    "1) reports/Step8D_mechanism_closure_QC_report.txt",
    "2) tables/03_strict_prioritized_mechanism_axes.tsv",
    "3) tables/06_focus_TF_ligand_closure_summary.tsv",
    "4) tables/07_strict_final_mechanism_model.tsv",
    "5) tables/08_final_decision.csv",
    "======================================================================"
  )

  writeLines(
    report_lines,
    con = REPORT_FILE,
    useBytes = TRUE
  )

  msg(
    "=============================================================="
  )
  msg(
    "STEP 8D MECHANISM CLOSURE AND QC COMPLETE"
  )
  msg(
    "Technical status: ",
    if (
      technical_pass
    ) {
      "PASS"
    } else {
      "FAIL"
    }
  )
  msg(
    "Strict Tier1 axes: ",
    n_strict_tier1
  )
  msg(
    "Strong direct TF-ligand closures: ",
    n_direct_strong
  )
  msg(
    "Report: ",
    safe_path(
      REPORT_FILE
    )
  )
  msg(
    "=============================================================="
  )

  invisible(
    result_object
  )
}

tryCatch(
  main(),
  error = function(e) {
    failed_call <- tryCatch(
      paste(
        deparse(
          sys.call(
            -1L
          )
        ),
        collapse = " "
      ),
      error = function(inner_error) {
        NA_character_
      }
    )

    error_lines <- c(
      "======================================================================",
      "GSE162610 | STEP 8D ERROR",
      "======================================================================",
      paste0(
        "Time: ",
        timestamp()
      ),
      paste0(
        "Message: ",
        conditionMessage(
          e
        )
      ),
      paste0(
        "Failed call: ",
        failed_call
      ),
      "",
      "Do not rerun Step1-Step8C.",
      "Correct only this Step8D file and rerun it.",
      "======================================================================"
    )

    writeLines(
      error_lines,
      con = ERROR_FILE,
      useBytes = TRUE
    )

    cat(
      paste0(
        error_lines,
        collapse = "\n"
      ),
      "\n"
    )

    stop(
      e
    )
  }
)
