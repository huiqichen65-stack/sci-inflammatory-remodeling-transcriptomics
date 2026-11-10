# ==============================================================================
# 47_step8C_GSE162610_recipient_program_integration_v2_VECTOR_FIX_ONE_FILE.R
# ==============================================================================
# GSE162610 adult SCI activated/DAM-like microglia project
#
# STEP 8C
#   Integrate the formally confirmed 7 dpi CellChat axes with transcriptional
#   response programs in recipient cells.
#
# INPUT EVIDENCE
#   - Step8A PASS manifest and full normalized RNA atlas.
#   - Step8B targeted full-cell nboot=100 confirmation.
#   - Only CONFIRMED_CORE axes are used as the primary input.
#
# ANALYSES
#   1) Score curated recipient-response programs in all eligible 1/3/7 dpi
#      recipient cells using deterministic gene-wise z scoring.
#   2) Summarize each program by biological sample and timepoint.
#   3) At 7 dpi, compare receptor-positive versus receptor-negative recipient
#      cells for each confirmed ligand-receptor axis.
#   4) Require directionally concordant support in both 7 dpi samples for the
#      strongest recipient-program evidence category.
#   5) Integrate branch, TF candidates, ligand, receptor, recipient and response
#      program into prioritized mechanism axes.
#
# INTERPRETATION BOUNDARY
#   - Receptor-positive versus receptor-negative comparisons are descriptive
#     single-cell associations, not biological-replicate hypothesis tests.
#   - CellChat plus recipient-program association does not prove direct signaling.
#   - Arm_1 and Arm_2 remain neutral names.
#   - Repair/scar terminology remains "biased" or "modulating" until external,
#     spatial and perturbation validation is complete.
#
# RSTUDIO
# source("<DOWNLOADS>/47_step8C_GSE162610_recipient_program_integration_v2_VECTOR_FIX_ONE_FILE.R")
# ==============================================================================

rm(list = ls())
gc()

options(
  stringsAsFactors = FALSE,
  width = 280,
  warn = 1,
  timeout = 86400,
  future.globals.maxSize = 32 * 1024^3
)

set.seed(20260803L)

cat(
  "R parser check PASS: the complete Step8C file was parsed by source().\n"
)

# ==============================================================================
# 0. SETTINGS
# ==============================================================================

PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
FORMAL_TIMEPOINTS <- c(
  "1dpi",
  "3dpi",
  "7dpi"
)

FORMAL_7DPI_SAMPLES_EXPECTED <- 2L

MIN_PROGRAM_GENES_PRESENT <- 4L
MIN_TARGET_CELLS_PER_SAMPLE <- 20L
MIN_RECEPTOR_POSITIVE_CELLS <- 10L
MIN_RECEPTOR_NEGATIVE_CELLS <- 10L
MIN_ABS_STANDARDIZED_EFFECT <- 0.10
MIN_ABS_SPEARMAN_RHO <- 0.05
MIN_RECEPTOR_POSITIVE_PERCENT <- 5

PRIMARY_CONFIRMATION_STATUS <- "CONFIRMED_CORE"

# ==============================================================================
# 1. INPUT PATHS
# ==============================================================================

ATLAS_RDS_FILE <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step1F_GSE162610_broad_annotation",
  "GSE162610_broad_annotation_audit.rds"
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

STEP8A_DECISION_FILE <- file.path(
  STEP8A_DIR,
  "tables",
  "08_final_decision.csv"
)

STEP8B_CONFIRM_DIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step8B_CONFIRM_GSE162610_7dpi_targeted_fullcell_nboot100"
)

STEP8B_CONFIRM_DECISION_FILE <- file.path(
  STEP8B_CONFIRM_DIR,
  "tables",
  "09_final_decision.csv"
)

STEP8B_CONFIRMED_AXES_FILE <- file.path(
  STEP8B_CONFIRM_DIR,
  "tables",
  "07_confirmed_core_axes.tsv"
)

STEP8B_PREVIEW_VS_FULL_FILE <- file.path(
  STEP8B_CONFIRM_DIR,
  "tables",
  "06_preview_vs_full_confirmation.tsv"
)

# ==============================================================================
# 2. OUTPUT PATHS
# ==============================================================================

OUTDIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step8C_GSE162610_recipient_program_integration"
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

AXIS_PROGRAM_MAP_FILE <- file.path(
  TABDIR,
  "01_confirmed_axis_to_recipient_program_map.tsv"
)

PROGRAM_GENE_COVERAGE_FILE <- file.path(
  TABDIR,
  "02_recipient_program_gene_coverage.tsv"
)

PROGRAM_TIMECOURSE_FILE <- file.path(
  TABDIR,
  "03_recipient_program_timecourse_by_sample.tsv"
)

PROGRAM_TIME_TREND_FILE <- file.path(
  TABDIR,
  "04_recipient_program_time_trend_summary.tsv"
)

RECEPTOR_AUDIT_LONG_FILE <- file.path(
  TABDIR,
  "05_7dpi_axis_program_receptor_audit_long.tsv"
)

AXIS_PROGRAM_EVIDENCE_FILE <- file.path(
  TABDIR,
  "06_axis_program_integrated_evidence.tsv"
)

PRIORITIZED_MECHANISM_AXES_FILE <- file.path(
  TABDIR,
  "07_prioritized_mechanism_axes.tsv"
)

TARGET_SUMMARY_FILE <- file.path(
  TABDIR,
  "08_recipient_target_summary.tsv"
)

FINAL_DECISION_FILE <- file.path(
  TABDIR,
  "09_final_decision.csv"
)

CELL_SCORE_RDS_FILE <- file.path(
  RDSDIR,
  "GSE162610_step8C_recipient_program_cell_scores.rds"
)

RESULT_RDS_FILE <- file.path(
  RDSDIR,
  "GSE162610_step8C_recipient_program_integration_results.rds"
)

REPORT_FILE <- file.path(
  REPORTDIR,
  "Step8C_recipient_program_integration_report.txt"
)

ERROR_FILE <- file.path(
  REPORTDIR,
  "Step8C_recipient_program_integration_ERROR.txt"
)

FIGURE_EFFECT_FILE <- file.path(
  FIGDIR,
  "Fig8C_01_axis_program_effects.pdf"
)

FIGURE_TIME_FILE <- file.path(
  FIGDIR,
  "Fig8C_02_recipient_program_timecourse.pdf"
)

FIGURE_PRIORITY_FILE <- file.path(
  FIGDIR,
  "Fig8C_03_prioritized_mechanism_axes.pdf"
)

# ==============================================================================
# 3. PACKAGES
# ==============================================================================

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "CellChat",
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
# 4. CURATED RECIPIENT PROGRAMS
# ==============================================================================
# These are directional, hypothesis-oriented gene sets. They are used only for
# response-program integration, not for redefining cell identities.

PROGRAM_GENE_SETS <- list(
  TNF_NFKB_RESPONSE = c(
    "Nfkbia",
    "Tnfaip3",
    "Relb",
    "Nfkb1",
    "Nfkb2",
    "Rela",
    "Birc3",
    "Traf1",
    "Socs3",
    "Ccl2",
    "Ccl7",
    "Cxcl2",
    "Cxcl10",
    "Icam1",
    "Vcam1",
    "Ptgs2",
    "Junb",
    "Fos",
    "Egr1",
    "Gadd45b"
  ),
  MYELOID_INFLAMMATION = c(
    "Il1b",
    "Tnf",
    "Ccl2",
    "Ccl3",
    "Ccl4",
    "Ccl7",
    "Cxcl2",
    "Cxcl10",
    "Nos2",
    "Ptgs2",
    "Nlrp3",
    "Stat1",
    "Irf1",
    "Irf7",
    "Socs3",
    "Cd14",
    "Lcn2",
    "Saa3",
    "Fcgr1",
    "Fcer1g"
  ),
  MYELOID_ADHESION_CHEMOTAXIS = c(
    "Cd44",
    "Lcp2",
    "Ptk2b",
    "Rac2",
    "Vav1",
    "Syk",
    "Lsp1",
    "Coro1a",
    "Cx3cr1",
    "Csf1r",
    "Tyrobp",
    "Fcer1g",
    "Was",
    "Dock2",
    "Rhoa",
    "Cdc42"
  ),
  EFFEROCYTOSIS_PHAGOCYTOSIS = c(
    "Mfge8",
    "Lrp1",
    "Cd36",
    "Marco",
    "Msr1",
    "Trem2",
    "Tyrobp",
    "Apoe",
    "Lpl",
    "Abca1",
    "Abcg1",
    "C1qa",
    "C1qb",
    "C1qc",
    "Gpnmb",
    "Fcgr1",
    "Fcgr3",
    "Lamp1",
    "Ctsb",
    "Ctsd",
    "Ctss",
    "Npc2"
  ),
  RESOLUTION_HOMEOSTASIS = c(
    "Il10",
    "Tgfb1",
    "Socs1",
    "Socs3",
    "Nfe2l2",
    "Hmox1",
    "Nqo1",
    "Gclc",
    "Gclm",
    "Abca1",
    "Apoe",
    "Pparg",
    "Nr4a1",
    "Dusp1",
    "Tnfaip3",
    "Klf4",
    "Maf",
    "Mafb"
  ),
  TAM_RECEPTOR_DOWNSTREAM = c(
    "Socs1",
    "Socs3",
    "Dusp1",
    "Dusp5",
    "Tnfaip3",
    "Ptpn6",
    "Ptpn11",
    "Cish",
    "Pik3r1",
    "Akt1",
    "Foxo1",
    "Abca1",
    "Apoe",
    "Lpl",
    "Bcl2l1",
    "Mapk1",
    "Mapk3"
  ),
  FIBROBLAST_ACTIVATION = c(
    "Pdgfra",
    "Pdgfrb",
    "Col1a1",
    "Col1a2",
    "Fn1",
    "Postn",
    "Ctgf",
    "Serpine1",
    "Thbs1",
    "Tnc",
    "Acta2",
    "Tagln",
    "Mmp2",
    "Mmp14",
    "Timp1"
  ),
  ECM_MATRIX_REMODELING = c(
    "Col1a1",
    "Col1a2",
    "Col3a1",
    "Col5a1",
    "Col5a2",
    "Fn1",
    "Postn",
    "Tnc",
    "Lox",
    "Loxl2",
    "Mmp2",
    "Mmp9",
    "Mmp14",
    "Timp1",
    "Timp2",
    "Serpine1",
    "Ctgf",
    "Thbs1",
    "Spp1",
    "Vim"
  ),
  COLLAGEN_FIBROSIS = c(
    "Col1a1",
    "Col1a2",
    "Col3a1",
    "Col5a1",
    "Col5a2",
    "Col6a1",
    "Col6a2",
    "Col6a3",
    "Pdgfra",
    "Pdgfrb",
    "Acta2",
    "Tagln",
    "Ctgf",
    "Tgfb1",
    "Lox",
    "Postn",
    "Dcn",
    "Lum",
    "Fbln1",
    "Fbln2"
  ),
  OSM_JAK_STAT_RESPONSE = c(
    "Socs3",
    "Stat3",
    "Jak1",
    "Jak2",
    "Cebpd",
    "Cebpb",
    "Junb",
    "Fos",
    "Serpine1",
    "Timp1",
    "Ccl2",
    "Icam1",
    "Lcn2",
    "Bcl3",
    "Pim1"
  ),
  VASCULAR_INFLAMMATION = c(
    "Icam1",
    "Vcam1",
    "Sele",
    "Selp",
    "Ccl2",
    "Cxcl1",
    "Cxcl2",
    "Cxcl10",
    "Nfkbia",
    "Tnfaip3",
    "Socs3",
    "Ptgs2",
    "Il6",
    "Edn1",
    "Serpine1"
  ),
  ENDOTHELIAL_ACTIVATION = c(
    "Icam1",
    "Vcam1",
    "Sele",
    "Selp",
    "Angpt2",
    "Plvap",
    "Edn1",
    "Serpine1",
    "Ccl2",
    "Cxcl10",
    "Ptgs2",
    "Nfkbia",
    "Tnfaip3",
    "Socs3"
  ),
  BARRIER_INTEGRITY = c(
    "Cldn5",
    "Ocln",
    "Tjp1",
    "Cdh5",
    "Klf2",
    "Klf4",
    "Nos3",
    "Mfsd2a",
    "Slc2a1",
    "Abcb1a",
    "Pecam1",
    "Tek"
  ),
  PERICYTE_ACTIVATION = c(
    "Rgs5",
    "Cspg4",
    "Pdgfrb",
    "Des",
    "Acta2",
    "Tagln",
    "Col4a1",
    "Col4a2",
    "Fn1",
    "Ctgf",
    "Serpine1",
    "Timp1",
    "Mmp14",
    "Vim"
  ),
  OPC_STRESS_RESPONSE = c(
    "Nfkbia",
    "Tnfaip3",
    "Socs3",
    "Stat1",
    "Irf1",
    "Junb",
    "Fos",
    "Ddit3",
    "Atf4",
    "Hspa1a",
    "Hspa1b",
    "Hmox1",
    "Gadd45b",
    "Btg2",
    "Ddit4"
  ),
  OLIGODENDROCYTE_DIFFERENTIATION = c(
    "Sox10",
    "Olig1",
    "Olig2",
    "Nkx2-2",
    "Myrf",
    "Tcf7l2",
    "Gpr17",
    "Enpp6",
    "Cnp",
    "Mbp",
    "Plp1",
    "Mag",
    "Mog",
    "Mal"
  ),
  T_CELL_CHEMOTAXIS = c(
    "Itgal",
    "Itgb2",
    "Sell",
    "Selplg",
    "Rac1",
    "Rac2",
    "Rhoa",
    "Cdc42",
    "Vav1",
    "Dock2",
    "Was",
    "Ptk2b",
    "Lcp2",
    "Syk",
    "Coro1a"
  ),
  T_CELL_ACTIVATION = c(
    "Cd3d",
    "Cd3e",
    "Cd247",
    "Lck",
    "Zap70",
    "Lat",
    "Trac",
    "Cd69",
    "Il2ra",
    "Icos",
    "Ctla4",
    "Pdcd1",
    "Ifng",
    "Tbx21",
    "Nfatc1",
    "Fos",
    "Jun",
    "Nr4a1"
  )
)

PROGRAM_ANNOTATION <- data.frame(
  response_program = names(
    PROGRAM_GENE_SETS
  ),
  program_class = c(
    "INFLAMMATORY",
    "INFLAMMATORY",
    "INFLAMMATORY",
    "MYELOID_CLEARANCE",
    "RESOLUTION_ADAPTATION",
    "TAM_SIGNALING",
    "SCAR_ECM",
    "SCAR_ECM",
    "SCAR_ECM",
    "OSM_RESPONSE",
    "VASCULAR",
    "VASCULAR",
    "VASCULAR_INTEGRITY",
    "VASCULAR",
    "OLIGODENDROGLIAL_STRESS",
    "OLIGODENDROGLIAL_DIFFERENTIATION",
    "T_CELL",
    "T_CELL"
  ),
  stringsAsFactors = FALSE
)

# ==============================================================================
# 5. HELPERS
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
      "GSE162610 | STEP 8C v2 ERROR",
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
    return(output)
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

load_cellchat_mouse_db <- function() {
  database <- tryCatch(
    get(
      "CellChatDB.mouse",
      envir = asNamespace(
        "CellChat"
      ),
      inherits = FALSE
    ),
    error = function(e) {
      NULL
    }
  )

  if (!is.null(database)) {
    return(database)
  }

  temporary_environment <- new.env(
    parent = emptyenv()
  )

  try(
    utils::data(
      "CellChatDB.mouse",
      package = "CellChat",
      envir = temporary_environment
    ),
    silent = TRUE
  )

  if (
    exists(
      "CellChatDB.mouse",
      envir = temporary_environment,
      inherits = FALSE
    )
  ) {
    return(
      get(
        "CellChatDB.mouse",
        envir = temporary_environment,
        inherits = FALSE
      )
    )
  }

  stop(
    "CellChatDB.mouse could not be loaded.",
    call. = FALSE
  )
}

entity_components <- function(
  entity,
  complex_table,
  expression_genes
) {
  components <- entity

  if (
    entity %in%
      rownames(
        complex_table
      )
  ) {
    components <- as.character(
      unlist(
        complex_table[
          entity,
          ,
          drop = FALSE
        ],
        use.names = FALSE
      )
    )

    components <- components[
      !is.na(
        components
      ) &
        nzchar(
          components
        )
    ]
  }

  intersect(
    unique(
      components
    ),
    expression_genes
  )
}

entity_cell_score <- function(
  entity,
  cells,
  expression_data,
  complex_table
) {
  cells <- intersect(
    cells,
    colnames(
      expression_data
    )
  )

  components <- entity_components(
    entity = entity,
    complex_table = complex_table,
    expression_genes = rownames(
      expression_data
    )
  )

  if (
    length(
      cells
    ) ==
      0L ||
    length(
      components
    ) ==
      0L
  ) {
    return(
      rep(
        NA_real_,
        length(
          cells
        )
      )
    )
  }

  current <- expression_data[
    components,
    cells,
    drop = FALSE
  ]

  if (
    length(
      components
    ) ==
      1L
  ) {
    return(
      as.numeric(
        current[
          1L,
          ,
          drop = TRUE
        ]
      )
    )
  }

  as.numeric(
    apply(
      as.matrix(
        current
      ),
      2L,
      min
    )
  )
}

score_program_within_target <- function(
  expression_data,
  cells,
  program_genes
) {
  cells <- intersect(
    cells,
    colnames(
      expression_data
    )
  )

  present_genes <- intersect(
    unique(
      program_genes
    ),
    rownames(
      expression_data
    )
  )

  if (
    length(
      cells
    ) ==
      0L ||
    length(
      present_genes
    ) <
      MIN_PROGRAM_GENES_PRESENT
  ) {
    return(
      list(
        score = rep(
          NA_real_,
          length(
            cells
          )
        ),
        cells = cells,
        present_genes = present_genes,
        variable_genes = character(0)
      )
    )
  }

  current <- as.matrix(
    expression_data[
      present_genes,
      cells,
      drop = FALSE
    ]
  )

  gene_means <- rowMeans(
    current
  )

  gene_sds <- apply(
    current,
    1L,
    stats::sd
  )

  variable <- is.finite(
    gene_sds
  ) &
    gene_sds >
      0

  current <- current[
    variable,
    ,
    drop = FALSE
  ]

  gene_means <- gene_means[
    variable
  ]

  gene_sds <- gene_sds[
    variable
  ]

  if (
    nrow(
      current
    ) <
      MIN_PROGRAM_GENES_PRESENT
  ) {
    return(
      list(
        score = rep(
          NA_real_,
          length(
            cells
          )
        ),
        cells = cells,
        present_genes = present_genes,
        variable_genes = rownames(
          current
        )
      )
    )
  }

  scaled <- sweep(
    current,
    1L,
    gene_means,
    "-"
  )

  scaled <- sweep(
    scaled,
    1L,
    gene_sds,
    "/"
  )

  score <- colMeans(
    scaled,
    na.rm = TRUE
  )

  list(
    score = as.numeric(
      score
    ),
    cells = cells,
    present_genes = present_genes,
    variable_genes = rownames(
      current
    )
  )
}

safe_wilcox_p <- function(
  positive_scores,
  negative_scores
) {
  positive_scores <- positive_scores[
    is.finite(
      positive_scores
    )
  ]

  negative_scores <- negative_scores[
    is.finite(
      negative_scores
    )
  ]

  if (
    length(
      positive_scores
    ) <
      2L ||
    length(
      negative_scores
    ) <
      2L
  ) {
    return(
      NA_real_
    )
  }

  tryCatch(
    stats::wilcox.test(
      positive_scores,
      negative_scores,
      exact = FALSE
    )$p.value,
    error = function(e) {
      NA_real_
    }
  )
}

safe_mean_finite <- function(x) {
  values <- x[
    is.finite(
      x
    )
  ]

  if (length(values) == 0L) {
    return(NA_real_)
  }

  mean(
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
    return(NA_real_)
  }

  max(
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
    return(NA_real_)
  }

  min(
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
      5L ||
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
        p_value = NA_real_
      )
    )
  }

  output <- tryCatch(
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

  if (is.null(output)) {
    return(
      c(
        rho = NA_real_,
        p_value = NA_real_
      )
    )
  }

  c(
    rho = unname(
      output$estimate
    ),
    p_value = output$p.value
  )
}

build_axis_program_map <- function(
  confirmed_axes
) {
  mapping_rows <- list()
  row_index <- 0L

  add_mapping <- function(
    axis_row,
    programs,
    expected_directions,
    interpretations
  ) {
    stopifnot(
      length(
        programs
      ) ==
        length(
          expected_directions
        ),
      length(
        programs
      ) ==
        length(
          interpretations
        )
    )

    output <- vector(
      "list",
      length(
        programs
      )
    )

    for (program_index in seq_along(programs)) {
      output[[program_index]] <- data.frame(
        target = axis_row$target,
        ligand = axis_row$ligand,
        receptor = axis_row$receptor,
        interaction_name = axis_row$interaction_name,
        pathway_name = axis_row$pathway_name,
        dominant_branch = axis_row$dominant_branch_full,
        mechanistic_priority_CellChat =
          axis_row$mechanistic_priority,
        Step7G_candidate_ligand = if (
          "Step7G_candidate_ligand_preview" %in%
            colnames(
              axis_row
            )
        ) {
          axis_row$Step7G_candidate_ligand_preview
        } else {
          NA
        },
        Step7G_candidate_TFs = if (
          "Step7G_candidate_TFs_preview" %in%
            colnames(
              axis_row
            )
        ) {
          axis_row$Step7G_candidate_TFs_preview
        } else {
          NA_character_
        },
        response_program = programs[
          program_index
        ],
        expected_direction = expected_directions[
          program_index
        ],
        expected_sign = ifelse(
          expected_directions[
            program_index
          ] ==
            "HIGHER_IN_RECEPTOR_POSITIVE",
          1,
          -1
        ),
        biological_interpretation =
          interpretations[
            program_index
          ],
        stringsAsFactors = FALSE
      )
    }

    output
  }

  for (
    axis_index in
    seq_len(
      nrow(
        confirmed_axes
      )
    )
  ) {
    axis_row <- confirmed_axes[
      axis_index,
      ,
      drop = FALSE
    ]

    target <- as.character(
      axis_row$target
    )

    pathway <- as.character(
      axis_row$pathway_name
    )

    programs <- character(0)
    directions <- character(0)
    interpretations <- character(0)

    if (
      target ==
        "Fibroblast" &&
      pathway ==
        "TNF"
    ) {
      programs <- c(
        "TNF_NFKB_RESPONSE",
        "FIBROBLAST_ACTIVATION",
        "ECM_MATRIX_REMODELING",
        "COLLAGEN_FIBROSIS"
      )
      directions <- rep(
        "HIGHER_IN_RECEPTOR_POSITIVE",
        4L
      )
      interpretations <- c(
        "TNF-responsive inflammatory activation",
        "fibroblast activation",
        "matrix remodeling",
        "collagen/fibrosis-associated program"
      )
    } else if (
      target ==
        "Fibroblast" &&
      pathway ==
        "GAS"
    ) {
      programs <- c(
        "TAM_RECEPTOR_DOWNSTREAM"
      )
      directions <- c(
        "HIGHER_IN_RECEPTOR_POSITIVE"
      )
      interpretations <- c(
        "GAS6-TAM downstream adaptation/homeostasis candidate"
      )
    } else if (
      target ==
        "Fibroblast" &&
      pathway ==
        "OSM"
    ) {
      programs <- c(
        "OSM_JAK_STAT_RESPONSE",
        "FIBROBLAST_ACTIVATION",
        "ECM_MATRIX_REMODELING"
      )
      directions <- rep(
        "HIGHER_IN_RECEPTOR_POSITIVE",
        3L
      )
      interpretations <- c(
        "OSM-JAK-STAT response",
        "fibroblast activation",
        "matrix remodeling"
      )
    } else if (
      target %in%
        c(
          "Microglia_other",
          "Other_myeloid"
        ) &&
      pathway ==
        "ICAM"
    ) {
      programs <- c(
        "MYELOID_ADHESION_CHEMOTAXIS",
        "MYELOID_INFLAMMATION"
      )
      directions <- rep(
        "HIGHER_IN_RECEPTOR_POSITIVE",
        2L
      )
      interpretations <- c(
        "myeloid adhesion/motility",
        "myeloid inflammatory activation"
      )
    } else if (
      target %in%
        c(
          "Microglia_other",
          "Other_myeloid"
        ) &&
      pathway ==
        "TNF"
    ) {
      programs <- c(
        "TNF_NFKB_RESPONSE",
        "MYELOID_INFLAMMATION"
      )
      directions <- rep(
        "HIGHER_IN_RECEPTOR_POSITIVE",
        2L
      )
      interpretations <- c(
        "TNF-responsive NF-kB activation",
        "myeloid inflammatory activation"
      )
    } else if (
      target ==
        "Microglia_other" &&
      pathway ==
        "GAS"
    ) {
      programs <- c(
        "EFFEROCYTOSIS_PHAGOCYTOSIS",
        "RESOLUTION_HOMEOSTASIS",
        "TAM_RECEPTOR_DOWNSTREAM"
      )
      directions <- rep(
        "HIGHER_IN_RECEPTOR_POSITIVE",
        3L
      )
      interpretations <- c(
        "efferocytosis/phagocytic clearance",
        "resolution and oxidative-stress adaptation",
        "GAS6-TAM downstream signaling"
      )
    } else if (
      target ==
        "Pericyte_VSMC" &&
      pathway ==
        "TNF"
    ) {
      programs <- c(
        "TNF_NFKB_RESPONSE",
        "VASCULAR_INFLAMMATION",
        "PERICYTE_ACTIVATION",
        "ECM_MATRIX_REMODELING"
      )
      directions <- rep(
        "HIGHER_IN_RECEPTOR_POSITIVE",
        4L
      )
      interpretations <- c(
        "TNF-responsive NF-kB activation",
        "vascular inflammatory activation",
        "pericyte/VSMC activation",
        "perivascular matrix remodeling"
      )
    } else if (
      target ==
        "Pericyte_VSMC" &&
      pathway ==
        "GAS"
    ) {
      programs <- c(
        "TAM_RECEPTOR_DOWNSTREAM",
        "RESOLUTION_HOMEOSTASIS"
      )
      directions <- rep(
        "HIGHER_IN_RECEPTOR_POSITIVE",
        2L
      )
      interpretations <- c(
        "GAS6-TAM downstream signaling",
        "homeostatic/adaptive response candidate"
      )
    } else if (
      target ==
        "Pericyte_VSMC" &&
      pathway ==
        "OSM"
    ) {
      programs <- c(
        "OSM_JAK_STAT_RESPONSE",
        "PERICYTE_ACTIVATION",
        "ECM_MATRIX_REMODELING"
      )
      directions <- rep(
        "HIGHER_IN_RECEPTOR_POSITIVE",
        3L
      )
      interpretations <- c(
        "OSM-JAK-STAT response",
        "pericyte/VSMC activation",
        "perivascular matrix remodeling"
      )
    } else if (
      target ==
        "Endothelial" &&
      pathway ==
        "TNF"
    ) {
      programs <- c(
        "TNF_NFKB_RESPONSE",
        "ENDOTHELIAL_ACTIVATION",
        "BARRIER_INTEGRITY"
      )
      directions <- c(
        "HIGHER_IN_RECEPTOR_POSITIVE",
        "HIGHER_IN_RECEPTOR_POSITIVE",
        "LOWER_IN_RECEPTOR_POSITIVE"
      )
      interpretations <- c(
        "TNF-responsive NF-kB activation",
        "endothelial inflammatory activation",
        "loss of vascular barrier-integrity program"
      )
    } else if (
      target ==
        "OPC" &&
      pathway ==
        "TNF"
    ) {
      programs <- c(
        "TNF_NFKB_RESPONSE",
        "OPC_STRESS_RESPONSE",
        "OLIGODENDROCYTE_DIFFERENTIATION"
      )
      directions <- c(
        "HIGHER_IN_RECEPTOR_POSITIVE",
        "HIGHER_IN_RECEPTOR_POSITIVE",
        "LOWER_IN_RECEPTOR_POSITIVE"
      )
      interpretations <- c(
        "TNF-responsive NF-kB activation",
        "OPC stress response",
        "reduced oligodendroglial differentiation program"
      )
    } else if (
      target ==
        "T_cell" &&
      pathway ==
        "TNF"
    ) {
      programs <- c(
        "TNF_NFKB_RESPONSE",
        "T_CELL_ACTIVATION"
      )
      directions <- rep(
        "HIGHER_IN_RECEPTOR_POSITIVE",
        2L
      )
      interpretations <- c(
        "TNF-responsive NF-kB activation",
        "T-cell activation"
      )
    } else if (
      target ==
        "T_cell" &&
      pathway ==
        "CCL"
    ) {
      programs <- c(
        "T_CELL_CHEMOTAXIS",
        "T_CELL_ACTIVATION"
      )
      directions <- rep(
        "HIGHER_IN_RECEPTOR_POSITIVE",
        2L
      )
      interpretations <- c(
        "chemotaxis/migration",
        "T-cell activation"
      )
    }

    if (length(programs) == 0L) {
      programs <- "TNF_NFKB_RESPONSE"
      directions <- "HIGHER_IN_RECEPTOR_POSITIVE"
      interpretations <- paste0(
        "generic response context for ",
        pathway,
        " signaling"
      )
    }

    new_rows <- add_mapping(
      axis_row = axis_row,
      programs = programs,
      expected_directions = directions,
      interpretations = interpretations
    )

    for (new_row in new_rows) {
      row_index <- row_index + 1L
      mapping_rows[[row_index]] <- new_row
    }
  }

  output <- do.call(
    rbind,
    mapping_rows
  )

  output <- merge(
    output,
    PROGRAM_ANNOTATION,
    by = "response_program",
    all.x = TRUE,
    sort = FALSE
  )

  output <- output[
    order(
      match(
        output$mechanistic_priority_CellChat,
        c(
          "TIER_1_CORE",
          "TIER_2_CONFIRMED_CONTEXT"
        )
      ),
      output$target,
      output$interaction_name,
      output$response_program
    ),
    ,
    drop = FALSE
  ]

  rownames(
    output
  ) <- NULL

  output
}

# ==============================================================================
# 6. MAIN
# ==============================================================================

main <- function() {
  start_time <- Sys.time()

  msg(
    "=============================================================="
  )
  msg(
    "GSE162610 Step 8C v2 | recipient-program integration"
  )
  msg(
    "=============================================================="
  )

  required_files <- c(
    ATLAS_RDS_FILE,
    STEP8A_MANIFEST_FILE,
    STEP8A_DECISION_FILE,
    STEP8B_CONFIRM_DECISION_FILE,
    STEP8B_CONFIRMED_AXES_FILE,
    STEP8B_PREVIEW_VS_FULL_FILE
  )

  input_audit <- data.frame(
    role = c(
      "Full normalized atlas",
      "Step8A input manifest",
      "Step8A final decision",
      "Step8B confirmation final decision",
      "Step8B confirmed core axes",
      "Step8B preview versus full table"
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

  step8a_decision <- data.table::fread(
    STEP8A_DECISION_FILE,
    data.table = FALSE,
    showProgress = FALSE
  )

  step8b_decision <- data.table::fread(
    STEP8B_CONFIRM_DECISION_FILE,
    data.table = FALSE,
    showProgress = FALSE
  )

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
      step8b_decision,
      "Overall Step 8B targeted confirmation"
    )
  ) {
    abort_step(
      "Step8B targeted confirmation PASS was not found."
    )
  }

  confirmed_axes <- data.table::fread(
    STEP8B_CONFIRMED_AXES_FILE,
    data.table = FALSE,
    showProgress = FALSE
  )

  required_axis_columns <- c(
    "target",
    "ligand",
    "receptor",
    "interaction_name",
    "pathway_name",
    "dominant_branch_full",
    "confirmation_status",
    "mechanistic_priority"
  )

  missing_axis_columns <- setdiff(
    required_axis_columns,
    colnames(
      confirmed_axes
    )
  )

  if (length(missing_axis_columns) > 0L) {
    abort_step(
      paste0(
        "Confirmed-axis table is missing: ",
        paste(
          missing_axis_columns,
          collapse = ", "
        )
      )
    )
  }

  confirmed_axes <- confirmed_axes[
    confirmed_axes$confirmation_status ==
      PRIMARY_CONFIRMATION_STATUS,
    ,
    drop = FALSE
  ]

  if (nrow(confirmed_axes) == 0L) {
    abort_step(
      "No CONFIRMED_CORE axes were available for Step8C."
    )
  }

  axis_program_map <- build_axis_program_map(
    confirmed_axes
  )

  write_tsv(
    axis_program_map,
    AXIS_PROGRAM_MAP_FILE
  )

  axis_key <- paste(
    confirmed_axes$target,
    confirmed_axes$interaction_name,
    sep = "||"
  )

  mapped_axis_key <- paste(
    axis_program_map$target,
    axis_program_map$interaction_name,
    sep = "||"
  )

  unmapped_axes <- setdiff(
    unique(
      axis_key
    ),
    unique(
      mapped_axis_key
    )
  )

  if (length(unmapped_axes) > 0L) {
    abort_step(
      paste0(
        "Recipient-program mapping is missing for: ",
        paste(
          unmapped_axes,
          collapse = ", "
        )
      )
    )
  }

  manifest <- readRDS(
    STEP8A_MANIFEST_FILE
  )

  required_manifest_fields <- c(
    "status",
    "cell_metadata",
    "eligible_groups_by_time"
  )

  missing_manifest_fields <- setdiff(
    required_manifest_fields,
    names(
      manifest
    )
  )

  if (length(missing_manifest_fields) > 0L) {
    abort_step(
      paste0(
        "Step8A manifest is missing: ",
        paste(
          missing_manifest_fields,
          collapse = ", "
        )
      )
    )
  }

  if (
    !identical(
      manifest$status,
      "PASS"
    )
  ) {
    abort_step(
      "Step8A manifest status is not PASS."
    )
  }

  metadata <- as.data.frame(
    manifest$cell_metadata,
    stringsAsFactors = FALSE
  )

  required_metadata_columns <- c(
    "cell",
    "sample_id",
    "timepoint_standard",
    "cellchat_group",
    "include_step8_timecourse"
  )

  missing_metadata_columns <- setdiff(
    required_metadata_columns,
    colnames(
      metadata
    )
  )

  if (length(missing_metadata_columns) > 0L) {
    abort_step(
      paste0(
        "Step8A metadata is missing: ",
        paste(
          missing_metadata_columns,
          collapse = ", "
        )
      )
    )
  }

  metadata$cell <- as.character(
    metadata$cell
  )
  metadata$sample_id <- as.character(
    metadata$sample_id
  )
  metadata$timepoint_standard <- as.character(
    metadata$timepoint_standard
  )
  metadata$cellchat_group <- as.character(
    metadata$cellchat_group
  )
  metadata$include_step8_timecourse <- as_logical_safe(
    metadata$include_step8_timecourse
  )

  target_groups <- sort(
    unique(
      as.character(
        confirmed_axes$target
      )
    )
  )

  metadata_target <- metadata[
    metadata$include_step8_timecourse %in%
      TRUE &
    metadata$timepoint_standard %in%
      FORMAL_TIMEPOINTS &
    metadata$cellchat_group %in%
      target_groups,
    ,
    drop = FALSE
  ]

  metadata_target <- metadata_target[
    order(
      metadata_target$cellchat_group,
      metadata_target$timepoint_standard,
      metadata_target$sample_id,
      metadata_target$cell
    ),
    ,
    drop = FALSE
  ]

  if (anyDuplicated(metadata_target$cell) > 0L) {
    abort_step(
      "Duplicated recipient cells were found."
    )
  }

  missing_target_groups <- setdiff(
    target_groups,
    unique(
      metadata_target$cellchat_group
    )
  )

  if (length(missing_target_groups) > 0L) {
    abort_step(
      paste0(
        "Recipient group(s) are absent from timecourse metadata: ",
        paste(
          missing_target_groups,
          collapse = ", "
        )
      )
    )
  }

  samples_7dpi <- sort(
    unique(
      metadata_target$sample_id[
        metadata_target$timepoint_standard ==
          "7dpi"
      ]
    )
  )

  if (
    length(
      samples_7dpi
    ) !=
      FORMAL_7DPI_SAMPLES_EXPECTED
  ) {
    abort_step(
      paste0(
        "Expected ",
        FORMAL_7DPI_SAMPLES_EXPECTED,
        " 7 dpi samples but found ",
        length(
          samples_7dpi
        ),
        "."
      )
    )
  }

  msg(
    "Loading normalized RNA data for ",
    nrow(
      metadata_target
    ),
    " recipient cells..."
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

  missing_cells <- setdiff(
    metadata_target$cell,
    colnames(
      expression_data
    )
  )

  if (length(missing_cells) > 0L) {
    abort_step(
      paste0(
        "Normalized RNA data are missing ",
        length(
          missing_cells
        ),
        " recipient cells."
      )
    )
  }

  all_program_genes <- sort(
    unique(
      unlist(
        PROGRAM_GENE_SETS,
        use.names = FALSE
      )
    )
  )

  CellChatDB_mouse <- load_cellchat_mouse_db()

  complex_table <- as.data.frame(
    CellChatDB_mouse$complex,
    stringsAsFactors = FALSE
  )

  all_receptor_components <- sort(
    unique(
      unlist(
        lapply(
          unique(
            confirmed_axes$receptor
          ),
          entity_components,
          complex_table = complex_table,
          expression_genes = rownames(
            expression_data
          )
        ),
        use.names = FALSE
      )
    )
  )

  required_genes <- intersect(
    unique(
      c(
        all_program_genes,
        all_receptor_components
      )
    ),
    rownames(
      expression_data
    )
  )

  expression_subset <- expression_data[
    required_genes,
    metadata_target$cell,
    drop = FALSE
  ]

  rm(
    atlas,
    expression_data
  )
  gc()

  # --------------------------------------------------------------------------
  # Program gene coverage
  # --------------------------------------------------------------------------

  coverage_rows <- vector(
    "list",
    length(
      PROGRAM_GENE_SETS
    )
  )

  for (
    program_index in
    seq_along(
      PROGRAM_GENE_SETS
    )
  ) {
    program_name <- names(
      PROGRAM_GENE_SETS
    )[
      program_index
    ]

    requested_genes <- unique(
      PROGRAM_GENE_SETS[[program_name]]
    )

    present_genes <- intersect(
      requested_genes,
      rownames(
        expression_subset
      )
    )

    coverage_rows[[program_index]] <- data.frame(
      response_program = program_name,
      n_requested = length(
        requested_genes
      ),
      n_present = length(
        present_genes
      ),
      coverage_fraction = length(
        present_genes
      ) /
        length(
          requested_genes
        ),
      sufficient_for_scoring =
        length(
          present_genes
        ) >=
          MIN_PROGRAM_GENES_PRESENT,
      present_genes = paste(
        present_genes,
        collapse = " | "
      ),
      missing_genes = paste(
        setdiff(
          requested_genes,
          present_genes
        ),
        collapse = " | "
      ),
      stringsAsFactors = FALSE
    )
  }

  program_gene_coverage <- do.call(
    rbind,
    coverage_rows
  )

  program_gene_coverage <- merge(
    program_gene_coverage,
    PROGRAM_ANNOTATION,
    by = "response_program",
    all.x = TRUE,
    sort = FALSE
  )

  write_tsv(
    program_gene_coverage,
    PROGRAM_GENE_COVERAGE_FILE
  )

  used_programs <- unique(
    axis_program_map$response_program
  )

  insufficient_used_programs <- program_gene_coverage$response_program[
    program_gene_coverage$response_program %in%
      used_programs &
    !program_gene_coverage$sufficient_for_scoring
  ]

  if (length(insufficient_used_programs) > 0L) {
    abort_step(
      paste0(
        "Insufficient gene coverage for used program(s): ",
        paste(
          insufficient_used_programs,
          collapse = ", "
        )
      )
    )
  }

  # --------------------------------------------------------------------------
  # Deterministic program scoring within each recipient cell type
  # --------------------------------------------------------------------------

  msg(
    "Scoring recipient programs within each target cell type..."
  )

  relevant_programs_by_target <- split(
    axis_program_map$response_program,
    axis_program_map$target
  )

  relevant_programs_by_target <- lapply(
    relevant_programs_by_target,
    unique
  )

  cell_score_parts <- list()
  score_index <- 0L

  for (target_group in names(relevant_programs_by_target)) {
    target_cells <- metadata_target$cell[
      metadata_target$cellchat_group ==
        target_group
    ]

    target_metadata <- metadata_target[
      match(
        target_cells,
        metadata_target$cell
      ),
      c(
        "cell",
        "sample_id",
        "timepoint_standard",
        "cellchat_group"
      ),
      drop = FALSE
    ]

    for (
      program_name in
      relevant_programs_by_target[[target_group]]
    ) {
      score_result <- score_program_within_target(
        expression_data = expression_subset,
        cells = target_cells,
        program_genes = PROGRAM_GENE_SETS[[program_name]]
      )

      if (
        length(
          score_result$score
        ) !=
          nrow(
            target_metadata
          )
      ) {
        abort_step(
          paste0(
            "Program score length mismatch for ",
            target_group,
            " / ",
            program_name,
            "."
          )
        )
      }

      score_index <- score_index + 1L

      cell_score_parts[[score_index]] <- data.frame(
        cell = target_metadata$cell,
        sample_id = target_metadata$sample_id,
        timepoint = target_metadata$timepoint_standard,
        target = target_group,
        response_program = program_name,
        program_score = score_result$score,
        n_present_genes = length(
          score_result$present_genes
        ),
        n_variable_genes = length(
          score_result$variable_genes
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  recipient_program_cell_scores <- do.call(
    rbind,
    cell_score_parts
  )

  rownames(
    recipient_program_cell_scores
  ) <- NULL

  saveRDS(
    recipient_program_cell_scores,
    CELL_SCORE_RDS_FILE,
    compress = FALSE
  )

  # --------------------------------------------------------------------------
  # Sample/timepoint summaries
  # --------------------------------------------------------------------------

  program_timecourse <- data.table::as.data.table(
    recipient_program_cell_scores
  )[
    is.finite(
      program_score
    ),
    .(
      n_cells = .N,
      mean_score = mean(
        program_score
      ),
      median_score = stats::median(
        program_score
      ),
      sd_score = stats::sd(
        program_score
      ),
      se_score = stats::sd(
        program_score
      ) /
        sqrt(
          .N
        )
    ),
    by = .(
      target,
      response_program,
      timepoint,
      sample_id
    )
  ]

  program_timecourse <- as.data.frame(
    program_timecourse,
    stringsAsFactors = FALSE
  )

  write_tsv(
    program_timecourse,
    PROGRAM_TIMECOURSE_FILE
  )

  timepoint_summary <- data.table::as.data.table(
    program_timecourse
  )[
    ,
    .(
      n_samples = .N,
      mean_of_sample_means = mean(
        mean_score
      ),
      sd_of_sample_means = stats::sd(
        mean_score
      )
    ),
    by = .(
      target,
      response_program,
      timepoint
    )
  ]

  timepoint_wide <- data.table::dcast(
    timepoint_summary,
    target +
      response_program ~
      timepoint,
    value.var = "mean_of_sample_means"
  )

  for (
    required_timepoint in
    FORMAL_TIMEPOINTS
  ) {
    if (
      !required_timepoint %in%
        colnames(
          timepoint_wide
        )
    ) {
      timepoint_wide[[required_timepoint]] <- NA_real_
    }
  }

  timepoint_wide$delta_7dpi_minus_1dpi <-
    timepoint_wide[["7dpi"]] -
    timepoint_wide[["1dpi"]]

  timepoint_wide$delta_7dpi_minus_3dpi <-
    timepoint_wide[["7dpi"]] -
    timepoint_wide[["3dpi"]]

  program_time_trend <- as.data.frame(
    timepoint_wide,
    stringsAsFactors = FALSE
  )

  write_tsv(
    program_time_trend,
    PROGRAM_TIME_TREND_FILE
  )

  # --------------------------------------------------------------------------
  # 7 dpi receptor-positive versus receptor-negative program audit
  # --------------------------------------------------------------------------

  msg(
    "Auditing receptor-program associations in both 7 dpi samples..."
  )

  seven_dpi_scores <- recipient_program_cell_scores[
    recipient_program_cell_scores$timepoint ==
      "7dpi",
    ,
    drop = FALSE
  ]

  receptor_score_cache <- new.env(
    parent = emptyenv()
  )

  get_receptor_score <- function(
    receptor,
    target,
    sample_id
  ) {
    cache_key <- paste(
      receptor,
      target,
      sample_id,
      sep = "||"
    )

    if (
      exists(
        cache_key,
        envir = receptor_score_cache,
        inherits = FALSE
      )
    ) {
      return(
        get(
          cache_key,
          envir = receptor_score_cache,
          inherits = FALSE
        )
      )
    }

    cells <- metadata_target$cell[
      metadata_target$timepoint_standard ==
        "7dpi" &
      metadata_target$sample_id ==
        sample_id &
      metadata_target$cellchat_group ==
        target
    ]

    scores <- entity_cell_score(
      entity = receptor,
      cells = cells,
      expression_data = expression_subset,
      complex_table = complex_table
    )

    output <- data.frame(
      cell = cells,
      receptor_score = scores,
      receptor_positive =
        is.finite(
          scores
        ) &
        scores >
          0,
      stringsAsFactors = FALSE
    )

    assign(
      cache_key,
      output,
      envir = receptor_score_cache
    )

    output
  }

  audit_rows <- list()
  audit_index <- 0L

  for (
    mapping_index in
    seq_len(
      nrow(
        axis_program_map
      )
    )
  ) {
    mapping_row <- axis_program_map[
      mapping_index,
      ,
      drop = FALSE
    ]

    for (sample_id in samples_7dpi) {
      current_scores <- seven_dpi_scores[
        seven_dpi_scores$sample_id ==
          sample_id &
        seven_dpi_scores$target ==
          mapping_row$target &
        seven_dpi_scores$response_program ==
          mapping_row$response_program,
        ,
        drop = FALSE
      ]

      receptor_scores <- get_receptor_score(
        receptor = mapping_row$receptor,
        target = mapping_row$target,
        sample_id = sample_id
      )

      current <- merge(
        current_scores,
        receptor_scores,
        by = "cell",
        all = FALSE,
        sort = FALSE
      )

      positive_scores <- current$program_score[
        current$receptor_positive
      ]

      negative_scores <- current$program_score[
        !current$receptor_positive
      ]

      positive_scores <- positive_scores[
        is.finite(
          positive_scores
        )
      ]

      negative_scores <- negative_scores[
        is.finite(
          negative_scores
        )
      ]

      effect <- if (
        length(
          positive_scores
        ) >
          0L &&
        length(
          negative_scores
        ) >
          0L
      ) {
        mean(
          positive_scores
        ) -
          mean(
            negative_scores
          )
      } else {
        NA_real_
      }

      correlation <- safe_spearman(
        current$receptor_score,
        current$program_score
      )

      n_target_cells <- nrow(
        current
      )

      n_receptor_positive <- sum(
        current$receptor_positive,
        na.rm = TRUE
      )

      n_receptor_negative <- sum(
        !current$receptor_positive,
        na.rm = TRUE
      )

      receptor_positive_percent <- if (
        n_target_cells >
          0L
      ) {
        n_receptor_positive /
          n_target_cells *
          100
      } else {
        NA_real_
      }

      evaluable <-
        n_target_cells >=
          MIN_TARGET_CELLS_PER_SAMPLE &&
        n_receptor_positive >=
          MIN_RECEPTOR_POSITIVE_CELLS &&
        n_receptor_negative >=
          MIN_RECEPTOR_NEGATIVE_CELLS &&
        is.finite(
          effect
        )

      expected_effect_match <-
        evaluable &&
        mapping_row$expected_sign *
          effect >
          0

      expected_correlation_match <-
        is.finite(
          correlation[["rho"]]
        ) &&
        mapping_row$expected_sign *
          correlation[["rho"]] >
          0

      sample_support <-
        evaluable &&
        expected_effect_match &&
        (
          abs(
            effect
          ) >=
            MIN_ABS_STANDARDIZED_EFFECT ||
          (
            expected_correlation_match &&
            abs(
              correlation[["rho"]]
            ) >=
              MIN_ABS_SPEARMAN_RHO
          )
        )

      audit_index <- audit_index + 1L

      audit_rows[[audit_index]] <- data.frame(
        sample_id = sample_id,
        target = mapping_row$target,
        ligand = mapping_row$ligand,
        receptor = mapping_row$receptor,
        interaction_name = mapping_row$interaction_name,
        pathway_name = mapping_row$pathway_name,
        dominant_branch = mapping_row$dominant_branch,
        response_program = mapping_row$response_program,
        program_class = mapping_row$program_class,
        expected_direction = mapping_row$expected_direction,
        expected_sign = mapping_row$expected_sign,
        n_target_cells = n_target_cells,
        n_receptor_positive = n_receptor_positive,
        n_receptor_negative = n_receptor_negative,
        receptor_positive_percent =
          receptor_positive_percent,
        receptor_positive_mean_program_score = if (
          length(
            positive_scores
          ) >
            0L
        ) {
          mean(
            positive_scores
          )
        } else {
          NA_real_
        },
        receptor_negative_mean_program_score = if (
          length(
            negative_scores
          ) >
            0L
        ) {
          mean(
            negative_scores
          )
        } else {
          NA_real_
        },
        effect_receptor_positive_minus_negative =
          effect,
        effect_in_expected_direction =
          mapping_row$expected_sign *
          effect,
        wilcox_cell_level_p = safe_wilcox_p(
          positive_scores,
          negative_scores
        ),
        spearman_rho_receptor_vs_program =
          correlation[["rho"]],
        spearman_cell_level_p =
          correlation[["p_value"]],
        evaluable = evaluable,
        expected_effect_match =
          expected_effect_match,
        expected_correlation_match =
          expected_correlation_match,
        sample_support = sample_support,
        stringsAsFactors = FALSE
      )
    }
  }

  receptor_audit_long <- do.call(
    rbind,
    audit_rows
  )

  receptor_audit_long$wilcox_cell_level_fdr <- stats::p.adjust(
    receptor_audit_long$wilcox_cell_level_p,
    method = "BH"
  )

  receptor_audit_long$spearman_cell_level_fdr <- stats::p.adjust(
    receptor_audit_long$spearman_cell_level_p,
    method = "BH"
  )

  write_tsv(
    receptor_audit_long,
    RECEPTOR_AUDIT_LONG_FILE
  )

  audit_summary <- data.table::as.data.table(
    receptor_audit_long
  )[
    ,
    .(
      n_samples = data.table::uniqueN(
        sample_id
      ),
      n_samples_evaluable = sum(
        evaluable,
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
      n_samples_supported = sum(
        sample_support,
        na.rm = TRUE
      ),
      mean_effect_receptor_positive_minus_negative =
        safe_mean_finite(
          effect_receptor_positive_minus_negative
        ),
      mean_effect_in_expected_direction =
        safe_mean_finite(
          effect_in_expected_direction
        ),
      mean_spearman_rho =
        safe_mean_finite(
          spearman_rho_receptor_vs_program
        ),
      minimum_receptor_positive_percent =
        safe_min_finite(
          receptor_positive_percent
        ),
      sample_effects = paste(
        paste0(
          sample_id,
          ":",
          round(
            effect_receptor_positive_minus_negative,
            4
          )
        ),
        collapse = " | "
      ),
      sample_rhos = paste(
        paste0(
          sample_id,
          ":",
          round(
            spearman_rho_receptor_vs_program,
            4
          )
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

  audit_summary <- as.data.frame(
    audit_summary,
    stringsAsFactors = FALSE
  )

  audit_summary$recipient_program_status <- ifelse(
    audit_summary$n_samples_evaluable ==
      FORMAL_7DPI_SAMPLES_EXPECTED &
    audit_summary$n_samples_supported ==
      FORMAL_7DPI_SAMPLES_EXPECTED,
    "SUPPORTED_IN_BOTH_7DPI_SAMPLES",
    ifelse(
      audit_summary$n_samples_supported >
        0L,
      "SUPPORTED_IN_ONE_7DPI_SAMPLE",
      ifelse(
        audit_summary$n_samples_effect_direction_match >
          0L,
        "DIRECTION_ONLY",
        ifelse(
          audit_summary$n_samples_evaluable ==
            0L,
          "NOT_EVALUABLE",
          "NOT_SUPPORTED"
        )
      )
    )
  )

  # --------------------------------------------------------------------------
  # Integrate CellChat, recipient audit and timecourse context
  # --------------------------------------------------------------------------

  axis_program_evidence <- merge(
    axis_program_map,
    audit_summary,
    by = c(
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
    ),
    all.x = TRUE,
    sort = FALSE
  )

  axis_program_evidence <- merge(
    axis_program_evidence,
    program_time_trend,
    by = c(
      "target",
      "response_program"
    ),
    all.x = TRUE,
    sort = FALSE
  )

  axis_program_evidence$temporal_7dpi_vs_1dpi_match <-
    is.finite(
      axis_program_evidence$delta_7dpi_minus_1dpi
    ) &
    axis_program_evidence$expected_sign *
      axis_program_evidence$delta_7dpi_minus_1dpi >
      0

  axis_program_evidence$temporal_7dpi_vs_3dpi_match <-
    is.finite(
      axis_program_evidence$delta_7dpi_minus_3dpi
    ) &
    axis_program_evidence$expected_sign *
      axis_program_evidence$delta_7dpi_minus_3dpi >
      0

  axis_program_evidence$program_evidence_tier <- ifelse(
    axis_program_evidence$recipient_program_status ==
      "SUPPORTED_IN_BOTH_7DPI_SAMPLES",
    "PROGRAM_TIER_1_TWO_SAMPLE",
    ifelse(
      axis_program_evidence$recipient_program_status ==
        "SUPPORTED_IN_ONE_7DPI_SAMPLE",
      "PROGRAM_TIER_2_ONE_SAMPLE",
      ifelse(
        axis_program_evidence$recipient_program_status ==
          "DIRECTION_ONLY",
        "PROGRAM_TIER_3_DIRECTION_ONLY",
        "PROGRAM_TIER_4_NO_SUPPORT"
      )
    )
  )

  axis_program_evidence <- axis_program_evidence[
    order(
      match(
        axis_program_evidence$program_evidence_tier,
        c(
          "PROGRAM_TIER_1_TWO_SAMPLE",
          "PROGRAM_TIER_2_ONE_SAMPLE",
          "PROGRAM_TIER_3_DIRECTION_ONLY",
          "PROGRAM_TIER_4_NO_SUPPORT"
        )
      ),
      -axis_program_evidence$mean_effect_in_expected_direction
    ),
    ,
    drop = FALSE
  ]

  write_tsv(
    axis_program_evidence,
    AXIS_PROGRAM_EVIDENCE_FILE
  )

  # --------------------------------------------------------------------------
  # Axis-level mechanism prioritization
  # --------------------------------------------------------------------------

  axis_level <- data.table::as.data.table(
    axis_program_evidence
  )[
    ,
    .(
      n_programs_tested = .N,
      n_programs_two_sample = sum(
        recipient_program_status ==
          "SUPPORTED_IN_BOTH_7DPI_SAMPLES",
        na.rm = TRUE
      ),
      n_programs_one_sample = sum(
        recipient_program_status ==
          "SUPPORTED_IN_ONE_7DPI_SAMPLE",
        na.rm = TRUE
      ),
      n_programs_direction_only = sum(
        recipient_program_status ==
          "DIRECTION_ONLY",
        na.rm = TRUE
      ),
      top_supported_program = {
        tier_order <- match(
          program_evidence_tier,
          c(
            "PROGRAM_TIER_1_TWO_SAMPLE",
            "PROGRAM_TIER_2_ONE_SAMPLE",
            "PROGRAM_TIER_3_DIRECTION_ONLY",
            "PROGRAM_TIER_4_NO_SUPPORT"
          )
        )

        effect_order <- -mean_effect_in_expected_direction

        best_index <- order(
          tier_order,
          effect_order,
          na.last = TRUE
        )[[1L]]

        response_program[
          best_index
        ]
      },
      top_program_class = {
        tier_order <- match(
          program_evidence_tier,
          c(
            "PROGRAM_TIER_1_TWO_SAMPLE",
            "PROGRAM_TIER_2_ONE_SAMPLE",
            "PROGRAM_TIER_3_DIRECTION_ONLY",
            "PROGRAM_TIER_4_NO_SUPPORT"
          )
        )

        effect_order <- -mean_effect_in_expected_direction

        best_index <- order(
          tier_order,
          effect_order,
          na.last = TRUE
        )[[1L]]

        program_class[
          best_index
        ]
      },
      top_program_status = {
        tier_order <- match(
          program_evidence_tier,
          c(
            "PROGRAM_TIER_1_TWO_SAMPLE",
            "PROGRAM_TIER_2_ONE_SAMPLE",
            "PROGRAM_TIER_3_DIRECTION_ONLY",
            "PROGRAM_TIER_4_NO_SUPPORT"
          )
        )

        effect_order <- -mean_effect_in_expected_direction

        best_index <- order(
          tier_order,
          effect_order,
          na.last = TRUE
        )[[1L]]

        recipient_program_status[
          best_index
        ]
      },
      top_mean_effect_in_expected_direction =
        safe_max_finite(
          mean_effect_in_expected_direction
        ),
      supported_programs_two_sample = collapse_unique(
        response_program[
          recipient_program_status ==
            "SUPPORTED_IN_BOTH_7DPI_SAMPLES"
        ]
      ),
      supported_programs_any = collapse_unique(
        response_program[
          recipient_program_status %in%
            c(
              "SUPPORTED_IN_BOTH_7DPI_SAMPLES",
              "SUPPORTED_IN_ONE_7DPI_SAMPLE"
            )
        ]
      ),
      temporal_support_7dpi_vs_1dpi = any(
        temporal_7dpi_vs_1dpi_match,
        na.rm = TRUE
      ),
      temporal_support_7dpi_vs_3dpi = any(
        temporal_7dpi_vs_3dpi_match,
        na.rm = TRUE
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

  axis_level <- as.data.frame(
    axis_level,
    stringsAsFactors = FALSE
  )

  axis_level$recipient_mechanism_status <- ifelse(
    axis_level$n_programs_two_sample >
      0L,
    "RECIPIENT_PROGRAM_CONFIRMED_TWO_SAMPLE",
    ifelse(
      axis_level$n_programs_one_sample >
        0L,
      "RECIPIENT_PROGRAM_PARTIAL",
      ifelse(
        axis_level$n_programs_direction_only >
          0L,
        "RECIPIENT_PROGRAM_DIRECTION_ONLY",
        "COMMUNICATION_ONLY"
      )
    )
  )

  axis_level$final_mechanism_tier <- ifelse(
    axis_level$mechanistic_priority_CellChat ==
      "TIER_1_CORE" &
    axis_level$recipient_mechanism_status ==
      "RECIPIENT_PROGRAM_CONFIRMED_TWO_SAMPLE",
    "MECHANISM_TIER_1",
    ifelse(
      axis_level$recipient_mechanism_status %in%
        c(
          "RECIPIENT_PROGRAM_CONFIRMED_TWO_SAMPLE",
          "RECIPIENT_PROGRAM_PARTIAL"
        ),
      "MECHANISM_TIER_2",
      ifelse(
        axis_level$recipient_mechanism_status ==
          "RECIPIENT_PROGRAM_DIRECTION_ONLY",
        "MECHANISM_TIER_3",
        "COMMUNICATION_ONLY"
      )
    )
  )

  axis_level$mechanism_chain <- paste0(
    axis_level$dominant_branch,
    " -> ",
    axis_level$ligand,
    " -> ",
    axis_level$receptor,
    " -> ",
    axis_level$target,
    " -> ",
    axis_level$top_supported_program
  )

  axis_level <- axis_level[
    order(
      match(
        axis_level$final_mechanism_tier,
        c(
          "MECHANISM_TIER_1",
          "MECHANISM_TIER_2",
          "MECHANISM_TIER_3",
          "COMMUNICATION_ONLY"
        )
      ),
      -axis_level$n_programs_two_sample,
      -axis_level$top_mean_effect_in_expected_direction
    ),
    ,
    drop = FALSE
  ]

  write_tsv(
    axis_level,
    PRIORITIZED_MECHANISM_AXES_FILE
  )

  target_summary <- data.table::as.data.table(
    axis_level
  )[
    ,
    .(
      n_confirmed_communication_axes = .N,
      n_mechanism_tier_1 = sum(
        final_mechanism_tier ==
          "MECHANISM_TIER_1"
      ),
      n_mechanism_tier_2 = sum(
        final_mechanism_tier ==
          "MECHANISM_TIER_2"
      ),
      dominant_branches = collapse_unique(
        dominant_branch
      ),
      pathways = collapse_unique(
        pathway_name
      ),
      ligands = collapse_unique(
        ligand
      ),
      supported_programs = collapse_unique(
        supported_programs_any
      )
    ),
    by = target
  ]

  target_summary <- target_summary[
    order(
      -n_mechanism_tier_1,
      -n_mechanism_tier_2,
      target
    )
  ]

  target_summary <- as.data.frame(
    target_summary,
    stringsAsFactors = FALSE
  )

  write_tsv(
    target_summary,
    TARGET_SUMMARY_FILE
  )

  # --------------------------------------------------------------------------
  # Figures
  # --------------------------------------------------------------------------

  effect_figure_data <- axis_program_evidence[
    is.finite(
      axis_program_evidence$mean_effect_in_expected_direction
    ),
    ,
    drop = FALSE
  ]

  if (nrow(effect_figure_data) > 0L) {
    effect_figure_data$axis_program_label <- paste0(
      effect_figure_data$interaction_name,
      " -> ",
      effect_figure_data$target,
      " | ",
      effect_figure_data$response_program
    )

    effect_figure_data <- utils::head(
      effect_figure_data[
        order(
          match(
            effect_figure_data$program_evidence_tier,
            c(
              "PROGRAM_TIER_1_TWO_SAMPLE",
              "PROGRAM_TIER_2_ONE_SAMPLE",
              "PROGRAM_TIER_3_DIRECTION_ONLY",
              "PROGRAM_TIER_4_NO_SUPPORT"
            )
          ),
          -effect_figure_data$mean_effect_in_expected_direction
        ),
        ,
        drop = FALSE
      ],
      50L
    )

    effect_figure_data$axis_program_label <- factor(
      effect_figure_data$axis_program_label,
      levels = rev(
        effect_figure_data$axis_program_label
      )
    )

    effect_plot <- ggplot2::ggplot(
      effect_figure_data,
      ggplot2::aes(
        x = mean_effect_in_expected_direction,
        y = axis_program_label
      )
    ) +
      ggplot2::geom_col() +
      ggplot2::geom_vline(
        xintercept = 0,
        linetype = 2
      ) +
      ggplot2::labs(
        x = "Mean receptor-program effect in expected direction",
        y = NULL,
        title = "Recipient-program support for confirmed CellChat axes"
      ) +
      ggplot2::theme_bw(
        base_size = 8
      )

    ggplot2::ggsave(
      filename = FIGURE_EFFECT_FILE,
      plot = effect_plot,
      width = 11,
      height = 12,
      units = "in"
    )
  }

  time_figure_data <- program_timecourse[
    program_timecourse$target %in%
      unique(
        confirmed_axes$target
      ),
    ,
    drop = FALSE
  ]

  if (nrow(time_figure_data) > 0L) {
    time_figure_data$timepoint <- factor(
      time_figure_data$timepoint,
      levels = FORMAL_TIMEPOINTS
    )

    time_plot <- ggplot2::ggplot(
      time_figure_data,
      ggplot2::aes(
        x = timepoint,
        y = mean_score,
        group = sample_id
      )
    ) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::facet_grid(
        target ~ response_program,
        scales = "free_y"
      ) +
      ggplot2::labs(
        x = "Time after SCI",
        y = "Within-target standardized program score",
        title = "Recipient response programs across 1, 3 and 7 dpi"
      ) +
      ggplot2::theme_bw(
        base_size = 7
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 45,
          hjust = 1
        )
      )

    ggplot2::ggsave(
      filename = FIGURE_TIME_FILE,
      plot = time_plot,
      width = 16,
      height = 12,
      units = "in",
      limitsize = FALSE
    )
  }

  priority_figure_data <- axis_level[
    axis_level$final_mechanism_tier !=
      "COMMUNICATION_ONLY",
    ,
    drop = FALSE
  ]

  if (nrow(priority_figure_data) > 0L) {
    priority_figure_data$mechanism_chain <- factor(
      priority_figure_data$mechanism_chain,
      levels = rev(
        priority_figure_data$mechanism_chain
      )
    )

    priority_plot <- ggplot2::ggplot(
      priority_figure_data,
      ggplot2::aes(
        x = n_programs_two_sample,
        y = mechanism_chain
      )
    ) +
      ggplot2::geom_col() +
      ggplot2::labs(
        x = "Recipient programs supported in both 7 dpi samples",
        y = NULL,
        title = "Prioritized microglia-to-recipient mechanism axes"
      ) +
      ggplot2::theme_bw(
        base_size = 8
      )

    ggplot2::ggsave(
      filename = FIGURE_PRIORITY_FILE,
      plot = priority_plot,
      width = 11,
      height = 8,
      units = "in"
    )
  }

  # --------------------------------------------------------------------------
  # Final decision
  # --------------------------------------------------------------------------

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
    "Step8A and Step8B formal PASS gates",
    TRUE,
    detail = "Both upstream formal decision gates passed"
  )

  decision_rows <- add_decision_row(
    decision_rows,
    "Confirmed CellChat axes loaded",
    nrow(
      confirmed_axes
    ) ==
      18L,
    detail = paste0(
      nrow(
        confirmed_axes
      ),
      " CONFIRMED_CORE axes"
    ),
    status_if_fail = "REVIEW_REQUIRED"
  )

  decision_rows <- add_decision_row(
    decision_rows,
    "Recipient-program mapping coverage",
    length(
      unmapped_axes
    ) ==
      0L,
    detail = paste0(
      length(
        unique(
          mapped_axis_key
        )
      ),
      "/",
      length(
        unique(
          axis_key
        )
      ),
      " communication axes mapped"
    )
  )

  decision_rows <- add_decision_row(
    decision_rows,
    "Program gene coverage",
    length(
      insufficient_used_programs
    ) ==
      0L,
    detail = paste0(
      length(
        used_programs
      ),
      " used programs passed >= ",
      MIN_PROGRAM_GENES_PRESENT,
      " genes"
    )
  )

  decision_rows <- add_decision_row(
    decision_rows,
    "Recipient timecourse cell scoring",
    nrow(
      recipient_program_cell_scores
    ) >
      0L &&
    all(
      FORMAL_TIMEPOINTS %in%
        recipient_program_cell_scores$timepoint
    ),
    detail = paste0(
      nrow(
        recipient_program_cell_scores
      ),
      " cell-program rows"
    )
  )

  decision_rows <- add_decision_row(
    decision_rows,
    "Two-sample 7 dpi receptor-program audit",
    nrow(
      receptor_audit_long
    ) >
      0L &&
    length(
      unique(
        receptor_audit_long$sample_id
      )
    ) ==
      FORMAL_7DPI_SAMPLES_EXPECTED,
    detail = paste0(
      nrow(
        receptor_audit_long
      ),
      " sample-axis-program rows across ",
      length(
        unique(
          receptor_audit_long$sample_id
        )
      ),
      " samples"
    )
  )

  n_tier1_mechanisms <- sum(
    axis_level$final_mechanism_tier ==
      "MECHANISM_TIER_1",
    na.rm = TRUE
  )

  decision_rows <- add_decision_row(
    decision_rows,
    "Recipient-program-supported mechanisms",
    n_tier1_mechanisms >
      0L,
    detail = paste0(
      n_tier1_mechanisms,
      " MECHANISM_TIER_1 axes; ",
      sum(
        axis_level$final_mechanism_tier ==
          "MECHANISM_TIER_2",
        na.rm = TRUE
      ),
      " MECHANISM_TIER_2 axes"
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
    "Overall Step 8C technical",
    technical_pass,
    detail = if (
      technical_pass
    ) {
      "PASS"
    } else {
      "At least one required technical domain failed"
    },
    status_if_fail = "FAIL"
  )

  biological_status <- if (
    n_tier1_mechanisms >
      0L
  ) {
    "RECIPIENT_PROGRAM_SUPPORT_IDENTIFIED"
  } else {
    "NO_TWO_SAMPLE_RECIPIENT_PROGRAM_SUPPORT"
  }

  decision_rows <- add_decision_row(
    decision_rows,
    "Overall Step 8C biological result",
    n_tier1_mechanisms >
      0L,
    detail = biological_status,
    status_if_fail = "REVIEW_REQUIRED"
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
    biological_status = biological_status,
    script_version = "Step8C_v2_vector_fix",
    generated = timestamp(),
    parameters = list(
      min_program_genes_present =
        MIN_PROGRAM_GENES_PRESENT,
      min_target_cells_per_sample =
        MIN_TARGET_CELLS_PER_SAMPLE,
      min_receptor_positive_cells =
        MIN_RECEPTOR_POSITIVE_CELLS,
      min_receptor_negative_cells =
        MIN_RECEPTOR_NEGATIVE_CELLS,
      min_abs_standardized_effect =
        MIN_ABS_STANDARDIZED_EFFECT,
      min_abs_spearman_rho =
        MIN_ABS_SPEARMAN_RHO
    ),
    confirmed_axes = confirmed_axes,
    axis_program_map = axis_program_map,
    program_gene_coverage =
      program_gene_coverage,
    program_timecourse =
      program_timecourse,
    program_time_trend =
      program_time_trend,
    receptor_audit_long =
      receptor_audit_long,
    axis_program_evidence =
      axis_program_evidence,
    prioritized_mechanism_axes =
      axis_level,
    target_summary =
      target_summary,
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

  top_mechanism_text <- if (
    nrow(
      axis_level
    ) >
      0L
  ) {
    paste(
      capture.output(
        print(
          axis_level[
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
                "recipient_mechanism_status",
                "n_programs_two_sample",
                "supported_programs_two_sample",
                "top_supported_program",
                "top_mean_effect_in_expected_direction",
                "Step7G_candidate_TFs",
                "mechanism_chain"
              ),
              colnames(
                axis_level
              )
            ),
            drop = FALSE
          ],
          row.names = FALSE
        )
      ),
      collapse = "\n"
    )
  } else {
    "No mechanism axes were produced."
  }

  report_lines <- c(
    "======================================================================",
    "GSE162610 | STEP 8C v2 RECIPIENT-PROGRAM INTEGRATION",
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
      "Biological result: ",
      biological_status
    ),
    paste0(
      "Elapsed time: ",
      round(
        elapsed_minutes,
        2
      ),
      " min"
    ),
    "",
    "INPUT:",
    paste0(
      nrow(
        confirmed_axes
      ),
      " formally confirmed 7 dpi CellChat axes."
    ),
    paste0(
      nrow(
        metadata_target
      ),
      " recipient cells across 1, 3 and 7 dpi."
    ),
    paste0(
      length(
        target_groups
      ),
      " recipient cell groups: ",
      paste(
        target_groups,
        collapse = " | "
      ),
      "."
    ),
    "",
    "RESULT SUMMARY:",
    paste0(
      "MECHANISM_TIER_1 axes: ",
      sum(
        axis_level$final_mechanism_tier ==
          "MECHANISM_TIER_1",
        na.rm = TRUE
      )
    ),
    paste0(
      "MECHANISM_TIER_2 axes: ",
      sum(
        axis_level$final_mechanism_tier ==
          "MECHANISM_TIER_2",
        na.rm = TRUE
      )
    ),
    paste0(
      "MECHANISM_TIER_3 axes: ",
      sum(
        axis_level$final_mechanism_tier ==
          "MECHANISM_TIER_3",
        na.rm = TRUE
      )
    ),
    paste0(
      "COMMUNICATION_ONLY axes: ",
      sum(
        axis_level$final_mechanism_tier ==
          "COMMUNICATION_ONLY",
        na.rm = TRUE
      )
    ),
    "",
    "PRIORITIZED MECHANISM AXES:",
    top_mechanism_text,
    "",
    "RECIPIENT TARGET SUMMARY:",
    paste(
      capture.output(
        print(
          target_summary,
          row.names = FALSE
        )
      ),
      collapse = "\n"
    ),
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
      "Recipient-program support means that receptor-positive cells showed ",
      "the expected program direction within the target cell type and that ",
      "the direction was evaluated separately in both 7 dpi samples."
    ),
    paste0(
      "Cell-level Wilcoxon and Spearman P values are descriptive and are not ",
      "biological-replicate inference."
    ),
    paste0(
      "The integrated chain remains associative until spatial proximity, ",
      "external replication or perturbation provides additional support."
    ),
    "",
    "NEXT ACTION:",
    paste0(
      "Use MECHANISM_TIER_1 axes to define the GSE162610 discovery mechanism. ",
      "Then proceed to external validation: first GSE142849 for long-term DAM/",
      "branch persistence, followed by GSE234774 spatial niche validation and ",
      "GSE196928 microglia-depletion testing."
    ),
    "",
    "UPLOAD THESE FILES:",
    "1) reports/Step8C_recipient_program_integration_report.txt",
    "2) tables/06_axis_program_integrated_evidence.tsv",
    "3) tables/07_prioritized_mechanism_axes.tsv",
    "4) tables/08_recipient_target_summary.tsv",
    "5) tables/09_final_decision.csv",
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
    "STEP 8C RECIPIENT-PROGRAM INTEGRATION COMPLETE"
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
    "MECHANISM_TIER_1 axes: ",
    n_tier1_mechanisms
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
      "GSE162610 | STEP 8C v2 ERROR",
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
      "Do not rerun Step1-Step8B.",
      "Correct only this Step8C file and rerun it.",
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
