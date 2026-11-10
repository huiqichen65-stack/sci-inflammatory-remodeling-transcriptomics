# ==============================================================================
# 70_step12G_GSE234774_recipient_neighborhood_CellChatDB_spatial_closure_v10_explicit_priority_program_merge_fix.R
#
# v2 repair
#   The first version calculated reference nUMI from only the candidate
#   ligand/receptor genes. Some valid nuclei, especially Fibroblast_VLMC,
#   have zero counts across that small candidate subset, which incorrectly
#   triggered "Invalid reference nUMI".
#
#   This version normalizes candidate-gene expression using the full-library
#   reference_nUMI saved by Step 12C. If that metadata field is unavailable,
#   it falls back to colSums() of the full Step 12C reference count matrix.
#   Candidate-gene zero expression remains valid and is represented by
#   detection_fraction = 0 and mean_logCP10K = 0.
#
# v3 speed repair
#   The v2 script processed each candidate and each lesion zone through
#   thousands of small R/data.table loops and recalculated a candidate-wise
#   Spearman correlation that is not used in downstream LR priority or animal
#   inference. This version vectorizes candidate x spot closure calculations
#   by recipient cell type and omits that unused diagnostic by default.
#
# v4 matrix/graph repair and acceleration
#   - Fixes dimension loss caused by pmax(0, matrix). That expression can
#     return a vector because pmax() preserves the attributes of its first
#     argument, producing the sweep() warning and "'dims' length cannot be zero".
#   - Uses exact RANN KNN search when RANN is installed, avoiding a dense
#     all-pairs distance matrix plus full sorting for every spot.
#   - Uses matrixStats row operations when available.
#   - Keeps the same CellChatDB candidates, expression thresholds, closure
#     formula, lesion zones, animal aggregation and priority rules.
#
#   No outcome-based candidate pruning is introduced.
#
# v5 repair
#   v4 preserved ligand/receptor matrix dimensions, but its fast zone summary
#   called base rowMeans on Matrix-class objects returned by sparse adjacency
#   multiplication. That function requires a base array and caused
#   "'x' must be an array of at least two dimensions" in every section.
#
#   v5 uses Matrix::rowMeans(), which supports both base matrices and Matrix
#   objects. It also stops immediately after the first failed section instead
#   of wasting time attempting all 36 sections with the same code-level error.
#
# v6 repair
#   The v5 section loop used:
#       spot_zones[spot_zones$section_id == section_id]
#   Inside data.table, the bare right-hand name section_id resolved to the
#   section_id column rather than the loop variable. The expression therefore
#   compared the column with itself and selected all 33,941 spots on every
#   iteration. score_section() then repeatedly used the first section ID, so
#   only 078_A_2months_3 had a matching checkpoint.
#
#   v6 uses a uniquely named loop variable (current_section_id), performs an
#   unambiguous subset, and enforces that each scoring call contains exactly
#   one section with the expected ID. All v5 speed and matrix fixes remain.
#
# v7 aggregation repair
#   The 36 v6 section checkpoints are valid. The failure occurred only after
#   checkpoint combination, in a data.table grouping expression that placed
#   the scalar constant `reference_zone = "remote"` inside by=.(). On the
#   user's data.table version, every grouping item must have one value per row.
#
#   v7:
#     1. creates reference_zone as a real full-length column before grouping;
#     2. uses explicit character vectors for every LR/program grouping stage;
#     3. validates all grouping columns before aggregation;
#     4. retains the v6 section-checkpoint signature so all 36 completed
#        sections are reused without recalculation.
#
# v8 priority/finalization repair
#   The v7 failure occurred after LR animal aggregation, when the target
#   priority table contained 75 rows. A vector-valued condition was coerced as
#   logical(1) inside the previous priority-building path.
#
#   v8:
#     1. rewrites build_target_priority() with explicit scalar checks and
#        vector-indexed tier assignment; no fifelse or vector-valued if();
#     2. rewrites grouping validation with an ordinary for-loop rather than
#        vapply(..., logical(1));
#     3. adds labeled Step 6 substage logs to identify any future failure
#        precisely;
#     4. hardens target-vs-shared and figure inputs before final output;
#     5. deliberately retains the v6 checkpoint signature, so all 36 valid
#        section checkpoints are reused.
#
# v9 parse repair
#   In v8, the data.table pronoun `..support_columns` was accidentally split
#   across two physical lines, producing a source-time syntax error before any
#   analysis ran. v9 replaces that construct with the more version-compatible:
#       DT[i, support_columns, with = FALSE]
#   The v6 checkpoint signature remains unchanged, so all 36 valid section
#   checkpoints are reused.
#
# v10 explicit-column repair
#   v9 reached Step 6D successfully, then build_target_priority() tried to
#   evaluate priority columns as free R objects outside a data.table j scope.
#   This caused:
#       object 'positive_animals_core_minus_remote' not found
#
#   v10:
#     1. builds core/border/perilesional priority columns by explicit merges,
#        avoiding multi-value dcast naming assumptions;
#     2. reads and writes every priority metric through DT[[column_name]];
#     3. builds target-versus-Shared paired animal rows by explicit merges,
#        avoiding another dcast/NSE compatibility point;
#     4. hardens the LR heatmap column selection;
#     5. retains the v6 checkpoint signature, so all 36 section checkpoints
#        are reused without recalculation.
#
# Step 12G
#   Recipient-cell neighborhood and CellChatDB-supported ligand-receptor
#   spatial-closure analysis for the frozen persistent inflammatory-remodeling
#   target program.
#
# Formal design
#   - Step 12F lesion zones are frozen and never redefined here.
#   - CellChatDB.mouse supplies curated ligand-receptor priors.
#   - Target-linked ligands come only from the frozen Step 12E target set.
#   - Shared_prebranch-linked ligands are retained as a parallel control.
#   - Source ligand support is checked in reference Microglia.
#   - Receptor support is checked independently in each recipient cell type.
#   - Spatial ligand signal is weighted by Microglia RCTD abundance.
#   - Neighbor receptor signal is weighted by recipient RCTD abundance.
#   - Animal, not spot, is the inferential unit.
#
# Why the uploaded full CellChat example is not copied literally
#   Visium spots are mixed-cell observations represented here by continuous
#   RCTD weights. Assigning each spot to one arbitrary discrete cell label
#   would discard deconvolution information. This script therefore uses
#   CellChatDB.mouse as the interaction prior and performs continuous,
#   deconvolution-aware spatial closure.
#
# Run:
# source(
#   "<PROJECT_ROOT>/70_step12G_GSE234774_recipient_neighborhood_CellChatDB_spatial_closure_v10_explicit_priority_program_merge_fix.R",
#   encoding = "UTF-8"
# )
# ==============================================================================

options(stringsAsFactors = FALSE, warn = 1, width = 200,
        future.globals.maxSize = 24 * 1024^3)

# ---- frozen controls ----------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
RECIPIENT_TYPES <- c(
  "BAM_macrophage", "Fibroblast_VLMC", "Endothelial", "Astrocyte", "Pericyte"
)
ALLOWED_DB_ANNOTATIONS <- c(
  "Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"
)
INJURED_TIMEPOINTS <- c("7days", "2months")

MIN_SOURCE_DETECTION <- 0.05
MIN_SOURCE_MEAN_LOGCP10K <- 0.05
MIN_RECEPTOR_DETECTION <- 0.05
MIN_RECEPTOR_MEAN_LOGCP10K <- 0.05

K_NEIGHBORS <- 12L
MIN_NEIGHBORS <- 6L
MAX_NEIGHBOR_DISTANCE_NN <- 3
LOGCP10K_SCALE <- 10000
MAX_CANDIDATE_RECEIVER_ROWS <- 8000L
MAX_TOP_LR_PLOT <- 40L

EXPECTED_SPOTS <- 33941L
EXPECTED_SECTIONS <- 36L
EXPECTED_INJURED_SECTIONS <- 24L
EXPECTED_ANIMALS_PER_TIMEPOINT <- 3L

RESUME_COMPLETED <- TRUE
FORCE_RERUN_SECTIONS <- character()

# This diagnostic is not used in any downstream animal-level inference or
# candidate-priority rule. FALSE substantially reduces runtime.
CALCULATE_UNUSED_PAIR_TARGET_SPEARMAN <- FALSE

# Optional accelerators. These are usually already installed with Seurat or
# Bioconductor workflows. Base-R fallbacks remain available.
USE_RANN_FAST_KNN_IF_AVAILABLE <- TRUE
USE_MATRIXSTATS_IF_AVAILABLE <- TRUE

# A code-level error usually affects every section. Stop after the first
# failure so it can be repaired without attempting all 36 sections.
STOP_ON_FIRST_SECTION_FAILURE <- TRUE

# ---- inputs ------------------------------------------------------------------
ROOT_RESULTS <- file.path(PROJECT_ROOT, "04_results")

STEP12F_DIR <- file.path(
  ROOT_RESULTS,
  "step12F_GSE234774_independent_lesion_zones_target_program_validation"
)
STEP12F_REPORT <- file.path(STEP12F_DIR, "10_STEP12F_FINAL_REPORT.txt")
STEP12F_SPOT_RDS <- file.path(
  STEP12F_DIR, "RDS",
  "GSE234774_step12F_spot_independent_lesion_zones_and_scores.rds"
)

STEP12E_DIR <- file.path(
  ROOT_RESULTS,
  "step12E_GSE234774_target_program_spatial_scoring_microglia_adjustment"
)
STEP12E_REPORT <- file.path(STEP12E_DIR, "10_STEP12E_FINAL_REPORT.txt")
STEP12E_HANDOFF_RDS <- file.path(
  STEP12E_DIR, "RDS",
  "GSE234774_step12E_target_program_spatial_score_handoff.rds"
)

STEP12D_DIR <- file.path(
  ROOT_RESULTS,
  "step12D_GSE234774_full_mode_RCTD_per_capture_area"
)
STEP12D_REPORT <- file.path(STEP12D_DIR, "10_STEP12D_FINAL_REPORT.txt")
STEP12D_WEIGHTS_RDS <- file.path(
  STEP12D_DIR, "RDS",
  "GSE234774_step12D_all_spots_normalized_RCTD_weights.rds"
)

STEP12C_DIR <- file.path(
  ROOT_RESULTS,
  "step12C_GSE234774_snRNA_reference_RCTD_input_audit"
)
STEP12C_REPORT <- file.path(STEP12C_DIR, "10_STEP12C_FINAL_REPORT.txt")
STEP12C_REFERENCE_RDS <- file.path(
  STEP12C_DIR, "RDS", "GSE234774_step12C_RCTD_reference.rds"
)
STEP12C_REFERENCE_META_RDS <- file.path(
  STEP12C_DIR, "RDS",
  "GSE234774_step12C_RCTD_reference_cell_metadata.rds"
)

STEP12B_DIR <- file.path(
  ROOT_RESULTS, "step12B_fast_BPCells_spatial_section_QC"
)
STEP12B_REPORT <- file.path(STEP12B_DIR, "10_STEP12B_FAST_FINAL_REPORT.txt")
SPATIAL_BPCELLS_DIR <- file.path(
  STEP12B_DIR, "BPCells", "GSE234774_spatial_2d_full_counts"
)

# ---- outputs -----------------------------------------------------------------
OUTPUT_DIR <- file.path(
  ROOT_RESULTS,
  "step12G_GSE234774_recipient_neighborhood_CellChatDB_spatial_closure"
)
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
RDS_DIR <- file.path(OUTPUT_DIR, "RDS")
CHECKPOINT_DIR <- file.path(OUTPUT_DIR, "section_checkpoints")

for (d in c(OUTPUT_DIR, TABLE_DIR, FIGURE_DIR, RDS_DIR, CHECKPOINT_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

RUN_LOG <- file.path(OUTPUT_DIR, "00_run_log.txt")
SESSION_INFO <- file.path(OUTPUT_DIR, "00_sessionInfo.txt")

OUT_INPUT_AUDIT <- file.path(
  TABLE_DIR, "01_step12G_input_database_parameter_audit.csv"
)
OUT_CANDIDATES <- file.path(
  TABLE_DIR, "02_CellChatDB_target_shared_LR_candidate_reference_support.csv"
)
OUT_REFERENCE_SUPPORT <- file.path(
  TABLE_DIR, "03_snRNA_reference_ligand_receptor_gene_support.csv.gz"
)
OUT_SPOT_NEIGHBORHOOD <- file.path(
  TABLE_DIR, "04_spot_recipient_neighborhood_scores.csv.gz"
)
OUT_SECTION_NEIGHBORHOOD <- file.path(
  TABLE_DIR, "05_section_recipient_neighborhood_zone_contrasts.csv"
)
OUT_ANIMAL_NEIGHBORHOOD <- file.path(
  TABLE_DIR, "06_animal_recipient_neighborhood_zone_contrasts.csv"
)
OUT_NEIGHBORHOOD_INFERENCE <- file.path(
  TABLE_DIR, "07_animal_recipient_neighborhood_inference.csv"
)
OUT_SECTION_LR <- file.path(
  TABLE_DIR, "08_section_LR_spatial_closure_zone_contrasts.csv.gz"
)
OUT_ANIMAL_LR <- file.path(
  TABLE_DIR, "09_animal_LR_spatial_closure_zone_contrasts.csv.gz"
)
OUT_LR_INFERENCE <- file.path(
  TABLE_DIR, "10_animal_LR_spatial_closure_inference.csv.gz"
)
OUT_TARGET_PRIORITY <- file.path(
  TABLE_DIR, "11_target_LR_priority_summary.csv"
)
OUT_PROGRAM_SECTION <- file.path(
  TABLE_DIR, "12_section_target_vs_shared_program_closure_contrasts.csv"
)
OUT_PROGRAM_ANIMAL <- file.path(
  TABLE_DIR, "13_animal_target_vs_shared_program_closure_contrasts.csv"
)
OUT_PROGRAM_DIFFERENTIAL <- file.path(
  TABLE_DIR, "14_animal_target_minus_shared_closure_differential.csv"
)
OUT_PROGRAM_INFERENCE <- file.path(
  TABLE_DIR, "15_target_minus_shared_closure_inference.csv"
)
OUT_CHECKS <- file.path(TABLE_DIR, "19_STEP12G_audit_checks.csv")

OUT_FIG_NEIGHBORHOOD <- file.path(
  FIGURE_DIR, "01_recipient_neighborhood_animal_contrasts.pdf"
)
OUT_FIG_PROGRAM <- file.path(
  FIGURE_DIR, "02_target_vs_Shared_prebranch_program_closure.pdf"
)
OUT_FIG_LR <- file.path(
  FIGURE_DIR, "03_top_target_LR_spatial_closure_heatmap.pdf"
)

OUT_FINAL_REPORT <- file.path(OUTPUT_DIR, "20_STEP12G_FINAL_REPORT.txt")
OUT_SPOT_RDS <- file.path(
  RDS_DIR, "GSE234774_step12G_spot_recipient_neighborhood_scores.rds"
)
OUT_PRIORITY_RDS <- file.path(
  RDS_DIR, "GSE234774_step12G_target_LR_priority_summary.rds"
)
OUT_HANDOFF_RDS <- file.path(
  RDS_DIR,
  "GSE234774_step12G_recipient_neighborhood_LR_spatial_closure_handoff.rds"
)

# ---- packages ----------------------------------------------------------------
required_packages <- c("data.table", "Matrix", "BPCells", "CellChat")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ",
       paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(BPCells)
  library(CellChat)
})

# ---- general helpers ----------------------------------------------------------
run_log <- character()

log_message <- function(...) {
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                 paste0(...))
  cat(line, "\n")
  run_log <<- c(run_log, line)
}

write_run_log <- function() writeLines(run_log, RUN_LOG, useBytes = TRUE)

elapsed_minutes <- function(start_time) {
  round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 3)
}

file_size_safe <- function(path) {
  out <- rep(NA_real_, length(path))
  ok <- file.exists(path)
  out[ok] <- as.numeric(file.info(path[ok])$size)
  out
}

required_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing) > 0L) {
    stop(label, " is missing column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
}

clean_genes <- function(x) {
  x <- trimws(as.character(x))
  unique(x[!is.na(x) & nzchar(x)])
}

resolve_features <- function(requested, available) {
  requested <- clean_genes(requested)
  upper <- toupper(available)
  keep <- !duplicated(upper)
  lookup <- setNames(available[keep], upper[keep])
  out <- unname(lookup[toupper(requested)])
  names(out) <- requested
  out
}

collapse_unique <- function(x, max_values = 40L) {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) return(NA_character_)
  if (length(x) > max_values) {
    x <- c(x[seq_len(max_values)],
           paste0("... +", length(x) - max_values, " more"))
  }
  paste(x, collapse = " | ")
}

robust_z <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (sum(ok) < 2L) {
    out[ok] <- 0
    return(out)
  }
  center <- median(x[ok])
  scale <- mad(x[ok], center = center, constant = 1.4826)
  if (!is.finite(scale) || scale <= 0) scale <- sd(x[ok])
  if (!is.finite(scale) || scale <= 0) {
    out[ok] <- 0
    return(out)
  }
  out[ok] <- pmax(-6, pmin(6, (x[ok] - center) / scale))
  out
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else mean(x)
}
safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else median(x)
}
safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else sd(x)
}
safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4L || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}
exact_sign_p <- function(x) {
  x <- x[is.finite(x) & x != 0]
  n <- length(x)
  if (n == 0L) return(NA_real_)
  k <- sum(x > 0)
  min(1, 2 * pbinom(min(k, n - k), size = n, prob = 0.5))
}
safe_wilcox_signed_p <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L || all(x == 0)) return(NA_real_)
  tryCatch(
    suppressWarnings(wilcox.test(
      x, mu = 0, alternative = "two.sided",
      exact = TRUE, paired = FALSE
    )$p.value),
    error = function(e) NA_real_
  )
}

to_dgCMatrix <- function(x, label) {
  if (inherits(x, "dgCMatrix")) return(x)
  out <- tryCatch(methods::as(x, "dgCMatrix"), error = function(e) e)
  if (inherits(out, "error")) {
    mid <- tryCatch(methods::as(x, "CsparseMatrix"), error = function(e) e)
    if (!inherits(mid, "error")) {
      out <- tryCatch(methods::as(mid, "dgCMatrix"), error = function(e) e)
    }
  }
  if (inherits(out, "error")) {
    stop("Could not convert ", label, " to dgCMatrix: ",
         conditionMessage(out), call. = FALSE)
  }
  out
}

get_reference_counts <- function(x) {
  if (base::isS4(x) && "counts" %in% methods::slotNames(x)) {
    return(methods::slot(x, "counts"))
  }
  if (is.list(x) && !is.null(x$counts)) return(x$counts)
  stop("Step 12C Reference has no accessible counts component.", call. = FALSE)
}

sanitize_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))

# ---- CellChatDB helpers --------------------------------------------------------
load_mouse_db <- function() {
  ns <- asNamespace("CellChat")
  if (exists("CellChatDB.mouse", envir = ns, inherits = FALSE)) {
    return(get("CellChatDB.mouse", envir = ns, inherits = FALSE))
  }
  env <- new.env(parent = emptyenv())
  suppressWarnings(utils::data(
    "CellChatDB.mouse", package = "CellChat", envir = env
  ))
  if (exists("CellChatDB.mouse", envir = env, inherits = FALSE)) {
    return(get("CellChatDB.mouse", envir = env, inherits = FALSE))
  }
  stop("CellChatDB.mouse could not be loaded.", call. = FALSE)
}

entity_subunits <- function(entity, complex_table) {
  entity <- trimws(as.character(entity))
  if (is.na(entity) || !nzchar(entity)) return(character())
  if (!is.null(complex_table) && nrow(complex_table) > 0L &&
      entity %in% rownames(complex_table)) {
    v <- trimws(as.character(unlist(
      complex_table[entity, , drop = FALSE], use.names = FALSE
    )))
    v <- unique(v[!is.na(v) & nzchar(v) & toupper(v) != "NA"])
    if (length(v) > 0L) return(v)
  }
  entity
}

standardize_db <- function(db, target_genes, shared_genes,
                           spatial_features, reference_features) {
  x <- as.data.table(db$interaction, keep.rownames = "db_row")

  find_col <- function(candidates) {
    hit <- intersect(candidates, names(x))
    if (length(hit) == 0L) NA_character_ else hit[[1L]]
  }

  ligand_col <- find_col(c("ligand", "ligand_name"))
  receptor_col <- find_col(c("receptor", "receptor_name"))
  annotation_col <- find_col(c("annotation", "classification"))
  pathway_col <- find_col(c("pathway_name", "pathway"))
  interaction_col <- find_col(c("interaction_name", "interaction_name_2"))

  if (anyNA(c(ligand_col, receptor_col, annotation_col, pathway_col))) {
    stop("Installed CellChatDB schema lacks required columns.", call. = FALSE)
  }

  out <- data.table(
    interaction_name = if (!is.na(interaction_col)) {
      as.character(x[[interaction_col]])
    } else {
      paste0(x[[ligand_col]], "_", x[[receptor_col]])
    },
    pathway_name = as.character(x[[pathway_col]]),
    annotation = as.character(x[[annotation_col]]),
    ligand_entity = as.character(x[[ligand_col]]),
    receptor_entity = as.character(x[[receptor_col]])
  )

  allowed_lower <- setNames(ALLOWED_DB_ANNOTATIONS,
                            tolower(ALLOWED_DB_ANNOTATIONS))
  out[, annotation := unname(allowed_lower[tolower(trimws(annotation))])]
  out <- out[!is.na(annotation)]
  out <- unique(out, by = c(
    "interaction_name", "pathway_name", "annotation",
    "ligand_entity", "receptor_entity"
  ))

  out[, ligand_requested := lapply(
    ligand_entity, entity_subunits, complex_table = db$complex
  )]
  out[, receptor_requested := lapply(
    receptor_entity, entity_subunits, complex_table = db$complex
  )]

  target_upper <- toupper(target_genes)
  shared_upper <- toupper(shared_genes)

  out[, target_hits := vapply(
    ligand_requested,
    function(v) sum(toupper(v) %in% target_upper),
    integer(1)
  )]
  out[, shared_hits := vapply(
    ligand_requested,
    function(v) sum(toupper(v) %in% shared_upper),
    integer(1)
  )]

  out[, program := fifelse(
    target_hits > 0L & shared_hits == 0L, "target_program",
    fifelse(
      shared_hits > 0L & target_hits == 0L,
      "Shared_prebranch", "EXCLUDE"
    )
  )]
  out <- out[program != "EXCLUDE"]

  map_list <- function(list_column, features) {
    lapply(list_column, function(v) unname(resolve_features(v, features)))
  }

  out[, ligand_spatial := map_list(ligand_requested, spatial_features)]
  out[, receptor_spatial := map_list(receptor_requested, spatial_features)]
  out[, ligand_reference := map_list(ligand_requested, reference_features)]
  out[, receptor_reference := map_list(receptor_requested, reference_features)]

  out[, mapping_pass := vapply(seq_len(.N), function(i) {
    all(
      length(ligand_requested[[i]]) > 0L,
      length(receptor_requested[[i]]) > 0L,
      !anyNA(ligand_spatial[[i]]),
      !anyNA(receptor_spatial[[i]]),
      !anyNA(ligand_reference[[i]]),
      !anyNA(receptor_reference[[i]])
    )
  }, logical(1))]

  out[, candidate_id := paste(
    interaction_name, ligand_entity, receptor_entity,
    pathway_name, annotation, sep = "||"
  )]

  setorder(out, program, annotation, pathway_name, interaction_name)
  out
}


# ---- reference support --------------------------------------------------------
build_reference_support <- function(counts, meta, genes, cell_types) {
  required_columns(meta, "rctd_cell_type", "Step 12C reference metadata")

  if (is.null(rownames(meta))) {
    stop("Reference metadata lacks nucleus row names.", call. = FALSE)
  }

  if (!setequal(colnames(counts), rownames(meta))) {
    stop("Reference count/metadata nucleus names do not match.", call. = FALSE)
  }

  meta <- meta[colnames(counts), , drop = FALSE]

  gene_index <- match(genes, rownames(counts))

  if (anyNA(gene_index)) {
    stop(
      "Mapped candidate subunit absent from Step 12C Reference.",
      call. = FALSE
    )
  }

  gene_counts <- counts[
    gene_index,
    ,
    drop = FALSE
  ]

  rownames(gene_counts) <- genes

  # IMPORTANT:
  # Library size must be calculated from the full reference transcriptome,
  # not from the small ligand/receptor candidate subset. A nucleus can
  # legitimately have zero counts for all candidate genes while still having
  # a valid positive total UMI count.
  if ("reference_nUMI" %in% names(meta)) {
    full_reference_nUMI <- suppressWarnings(
      as.numeric(
        meta$reference_nUMI
      )
    )

    nUMI_source <- "Step12C_metadata_reference_nUMI"
  } else {
    full_reference_nUMI <- as.numeric(
      Matrix::colSums(
        counts
      )
    )

    nUMI_source <- "full_reference_count_colSums_fallback"
  }

  names(full_reference_nUMI) <- colnames(counts)

  invalid_full_nUMI <- !is.finite(
    full_reference_nUMI
  ) |
    full_reference_nUMI <= 0

  if (any(invalid_full_nUMI)) {
    stop(
      "Invalid full-library reference nUMI for ",
      sum(invalid_full_nUMI),
      " nuclei. Source: ",
      nUMI_source,
      ".",
      call. = FALSE
    )
  }

  result <- lapply(
    cell_types,
    function(cell_type) {
      cells <- rownames(meta)[
        as.character(
          meta$rctd_cell_type
        ) ==
          cell_type
      ]

      if (length(cells) == 0L) {
        stop(
          "No reference nuclei for ",
          cell_type,
          ".",
          call. = FALSE
        )
      }

      current <- gene_counts[
        ,
        cells,
        drop = FALSE
      ]

      nUMI <- full_reference_nUMI[
        cells
      ]

      if (
        length(nUMI) != length(cells) ||
          anyNA(nUMI) ||
          any(
            !is.finite(nUMI) |
              nUMI <= 0
          )
      ) {
        stop(
          "Invalid aligned full-library reference nUMI in ",
          cell_type,
          ". Source: ",
          nUMI_source,
          ".",
          call. = FALSE
        )
      }

      normalized <- current %*%
        Matrix::Diagonal(
          x = LOGCP10K_SCALE /
            nUMI
        )

      # Sparse zeros remain zero. Only stored nonzero values are transformed.
      normalized@x <- log1p(
        normalized@x
      )

      data.table(
        gene = genes,
        rctd_cell_type = cell_type,
        n_reference_nuclei = length(cells),
        detection_fraction = as.numeric(
          Matrix::rowMeans(
            current > 0
          )
        ),
        mean_logCP10K = as.numeric(
          Matrix::rowMeans(
            normalized
          )
        ),
        reference_nUMI_source = nUMI_source
      )
    }
  )

  rbindlist(
    result,
    fill = TRUE
  )
}

entity_support <- function(subunits, cell_type, support_table) {
  rows <- support_table[
    rctd_cell_type == cell_type & gene %in% subunits
  ]
  if (nrow(rows) != length(unique(subunits))) {
    return(c(detection = NA_real_, mean_logCP10K = NA_real_))
  }
  c(
    detection = min(rows$detection_fraction),
    mean_logCP10K = min(rows$mean_logCP10K)
  )
}

expand_receivers <- function(interactions, support_table) {
  rows <- vector("list", nrow(interactions) * length(RECIPIENT_TYPES))
  k <- 0L

  for (i in seq_len(nrow(interactions))) {
    source <- entity_support(
      interactions$ligand_reference[[i]],
      "Microglia", support_table
    )

    for (recipient in RECIPIENT_TYPES) {
      receptor <- entity_support(
        interactions$receptor_reference[[i]],
        recipient, support_table
      )
      k <- k + 1L

      source_pass <- is.finite(source[["detection"]]) &&
        is.finite(source[["mean_logCP10K"]]) &&
        source[["detection"]] >= MIN_SOURCE_DETECTION &&
        source[["mean_logCP10K"]] >= MIN_SOURCE_MEAN_LOGCP10K

      receptor_pass <- is.finite(receptor[["detection"]]) &&
        is.finite(receptor[["mean_logCP10K"]]) &&
        receptor[["detection"]] >= MIN_RECEPTOR_DETECTION &&
        receptor[["mean_logCP10K"]] >= MIN_RECEPTOR_MEAN_LOGCP10K

      rows[[k]] <- data.table(
        candidate_id = interactions$candidate_id[i],
        candidate_receiver_id = paste(
          interactions$candidate_id[i], recipient, sep = "||"
        ),
        program = interactions$program[i],
        annotation = interactions$annotation[i],
        pathway_name = interactions$pathway_name[i],
        interaction_name = interactions$interaction_name[i],
        ligand_entity = interactions$ligand_entity[i],
        receptor_entity = interactions$receptor_entity[i],
        ligand_requested = list(interactions$ligand_requested[[i]]),
        receptor_requested = list(interactions$receptor_requested[[i]]),
        ligand_spatial = list(interactions$ligand_spatial[[i]]),
        receptor_spatial = list(interactions$receptor_spatial[[i]]),
        ligand_reference = list(interactions$ligand_reference[[i]]),
        receptor_reference = list(interactions$receptor_reference[[i]]),
        mapping_pass = interactions$mapping_pass[i],
        recipient_type = recipient,
        source_Microglia_detection_fraction = source[["detection"]],
        source_Microglia_mean_logCP10K = source[["mean_logCP10K"]],
        recipient_receptor_detection_fraction = receptor[["detection"]],
        recipient_receptor_mean_logCP10K = receptor[["mean_logCP10K"]],
        source_support_pass = source_pass,
        receptor_support_pass = receptor_pass,
        candidate_pass = interactions$mapping_pass[i] &&
          source_pass && receptor_pass
      )
    }
  }

  out <- rbindlist(rows)
  out[, candidate_status := fifelse(
    !mapping_pass, "FAIL_MAPPING",
    fifelse(
      !source_support_pass, "FAIL_MICROGLIA_LIGAND_SUPPORT",
      fifelse(!receptor_support_pass, "FAIL_RECIPIENT_RECEPTOR_SUPPORT", "PASS")
    )
  )]
  setorder(out, -candidate_pass, program, recipient_type,
           annotation, pathway_name, interaction_name)
  out
}

candidate_to_csv <- function(x) {
  out <- copy(x)
  list_columns <- c(
    "ligand_requested", "receptor_requested",
    "ligand_spatial", "receptor_spatial",
    "ligand_reference", "receptor_reference"
  )
  for (column in list_columns) {
    out[, (column) := vapply(
      get(column), paste, collapse = " | ", FUN.VALUE = character(1)
    )]
  }
  out
}

build_entity_maps <- function(candidates) {
  ligands <- unique(candidates$ligand_entity)
  receptors <- unique(candidates$receptor_entity)

  ligand_map <- setNames(lapply(ligands, function(entity) {
    candidates[ligand_entity == entity, ligand_spatial[[1L]]]
  }), ligands)

  receptor_map <- setNames(lapply(receptors, function(entity) {
    candidates[receptor_entity == entity, receptor_spatial[[1L]]]
  }), receptors)

  list(ligand = ligand_map, receptor = receptor_map)
}

# ---- spatial graph and entity expression --------------------------------------
build_adjacency <- function(x, y) {
  coordinates <- cbind(
    as.numeric(x),
    as.numeric(y)
  )

  n_spots <- nrow(
    coordinates
  )

  if (
    n_spots < 2L ||
      any(
        !is.finite(
          coordinates
        )
      )
  ) {
    stop(
      "Invalid section coordinates.",
      call. = FALSE
    )
  }

  use_rann <- isTRUE(
    USE_RANN_FAST_KNN_IF_AVAILABLE
  ) &&
    requireNamespace(
      "RANN",
      quietly = TRUE
    )

  selected_indices <- vector(
    "list",
    n_spots
  )

  selected_distances <- vector(
    "list",
    n_spots
  )

  neighbor_n <- integer(
    n_spots
  )

  if (use_rann) {
    k_query <- min(
      n_spots,
      K_NEIGHBORS +
        2L
    )

    knn <- RANN::nn2(
      data = coordinates,
      query = coordinates,
      k = k_query,
      eps = 0
    )

    nearest_positive <- numeric(
      n_spots
    )

    for (i in seq_len(
      n_spots
    )) {
      candidate_indices <- as.integer(
        knn$nn.idx[
          i,
          ,
          drop = TRUE
        ]
      )

      candidate_distances <- as.numeric(
        knn$nn.dists[
          i,
          ,
          drop = TRUE
        ]
      )

      valid <- candidate_indices > 0L &
        candidate_indices != i &
        is.finite(
          candidate_distances
        )

      candidate_indices <- candidate_indices[
        valid
      ]

      candidate_distances <- candidate_distances[
        valid
      ]

      positive <- candidate_distances[
        candidate_distances > 0
      ]

      nearest_positive[
        i
      ] <- if (
        length(
          positive
        ) > 0L
      ) {
        min(
          positive
        )
      } else {
        NA_real_
      }

      selected_indices[[
        i
      ]] <- candidate_indices

      selected_distances[[
        i
      ]] <- candidate_distances
    }

    nn_spacing <- stats::median(
      nearest_positive[
        is.finite(
          nearest_positive
        ) &
          nearest_positive > 0
      ]
    )

    graph_engine <- "RANN_exact_KNN"
  } else {
    distances <- as.matrix(
      stats::dist(
        coordinates,
        diag = TRUE,
        upper = TRUE
      )
    )

    diag(
      distances
    ) <- Inf

    nearest <- apply(
      distances,
      1L,
      min
    )

    nn_spacing <- stats::median(
      nearest[
        is.finite(
          nearest
        ) &
          nearest > 0
      ]
    )

    for (i in seq_len(
      n_spots
    )) {
      ordered <- order(
        distances[
          i,
        ],
        na.last = NA
      )

      ordered <- ordered[
        ordered != i
      ]

      selected_indices[[
        i
      ]] <- ordered

      selected_distances[[
        i
      ]] <- distances[
        i,
        ordered
      ]
    }

    graph_engine <- "base_exact_all_pairs_fallback"
  }

  if (
    !is.finite(
      nn_spacing
    ) ||
      nn_spacing <= 0
  ) {
    stop(
      "Could not estimate positive local NN spacing.",
      call. = FALSE
    )
  }

  for (i in seq_len(
    n_spots
  )) {
    current_indices <- selected_indices[[
      i
    ]]

    current_distances <- selected_distances[[
      i
    ]]

    valid <- current_indices > 0L &
      current_indices != i &
      is.finite(
        current_distances
      )

    current_indices <- current_indices[
      valid
    ]

    current_distances <- current_distances[
      valid
    ]

    within <- which(
      current_distances <=
        MAX_NEIGHBOR_DISTANCE_NN *
        nn_spacing
    )

    selected_positions <- head(
      within,
      K_NEIGHBORS
    )

    if (length(
      selected_positions
    ) < MIN_NEIGHBORS) {
      selected_positions <- seq_len(
        min(
          K_NEIGHBORS,
          length(
            current_indices
          )
        )
      )
    }

    selected_indices[[
      i
    ]] <- current_indices[
      selected_positions
    ]

    selected_distances[[
      i
    ]] <- current_distances[
      selected_positions
    ]

    neighbor_n[
      i
    ] <- length(
      selected_positions
    )
  }

  if (any(
    neighbor_n <= 0L
  )) {
    stop(
      "At least one spot has no spatial neighbors.",
      call. = FALSE
    )
  }

  ii <- rep(
    seq_len(
      n_spots
    ),
    lengths(
      selected_indices
    )
  )

  jj <- unlist(
    selected_indices,
    use.names = FALSE
  )

  distance_values <- unlist(
    selected_distances,
    use.names = FALSE
  )

  weight_values <- exp(
    -distance_values /
      nn_spacing
  )

  row_totals <- ave(
    weight_values,
    ii,
    FUN = sum
  )

  if (any(
    !is.finite(
      row_totals
    ) |
      row_totals <= 0
  )) {
    stop(
      "Invalid spatial-neighbor weight totals.",
      call. = FALSE
    )
  }

  weight_values <- weight_values /
    row_totals

  adjacency <- Matrix::sparseMatrix(
    i = ii,
    j = jj,
    x = weight_values,
    dims = c(
      n_spots,
      n_spots
    ),
    giveCsparse = TRUE
  )

  maximum_row_sum_deviation <- max(
    abs(
      as.numeric(
        Matrix::rowSums(
          adjacency
        )
      ) -
        1
    )
  )

  if (
    !is.finite(
      maximum_row_sum_deviation
    ) ||
      maximum_row_sum_deviation >
        1e-8
  ) {
    stop(
      "Spatial adjacency rows do not sum to 1.",
      call. = FALSE
    )
  }

  list(
    adjacency =
      adjacency,
    nn_spacing =
      nn_spacing,
    neighbor_n =
      neighbor_n,
    graph_engine =
      graph_engine
  )
}

entity_expression <- function(log_expression, entity_map) {
  entities <- names(entity_map)
  out <- matrix(
    NA_real_, nrow = length(entities), ncol = ncol(log_expression),
    dimnames = list(entities, colnames(log_expression))
  )
  for (i in seq_along(entities)) {
    subunits <- entity_map[[entities[i]]]
    index <- match(subunits, rownames(log_expression))
    if (anyNA(index)) {
      stop("Spatial subunit missing for entity ", entities[i], ".",
           call. = FALSE)
    }
    current <- log_expression[index, , drop = FALSE]
    out[i, ] <- if (length(index) == 1L) {
      as.numeric(current[1L, ])
    } else {
      apply(as.matrix(current), 2L, min)
    }
  }
  out
}

row_robust_z <- function(x) {
  x <- as.matrix(
    x
  )

  if (
    length(
      dim(
        x
      )
    ) != 2L
  ) {
    stop(
      "row_robust_z requires a two-dimensional matrix.",
      call. = FALSE
    )
  }

  use_matrixstats <- isTRUE(
    USE_MATRIXSTATS_IF_AVAILABLE
  ) &&
    requireNamespace(
      "matrixStats",
      quietly = TRUE
    ) &&
    all(
      is.finite(
        x
      )
    )

  if (use_matrixstats) {
    centers <- matrixStats::rowMedians(
      x,
      na.rm = FALSE
    )

    scales <- matrixStats::rowMads(
      x,
      center = centers,
      constant = 1.4826,
      na.rm = FALSE
    )

    invalid_scale <- !is.finite(
      scales
    ) |
      scales <= 0

    if (any(
      invalid_scale
    )) {
      scales[
        invalid_scale
      ] <- matrixStats::rowSds(
        x[
          invalid_scale,
          ,
          drop = FALSE
        ],
        na.rm = FALSE
      )
    }

    output <- sweep(
      x,
      1L,
      centers,
      "-"
    )

    valid_scale <- is.finite(
      scales
    ) &
      scales > 0

    if (any(
      valid_scale
    )) {
      output[
        valid_scale,
      ] <- sweep(
        output[
          valid_scale,
          ,
          drop = FALSE
        ],
        1L,
        scales[
          valid_scale
        ],
        "/"
      )
    }

    if (any(
      !valid_scale
    )) {
      output[
        !valid_scale,
      ] <- 0
    }

    output[
      output > 6
    ] <- 6

    output[
      output < -6
    ] <- -6

    dimnames(
      output
    ) <- dimnames(
      x
    )

    return(
      output
    )
  }

  output <- matrix(
    NA_real_,
    nrow = nrow(
      x
    ),
    ncol = ncol(
      x
    ),
    dimnames = dimnames(
      x
    )
  )

  for (i in seq_len(
    nrow(
      x
    )
  )) {
    output[
      i,
    ] <- robust_z(
      x[
        i,
      ]
    )
  }

  output
}

fast_row_medians <- function(x) {
  x <- as.matrix(
    x
  )

  if (
    isTRUE(
      USE_MATRIXSTATS_IF_AVAILABLE
    ) &&
      requireNamespace(
        "matrixStats",
        quietly = TRUE
      )
  ) {
    return(
      matrixStats::rowMedians(
        x,
        na.rm = FALSE
      )
    )
  }

  apply(
    x,
    1L,
    stats::median
  )
}

# ---- section checkpoints ------------------------------------------------------
section_files <- function(section_id) {
  directory <- file.path(CHECKPOINT_DIR, sanitize_name(section_id))
  list(
    directory = directory,
    result = file.path(directory, "section_result.rds"),
    flag = file.path(directory, "SUCCESS.flag"),
    failure = file.path(directory, "FAILURE.txt")
  )
}

section_complete <- function(section_id, signature) {
  f <- section_files(section_id)
  if (!file.exists(f$result) || !file.exists(f$flag) ||
      file_size_safe(f$result) <= 0) return(FALSE)
  lines <- readLines(f$flag, warn = FALSE)
  any(lines == paste0("candidate_signature=", signature)) &&
    any(lines == "status=PASS")
}

write_section_failure <- function(section_id, reason) {
  f <- section_files(section_id)
  dir.create(f$directory, recursive = TRUE, showWarnings = FALSE)
  writeLines(c(
    paste0("section_id=", section_id),
    "status=FAIL",
    paste0("error=", reason),
    "Completed sections remain reusable."
  ), f$failure)
  if (file.exists(f$flag)) unlink(f$flag, force = TRUE)
}


# ---- score one section --------------------------------------------------------
score_section <- function(
    section_data, weights, spatial_counts, candidates,
    ligand_map, receptor_map, signature,
    expected_section_id = NULL) {

  if (
    nrow(
      section_data
    ) == 0L
  ) {
    stop(
      "score_section received zero spots.",
      call. = FALSE
    )
  }

  observed_section_ids <- unique(
    as.character(
      section_data$section_id
    )
  )

  if (
    length(
      observed_section_ids
    ) != 1L
  ) {
    stop(
      "score_section must receive exactly one section; observed ",
      length(
        observed_section_ids
      ),
      ": ",
      paste(
        observed_section_ids,
        collapse = " | "
      ),
      call. = FALSE
    )
  }

  section_id <- observed_section_ids[[1L]]

  if (
    !is.null(
      expected_section_id
    ) &&
      !identical(
        section_id,
        as.character(
          expected_section_id
        )
      )
  ) {
    stop(
      "Section-subsetting mismatch: expected ",
      expected_section_id,
      " but score_section received ",
      section_id,
      ".",
      call. = FALSE
    )
  }

  timepoint_values <- unique(
    as.character(
      section_data$timepoint
    )
  )

  if (
    length(
      timepoint_values
    ) != 1L
  ) {
    stop(
      "Section ",
      section_id,
      " contains multiple timepoints: ",
      paste(
        timepoint_values,
        collapse = " | "
      ),
      call. = FALSE
    )
  }

  timepoint <- timepoint_values[[1L]]
  f <- section_files(section_id)
  dir.create(f$directory, recursive = TRUE, showWarnings = FALSE)

  if (section_id %in% FORCE_RERUN_SECTIONS) {
    unlink(f$directory, recursive = TRUE, force = TRUE)
    dir.create(f$directory, recursive = TRUE, showWarnings = FALSE)
  }

  barcodes <- as.character(section_data$barcode)
  spot_index <- match(barcodes, colnames(spatial_counts))
  if (anyNA(spot_index)) {
    stop("Spatial matrix is missing barcodes for ", section_id, ".",
         call. = FALSE)
  }
  if (!all(barcodes %in% rownames(weights))) {
    stop("RCTD weights are missing barcodes for ", section_id, ".",
         call. = FALSE)
  }
  section_weights <- weights[barcodes, , drop = FALSE]

  graph <- build_adjacency(
    section_data$coordinate_x_fullres,
    section_data$coordinate_y_fullres
  )
  adjacency <- graph$adjacency

  target_positive <- pmax(
    0, as.numeric(section_data$target_program_adjusted_z_within_section)
  )
  shared_positive <- pmax(
    0, as.numeric(section_data$shared_prebranch_adjusted_z_within_section)
  )

  neighborhood <- data.table(
    barcode = barcodes,
    capture_area_id = as.character(section_data$capture_area_id),
    section_id = section_id,
    timepoint = timepoint,
    biological_replicate_id =
      as.character(section_data$biological_replicate_id),
    lesion_zone = as.character(section_data$lesion_zone),
    target_program_adjusted_z_within_section =
      as.numeric(section_data$target_program_adjusted_z_within_section),
    shared_prebranch_adjusted_z_within_section =
      as.numeric(section_data$shared_prebranch_adjusted_z_within_section),
    local_nn_spacing_fullres = graph$nn_spacing,
    neighbor_count = graph$neighbor_n
  )

  for (recipient in RECIPIENT_TYPES) {
    neighbor_weight <- as.numeric(
      adjacency %*% as.numeric(section_weights[, recipient])
    )
    neighborhood[, (paste0("neighbor_", recipient, "_weight")) :=
                   neighbor_weight]
    neighborhood[, (paste0("target_neighbor_", recipient,
                            "_colocalization")) :=
                   target_positive * neighbor_weight]
    neighborhood[, (paste0("shared_neighbor_", recipient,
                            "_colocalization")) :=
                   shared_positive * neighbor_weight]
  }

  # Uninjured sections contribute recipient-neighborhood reference values only.
  if (!timepoint %in% INJURED_TIMEPOINTS) {
    result <- list(
      neighborhood = neighborhood,
      pair_contrasts = data.table(),
      program_contrasts = data.table(),
      audit = data.table(
        section_id = section_id, timepoint = timepoint,
        n_spots = nrow(section_data),
        n_candidates = nrow(candidates),
        n_pair_contrasts = 0L,
        nn_spacing = graph$nn_spacing,
        graph_engine = graph$graph_engine,
        min_neighbors = min(graph$neighbor_n),
        median_neighbors = median(graph$neighbor_n),
        max_neighbors = max(graph$neighbor_n),
        status = "PASS_UNINJURED"
      )
    )
    saveRDS(result, f$result, compress = FALSE)
    writeLines(c(
      paste0("candidate_signature=", signature),
      "status=PASS", paste0("section_id=", section_id),
      "analysis=UNINJURED_NEIGHBORHOOD_ONLY"
    ), f$flag)
    return(result)
  }

  required_genes <- unique(c(
    unlist(ligand_map, use.names = FALSE),
    unlist(receptor_map, use.names = FALSE)
  ))
  gene_index <- match(required_genes, rownames(spatial_counts))
  if (anyNA(gene_index)) {
    stop("A candidate LR subunit is absent from the spatial matrix.",
         call. = FALSE)
  }

  counts <- to_dgCMatrix(
    spatial_counts[gene_index, spot_index, drop = FALSE],
    paste0(section_id, " LR counts")
  )
  rownames(counts) <- required_genes
  colnames(counts) <- barcodes

  totals <- as.numeric(section_data$total_spatial_counts)
  if (any(!is.finite(totals) | totals <= 0)) {
    stop("Invalid total_spatial_counts in ", section_id, ".", call. = FALSE)
  }

  log_expression <- counts %*% Matrix::Diagonal(
    x = LOGCP10K_SCALE / totals
  )
  log_expression@x <- log1p(log_expression@x)

  ligand_expression <- entity_expression(log_expression, ligand_map)
  receptor_expression <- entity_expression(log_expression, receptor_map)

  ligand_positive_z <- row_robust_z(
    ligand_expression
  )

  receptor_positive_z <- row_robust_z(
    receptor_expression
  )

  # pmax(0, matrix) can drop dim attributes because the first argument is a
  # scalar. Truncate in-place to preserve candidate x spot dimensions.
  ligand_positive_z[
    !is.finite(
      ligand_positive_z
    ) |
      ligand_positive_z < 0
  ] <- 0

  receptor_positive_z[
    !is.finite(
      receptor_positive_z
    ) |
      receptor_positive_z < 0
  ] <- 0

  expected_spots <- nrow(
    section_data
  )

  matrix_dimensions_ok <-
    length(
      dim(
        ligand_positive_z
      )
    ) ==
      2L &&
    length(
      dim(
        receptor_positive_z
      )
    ) ==
      2L &&
    ncol(
      ligand_positive_z
    ) ==
      expected_spots &&
    ncol(
      receptor_positive_z
    ) ==
      expected_spots

  if (!matrix_dimensions_ok) {
    stop(
      "Ligand/receptor robust-Z matrices lost their expected dimensions in ",
      section_id,
      ".",
      call. = FALSE
    )
  }

  microglia_weight <- as.numeric(
    section_weights[
      ,
      "Microglia"
    ]
  )

  if (
    length(
      microglia_weight
    ) !=
      ncol(
        ligand_positive_z
      )
  ) {
    stop(
      "Microglia weight length does not match ligand matrix columns in ",
      section_id,
      ".",
      call. = FALSE
    )
  }

  source_signal <- sweep(
    ligand_positive_z,
    2L,
    microglia_weight,
    "*"
  )

  neighbor_receptor <- setNames(
    vector(
      "list",
      length(
        RECIPIENT_TYPES
      )
    ),
    RECIPIENT_TYPES
  )

  for (recipient in RECIPIENT_TYPES) {
    recipient_weight <- as.numeric(
      section_weights[
        ,
        recipient
      ]
    )

    if (
      length(
        recipient_weight
      ) !=
        ncol(
          receptor_positive_z
        )
    ) {
      stop(
        "Recipient weight length does not match receptor matrix columns for ",
        recipient,
        " in ",
        section_id,
        ".",
        call. = FALSE
      )
    }

    receiver_signal <- sweep(
      receptor_positive_z,
      2L,
      recipient_weight,
      "*"
    )

    neighbor_receptor[[
      recipient
    ]] <- t(
      adjacency %*%
        t(
          receiver_signal
        )
    )
  }

  zones <- c("core", "border", "perilesional", "remote")

  # Vectorized candidate processing: one candidate x spot closure matrix per
  # recipient cell type, followed by matrix row summaries for each lesion zone.
  pair_zone_rows <- list()
  out_i <- 0L

  target_vector <- as.numeric(
    section_data$target_program_adjusted_z_within_section
  )

  for (recipient in RECIPIENT_TYPES) {
    candidate_indices <- which(
      candidates$recipient_type ==
        recipient
    )

    if (length(candidate_indices) == 0L) {
      next
    }

    current_candidates <- candidates[
      candidate_indices
    ]

    ligand_matrix <- source_signal[
      current_candidates$ligand_entity,
      ,
      drop = FALSE
    ]

    receptor_matrix <- neighbor_receptor[[
      recipient
    ]][
      current_candidates$receptor_entity,
      ,
      drop = FALSE
    ]

    closure_matrix <- ligand_matrix *
      receptor_matrix

    if (
      isTRUE(
        CALCULATE_UNUSED_PAIR_TARGET_SPEARMAN
      )
    ) {
      target_correlations <- apply(
        closure_matrix,
        1L,
        function(closure) {
          safe_spearman(
            closure,
            target_vector
          )
        }
      )
    } else {
      target_correlations <- rep(
        NA_real_,
        nrow(
          closure_matrix
        )
      )
    }

    candidate_metadata <- current_candidates[
      ,
      .(
        candidate_id,
        candidate_receiver_id,
        program,
        annotation,
        pathway_name,
        interaction_name,
        ligand_entity,
        receptor_entity,
        recipient_type
      )
    ]

    for (zone in zones) {
      zone_indices <- which(
        as.character(
          section_data$lesion_zone
        ) ==
          zone
      )

      if (length(zone_indices) == 0L) {
        stop(
          "No spots found in lesion zone ",
          zone,
          " for ",
          section_id,
          ".",
          call. = FALSE
        )
      }

      zone_closure <- closure_matrix[
        ,
        zone_indices,
        drop = FALSE
      ]

      zone_ligand <- ligand_matrix[
        ,
        zone_indices,
        drop = FALSE
      ]

      zone_receptor <- receptor_matrix[
        ,
        zone_indices,
        drop = FALSE
      ]

      out_i <- out_i +
        1L

      pair_zone_rows[[
        out_i
      ]] <- cbind(
        data.table(
          capture_area_id =
            as.character(
              section_data$capture_area_id[
                1L
              ]
            ),
          section_id =
            section_id,
          timepoint =
            timepoint,
          biological_replicate_id =
            as.character(
              section_data$biological_replicate_id[
                1L
              ]
            )
        ),
        candidate_metadata,
        data.table(
          lesion_zone =
            zone,
          n_spots =
            length(
              zone_indices
            ),
          closure_mean =
            as.numeric(
              Matrix::rowMeans(
                zone_closure
              )
            ),
          closure_median =
            fast_row_medians(
              zone_closure
            ),
          source_signal_mean =
            as.numeric(
              Matrix::rowMeans(
                zone_ligand
              )
            ),
          neighbor_receptor_signal_mean =
            as.numeric(
              Matrix::rowMeans(
                zone_receptor
              )
            ),
          closure_nonzero_fraction =
            as.numeric(
              Matrix::rowMeans(
                zone_closure >
                  0
              )
            ),
          closure_vs_target_spearman =
            target_correlations
        )
      )
    }

    rm(
      ligand_matrix,
      receptor_matrix,
      closure_matrix
    )

    invisible(
      gc()
    )
  }

  pair_zone <- rbindlist(
    pair_zone_rows,
    fill = TRUE
  )

  remote <- pair_zone[
    lesion_zone == "remote",
    .(
      candidate_receiver_id,
      remote_n_spots = n_spots,
      remote_closure_mean = closure_mean
    )
  ]

  pair_contrasts <- merge(
    pair_zone[lesion_zone != "remote"],
    remote,
    by = "candidate_receiver_id",
    all.x = TRUE,
    sort = FALSE
  )
  pair_contrasts[, contrast := paste0(lesion_zone, "_minus_remote")]
  pair_contrasts[, closure_zone_minus_remote :=
                   closure_mean - remote_closure_mean]

  program_by_annotation <- pair_zone[, .(
    n_candidate_pairs = uniqueN(candidate_receiver_id),
    mean_pair_closure = safe_mean(closure_mean),
    median_pair_closure = safe_median(closure_mean)
  ), by = .(
    capture_area_id, section_id, timepoint, biological_replicate_id,
    program, recipient_type, annotation_group = annotation, lesion_zone
  )]

  program_overall <- pair_zone[, .(
    n_candidate_pairs = uniqueN(candidate_receiver_id),
    mean_pair_closure = safe_mean(closure_mean),
    median_pair_closure = safe_median(closure_mean)
  ), by = .(
    capture_area_id, section_id, timepoint, biological_replicate_id,
    program, recipient_type, lesion_zone
  )]
  program_overall[, annotation_group := "ALL_PROTEIN"]

  program_zone <- rbind(program_by_annotation, program_overall, fill = TRUE)
  program_remote <- program_zone[
    lesion_zone == "remote",
    .(
      section_id, program, recipient_type, annotation_group,
      remote_n_candidate_pairs = n_candidate_pairs,
      remote_mean_pair_closure = mean_pair_closure
    )
  ]
  program_contrasts <- merge(
    program_zone[lesion_zone != "remote"],
    program_remote,
    by = c("section_id", "program", "recipient_type", "annotation_group"),
    all.x = TRUE,
    sort = FALSE
  )
  program_contrasts[, contrast := paste0(lesion_zone, "_minus_remote")]
  program_contrasts[, program_closure_zone_minus_remote :=
                      mean_pair_closure - remote_mean_pair_closure]

  result <- list(
    neighborhood = neighborhood,
    pair_contrasts = pair_contrasts,
    program_contrasts = program_contrasts,
    audit = data.table(
      section_id = section_id, timepoint = timepoint,
      n_spots = nrow(section_data),
      n_candidates = nrow(candidates),
      n_pair_contrasts = nrow(pair_contrasts),
      n_program_contrasts = nrow(program_contrasts),
      nn_spacing = graph$nn_spacing,
      graph_engine = graph$graph_engine,
      min_neighbors = min(graph$neighbor_n),
      median_neighbors = median(graph$neighbor_n),
      max_neighbors = max(graph$neighbor_n),
      closure_min = min(pair_zone$closure_mean),
      closure_max = max(pair_zone$closure_mean),
      status = "PASS_INJURED"
    )
  )

  saveRDS(result, f$result, compress = FALSE)
  writeLines(c(
    paste0("candidate_signature=", signature),
    "status=PASS", paste0("section_id=", section_id),
    "analysis=INJURED_NEIGHBORHOOD_AND_LR_CLOSURE"
  ), f$flag)

  result
}

# ---- contrast helpers ---------------------------------------------------------
build_neighborhood_section_contrasts <- function(spot_table) {
  injured <- spot_table[timepoint %in% INJURED_TIMEPOINTS]
  rows <- list()
  k <- 0L

  for (current_section in unique(as.character(injured$section_id))) {
    current <- injured[section_id == current_section]
    remote <- current[lesion_zone == "remote"]

    for (recipient in RECIPIENT_TYPES) {
      metrics <- c(
        neighbor_weight = paste0("neighbor_", recipient, "_weight"),
        target_neighbor_colocalization = paste0(
          "target_neighbor_", recipient, "_colocalization"
        ),
        shared_neighbor_colocalization = paste0(
          "shared_neighbor_", recipient, "_colocalization"
        )
      )

      for (metric_name in names(metrics)) {
        metric_column <- unname(metrics[[metric_name]])
        remote_mean <- safe_mean(remote[[metric_column]])

        for (zone in c("core", "border", "perilesional")) {
          zone_data <- current[lesion_zone == zone]
          zone_mean <- safe_mean(zone_data[[metric_column]])
          k <- k + 1L

          rows[[k]] <- data.table(
            capture_area_id = as.character(current$capture_area_id[1L]),
            section_id = current_section,
            timepoint = as.character(current$timepoint[1L]),
            biological_replicate_id =
              as.character(current$biological_replicate_id[1L]),
            recipient_type = recipient,
            metric = metric_name,
            contrast = paste0(zone, "_minus_remote"),
            zone = zone,
            reference_zone = "remote",
            n_zone_spots = nrow(zone_data),
            n_remote_spots = nrow(remote),
            zone_mean = zone_mean,
            remote_mean = remote_mean,
            zone_minus_remote = zone_mean - remote_mean
          )
        }
      }
    }
  }
  rbindlist(rows)
}

validate_grouping_columns <- function(
    table,
    columns,
    label) {

  if (!data.table::is.data.table(
    table
  )) {
    table <- data.table::as.data.table(
      table
    )
  }

  columns <- as.character(
    columns
  )

  required_columns(
    table,
    columns,
    label
  )

  invalid_details <- character()

  for (column_name in columns) {
    value <- table[[
      column_name
    ]]

    column_invalid <-
      length(
        value
      ) !=
        nrow(
          table
        ) ||
      !is.null(
        dim(
          value
        )
      ) ||
      is.list(
        value
      )

    if (isTRUE(
      column_invalid
    )) {
      invalid_details <- c(
        invalid_details,
        paste0(
          column_name,
          "{class=",
          paste(
            class(
              value
            ),
            collapse = "/"
          ),
          ", length=",
          length(
            value
          ),
          ", nrow=",
          nrow(
            table
          ),
          ", dim=",
          paste(
            dim(
              value
            ),
            collapse = "x"
          ),
          "}"
        )
      )
    }
  }

  if (length(
    invalid_details
  ) > 0L) {
    stop(
      label,
      " has invalid grouping column(s): ",
      paste(
        invalid_details,
        collapse = " | "
      ),
      call. = FALSE
    )
  }

  invisible(
    TRUE
  )
}

aggregate_to_animals <- function(section_table, group_columns, value_column) {
  section_table[, .(
    n_sections = uniqueN(section_id),
    animal_mean_contrast = safe_mean(get(value_column)),
    animal_median_contrast = safe_median(get(value_column)),
    positive_sections = sum(get(value_column) > 0),
    negative_sections = sum(get(value_column) < 0),
    minimum_section_contrast = min(get(value_column)),
    maximum_section_contrast = max(get(value_column))
  ), by = c("timepoint", "biological_replicate_id", group_columns)]
}

animal_inference <- function(animal_table, group_columns,
                             value_column = "animal_mean_contrast") {
  animal_table[, .(
    n_animals = .N,
    mean_animal_contrast = safe_mean(get(value_column)),
    median_animal_contrast = safe_median(get(value_column)),
    sd_animal_contrast = safe_sd(get(value_column)),
    minimum_animal_contrast = min(get(value_column)),
    maximum_animal_contrast = max(get(value_column)),
    positive_animals = sum(get(value_column) > 0),
    negative_animals = sum(get(value_column) < 0),
    exact_sign_test_p = exact_sign_p(get(value_column)),
    wilcoxon_signed_rank_p = safe_wilcox_signed_p(get(value_column)),
    inference_unit = "biological_animal",
    spot_level_p_values_used = FALSE
  ), by = c("timepoint", group_columns)]
}

build_target_priority <- function(
    pair_animals,
    candidate_support) {

  pair_animals <- data.table::as.data.table(
    pair_animals
  )

  candidate_support <- data.table::as.data.table(
    candidate_support
  )

  required_columns(
    pair_animals,
    c(
      "program",
      "timepoint",
      "biological_replicate_id",
      "candidate_receiver_id",
      "candidate_id",
      "annotation",
      "pathway_name",
      "interaction_name",
      "ligand_entity",
      "receptor_entity",
      "recipient_type",
      "contrast",
      "animal_mean_contrast"
    ),
    "Pair-level animal table for target priority"
  )

  target <- pair_animals[
    program ==
      "target_program"
  ]

  if (nrow(
    target
  ) == 0L) {
    stop(
      "No target-program pair-animal rows are available for priority ranking.",
      call. = FALSE
    )
  }

  identifier_columns <- c(
    "candidate_receiver_id",
    "candidate_id",
    "annotation",
    "pathway_name",
    "interaction_name",
    "ligand_entity",
    "receptor_entity",
    "recipient_type"
  )

  summary_group_columns <- c(
    identifier_columns,
    "contrast"
  )

  validate_grouping_columns(
    target,
    summary_group_columns,
    "Target priority input"
  )

  summary <- target[
    ,
    .(
      n_animals =
        .N,
      mean_animal_contrast =
        safe_mean(
          animal_mean_contrast
        ),
      minimum_animal_contrast =
        min(
          animal_mean_contrast
        ),
      positive_animals =
        sum(
          animal_mean_contrast >
            0
        ),
      mean_7days =
        safe_mean(
          animal_mean_contrast[
            timepoint ==
              "7days"
          ]
        ),
      mean_2months =
        safe_mean(
          animal_mean_contrast[
            timepoint ==
              "2months"
          ]
        )
    ),
    by =
      summary_group_columns
  ]

  required_contrasts <- c(
    "core_minus_remote",
    "border_minus_remote",
    "perilesional_minus_remote"
  )

  observed_contrasts <- unique(
    as.character(
      summary$contrast
    )
  )

  missing_contrasts <- setdiff(
    required_contrasts,
    observed_contrasts
  )

  if (length(
    missing_contrasts
  ) > 0L) {
    stop(
      "Target priority summary is missing contrast(s): ",
      paste(
        missing_contrasts,
        collapse = " | "
      ),
      ". Observed: ",
      paste(
        observed_contrasts,
        collapse = " | "
      ),
      call. = FALSE
    )
  }

  metric_columns <- c(
    "n_animals",
    "mean_animal_contrast",
    "minimum_animal_contrast",
    "positive_animals",
    "mean_7days",
    "mean_2months"
  )

  wide <- unique(
    summary[
      ,
      identifier_columns,
      with = FALSE
    ],
    by =
      identifier_columns
  )

  for (current_contrast in required_contrasts) {
    current <- summary[
      contrast ==
        current_contrast,
      c(
        identifier_columns,
        metric_columns
      ),
      with = FALSE
    ]

    duplicate_keys <- current[
      ,
      .N,
      by =
        identifier_columns
    ][
      N >
        1L
    ]

    if (nrow(
      duplicate_keys
    ) > 0L) {
      stop(
        "Duplicate target-priority rows for contrast ",
        current_contrast,
        ": ",
        nrow(
          duplicate_keys
        ),
        " duplicated candidate keys.",
        call. = FALSE
      )
    }

    renamed_metrics <- paste0(
      metric_columns,
      "_",
      current_contrast
    )

    data.table::setnames(
      current,
      metric_columns,
      renamed_metrics
    )

    wide <- merge(
      wide,
      current,
      by =
        identifier_columns,
      all.x =
        TRUE,
      sort =
        FALSE
    )
  }

  required_priority_columns <- c(
    "mean_animal_contrast_core_minus_remote",
    "mean_animal_contrast_border_minus_remote",
    "mean_animal_contrast_perilesional_minus_remote",
    "positive_animals_core_minus_remote",
    "positive_animals_border_minus_remote",
    "mean_7days_core_minus_remote",
    "mean_2months_core_minus_remote"
  )

  required_columns(
    wide,
    required_priority_columns,
    "Merged target priority table"
  )

  missing_priority_values <- vapply(
    required_priority_columns,
    function(column_name) {
      sum(
        !is.finite(
          as.numeric(
            wide[[
              column_name
            ]]
          )
        )
      )
    },
    integer(1)
  )

  if (any(
    missing_priority_values >
      0L
  )) {
    stop(
      "Target priority contains non-finite required metrics: ",
      paste(
        paste0(
          names(
            missing_priority_values
          ),
          "=",
          missing_priority_values
        ),
        collapse = " | "
      ),
      call. = FALSE
    )
  }

  core_mean <- as.numeric(
    wide[[
      "mean_animal_contrast_core_minus_remote"
    ]]
  )

  border_mean <- as.numeric(
    wide[[
      "mean_animal_contrast_border_minus_remote"
    ]]
  )

  perilesional_mean <- as.numeric(
    wide[[
      "mean_animal_contrast_perilesional_minus_remote"
    ]]
  )

  core_positive <- as.numeric(
    wide[[
      "positive_animals_core_minus_remote"
    ]]
  )

  border_positive <- as.numeric(
    wide[[
      "positive_animals_border_minus_remote"
    ]]
  )

  core_7days <- as.numeric(
    wide[[
      "mean_7days_core_minus_remote"
    ]]
  )

  core_2months <- as.numeric(
    wide[[
      "mean_2months_core_minus_remote"
    ]]
  )

  wide[[
    "priority_score"
  ]] <- core_mean +
    0.5 *
      border_mean +
    0.25 *
      perilesional_mean

  wide[[
    "priority_tier"
  ]] <- rep(
    "Tier_C_other",
    nrow(
      wide
    )
  )

  tier_B_index <-
    is.finite(
      core_positive
    ) &
    core_positive >=
      5L &
    is.finite(
      core_mean
    ) &
    core_mean >
      0

  wide[[
    "priority_tier"
  ]][
    tier_B_index
  ] <- "Tier_B_consistent_core"

  tier_A_index <-
    is.finite(
      core_positive
    ) &
    core_positive ==
      6L &
    is.finite(
      border_positive
    ) &
    border_positive >=
      5L &
    is.finite(
      core_7days
    ) &
    core_7days >
      0 &
    is.finite(
      core_2months
    ) &
    core_2months >
      0

  wide[[
    "priority_tier"
  ]][
    tier_A_index
  ] <- "Tier_A_persistent_core_border"

  support_columns <- c(
    "candidate_receiver_id",
    "source_Microglia_detection_fraction",
    "source_Microglia_mean_logCP10K",
    "recipient_receptor_detection_fraction",
    "recipient_receptor_mean_logCP10K"
  )

  required_columns(
    candidate_support,
    c(
      "candidate_pass",
      support_columns
    ),
    "Candidate reference-support table"
  )

  support <- candidate_support[
    candidate_pass ==
      TRUE,
    support_columns,
    with = FALSE
  ]

  support_duplicate_keys <- support[
    ,
    .N,
    by =
      "candidate_receiver_id"
  ][
    N >
      1L
  ]

  if (nrow(
    support_duplicate_keys
  ) > 0L) {
    support <- support[
      ,
      lapply(
        .SD,
        function(values) {
          values[
            1L
          ]
        }
      ),
      by =
        "candidate_receiver_id"
    ]
  }

  output <- merge(
    wide,
    support,
    by =
      "candidate_receiver_id",
    all.x =
      TRUE,
    sort =
      FALSE
  )

  if (nrow(
    output
  ) !=
      nrow(
        wide
      )) {
    stop(
      "Target priority support merge changed row count: before=",
      nrow(
        wide
      ),
      "; after=",
      nrow(
        output
      ),
      ".",
      call. = FALSE
    )
  }

  output[[
    "tier_order__"
  ]] <- match(
    output[[
      "priority_tier"
    ]],
    c(
      "Tier_A_persistent_core_border",
      "Tier_B_consistent_core",
      "Tier_C_other"
    )
  )

  data.table::setorderv(
    output,
    cols = c(
      "tier_order__",
      "priority_score",
      "mean_animal_contrast_core_minus_remote"
    ),
    order = c(
      1L,
      -1L,
      -1L
    ),
    na.last = TRUE
  )

  output[
    ,
    tier_order__ :=
      NULL
  ]

  output
}


# ---- plotting ----------------------------------------------------------------
draw_boxplot <- function(groups, labels, ylab, title, output_path) {
  if (
    length(
      groups
    ) == 0L ||
      length(
        labels
      ) !=
        length(
          groups
        )
  ) {
    stop(
      "Invalid boxplot groups/labels for ",
      basename(
        output_path
      ),
      ".",
      call. = FALSE
    )
  }

  groups <- lapply(
    groups,
    function(values) {
      values <- as.numeric(
        values
      )

      values[
        is.finite(
          values
        )
      ]
    }
  )

  if (any(
    lengths(
      groups
    ) == 0L
  )) {
    stop(
      "At least one animal-level boxplot group is empty in ",
      basename(
        output_path
      ),
      ".",
      call. = FALSE
    )
  }

  pdf(
    output_path,
    width = max(9, 0.55 * length(groups)),
    height = 6.5,
    useDingbats = FALSE
  )
  par(mar = c(8, 4.5, 2.8, 1))
  boxplot(
    groups, names = labels, ylab = ylab, main = title,
    las = 2, outline = FALSE
  )
  abline(h = 0, lty = 2)
  for (i in seq_along(groups)) {
    points(rep(i, length(groups[[i]])), groups[[i]], pch = 16, cex = 0.85)
  }
  dev.off()
}

draw_lr_heatmap <- function(priority, pair_inference, output_path) {
  priority <- data.table::as.data.table(
    priority
  )

  pair_inference <- data.table::as.data.table(
    pair_inference
  )

  required_columns(
    priority,
    c(
      "candidate_receiver_id",
      "priority_score",
      "pathway_name",
      "ligand_entity",
      "receptor_entity",
      "recipient_type"
    ),
    "LR priority table for heatmap"
  )

  top <- head(
    priority[
      is.finite(
        priority_score
      )
    ],
    MAX_TOP_LR_PLOT
  )

  if (nrow(
    top
  ) == 0L) {
    stop("No finite target LR priority row for heatmap.", call. = FALSE)
  }

  values <- pair_inference[
    program == "target_program" &
      candidate_receiver_id %in% top$candidate_receiver_id
  ]
  values[, heatmap_column := paste(timepoint, contrast, sep = "|")]

  wide <- data.table::dcast(
    values,
    candidate_receiver_id ~
      heatmap_column,
    value.var =
      "mean_animal_contrast"
  )

  wide <- wide[
    match(
      top$candidate_receiver_id,
      candidate_receiver_id
    )
  ]

  value_columns <- setdiff(
    names(
      wide
    ),
    "candidate_receiver_id"
  )

  if (length(
    value_columns
  ) == 0L) {
    stop(
      "No LR heatmap value columns were produced.",
      call. = FALSE
    )
  }

  matrix_values <- as.matrix(
    wide[
      ,
      value_columns,
      with = FALSE
    ]
  )
  rownames(matrix_values) <- paste0(
    top$pathway_name, ": ", top$ligand_entity, "→",
    top$receptor_entity, " | ", top$recipient_type
  )

  finite <- matrix_values[is.finite(matrix_values)]
  limit <- if (length(finite) > 0L) max(abs(finite)) else 1
  if (!is.finite(limit) || limit <= 0) limit <- 1

  pdf(
    output_path,
    width = 12,
    height = max(8, 0.28 * nrow(matrix_values) + 3),
    useDingbats = FALSE
  )
  par(mar = c(9, 16, 3, 2))
  image(
    seq_len(ncol(matrix_values)),
    seq_len(nrow(matrix_values)),
    t(matrix_values),
    col = hcl.colors(101, "Blue-Red 3"),
    zlim = c(-limit, limit),
    axes = FALSE, xlab = "", ylab = "",
    main = "Top target-linked ligand–receptor spatial-closure contrasts"
  )
  axis(1, at = seq_len(ncol(matrix_values)),
       labels = colnames(matrix_values), las = 2, cex.axis = 0.72)
  axis(2, at = seq_len(nrow(matrix_values)),
       labels = rownames(matrix_values), las = 2, cex.axis = 0.58)
  box()
  dev.off()
}

# ---- checks and fatal report --------------------------------------------------
checks <- data.table(
  check_id = character(),
  category = character(),
  description = character(),
  status = character(),
  observed = character(),
  required = character(),
  critical = logical(),
  action_if_not_pass = character()
)

add_check <- function(check_id, category, description, passed,
                      observed, required, critical = TRUE,
                      action_if_not_pass = "") {
  checks <<- rbind(
    checks,
    data.table(
      check_id = check_id,
      category = category,
      description = description,
      status = if (isTRUE(passed)) "PASS" else "FAIL",
      observed = as.character(observed),
      required = as.character(required),
      critical = isTRUE(critical),
      action_if_not_pass = action_if_not_pass
    ),
    fill = TRUE
  )
}

add_review <- function(check_id, category, description,
                       observed, required, action_if_not_pass) {
  checks <<- rbind(
    checks,
    data.table(
      check_id = check_id,
      category = category,
      description = description,
      status = "REVIEW",
      observed = as.character(observed),
      required = as.character(required),
      critical = FALSE,
      action_if_not_pass = action_if_not_pass
    ),
    fill = TRUE
  )
}

write_fatal_report <- function(reason) {
  writeLines(c(
    "============================================================",
    "GSE234774 Step 12G recipient-neighborhood / LR spatial closure",
    "============================================================",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "FINAL_STATUS: STOP",
    "",
    paste0("FATAL ERROR: ", reason),
    "",
    "- Existing successful section checkpoints were preserved.",
    "- No spot-level inferential P values were used.",
    "- The program was not renamed scar-associated.",
    "============================================================"
  ), OUT_FINAL_REPORT, useBytes = TRUE)
  log_message("FATAL ERROR: ", reason)
  write_run_log()
}

# ---- main --------------------------------------------------------------------
main <- function() {
  started <- Sys.time()

  log_message("============================================================")
  log_message("GSE234774 Step 12G started.")
  log_message("============================================================")

  # 1. Frozen input and release audit.
  required_inputs <- c(
    STEP12F_REPORT, STEP12F_SPOT_RDS,
    STEP12E_REPORT, STEP12E_HANDOFF_RDS,
    STEP12D_REPORT, STEP12D_WEIGHTS_RDS,
    STEP12C_REPORT, STEP12C_REFERENCE_RDS,
    STEP12C_REFERENCE_META_RDS, STEP12B_REPORT
  )
  input_ok <- all(file.exists(required_inputs)) &&
    dir.exists(SPATIAL_BPCELLS_DIR)

  add_check(
    "C01", "input",
    "All frozen upstream inputs and the spatial BPCells directory exist.",
    input_ok,
    paste0("files=", sum(file.exists(required_inputs)), "/",
           length(required_inputs),
           "; spatial_BPCells=", dir.exists(SPATIAL_BPCELLS_DIR)),
    "All files and BPCells directory present",
    TRUE,
    "Restore upstream outputs; do not reconstruct them."
  )
  if (!input_ok) stop("Required Step 12G input is missing.", call. = FALSE)

  get_status <- function(path) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    lines[grepl("^FINAL_STATUS:", lines)]
  }
  status <- c(
    Step12B = get_status(STEP12B_REPORT),
    Step12C = get_status(STEP12C_REPORT),
    Step12D = get_status(STEP12D_REPORT),
    Step12E = get_status(STEP12E_REPORT),
    Step12F = get_status(STEP12F_REPORT)
  )
  upstream_ok <- length(status) == 5L &&
    all(grepl("^FINAL_STATUS: (PASS|PASS_WITH_REVIEW)$", status))

  add_check(
    "C02", "upstream",
    "Step 12B through Step 12F release Step 12G.",
    upstream_ok,
    paste(names(status), status, sep = "=", collapse = "; "),
    "All PASS/PASS_WITH_REVIEW",
    TRUE,
    "Do not continue before upstream release."
  )
  if (!upstream_ok) stop("Upstream stage does not release Step 12G.",
                         call. = FALSE)

  log_message("[1/8] Loading zones, scores, RCTD weights, reference and genes.")

  spot_zones <- as.data.table(readRDS(STEP12F_SPOT_RDS))
  weights <- as.matrix(readRDS(STEP12D_WEIGHTS_RDS))
  spatial_counts <- BPCells::open_matrix_dir(SPATIAL_BPCELLS_DIR)
  handoff <- readRDS(STEP12E_HANDOFF_RDS)
  reference_object <- readRDS(STEP12C_REFERENCE_RDS)
  reference_meta <- as.data.frame(readRDS(STEP12C_REFERENCE_META_RDS))
  reference_counts <- get_reference_counts(reference_object)

  required_columns(
    spot_zones,
    c(
      "barcode", "capture_area_id", "section_id", "timepoint",
      "biological_replicate_id", "lesion_zone",
      "coordinate_x_fullres", "coordinate_y_fullres",
      "target_program_adjusted_z_within_section",
      "shared_prebranch_adjusted_z_within_section",
      "total_spatial_counts"
    ),
    "Step 12F spot table"
  )

  if (!all(c("target_genes_spatial", "shared_genes_spatial") %in%
           names(handoff))) {
    stop("Step 12E handoff lacks frozen target/shared genes.", call. = FALSE)
  }
  target_genes <- clean_genes(handoff$target_genes_spatial)
  shared_genes <- clean_genes(handoff$shared_genes_spatial)

  alignment_ok <- nrow(spot_zones) == EXPECTED_SPOTS &&
    nrow(weights) == EXPECTED_SPOTS &&
    ncol(spatial_counts) == EXPECTED_SPOTS &&
    !anyDuplicated(spot_zones$barcode) &&
    setequal(spot_zones$barcode, rownames(weights)) &&
    setequal(spot_zones$barcode, colnames(spatial_counts)) &&
    all(c("Microglia", RECIPIENT_TYPES) %in% colnames(weights)) &&
    length(intersect(toupper(target_genes), toupper(shared_genes))) == 0L

  add_check(
    "C03", "alignment",
    "Spots, weights and frozen target/shared programs are aligned.",
    alignment_ok,
    paste0("spots=", nrow(spot_zones),
           "; weight rows=", nrow(weights),
           "; target genes=", length(target_genes),
           "; shared genes=", length(shared_genes)),
    "33,941 aligned spots and disjoint frozen programs",
    TRUE,
    "Do not use misaligned spatial or program inputs."
  )
  if (!alignment_ok) stop("Frozen input alignment failed.", call. = FALSE)

  spot_zones <- spot_zones[match(rownames(weights), barcode)]
  if (!identical(spot_zones$barcode, rownames(weights))) {
    stop("Could not reorder spot zones to RCTD weights.", call. = FALSE)
  }

  animals_by_timepoint <- vapply(
    INJURED_TIMEPOINTS,
    function(current_timepoint) {
      uniqueN(spot_zones[
        timepoint == current_timepoint,
        biological_replicate_id
      ])
    },
    integer(1)
  )
  design_ok <- uniqueN(spot_zones$section_id) == EXPECTED_SECTIONS &&
    uniqueN(spot_zones[
      timepoint %in% INJURED_TIMEPOINTS, section_id
    ]) == EXPECTED_INJURED_SECTIONS &&
    all(animals_by_timepoint == EXPECTED_ANIMALS_PER_TIMEPOINT)

  add_check(
    "C04", "design",
    "The design contains 36 sections, 24 injured sections and 3 animals/timepoint.",
    design_ok,
    paste0("sections=", uniqueN(spot_zones$section_id),
           "; injured sections=",
           uniqueN(spot_zones[
             timepoint %in% INJURED_TIMEPOINTS, section_id
           ]),
           "; animals=", paste(
             names(animals_by_timepoint), animals_by_timepoint,
             sep = ":", collapse = " | "
           )),
    "36 sections, 24 injured, 3 animals/timepoint",
    TRUE,
    "Animal is the inferential unit."
  )
  if (!design_ok) stop("Section/animal design audit failed.", call. = FALSE)

  # 2. CellChatDB candidate generation and independent reference support.
  log_message("[2/8] Building CellChatDB candidates and snRNA support.")

  db <- load_mouse_db()
  interactions <- standardize_db(
    db, target_genes, shared_genes,
    rownames(spatial_counts), rownames(reference_counts)
  )

  mapped <- interactions[mapping_pass == TRUE]
  required_reference_genes <- unique(c(
    unlist(mapped$ligand_reference, use.names = FALSE),
    unlist(mapped$receptor_reference, use.names = FALSE)
  ))
  if (length(required_reference_genes) == 0L) {
    stop("No mapped target/shared CellChatDB interaction remains.",
         call. = FALSE)
  }

  reference_support <- build_reference_support(
    reference_counts, reference_meta,
    required_reference_genes,
    c("Microglia", RECIPIENT_TYPES)
  )
  candidates_all <- expand_receivers(interactions, reference_support)
  candidates <- candidates_all[candidate_pass == TRUE]

  candidate_counts <- candidates[, .(
    n_rows = .N,
    n_interactions = uniqueN(candidate_id),
    n_ligands = uniqueN(ligand_entity),
    n_receptors = uniqueN(receptor_entity)
  ), by = .(program, recipient_type)]

  target_all_recipients <- all(
    RECIPIENT_TYPES %in%
      candidate_counts[program == "target_program", recipient_type]
  )
  shared_available <- nrow(candidates[program == "Shared_prebranch"]) > 0L
  candidate_ok <- nrow(candidates) > 0L &&
    nrow(candidates) <= MAX_CANDIDATE_RECEIVER_ROWS &&
    target_all_recipients && shared_available

  add_check(
    "C05", "CellChatDB_candidates",
    "Reference-supported target candidates cover all recipients and a shared control exists.",
    candidate_ok,
    paste0("passed rows=", nrow(candidates),
           "; max=", MAX_CANDIDATE_RECEIVER_ROWS,
           "; target all recipients=", target_all_recipients,
           "; shared available=", shared_available),
    "Nonzero target/shared candidates and target coverage of five recipients",
    TRUE,
    "Review database schema, mapping or reference thresholds."
  )
  if (!candidate_ok) stop("Candidate-generation boundary failed.",
                          call. = FALSE)

  fwrite(candidate_to_csv(candidates_all), OUT_CANDIDATES, bom = TRUE)
  fwrite(
    reference_support, OUT_REFERENCE_SUPPORT,
    bom = TRUE, compress = "gzip"
  )
  signature <- paste0(
    unname(
      tools::md5sum(
        OUT_CANDIDATES
      )
    ),
    "|algorithm=v6_section_subsetting_fix",
    "|unused_pair_target_spearman=",
    CALCULATE_UNUSED_PAIR_TARGET_SPEARMAN,
    "|RANN_available=",
    requireNamespace(
      "RANN",
      quietly = TRUE
    ),
    "|matrixStats_available=",
    requireNamespace(
      "matrixStats",
      quietly = TRUE
    ),
    "|stop_on_first_section_failure=",
    STOP_ON_FIRST_SECTION_FAILURE
  )

  entity_maps <- build_entity_maps(candidates)

  input_audit <- rbindlist(list(
    data.table(
      item_type = "package",
      item = required_packages,
      value = vapply(
        required_packages,
        function(p) as.character(packageVersion(p)),
        character(1)
      ),
      note = "Version used in Step 12G"
    ),
    data.table(
      item_type = "parameter",
      item = c(
        "CellChatDB_annotations", "recipient_types",
        "MIN_SOURCE_DETECTION", "MIN_SOURCE_MEAN_LOGCP10K",
        "MIN_RECEPTOR_DETECTION", "MIN_RECEPTOR_MEAN_LOGCP10K",
        "K_NEIGHBORS", "MIN_NEIGHBORS", "MAX_NEIGHBOR_DISTANCE_NN",
        "target_genes", "shared_genes", "candidate_signature",
        "candidate_receiver_rows", "inferential_unit"
      ),
      value = c(
        paste(ALLOWED_DB_ANNOTATIONS, collapse = " | "),
        paste(RECIPIENT_TYPES, collapse = " | "),
        MIN_SOURCE_DETECTION, MIN_SOURCE_MEAN_LOGCP10K,
        MIN_RECEPTOR_DETECTION, MIN_RECEPTOR_MEAN_LOGCP10K,
        K_NEIGHBORS, MIN_NEIGHBORS, MAX_NEIGHBOR_DISTANCE_NN,
        length(target_genes), length(shared_genes), signature,
        nrow(candidates), "biological_animal"
      ),
      note = c(
        "Non-protein signaling excluded",
        "Prespecified recipient types",
        "Minimum across ligand complex subunits",
        "Minimum across ligand complex subunits",
        "Minimum across receptor complex subunits",
        "Minimum across receptor complex subunits",
        "Maximum graph neighbors",
        "Fallback minimum neighbors",
        "Radius in local NN units",
        "Frozen Step 12E target",
        "Frozen Step 12E parallel control",
        "MD5 of candidate audit",
        "No spatial-outcome filtering",
        "Four sections averaged within animal"
      )
    ),
    candidate_counts[, .(
      item_type = "candidate_count",
      item = paste(program, recipient_type, sep = "|"),
      value = n_rows,
      note = paste0(
        "interactions=", n_interactions,
        "; ligands=", n_ligands,
        "; receptors=", n_receptors
      )
    )]
  ), fill = TRUE)
  fwrite(input_audit, OUT_INPUT_AUDIT, bom = TRUE)

  missing_shared <- setdiff(
    RECIPIENT_TYPES,
    candidate_counts[program == "Shared_prebranch", recipient_type]
  )
  if (length(missing_shared) > 0L) {
    add_review(
      "R01", "parallel_control",
      "Shared_prebranch lacks a supported candidate for some recipients.",
      paste(missing_shared, collapse = " | "),
      "Use target-vs-shared only where both programs are available.",
      "Do not add unrelated control ligands."
    )
  }

  # 3. Run/reuse 36 section checkpoints.
  log_message("[3/8] Running or reusing 36 section checkpoints.")

  section_ids <- unique(
    as.character(
      spot_zones$section_id
    )
  )

  for (index in seq_along(
    section_ids
  )) {
    current_section_id <- section_ids[
      index
    ]

    log_message(
      "[",
      index,
      "/",
      length(
        section_ids
      ),
      "] ",
      current_section_id
    )

    if (
      RESUME_COMPLETED &&
        !current_section_id %in%
          FORCE_RERUN_SECTIONS &&
        section_complete(
          current_section_id,
          signature
        )
    ) {
      log_message(
        current_section_id,
        ": matching checkpoint found; skipping."
      )

      next
    }

    # Use a loop-variable name that cannot be confused with the data.table
    # column. This must yield one and only one section.
    current <- spot_zones[
      section_id ==
        current_section_id
    ]

    current_subset_ok <-
      nrow(
        current
      ) > 0L &&
      data.table::uniqueN(
        current$section_id
      ) ==
        1L &&
      identical(
        as.character(
          current$section_id[
            1L
          ]
        ),
        current_section_id
      )

    if (!current_subset_ok) {
      stop(
        "Section subset validation failed for ",
        current_section_id,
        ": rows=",
        nrow(
          current
        ),
        "; observed IDs=",
        paste(
          unique(
            as.character(
              current$section_id
            )
          ),
          collapse = " | "
        ),
        ".",
        call. = FALSE
      )
    }

    section_passed <- tryCatch(
      {
        score_section(
          current,
          weights,
          spatial_counts,
          candidates,
          entity_maps$ligand,
          entity_maps$receptor,
          signature,
          expected_section_id =
            current_section_id
        )

        TRUE
      },
      error = function(e) {
        write_section_failure(
          current_section_id,
          conditionMessage(
            e
          )
        )

        log_message(
          current_section_id,
          " failed: ",
          conditionMessage(
            e
          )
        )

        FALSE
      }
    )

    invisible(
      gc()
    )

    if (
      !section_passed &&
        isTRUE(
          STOP_ON_FIRST_SECTION_FAILURE
        )
    ) {
      stop(
        "Stopped after first failed section ",
        current_section_id,
        ". See its FAILURE.txt and the console error above.",
        call. = FALSE
      )
    }
  }

  completed <- section_ids[
    vapply(
      section_ids, section_complete, logical(1),
      signature = signature
    )
  ]
  incomplete <- setdiff(section_ids, completed)
  checkpoint_ok <- length(incomplete) == 0L

  add_check(
    "C06", "section_checkpoints",
    "All 36 sections have successful matching checkpoints.",
    checkpoint_ok,
    paste0("complete=", length(completed), "/36; incomplete=",
           paste(incomplete, collapse = " | ")),
    "36/36",
    TRUE,
    "Rerun; matching successful sections will be reused."
  )
  if (!checkpoint_ok) {
    fwrite(checks, OUT_CHECKS, bom = TRUE)
    stop("Incomplete sections: ", paste(incomplete, collapse = " | "),
         call. = FALSE)
  }


  # 4. Combine checkpoints.
  log_message("[4/8] Combining section outputs.")

  section_results <- lapply(
    section_ids,
    function(section_id) readRDS(section_files(section_id)$result)
  )

  spot_neighborhood <- rbindlist(
    lapply(section_results, `[[`, "neighborhood"),
    fill = TRUE
  )
  spot_neighborhood <- spot_neighborhood[
    match(spot_zones$barcode, barcode)
  ]

  pair_section <- rbindlist(
    lapply(section_results, `[[`, "pair_contrasts"),
    fill = TRUE
  )
  program_section <- rbindlist(
    lapply(section_results, `[[`, "program_contrasts"),
    fill = TRUE
  )
  section_audit <- rbindlist(
    lapply(section_results, `[[`, "audit"),
    fill = TRUE
  )

  score_columns <- grep(
    "^(neighbor_|target_neighbor_|shared_neighbor_)",
    names(spot_neighborhood),
    value = TRUE
  )
  spot_ok <- nrow(spot_neighborhood) == EXPECTED_SPOTS &&
    !anyNA(spot_neighborhood$barcode) &&
    !anyDuplicated(spot_neighborhood$barcode) &&
    identical(spot_neighborhood$barcode, spot_zones$barcode) &&
    all(is.finite(as.matrix(spot_neighborhood[, ..score_columns])))

  add_check(
    "C07", "neighborhood_output",
    "Recipient-neighborhood scores cover all spots exactly once and are finite.",
    spot_ok,
    paste0("rows=", nrow(spot_neighborhood),
           "; duplicated=", anyDuplicated(spot_neighborhood$barcode),
           "; sections=", uniqueN(spot_neighborhood$section_id)),
    "33,941 unique finite spot rows",
    TRUE,
    "Do not aggregate incomplete neighborhood outputs."
  )
  if (!spot_ok) stop("Combined spot-neighborhood output failed.",
                     call. = FALSE)

  fwrite(
    spot_neighborhood, OUT_SPOT_NEIGHBORHOOD,
    bom = TRUE, compress = "gzip"
  )
  saveRDS(spot_neighborhood, OUT_SPOT_RDS, compress = FALSE)

  # 5. Recipient-neighborhood section and animal contrasts.
  log_message("[5/8] Building recipient-neighborhood animal summaries.")

  neighborhood_section <- build_neighborhood_section_contrasts(
    spot_neighborhood
  )
  neighborhood_animals <- aggregate_to_animals(
    neighborhood_section,
    c("recipient_type", "metric", "contrast", "zone", "reference_zone"),
    "zone_minus_remote"
  )
  neighborhood_inference <- animal_inference(
    neighborhood_animals,
    c("recipient_type", "metric", "contrast", "zone", "reference_zone")
  )

  neighborhood_counts <- neighborhood_animals[, uniqueN(
    biological_replicate_id
  ), by = .(timepoint, recipient_type, metric, contrast)]
  neighborhood_unit_ok <- all(
    neighborhood_counts$V1 == EXPECTED_ANIMALS_PER_TIMEPOINT
  )

  add_check(
    "C08", "recipient_neighborhood",
    "Every recipient-neighborhood contrast has three animals/timepoint.",
    neighborhood_unit_ok,
    paste0("minimum animals=", min(neighborhood_counts$V1),
           "; maximum animals=", max(neighborhood_counts$V1)),
    "Exactly 3 animals/group",
    TRUE,
    "Do not use spots as independent replicates."
  )
  if (!neighborhood_unit_ok) {
    stop("Recipient-neighborhood animal-unit audit failed.", call. = FALSE)
  }

  fwrite(neighborhood_section, OUT_SECTION_NEIGHBORHOOD, bom = TRUE)
  fwrite(neighborhood_animals, OUT_ANIMAL_NEIGHBORHOOD, bom = TRUE)
  fwrite(neighborhood_inference, OUT_NEIGHBORHOOD_INFERENCE, bom = TRUE)

  # 6. Pair-level animal summaries and priorities.
  log_message("[6/8] Building LR and target-vs-shared animal summaries.")

  required_columns(
    pair_section,
    c(
      "section_id",
      "timepoint",
      "biological_replicate_id",
      "candidate_receiver_id",
      "candidate_id",
      "program",
      "annotation",
      "pathway_name",
      "interaction_name",
      "ligand_entity",
      "receptor_entity",
      "recipient_type",
      "contrast",
      "lesion_zone",
      "closure_zone_minus_remote"
    ),
    "Pair-level section contrasts"
  )

  pair_ok <- nrow(pair_section) > 0L &&
    all(is.finite(pair_section$closure_zone_minus_remote)) &&
    uniqueN(pair_section$section_id) == EXPECTED_INJURED_SECTIONS

  add_check(
    "C09", "LR_section_output",
    "Pair-level closure contrasts are finite for all injured sections.",
    pair_ok,
    paste0("rows=", nrow(pair_section),
           "; candidates=", uniqueN(pair_section$candidate_receiver_id),
           "; sections=", uniqueN(pair_section$section_id)),
    "Finite rows and 24 injured sections",
    TRUE,
    "Do not prioritize incomplete LR outputs."
  )
  if (!pair_ok) stop("Pair-level section closure failed.", call. = FALSE)

  # data.table version compatibility:
  # materialize constants as full-length columns before grouping.
  log_message("[6A/8] Materializing LR contrast grouping columns.")

  pair_section[
    ,
    reference_zone :=
      "remote"
  ]

  pair_group_columns <- c(
    "timepoint",
    "biological_replicate_id",
    "candidate_receiver_id",
    "candidate_id",
    "program",
    "annotation",
    "pathway_name",
    "interaction_name",
    "ligand_entity",
    "receptor_entity",
    "recipient_type",
    "contrast",
    "lesion_zone",
    "reference_zone"
  )

  validate_grouping_columns(
    pair_section,
    pair_group_columns,
    "Pair-level section contrasts"
  )

  log_message("[6B/8] Aggregating LR section contrasts to animals.")

  pair_animals <- pair_section[
    ,
    .(
      n_sections =
        uniqueN(
          section_id
        ),
      animal_mean_contrast =
        safe_mean(
          closure_zone_minus_remote
        ),
      animal_median_contrast =
        safe_median(
          closure_zone_minus_remote
        ),
      positive_sections =
        sum(
          closure_zone_minus_remote >
            0
        ),
      negative_sections =
        sum(
          closure_zone_minus_remote <
            0
        ),
      minimum_section_contrast =
        min(
          closure_zone_minus_remote
        ),
      maximum_section_contrast =
        max(
          closure_zone_minus_remote
        )
    ),
    by =
      pair_group_columns
  ]

  pair_inference_group_columns <- c(
    "timepoint",
    "candidate_receiver_id",
    "candidate_id",
    "program",
    "annotation",
    "pathway_name",
    "interaction_name",
    "ligand_entity",
    "receptor_entity",
    "recipient_type",
    "contrast",
    "lesion_zone",
    "reference_zone"
  )

  validate_grouping_columns(
    pair_animals,
    pair_inference_group_columns,
    "Pair-level animal contrasts"
  )

  log_message("[6C/8] Building LR animal-level inference.")

  pair_inference <- pair_animals[
    ,
    .(
      n_animals =
        .N,
      mean_animal_contrast =
        safe_mean(
          animal_mean_contrast
        ),
      median_animal_contrast =
        safe_median(
          animal_mean_contrast
        ),
      sd_animal_contrast =
        safe_sd(
          animal_mean_contrast
        ),
      minimum_animal_contrast =
        min(
          animal_mean_contrast
        ),
      maximum_animal_contrast =
        max(
          animal_mean_contrast
        ),
      positive_animals =
        sum(
          animal_mean_contrast >
            0
        ),
      negative_animals =
        sum(
          animal_mean_contrast <
            0
        ),
      exact_sign_test_p =
        exact_sign_p(
          animal_mean_contrast
        ),
      wilcoxon_signed_rank_p =
        safe_wilcox_signed_p(
          animal_mean_contrast
        ),
      inference_unit =
        "biological_animal",
      spot_level_p_values_used =
        FALSE
    ),
    by =
      pair_inference_group_columns
  ]

  pair_count_group_columns <- c(
    "timepoint",
    "candidate_receiver_id",
    "contrast"
  )

  validate_grouping_columns(
    pair_animals,
    pair_count_group_columns,
    "Pair-level animal-unit audit"
  )

  pair_counts <- pair_animals[
    ,
    .(
      n_animals =
        uniqueN(
          biological_replicate_id
        )
    ),
    by =
      pair_count_group_columns
  ]

  pair_unit_ok <- all(
    pair_counts$n_animals ==
      EXPECTED_ANIMALS_PER_TIMEPOINT
  )

  add_check(
    "C10",
    "LR_animal_unit",
    "Every LR contrast has three animals/timepoint.",
    pair_unit_ok,
    paste0(
      "minimum animals=",
      min(
        pair_counts$n_animals
      ),
      "; maximum animals=",
      max(
        pair_counts$n_animals
      )
    ),
    "Exactly 3 animals/group",
    TRUE,
    "Do not infer LR closure from spot-level replication."
  )

  if (!pair_unit_ok) {
    stop(
      "LR animal-unit audit failed.",
      call. = FALSE
    )
  }

  fwrite(
    pair_section,
    OUT_SECTION_LR,
    bom = TRUE,
    compress = "gzip"
  )

  fwrite(
    pair_animals,
    OUT_ANIMAL_LR,
    bom = TRUE,
    compress = "gzip"
  )

  fwrite(
    pair_inference,
    OUT_LR_INFERENCE,
    bom = TRUE,
    compress = "gzip"
  )

  log_message("[6D/8] Ranking target-linked LR candidates.")

  priority <- build_target_priority(
    pair_animals,
    candidates_all
  )

  if (
    nrow(
      priority
    ) == 0L
  ) {
    stop(
      "Target LR priority table is empty.",
      call. = FALSE
    )
  }

  log_message(
    "[6D/8] Target priority completed. Rows: ",
    nrow(
      priority
    ),
    "; finite scores: ",
    sum(
      is.finite(
        priority$priority_score
      )
    )
  )

  fwrite(
    priority,
    OUT_TARGET_PRIORITY,
    bom = TRUE
  )

  saveRDS(
    priority,
    OUT_PRIORITY_RDS,
    compress = FALSE
  )

  log_message("[6E/8] Building target-versus-Shared program summaries.")

  # Program-level target versus Shared_prebranch closure.
  program_section[
    ,
    reference_zone :=
      "remote"
  ]

  program_group_columns <- c(
    "timepoint",
    "biological_replicate_id",
    "program",
    "recipient_type",
    "annotation_group",
    "contrast",
    "lesion_zone",
    "reference_zone"
  )

  validate_grouping_columns(
    program_section,
    program_group_columns,
    "Program-level section contrasts"
  )

  program_animals <- program_section[
    ,
    .(
      n_sections =
        uniqueN(
          section_id
        ),
      animal_mean_contrast =
        safe_mean(
          program_closure_zone_minus_remote
        ),
      animal_median_contrast =
        safe_median(
          program_closure_zone_minus_remote
        ),
      positive_sections =
        sum(
          program_closure_zone_minus_remote >
            0
        ),
      negative_sections =
        sum(
          program_closure_zone_minus_remote <
            0
        )
    ),
    by =
      program_group_columns
  ]

  log_message("[6F/8] Merging paired target and Shared program columns.")

  program_pair_keys <- c(
    "timepoint",
    "biological_replicate_id",
    "recipient_type",
    "annotation_group",
    "contrast",
    "lesion_zone",
    "reference_zone"
  )

  program_all_protein <- program_animals[
    annotation_group ==
      "ALL_PROTEIN"
  ]

  validate_grouping_columns(
    program_all_protein,
    c(
      program_pair_keys,
      "program"
    ),
    "ALL_PROTEIN program animal rows"
  )

  target_program_rows <- program_all_protein[
    program ==
      "target_program",
    c(
      program_pair_keys,
      "animal_mean_contrast"
    ),
    with = FALSE
  ]

  shared_program_rows <- program_all_protein[
    program ==
      "Shared_prebranch",
    c(
      program_pair_keys,
      "animal_mean_contrast"
    ),
    with = FALSE
  ]

  data.table::setnames(
    target_program_rows,
    "animal_mean_contrast",
    "target_program"
  )

  data.table::setnames(
    shared_program_rows,
    "animal_mean_contrast",
    "Shared_prebranch"
  )

  target_duplicates <- target_program_rows[
    ,
    .N,
    by =
      program_pair_keys
  ][
    N >
      1L
  ]

  shared_duplicates <- shared_program_rows[
    ,
    .N,
    by =
      program_pair_keys
  ][
    N >
      1L
  ]

  if (
    nrow(
      target_duplicates
    ) >
      0L ||
      nrow(
        shared_duplicates
      ) >
        0L
  ) {
    stop(
      "Duplicate paired program rows: target=",
      nrow(
        target_duplicates
      ),
      "; shared=",
      nrow(
        shared_duplicates
      ),
      ".",
      call. = FALSE
    )
  }

  program_wide <- merge(
    target_program_rows,
    shared_program_rows,
    by =
      program_pair_keys,
    all =
      TRUE,
    sort =
      FALSE
  )

  program_columns_present <- isTRUE(
    all(
      c(
        "target_program",
        "Shared_prebranch"
      ) %in%
        names(
          program_wide
        )
    )
  )

  if (isTRUE(
    program_columns_present
  )) {
    program_wide[[
      "target_minus_shared_closure"
    ]] <- as.numeric(
      program_wide[[
        "target_program"
      ]]
    ) -
      as.numeric(
        program_wide[[
          "Shared_prebranch"
        ]]
      )
  } else {
    program_wide[[
      "target_minus_shared_closure"
    ]] <- rep(
      NA_real_,
      nrow(
        program_wide
      )
    )
  }

  incomplete_program_rows <- program_wide[
    !is.finite(
      target_minus_shared_closure
    )
  ]

  program_wide_paired <- program_wide[
    is.finite(
      target_minus_shared_closure
    )
  ]

  program_ok <- isTRUE(
    program_columns_present
  ) &&
    nrow(
      program_wide_paired
    ) >
      0L

  add_check(
    "C11",
    "target_vs_shared",
    "Target and Shared_prebranch closure have finite paired animal rows.",
    program_ok,
    paste0(
      "all rows=",
      nrow(
        program_wide
      ),
      "; paired rows=",
      nrow(
        program_wide_paired
      ),
      "; incomplete rows=",
      nrow(
        incomplete_program_rows
      ),
      "; target column=",
      "target_program" %in%
        names(
          program_wide
        ),
      "; shared column=",
      "Shared_prebranch" %in%
        names(
          program_wide
        )
    ),
    "At least one finite paired animal comparison",
    TRUE,
    "Do not claim target specificity without the frozen shared control."
  )

  if (!isTRUE(
    program_ok
  )) {
    stop(
      "Target-vs-Sharedprebranch comparison failed.",
      call. = FALSE
    )
  }

  if (
    nrow(
      incomplete_program_rows
    ) >
      0L
  ) {
    add_review(
      "R05",
      "target_vs_shared_coverage",
      "Some recipient/timepoint/contrast rows lack one frozen program.",
      paste0(
        "incomplete rows=",
        nrow(
          incomplete_program_rows
        )
      ),
      "Restrict target-minus-shared interpretation to finite paired rows.",
      "Do not impute missing Shared_prebranch or target closure."
    )
  }

  program_inference_group_columns <- c(
    "timepoint",
    "recipient_type",
    "annotation_group",
    "contrast",
    "lesion_zone",
    "reference_zone"
  )

  validate_grouping_columns(
    program_wide_paired,
    program_inference_group_columns,
    "Target-minus-shared paired animal rows"
  )

  log_message(
    "[6F/8] Paired program rows: ",
    nrow(
      program_wide
    ),
    "; finite pairs: ",
    sum(
      is.finite(
        program_wide$target_minus_shared_closure
      )
    )
  )

  log_message("[6G/8] Building target-minus-Shared animal inference.")

  program_inference <- program_wide_paired[
    ,
    .(
      n_animals =
        .N,
      mean_target_minus_shared =
        safe_mean(
          target_minus_shared_closure
        ),
      median_target_minus_shared =
        safe_median(
          target_minus_shared_closure
        ),
      sd_target_minus_shared =
        safe_sd(
          target_minus_shared_closure
        ),
      positive_animals =
        sum(
          target_minus_shared_closure >
            0
        ),
      negative_animals =
        sum(
          target_minus_shared_closure <
            0
        ),
      exact_sign_test_p =
        exact_sign_p(
          target_minus_shared_closure
        ),
      wilcoxon_signed_rank_p =
        safe_wilcox_signed_p(
          target_minus_shared_closure
        ),
      inference_unit =
        "biological_animal",
      spot_level_p_values_used =
        FALSE
    ),
    by =
      program_inference_group_columns
  ]

  fwrite(
    program_section,
    OUT_PROGRAM_SECTION,
    bom = TRUE
  )

  fwrite(
    program_animals,
    OUT_PROGRAM_ANIMAL,
    bom = TRUE
  )

  fwrite(
    program_wide,
    OUT_PROGRAM_DIFFERENTIAL,
    bom = TRUE
  )

  fwrite(
    program_inference,
    OUT_PROGRAM_INFERENCE,
    bom = TRUE
  )

  # 7. Figures.
  log_message("[7/8] Drawing summary figures.")

  neighborhood_plot <- neighborhood_animals[
    metric == "target_neighbor_colocalization" &
      contrast == "core_minus_remote"
  ]
  neighborhood_grid <- CJ(
    timepoint = INJURED_TIMEPOINTS,
    recipient_type = RECIPIENT_TYPES,
    sorted = FALSE
  )
  neighborhood_groups <- lapply(seq_len(nrow(neighborhood_grid)), function(i) {
    neighborhood_plot[
      timepoint == neighborhood_grid$timepoint[i] &
        recipient_type == neighborhood_grid$recipient_type[i],
      animal_mean_contrast
    ]
  })
  draw_boxplot(
    neighborhood_groups,
    paste(
      neighborhood_grid$timepoint,
      neighborhood_grid$recipient_type,
      sep = "\n"
    ),
    "Animal mean: target–recipient neighborhood core − remote",
    "Recipient-cell neighborhoods around adjusted target activity",
    OUT_FIG_NEIGHBORHOOD
  )

  program_plot <- program_animals[
    annotation_group == "ALL_PROTEIN" &
      contrast == "core_minus_remote"
  ]
  program_grid <- unique(
    program_plot[
      ,
      .(
        timepoint,
        recipient_type,
        program
      )
    ]
  )
  program_grid[
    ,
    timepoint_order__ :=
      match(
        timepoint,
        INJURED_TIMEPOINTS
      )
  ]
  program_grid[
    ,
    recipient_order__ :=
      match(
        recipient_type,
        RECIPIENT_TYPES
      )
  ]
  program_grid[
    ,
    program_order__ :=
      match(
        program,
        c(
          "target_program",
          "Shared_prebranch"
        )
      )
  ]
  setorder(
    program_grid,
    timepoint_order__,
    recipient_order__,
    program_order__
  )
  program_grid[
    ,
    c(
      "timepoint_order__",
      "recipient_order__",
      "program_order__"
    ) := NULL
  ]

  program_groups <- lapply(seq_len(nrow(program_grid)), function(i) {
    program_plot[
      timepoint == program_grid$timepoint[i] &
        recipient_type == program_grid$recipient_type[i] &
        program == program_grid$program[i],
      animal_mean_contrast
    ]
  })
  draw_boxplot(
    program_groups,
    paste(
      program_grid$timepoint,
      program_grid$recipient_type,
      program_grid$program,
      sep = "\n"
    ),
    "Animal mean program closure: core − remote",
    "Target-program versus Shared_prebranch LR closure",
    OUT_FIG_PROGRAM
  )

  draw_lr_heatmap(priority, pair_inference, OUT_FIG_LR)

  figure_paths <- c(
    OUT_FIG_NEIGHBORHOOD, OUT_FIG_PROGRAM, OUT_FIG_LR
  )
  figures_ok <- all(file.exists(figure_paths)) &&
    all(file_size_safe(figure_paths) > 0)

  add_check(
    "C12", "figures",
    "All three animal-level summary PDFs were written.",
    figures_ok,
    paste(
      basename(figure_paths), file_size_safe(figure_paths),
      sep = "=", collapse = "; "
    ),
    "Three nonempty PDFs",
    TRUE,
    "Review animal points and top-candidate heatmap."
  )
  if (!figures_ok) stop("Figure generation failed.", call. = FALSE)

  # 8. Review flags, handoff and final report.
  log_message("[8/8] Writing final report.")

  tier_A_count <- nrow(priority[
    priority_tier == "Tier_A_persistent_core_border"
  ])
  if (tier_A_count == 0L) {
    add_review(
      "R02", "target_LR_priority",
      "No target-linked LR candidate met the strict Tier A rule.",
      "Tier A count=0",
      "Retain broader lesion-proximal program interpretation.",
      "Do not promote isolated single-timepoint pairs."
    )
  }

  target_specific_core <- program_inference[
    contrast == "core_minus_remote"
  ]
  positive_fraction <- mean(
    target_specific_core$mean_target_minus_shared > 0
  )
  if (is.finite(positive_fraction) && positive_fraction < 0.60) {
    add_review(
      "R03", "target_specificity",
      "Target-minus-shared core closure is not positive in most summaries.",
      paste0("positive fraction=", round(positive_fraction, 3)),
      "Use recipient-specific rather than global target-specific claims.",
      "Do not omit Shared_prebranch."
    )
  }

  neighborhood_primary <- neighborhood_inference[
    metric == "target_neighbor_colocalization" &
      contrast == "core_minus_remote"
  ]
  if (!isTRUE(
    all(
      neighborhood_primary$mean_animal_contrast >
        0
    )
  )) {
    add_review(
      "R04", "recipient_neighborhood",
      "Not all recipient types have positive target-neighborhood core contrasts.",
      paste(
        paste0(
          neighborhood_primary$timepoint, ":",
          neighborhood_primary$recipient_type, "=",
          round(neighborhood_primary$mean_animal_contrast, 5)
        ),
        collapse = " | "
      ),
      "Report recipient-specific heterogeneity.",
      "Do not generalize one recipient niche to all recipients."
    )
  }

  critical_failures <- checks[critical == TRUE & status == "FAIL"]
  reviews <- checks[status == "REVIEW"]
  final_status <- if (nrow(critical_failures) > 0L) {
    "STOP"
  } else if (nrow(reviews) > 0L) {
    "PASS_WITH_REVIEW"
  } else {
    "PASS"
  }

  fwrite(checks, OUT_CHECKS, bom = TRUE)

  top_candidates <- head(
    priority,
    10L
  )[, .(
    priority_tier, pathway_name, ligand_entity,
    receptor_entity, recipient_type, priority_score,
    mean_animal_contrast_core_minus_remote,
    positive_animals_core_minus_remote
  )]

  saveRDS(list(
    step = "Step12G",
    completed_at = as.character(Sys.time()),
    final_status = final_status,
    CellChat_version = as.character(packageVersion("CellChat")),
    candidate_signature = signature,
    recipient_types = RECIPIENT_TYPES,
    candidate_receiver_rows_passed = nrow(candidates),
    target_Tier_A_candidates = tier_A_count,
    direct_binding_proven = FALSE,
    spot_level_p_values_used = FALSE,
    program_name_upgraded_to_scar_associated = FALSE,
    priority_CSV = OUT_TARGET_PRIORITY,
    target_vs_shared_CSV = OUT_PROGRAM_INFERENCE,
    next_step = paste(
      "Review recipient specificity and Tier A LR axes.",
      "If consistent, validate temporal expression in snRNA data",
      "and assemble manuscript-ready multi-layer evidence."
    )
  ), OUT_HANDOFF_RDS, compress = FALSE)

  neighborhood_lines <- paste0(
    "- ", neighborhood_primary$timepoint, " | ",
    neighborhood_primary$recipient_type,
    ": mean core-minus-remote=",
    round(neighborhood_primary$mean_animal_contrast, 5),
    "; positive animals=", neighborhood_primary$positive_animals,
    "/", neighborhood_primary$n_animals,
    "; exact sign P=",
    format(neighborhood_primary$exact_sign_test_p, digits = 4)
  )

  specificity_lines <- paste0(
    "- ", target_specific_core$timepoint, " | ",
    target_specific_core$recipient_type,
    ": target-minus-shared=",
    round(target_specific_core$mean_target_minus_shared, 6),
    "; positive animals=", target_specific_core$positive_animals,
    "/", target_specific_core$n_animals
  )

  top_lines <- if (nrow(top_candidates) > 0L) {
    paste0(
      "- ", top_candidates$priority_tier, " | ",
      top_candidates$pathway_name, " | ",
      top_candidates$ligand_entity, "→",
      top_candidates$receptor_entity, " | ",
      top_candidates$recipient_type, " | priority=",
      round(top_candidates$priority_score, 6),
      " | core positive animals=",
      top_candidates$positive_animals_core_minus_remote, "/6"
    )
  } else {
    "- No finite target candidate priority row."
  }

  report <- c(
    "============================================================",
    "GSE234774 Step 12G recipient-neighborhood / LR spatial closure",
    "============================================================",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("Output directory: ", OUTPUT_DIR),
    "",
    paste0("FINAL_STATUS: ", final_status),
    "",
    "METHOD BOUNDARY",
    "- CellChatDB.mouse supplied literature-supported LR priors.",
    "- Full discrete-label CellChat probabilities were not run because Visium spots are mixed-cell observations represented by continuous RCTD weights.",
    "- Ligand expression was Microglia-weighted and neighboring receptor expression was recipient-weighted.",
    "- Target and Shared_prebranch remained frozen and separate.",
    "- Step 12F lesion zones were not changed.",
    "- Spatial closure supports proximity and compatibility, not direct molecular binding.",
    "",
    "CANDIDATE COVERAGE",
    paste0("- CellChat version: ", packageVersion("CellChat"), "."),
    paste0("- Passed candidate-recipient rows: ", nrow(candidates), "."),
    paste0("- Target Tier A persistent candidates: ", tier_A_count, "."),
    paste0("- Successful section checkpoints: ", length(completed), "/36."),
    "",
    "RECIPIENT-NEIGHBORHOOD CORE RESULTS",
    neighborhood_lines,
    "",
    "TARGET MINUS SHAREDPREBRANCH CORE-CLOSURE RESULTS",
    specificity_lines,
    "",
    "TOP TARGET-LINKED LR CANDIDATES",
    top_lines,
    "",
    "OUTPUT COVERAGE",
    paste0("- Spots: ", nrow(spot_neighborhood), "/33,941."),
    paste0("- Injured sections: ", uniqueN(pair_section$section_id), "/24."),
    paste0(
      "- Independent injured animals: ",
      uniqueN(pair_animals[, paste(
        timepoint, biological_replicate_id, sep = "|"
      )]),
      "/6."
    ),
    "",
    "CRITICAL CHECKS",
    paste0(
      "- PASS: ",
      nrow(checks[critical == TRUE & status == "PASS"]),
      "; FAIL: ", nrow(critical_failures), "."
    )
  )

  if (nrow(reviews) > 0L) {
    report <- c(
      report, "", "REVIEW ITEMS",
      paste0(
        "- ", reviews$check_id, ": ", reviews$description,
        " Observed: ", reviews$observed,
        " Action: ", reviews$action_if_not_pass
      )
    )
  }

  report <- c(
    report,
    "",
    "INTERPRETATION",
    if (final_status %in% c("PASS", "PASS_WITH_REVIEW")) {
      "- Step 12G is technically complete. Recipient-neighborhood, CellChatDB-supported LR closure, and Shared_prebranch-control results are available at section and animal levels."
    } else {
      "- Step 12G is not released."
    },
    "- Strong evidence requires animal-level consistency, recipient specificity, and positive target-minus-Sharedprebranch closure.",
    "- With three animals/timepoint, exact P values are low-resolution; effect direction and replication are emphasized.",
    "- The program is not automatically renamed scar-associated.",
    "",
    "NEXT ACTION",
    "- Review the three PDFs, candidate support table, Tier A list, and target-minus-shared table before temporal snRNA validation.",
    "",
    "KEY OUTPUTS TO UPLOAD",
    "- tables/01_step12G_input_database_parameter_audit.csv",
    "- tables/02_CellChatDB_target_shared_LR_candidate_reference_support.csv",
    "- tables/04_spot_recipient_neighborhood_scores.csv.gz",
    "- tables/05_section_recipient_neighborhood_zone_contrasts.csv",
    "- tables/06_animal_recipient_neighborhood_zone_contrasts.csv",
    "- tables/07_animal_recipient_neighborhood_inference.csv",
    "- tables/08_section_LR_spatial_closure_zone_contrasts.csv.gz",
    "- tables/09_animal_LR_spatial_closure_zone_contrasts.csv.gz",
    "- tables/10_animal_LR_spatial_closure_inference.csv.gz",
    "- tables/11_target_LR_priority_summary.csv",
    "- tables/13_animal_target_vs_shared_program_closure_contrasts.csv",
    "- tables/15_target_minus_shared_closure_inference.csv",
    "- figures/01_recipient_neighborhood_animal_contrasts.pdf",
    "- figures/02_target_vs_Shared_prebranch_program_closure.pdf",
    "- figures/03_top_target_LR_spatial_closure_heatmap.pdf",
    "- tables/19_STEP12G_audit_checks.csv",
    "- 20_STEP12G_FINAL_REPORT.txt",
    "",
    "DO NOT UPLOAD UNLESS REQUESTED",
    paste0("- ", OUT_REFERENCE_SUPPORT),
    paste0("- ", OUT_SPOT_RDS),
    paste0("- ", OUT_PRIORITY_RDS),
    paste0("- ", OUT_HANDOFF_RDS),
    "- section_checkpoints/",
    "",
    paste0("TOTAL ELAPSED TIME: ", elapsed_minutes(started), " minutes"),
    "============================================================"
  )

  writeLines(report, OUT_FINAL_REPORT, useBytes = TRUE)
  writeLines(capture.output(sessionInfo()), SESSION_INFO, useBytes = TRUE)

  log_message("FINAL_STATUS: ", final_status)
  log_message("Passed candidate-recipient rows: ", nrow(candidates))
  log_message("Tier A candidates: ", tier_A_count)
  log_message("Section checkpoints: ", length(completed), "/36")
  log_message("Open first: ", OUT_FINAL_REPORT)
  write_run_log()

  cat("\n============================================================\n")
  cat("Step 12G completed: ", final_status, "\n", sep = "")
  cat("Candidate-recipient rows: ", nrow(candidates), "\n", sep = "")
  cat("Tier A target candidates: ", tier_A_count, "\n", sep = "")
  cat("Section checkpoints: ", length(completed), "/36\n", sep = "")
  cat("Open first: ", OUT_FINAL_REPORT, "\n", sep = "")
  cat("============================================================\n")

  invisible(list(
    final_status = final_status,
    target_priority = priority,
    target_vs_shared = program_inference,
    neighborhood_inference = neighborhood_inference,
    checks = checks
  ))
}

tryCatch(
  main(),
  error = function(e) {
    write_fatal_report(conditionMessage(e))
    message("Step 12G stopped: ", conditionMessage(e))
  }
)
