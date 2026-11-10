# ============================================================
# GLIA Reviewer 2 - Build public GSE162610 cell-level annotations (v2)
#
# GOAL
#   Create one public master annotation table for GSE162610 containing
#   cell barcodes, lineage/state information, transition/trajectory fields,
#   pseudotime, and branch labels used in the revised GLIA manuscript.
#
# IMPORTANT
#   This script reads EXACTLY TWO RDS files:
#     1) GSE162610_step2_full_myeloid_with_final_lineage.rds
#     2) GSE162610_step3B_clean_state_scored.rds
#
#   readRDS() necessarily deserializes the full Seurat object into memory,
#   including its assays/expression data. HOWEVER:
#     - no expression values are accessed;
#     - no Seurat analysis function is called;
#     - no normalization, PCA, clustering, DE, scoring, trajectory,
#       CellChat, SCENIC, or statistical analysis is run;
#     - only the meta.data slot is extracted;
#     - each RDS object is removed immediately after metadata extraction.
#
#   All remaining information is read only from existing metadata/result
#   CSV/TSV/XLSX files.
#
# OUTPUT
#   Desktop/GLIA_public_metadata_export_<timestamp>/
#     GSE162610_cell_annotations.csv.gz
#     GSE162610_cell_annotations_preview100.csv
#     00_source_file_manifest.csv
#     01_annotation_field_dictionary.csv
#     02_export_audit.csv
#     03_value_dictionary_branch_and_state.csv
#     README.txt
#
# This is an EXPORT/ASSEMBLY script only, not an analysis script.
# ============================================================

options(stringsAsFactors = FALSE)

# --------------------------- 1. Paths ---------------------------
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
DESKTOP      <- Sys.getenv("GLIA_METADATA_EXPORT_DIR", unset = PROJECT_ROOT)

if (!dir.exists(PROJECT_ROOT)) stop("Project root not found: ", PROJECT_ROOT)
if (!dir.exists(DESKTOP)) stop("Desktop not found: ", DESKTOP)

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTDIR <- file.path(DESKTOP, paste0("GLIA_public_metadata_export_", stamp))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# --------------------------- 2. Frozen expected counts ---------------------------
# These counts correspond to the revised analysis chain already audited.
EXPECTED <- list(
  atlas_cells       = 62842L,
  step2_myeloid     = 32405L,
  clean_microglia   = 19958L,
  trajectory_cells  = 14005L,
  formal_branch     = 4203L
)

# --------------------------- 3. Helper functions ---------------------------
normalize_slashes <- function(x) gsub("\\\\", "/", x)

find_exact_one <- function(root, filename) {
  hits <- list.files(
    root,
    recursive = TRUE,
    full.names = TRUE,
    all.files = FALSE,
    no.. = TRUE
  )
  hits <- hits[basename(hits) == filename]

  if (length(hits) == 0L) {
    stop("Required file not found:\n  ", filename)
  }
  if (length(hits) > 1L) {
    stop(
      "More than one file with this exact name was found:\n  ",
      filename, "\n\nCandidates:\n",
      paste0("  - ", normalize_slashes(hits), collapse = "\n"),
      "\n\nPlease keep only the intended project copy or set the path manually."
    )
  }
  normalizePath(hits, winslash = "/", mustWork = TRUE)
}

pick_col <- function(df, candidates, required = TRUE, label = NULL) {
  # Candidates are listed in PRIORITY ORDER.
  # If more than one candidate exists, use the first one and report the choice.
  found <- candidates[candidates %in% names(df)]

  if (length(found) >= 1L) {
    chosen <- found[1L]

    if (length(found) > 1L) {
      message(
        "Multiple candidate columns found for ",
        ifelse(is.null(label), "field", label),
        ": ", paste(found, collapse = ", "),
        ". Using highest-priority column: ", chosen
      )
    }

    return(chosen)
  }

  if (required) {
    stop(
      "Could not find required column for ",
      ifelse(is.null(label), "field", label),
      ". Tried:\n",
      paste0("  - ", candidates, collapse = "\n"),
      "\nAvailable columns:\n",
      paste(names(df), collapse = ", ")
    )
  }

  NA_character_
}

extract_meta_only <- function(rds_path, object_label) {
  message("Reading RDS ONLY to export meta.data: ", basename(rds_path))

  obj <- readRDS(rds_path)

  md <- NULL

  if (isS4(obj) && "meta.data" %in% methods::slotNames(obj)) {
    md <- methods::slot(obj, "meta.data")
  } else if (is.list(obj) && !is.null(obj$meta.data)) {
    md <- obj$meta.data
  } else {
    rm(obj)
    gc(verbose = FALSE)
    stop(object_label, ": object does not expose a meta.data slot/component.")
  }

  md <- as.data.frame(md, stringsAsFactors = FALSE)

  # Preserve the cell identifier exactly as stored in metadata row names,
  # unless an explicit unambiguous cell/barcode column already exists.
  explicit_cell_cols <- c("cell", "cell_id", "barcode")
  explicit_found <- explicit_cell_cols[explicit_cell_cols %in% names(md)]

  if (length(explicit_found) == 1L) {
    md$.__cell_id__ <- as.character(md[[explicit_found]])
  } else if (length(explicit_found) > 1L) {
    vals <- lapply(explicit_found, function(z) as.character(md[[z]]))
    same <- all(vapply(vals[-1], function(v) identical(v, vals[[1]]), logical(1)))
    if (!same) {
      rm(obj)
      gc(verbose = FALSE)
      stop(
        object_label,
        ": multiple non-identical cell identifier columns detected: ",
        paste(explicit_found, collapse = ", ")
      )
    }
    md$.__cell_id__ <- vals[[1]]
  } else {
    md$.__cell_id__ <- rownames(md)
  }

  if (anyNA(md$.__cell_id__) || any(md$.__cell_id__ == "")) {
    rm(obj)
    gc(verbose = FALSE)
    stop(object_label, ": missing cell identifiers in metadata.")
  }

  if (anyDuplicated(md$.__cell_id__)) {
    rm(obj)
    gc(verbose = FALSE)
    stop(object_label, ": duplicated cell identifiers in metadata.")
  }

  # Remove the full object immediately after metadata extraction.
  rm(obj)
  invisible(gc(verbose = FALSE))

  md
}

read_xlsx_sheet <- function(path, sheet) {
  if (requireNamespace("readxl", quietly = TRUE)) {
    return(as.data.frame(
      readxl::read_excel(path, sheet = sheet),
      stringsAsFactors = FALSE
    ))
  }

  if (requireNamespace("openxlsx", quietly = TRUE)) {
    return(as.data.frame(
      openxlsx::read.xlsx(path, sheet = sheet),
      stringsAsFactors = FALSE
    ))
  }

  stop(
    "Reading the trajectory XLSX requires either the 'readxl' or 'openxlsx' ",
    "package. Install ONE of them, then rerun this script."
  )
}

assert_n <- function(x, expected, label) {
  n <- nrow(x)
  if (!identical(as.integer(n), as.integer(expected))) {
    stop(
      label, " row count mismatch.\n",
      "Expected: ", expected, "\n",
      "Observed: ", n, "\n",
      "The script stops rather than silently exporting a non-frozen dataset."
    )
  }
}

assert_unique_cells <- function(x, cell_col, label) {
  ids <- as.character(x[[cell_col]])
  if (anyNA(ids) || any(ids == "")) stop(label, ": missing cell IDs.")
  if (anyDuplicated(ids)) stop(label, ": duplicated cell IDs.")
}

add_by_match <- function(master, source, master_id, source_id, mapping) {
  idx <- match(master[[master_id]], source[[source_id]])

  for (dest in names(mapping)) {
    src <- mapping[[dest]]
    if (is.na(src) || !(src %in% names(source))) {
      master[[dest]] <- NA
    } else {
      master[[dest]] <- source[[src]][idx]
    }
  }

  master
}

yes_no <- function(x) ifelse(x, TRUE, FALSE)

# --------------------------- 4. Locate required files ---------------------------
RDS_STEP2 <- find_exact_one(
  PROJECT_ROOT,
  "GSE162610_step2_full_myeloid_with_final_lineage.rds"
)

RDS_STEP3 <- find_exact_one(
  PROJECT_ROOT,
  "GSE162610_step3B_clean_state_scored.rds"
)

ATLAS_TSV <- find_exact_one(
  PROJECT_ROOT,
  "02_CellChat_cell_metadata.tsv.gz"
)

TRAJ_XLSX <- find_exact_one(
  PROJECT_ROOT,
  "Table5BR2_single_root_branch_audit.xlsx"
)

BRANCH_CSV <- find_exact_one(
  PROJECT_ROOT,
  "GSE162610_step7A_frozen_4203_cell_metadata.csv"
)

source_paths <- c(
  RDS_STEP2,
  RDS_STEP3,
  ATLAS_TSV,
  TRAJ_XLSX,
  BRANCH_CSV
)

# --------------------------- 5. Read base 62,842-cell atlas metadata ---------------------------
message("Reading existing atlas metadata table...")
atlas <- utils::read.delim(
  gzfile(ATLAS_TSV, open = "rt"),
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

assert_n(atlas, EXPECTED$atlas_cells, "Step8A atlas metadata")

atlas_cell_col <- pick_col(
  atlas,
  c("cell", "cell_id", "barcode"),
  TRUE,
  "atlas cell identifier"
)
assert_unique_cells(atlas, atlas_cell_col, "Step8A atlas metadata")

atlas_sample_col <- pick_col(
  atlas,
  c("sample_id", "sample", "orig.ident"),
  TRUE,
  "sample ID"
)

atlas_time_col <- pick_col(
  atlas,
  c("timepoint_standard", "timepoint", "timepoint_raw"),
  TRUE,
  "timepoint"
)

atlas_broad_col <- pick_col(
  atlas,
  c("broad_type_reconstructed", "broad_cell_type", "broad_type"),
  TRUE,
  "broad cell type"
)

atlas_clean_col <- pick_col(
  atlas,
  c("is_clean_microglia", "clean_microglia"),
  FALSE,
  "clean microglia flag"
)

atlas_frozen_branch_col <- pick_col(
  atlas,
  c("frozen_branch_state", "branch_state"),
  FALSE,
  "frozen branch state"
)

atlas_cellchat_group_col <- pick_col(
  atlas,
  c("cellchat_group"),
  FALSE,
  "CellChat group"
)

master <- data.frame(
  cell_id = as.character(atlas[[atlas_cell_col]]),
  sample_id = as.character(atlas[[atlas_sample_col]]),
  timepoint = as.character(atlas[[atlas_time_col]]),
  broad_cell_type = as.character(atlas[[atlas_broad_col]]),
  stringsAsFactors = FALSE
)

if (!is.na(atlas_clean_col)) {
  master$clean_microglia_step8A <- atlas[[atlas_clean_col]]
} else {
  master$clean_microglia_step8A <- NA
}

if (!is.na(atlas_frozen_branch_col)) {
  master$frozen_branch_state_step8A <- as.character(atlas[[atlas_frozen_branch_col]])
} else {
  master$frozen_branch_state_step8A <- NA_character_
}

if (!is.na(atlas_cellchat_group_col)) {
  master$cellchat_group <- as.character(atlas[[atlas_cellchat_group_col]])
} else {
  master$cellchat_group <- NA_character_
}

# --------------------------- 6. Read Step2C RDS metadata ONLY ---------------------------
step2 <- extract_meta_only(RDS_STEP2, "Step2C final myeloid object")
assert_n(step2, EXPECTED$step2_myeloid, "Step2C final myeloid metadata")
assert_unique_cells(step2, ".__cell_id__", "Step2C final myeloid metadata")

step2_lineage_col <- pick_col(
  step2,
  c(
    "final_step2_lineage",
    "final_lineage_class",
    "final_myeloid_lineage",
    "final_lineage"
  ),
  TRUE,
  "final Step2 lineage"
)

step2_cluster_col <- pick_col(
  step2,
  c(
    "myeloid_cluster_r04",
    "step2_final_cluster",
    "seurat_clusters"
  ),
  FALSE,
  "Step2 myeloid cluster"
)

# Keep only the minimal metadata required for public annotation.
step2_min <- data.frame(
  cell_id = as.character(step2$.__cell_id__),
  myeloid_lineage = as.character(step2[[step2_lineage_col]]),
  stringsAsFactors = FALSE
)

if (!is.na(step2_cluster_col)) {
  step2_min$myeloid_cluster <- as.character(step2[[step2_cluster_col]])
} else {
  step2_min$myeloid_cluster <- NA_character_
}

rm(step2)
invisible(gc(verbose = FALSE))

master$myeloid_included <- master$cell_id %in% step2_min$cell_id
master <- add_by_match(
  master, step2_min,
  "cell_id", "cell_id",
  c(
    myeloid_lineage = "myeloid_lineage",
    myeloid_cluster = "myeloid_cluster"
  )
)

# --------------------------- 7. Read Step3B RDS metadata ONLY ---------------------------
step3 <- extract_meta_only(RDS_STEP3, "Step3B clean resident microglia object")
assert_n(step3, EXPECTED$clean_microglia, "Step3B clean resident metadata")
assert_unique_cells(step3, ".__cell_id__", "Step3B clean resident metadata")

resident_cluster_col <- pick_col(
  step3,
  c(
    "step2_final_cluster",
    "step3A_cluster_r04",
    "seurat_clusters"
  ),
  TRUE,
  "resident-state cluster"
)

# Preserve an existing explicit final state field if one exists AND is populated.
state_candidates <- c(
  "resident_state",
  "final_resident_state",
  "final_state",
  "state_label",
  "state_annotation",
  "final_state_annotation"
)
state_found <- state_candidates[state_candidates %in% names(step3)]

existing_state_col <- NA_character_
if (length(state_found)) {
  nonempty <- state_found[vapply(
    state_found,
    function(z) {
      v <- as.character(step3[[z]])
      any(!is.na(v) & nzchar(v))
    },
    logical(1)
  )]

  if (length(nonempty) == 1L) {
    existing_state_col <- nonempty
  } else if (length(nonempty) > 1L) {
    stop(
      "More than one populated resident-state metadata field exists in Step3B: ",
      paste(nonempty, collapse = ", "),
      ". Please resolve explicitly rather than guessing."
    )
  }
}

step3_min <- data.frame(
  cell_id = as.character(step3$.__cell_id__),
  resident_cluster = as.character(step3[[resident_cluster_col]]),
  stringsAsFactors = FALSE
)

if (!is.na(existing_state_col)) {
  step3_min$resident_state <- as.character(step3[[existing_state_col]])
  step3_min$resident_state_source <- paste0(
    "existing Step3B metadata field: ", existing_state_col
  )
} else {
  # Frozen Figure 1 / Figure S2 cluster-to-state mapping used in the revised
  # manuscript. No clustering or state inference is rerun here.
  cl <- as.character(step3_min$resident_cluster)

  state <- rep(NA_character_, length(cl))
  state[cl %in% c("0")]     <- "Canonical resident"
  state[cl %in% c("3", "4")] <- "Activated resident"
  state[cl %in% c("9")]     <- "IFN-responsive resident"
  state[cl %in% c("5", "7")] <- "Cycling resident"

  step3_min$resident_state <- state
  step3_min$resident_state_source <- ifelse(
    is.na(state),
    NA_character_,
    "frozen revised-manuscript Figure 1 cluster-to-state mapping"
  )
}

rm(step3)
invisible(gc(verbose = FALSE))

master$clean_resident_microglia <- master$cell_id %in% step3_min$cell_id
master <- add_by_match(
  master, step3_min,
  "cell_id", "cell_id",
  c(
    resident_cluster = "resident_cluster",
    resident_state = "resident_state",
    resident_state_source = "resident_state_source"
  )
)

# --------------------------- 8. Read 14,005-cell trajectory metadata ---------------------------
message("Reading trajectory Cell_results sheet...")
traj <- read_xlsx_sheet(TRAJ_XLSX, "Cell_results")
assert_n(traj, EXPECTED$trajectory_cells, "Trajectory Cell_results")

traj_cell_col <- pick_col(
  traj,
  c("cell", "cell_id", "barcode"),
  TRUE,
  "trajectory cell identifier"
)
assert_unique_cells(traj, traj_cell_col, "Trajectory Cell_results")

traj_transition_col <- pick_col(
  traj,
  c("transition_status"),
  FALSE,
  "transition status"
)
traj_highconf_col <- pick_col(
  traj,
  c("shared_transition_3dpi_highconf"),
  FALSE,
  "3-dpi high-confidence transition flag"
)
traj_pt_col <- pick_col(
  traj,
  c("pseudotime_single_root", "pseudotime"),
  TRUE,
  "final pseudotime"
)
traj_vertex_col <- pick_col(
  traj,
  c("principal_vertex"),
  FALSE,
  "principal vertex"
)
traj_branch_arm_col <- pick_col(
  traj,
  c("branch_arm"),
  TRUE,
  "branch arm"
)
traj_component_col <- pick_col(
  traj,
  c("branch_component"),
  FALSE,
  "branch component"
)

traj_min <- data.frame(
  cell_id = as.character(traj[[traj_cell_col]]),
  pseudotime = traj[[traj_pt_col]],
  trajectory_branch_arm_raw = as.character(traj[[traj_branch_arm_col]]),
  stringsAsFactors = FALSE
)

traj_min$transition_status <- if (!is.na(traj_transition_col)) {
  as.character(traj[[traj_transition_col]])
} else NA_character_

traj_min$shared_transition_3dpi_highconf <- if (!is.na(traj_highconf_col)) {
  traj[[traj_highconf_col]]
} else NA

traj_min$principal_vertex <- if (!is.na(traj_vertex_col)) {
  as.character(traj[[traj_vertex_col]])
} else NA_character_

traj_min$branch_component <- if (!is.na(traj_component_col)) {
  as.character(traj[[traj_component_col]])
} else NA_character_

master$trajectory_included <- master$cell_id %in% traj_min$cell_id
master <- add_by_match(
  master, traj_min,
  "cell_id", "cell_id",
  c(
    transition_status = "transition_status",
    shared_transition_3dpi_highconf = "shared_transition_3dpi_highconf",
    pseudotime = "pseudotime",
    principal_vertex = "principal_vertex",
    trajectory_branch_arm_raw = "trajectory_branch_arm_raw",
    branch_component = "branch_component"
  )
)

# --------------------------- 9. Read frozen 4,203 formal-branch metadata ---------------------------
message("Reading frozen 4,203-cell formal branch metadata...")
branch <- utils::read.csv(
  BRANCH_CSV,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
assert_n(branch, EXPECTED$formal_branch, "Step7A frozen branch metadata")

branch_cell_col <- pick_col(
  branch,
  c("cell", "cell_id", "barcode"),
  TRUE,
  "formal branch cell identifier"
)
assert_unique_cells(branch, branch_cell_col, "Step7A frozen branch metadata")

branch_state_col <- pick_col(
  branch,
  c(
    "step5CR_state",
    "frozen_branch_state",
    "branch_state",
    "branch_label"
  ),
  FALSE,
  "formal branch state"
)

branch_min <- data.frame(
  cell_id = as.character(branch[[branch_cell_col]]),
  stringsAsFactors = FALSE
)

branch_min$formal_branch_state_raw <- if (!is.na(branch_state_col)) {
  as.character(branch[[branch_state_col]])
} else {
  NA_character_
}

master$formal_branch_included <- master$cell_id %in% branch_min$cell_id
master <- add_by_match(
  master, branch_min,
  "cell_id", "cell_id",
  c(formal_branch_state_raw = "formal_branch_state_raw")
)

# --------------------------- 10. Public branch labels ---------------------------
# Preserve raw source labels AND create a reader-friendly public label.
#
# Previously audited frozen correspondence:
#   Arm_1 = Alternative branch
#   Arm_2 = Target-program branch
#   Pre_branch_root_side / Branchpoint_vertex = Shared prebranch
#
# Only formal_branch_included cells receive a formal public branch label.

master$branch_label <- NA_character_

formal <- master$formal_branch_included %in% TRUE
arm <- as.character(master$trajectory_branch_arm_raw)

master$branch_label[formal & arm %in% c("Arm_1")] <- "Alternative"
master$branch_label[formal & arm %in% c("Arm_2")] <- "Target"
master$branch_label[
  formal & arm %in% c("Pre_branch_root_side", "Branchpoint_vertex")
] <- "Shared_prebranch"

# If the frozen Step7A state already contains recognizable public labels,
# use it only for still-unassigned formal cells; never overwrite Arm_1/Arm_2.
raw_state <- as.character(master$formal_branch_state_raw)
fill <- formal & is.na(master$branch_label) & !is.na(raw_state)

master$branch_label[fill & grepl("shared|prebranch", raw_state, ignore.case = TRUE)] <-
  "Shared_prebranch"
master$branch_label[fill & grepl("target", raw_state, ignore.case = TRUE)] <-
  "Target"
master$branch_label[fill & grepl("alternative|alt", raw_state, ignore.case = TRUE)] <-
  "Alternative"

master$target_program_branch <- ifelse(
  is.na(master$branch_label),
  NA,
  master$branch_label == "Target"
)

# --------------------------- 11. Consistency checks ---------------------------
# 11A. Every source cell must map to the 62,842-cell atlas.
unmatched_step2 <- setdiff(step2_min$cell_id, master$cell_id)
unmatched_step3 <- setdiff(step3_min$cell_id, master$cell_id)
unmatched_traj  <- setdiff(traj_min$cell_id, master$cell_id)
unmatched_branch <- setdiff(branch_min$cell_id, master$cell_id)

if (length(unmatched_step2) ||
    length(unmatched_step3) ||
    length(unmatched_traj) ||
    length(unmatched_branch)) {
  stop(
    "At least one downstream cell ID does not map to the 62,842-cell atlas.\n",
    "Unmatched Step2: ", length(unmatched_step2), "\n",
    "Unmatched Step3: ", length(unmatched_step3), "\n",
    "Unmatched trajectory: ", length(unmatched_traj), "\n",
    "Unmatched formal branch: ", length(unmatched_branch), "\n",
    "Export stopped to prevent an inconsistent public annotation file."
  )
}

# 11B. Verify membership counts after matching.
membership_counts <- c(
  atlas_cells = nrow(master),
  step2_myeloid = sum(master$myeloid_included %in% TRUE, na.rm = TRUE),
  clean_microglia = sum(master$clean_resident_microglia %in% TRUE, na.rm = TRUE),
  trajectory_cells = sum(master$trajectory_included %in% TRUE, na.rm = TRUE),
  formal_branch = sum(master$formal_branch_included %in% TRUE, na.rm = TRUE)
)

expected_vec <- unlist(EXPECTED)
if (!all(as.integer(membership_counts) == as.integer(expected_vec[names(membership_counts)]))) {
  stop(
    "Membership count audit failed.\nObserved:\n",
    paste0("  ", names(membership_counts), " = ", membership_counts, collapse = "\n")
  )
}

# 11C. Formal branch cells should be a subset of trajectory cells.
bad_formal <- master$formal_branch_included %in% TRUE &
  !(master$trajectory_included %in% TRUE)

if (any(bad_formal)) {
  stop("Formal branch membership is not fully contained within trajectory cells.")
}

# 11D. Preserve missing state labels rather than guessing.
n_clean_missing_state <- sum(
  master$clean_resident_microglia %in% TRUE &
    (is.na(master$resident_state) | master$resident_state == "")
)

# --------------------------- 12. Reorder public columns ---------------------------
public_cols <- c(
  "cell_id",
  "sample_id",
  "timepoint",
  "broad_cell_type",

  "myeloid_included",
  "myeloid_cluster",
  "myeloid_lineage",

  "clean_resident_microglia",
  "resident_cluster",
  "resident_state",

  "transition_status",
  "shared_transition_3dpi_highconf",

  "trajectory_included",
  "pseudotime",
  "principal_vertex",
  "trajectory_branch_arm_raw",
  "branch_component",

  "formal_branch_included",
  "formal_branch_state_raw",
  "branch_label",
  "target_program_branch",

  "cellchat_group"
)

master_public <- master[, public_cols, drop = FALSE]

# --------------------------- 13. Write master annotation ---------------------------
out_gz <- file.path(OUTDIR, "GSE162610_cell_annotations.csv.gz")

con <- gzfile(out_gz, open = "wt")
utils::write.csv(
  master_public,
  con,
  row.names = FALSE,
  na = ""
)
close(con)

utils::write.csv(
  utils::head(master_public, 100),
  file.path(OUTDIR, "GSE162610_cell_annotations_preview100.csv"),
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)

# --------------------------- 14. Source manifest ---------------------------
fi <- file.info(source_paths)

source_manifest <- data.frame(
  role = c(
    "Step2C final myeloid lineage metadata",
    "Step3B clean resident microglia metadata",
    "Step8A 62,842-cell atlas/CellChat metadata",
    "Step5BR2 14,005-cell trajectory metadata",
    "Step7A frozen 4,203-cell formal branch metadata"
  ),
  source_path = normalize_slashes(source_paths),
  size_bytes = fi$size,
  modified_time = as.character(fi$mtime),
  md5 = unname(tools::md5sum(source_paths)),
  stringsAsFactors = FALSE
)

utils::write.csv(
  source_manifest,
  file.path(OUTDIR, "00_source_file_manifest.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# --------------------------- 15. Field dictionary ---------------------------
field_dictionary <- data.frame(
  field = public_cols,
  description = c(
    "Unique cell identifier/barcode used for joins across the analysis.",
    "Sample identifier.",
    "Standardized experimental timepoint.",
    "Broad atlas cell-type/lineage label reconstructed for Step8A.",

    "TRUE if the cell belongs to the final Step2C myeloid analysis set.",
    "Final Step2 myeloid cluster identifier when available.",
    "Final Step2C myeloid lineage annotation.",

    "TRUE if the cell belongs to the Step3B clean resident microglia set.",
    "Resident-state cluster used for the revised-manuscript state display.",
    "Final resident-state label used in the revised manuscript.",

    "Existing trajectory/transition status from the frozen trajectory audit.",
    "Existing high-confidence 3-dpi shared-transition flag.",

    "TRUE if the cell belongs to the 14,005-cell trajectory set.",
    "Final single-root pseudotime from the frozen trajectory audit.",
    "Principal-graph vertex from the frozen trajectory audit.",
    "Raw branch-arm label retained from the trajectory audit.",
    "Raw branch-component label retained from the trajectory audit.",

    "TRUE if the cell belongs to the final 4,203-cell branch analysis set.",
    "Raw frozen formal-branch state from Step7A, when available.",
    "Reader-facing branch label: Shared_prebranch, Alternative, or Target.",
    "TRUE only for formal branch cells assigned to the Target branch.",

    "Frozen Step8A CellChat grouping, when available."
  ),
  source = c(
    "Step8A atlas metadata",
    "Step8A atlas metadata",
    "Step8A atlas metadata",
    "Step8A atlas metadata",

    "membership in Step2C metadata",
    "Step2C metadata",
    "Step2C metadata",

    "membership in Step3B metadata",
    "Step3B metadata",
    "Step3B metadata or frozen Figure 1 cluster-to-state mapping",

    "Step5BR2 Cell_results",
    "Step5BR2 Cell_results",

    "membership in Step5BR2 Cell_results",
    "Step5BR2 Cell_results",
    "Step5BR2 Cell_results",
    "Step5BR2 Cell_results",
    "Step5BR2 Cell_results",

    "membership in Step7A frozen 4,203-cell table",
    "Step7A frozen 4,203-cell table",
    "derived only from frozen branch-arm/state labels",
    "branch_label == Target",

    "Step8A atlas metadata"
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  field_dictionary,
  file.path(OUTDIR, "01_annotation_field_dictionary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# --------------------------- 16. Export audit ---------------------------
audit <- data.frame(
  check = c(
    "atlas_cells",
    "step2_myeloid_cells",
    "clean_resident_microglia_cells",
    "trajectory_cells",
    "formal_branch_cells",
    "unmatched_step2_cells",
    "unmatched_step3_cells",
    "unmatched_trajectory_cells",
    "unmatched_formal_branch_cells",
    "clean_resident_cells_missing_resident_state",
    "formal_branch_cells_missing_public_branch_label"
  ),
  value = c(
    nrow(master_public),
    sum(master$myeloid_included %in% TRUE, na.rm = TRUE),
    sum(master$clean_resident_microglia %in% TRUE, na.rm = TRUE),
    sum(master$trajectory_included %in% TRUE, na.rm = TRUE),
    sum(master$formal_branch_included %in% TRUE, na.rm = TRUE),
    length(unmatched_step2),
    length(unmatched_step3),
    length(unmatched_traj),
    length(unmatched_branch),
    n_clean_missing_state,
    sum(
      master$formal_branch_included %in% TRUE &
        is.na(master$branch_label)
    )
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  audit,
  file.path(OUTDIR, "02_export_audit.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# --------------------------- 17. Value dictionary ---------------------------
value_dictionary <- data.frame(
  field = c(
    "resident_state",
    "resident_state",
    "resident_state",
    "resident_state",
    "branch_label",
    "branch_label",
    "branch_label"
  ),
  value = c(
    "Canonical resident",
    "Activated resident",
    "IFN-responsive resident",
    "Cycling resident",
    "Shared_prebranch",
    "Alternative",
    "Target"
  ),
  meaning = c(
    "Resident microglial state displayed in the revised manuscript.",
    "Resident microglial state displayed in the revised manuscript.",
    "Resident microglial state displayed in the revised manuscript.",
    "Resident microglial state displayed in the revised manuscript.",
    "Formal shared/prebranch portion of the frozen branch analysis.",
    "Formal alternative branch (frozen Arm_1 correspondence).",
    "Formal target-program branch (frozen Arm_2 correspondence)."
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  value_dictionary,
  file.path(OUTDIR, "03_value_dictionary_branch_and_state.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# --------------------------- 18. README ---------------------------
readme <- c(
  "GSE162610 public cell-level annotation export",
  paste0("Created: ", Sys.time()),
  "",
  "Purpose:",
  "This folder contains the cell-level annotation table requested for public",
  "reproducibility of the revised GLIA manuscript.",
  "",
  "IMPORTANT RDS handling:",
  "Exactly two RDS files were read solely to export existing Seurat meta.data:",
  paste0("1. ", basename(RDS_STEP2)),
  paste0("2. ", basename(RDS_STEP3)),
  "readRDS() necessarily deserializes the full object, including assays, into memory.",
  "No expression values were accessed and no statistical or bioinformatic analysis",
  "was rerun. Each object was removed immediately after meta.data extraction.",
  "",
  "No annotation was re-inferred.",
  "The export assembles already-frozen metadata and preserves raw trajectory/branch",
  "labels alongside reader-facing labels.",
  "",
  "Expected frozen memberships:",
  paste0("- Atlas cells: ", EXPECTED$atlas_cells),
  paste0("- Final Step2 myeloid cells: ", EXPECTED$step2_myeloid),
  paste0("- Clean resident microglia: ", EXPECTED$clean_microglia),
  paste0("- Trajectory cells: ", EXPECTED$trajectory_cells),
  paste0("- Formal branch cells: ", EXPECTED$formal_branch),
  "",
  paste0(
    "Clean resident cells without a public resident_state after export: ",
    n_clean_missing_state
  ),
  "",
  "Primary public file:",
  "GSE162610_cell_annotations.csv.gz",
  "",
  "Supporting provenance files:",
  "- 00_source_file_manifest.csv",
  "- 01_annotation_field_dictionary.csv",
  "- 02_export_audit.csv",
  "- 03_value_dictionary_branch_and_state.csv",
  "- GSE162610_cell_annotations_preview100.csv"
)

writeLines(
  readme,
  file.path(OUTDIR, "README.txt"),
  useBytes = TRUE
)

# --------------------------- 19. Final report ---------------------------
cat("\n============================================================\n")
cat("PUBLIC CELL-ANNOTATION EXPORT COMPLETE\n")
cat("============================================================\n")
cat("Output directory:\n", normalize_slashes(OUTDIR), "\n\n", sep = "")
cat("Master annotation:\n", normalize_slashes(out_gz), "\n\n", sep = "")
cat("Rows: ", nrow(master_public), "\n", sep = "")
cat("Columns: ", ncol(master_public), "\n", sep = "")
cat("Clean resident cells missing resident_state: ",
    n_clean_missing_state, "\n", sep = "")
cat("Formal branch cells missing public branch_label: ",
    sum(master$formal_branch_included %in% TRUE & is.na(master$branch_label)),
    "\n", sep = "")
cat("============================================================\n")
