options(stringsAsFactors = FALSE, warn = 1)
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
STEP <- "step11C_GSE234774_author_myeloid_target_program_comparison"
OUT_DIR <- file.path(PROJECT_ROOT, "04_results", STEP)
TABLE_DIR <- file.path(OUT_DIR, "tables")
DESKTOP <- dirname(PROJECT_ROOT)
TRANSFER_ZIP <- file.path(DESKTOP, paste0(STEP, ".zip"))
ZIP_LIMIT <- 512 * 1024^2
INPUT <- c( atlas = file.path(PROJECT_ROOT, "04_results", "step10D_GSE234774_marker_first_major_lineage_annotation", "RDS", "GSE234774_step10D_major_lineage_annotation_evidence_BPCells.rds"), full_myeloid = file.path(PROJECT_ROOT, "04_results", "step10G_GSE234774_frozen_myeloid_annotations_resident_objects", "RDS", "GSE234774_step10G_full_myeloid_final_annotations_BPCells.rds"), resident_scored = file.path(PROJECT_ROOT, "04_results", "step11B_GSE234774_refined_main_Arm_specific_shared_AUCell_ssGSEA", "RDS", "GSE234774_step11B_refined_main_Arm_scored_BPCells.rds"), frozen_gene_sets = file.path(PROJECT_ROOT, "04_results", "step11B_GSE234774_refined_main_Arm_specific_shared_AUCell_ssGSEA", "tables", "04_refined_Arm_specific_shared_signature_audit.csv") )
if (!dir.exists(PROJECT_ROOT)) stop("Project root not found: ", PROJECT_ROOT)
missing <- INPUT[!file.exists(INPUT)]
if (length(missing)) stop("Missing ORIGINAL inputs:\n", paste(missing, collapse = "\n"))
if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  stop("SeuratObject is required; install it in the existing project R environment.")
}
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
logfile <- file.path(OUT_DIR, "00_run_log.txt")
log <- function(...) {
  msg <- paste0(..., collapse = "")
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", msg, "\n", sep = "")
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", msg, "\n", file = logfile, append = TRUE, sep = "")
}
log("BEGIN ", STEP)
input_audit <- data.frame( role = names(INPUT), absolute_path = unname(INPUT), size_MB = round(as.numeric(file.info(INPUT)$size) / 1024^2, 3), modified_at = as.character(file.info(INPUT)$mtime), stringsAsFactors = FALSE )
write.csv(input_audit, file.path(TABLE_DIR, "01_original_input_file_audit.csv"), row.names = FALSE)
write_table <- function(x, filename) {
  destination <- file.path(TABLE_DIR, filename)
  if (grepl("\\.gz$", destination)) {
    con <- gzfile(destination, open = "wt", compression = 6)
    on.exit(close(con), add = TRUE)
    utils::write.csv(x, file = con, row.names = FALSE, na = "")
  } else {
    utils::write.csv(x, file = destination, row.names = FALSE, na = "")
  }
  log("OUTPUT ", filename, " rows=", nrow(x), " size_MB=", round(file.info(destination)$size / 1024^2, 3))
  invisible(destination)
}
meta_from_seurat <- function(object, role) {
  if (!inherits(object, "Seurat")) stop(role, " is not a Seurat object.")
  m <- object@meta.data
  cells <- colnames(object)
  if (anyDuplicated(cells) || anyDuplicated(rownames(m)) || !all(cells %in% rownames(m))) {
    stop("Invalid or unmatched barcodes in ", role)
  }
  m[cells, , drop = FALSE]
}
require_columns <- function(m, required, role) {
  lacking <- setdiff(required, colnames(m))
  if (length(lacking)) {
    stop("Missing required metadata columns in ", role, ": ", paste(lacking, collapse = ", "), ". Do not silently substitute another annotation field.")
  }
}
make_metadata <- function(m, fields) {
  out <- data.frame(cell_barcode = rownames(m), stringsAsFactors = FALSE)
  for (field in fields) {
    out[[field]] <- if (field %in% colnames(m)) as.character(m[[field]])
                    else rep(NA_character_, nrow(m))
  }
  out
}
summarise_scores <- function(df, by, scores) {
  require_columns(df, c(by, scores), "score-summary table")
  if (!nrow(df)) return(df[FALSE, c(by, scores), drop = FALSE])
  keys <- do.call(paste, c(lapply(df[by], function(z) {
    z <- as.character(z)
    z[is.na(z)] <- "<missing>"
    z
  }), sep = "\034"))
  groups <- split(seq_len(nrow(df)), keys, drop = TRUE)
  result <- lapply(groups, function(index) {
    row <- df[index[1], by, drop = FALSE]
    row$n_cells <- length(index)
    for (score in scores) {
      values <- suppressWarnings(as.numeric(df[[score]][index]))
      values <- values[is.finite(values)]
      row[[paste0(score, "_median")]] <- if (length(values)) median(values) else NA_real_
      row[[paste0(score, "_mean")]] <- if (length(values)) mean(values) else NA_real_
      row[[paste0(score, "_n_scored")]] <- length(values)
    }
    row
  })
  result <- do.call(rbind, result)
  rownames(result) <- NULL
  result$descriptive_n_ge_10 <- result$n_cells >= 10
  result
}
get_counts <- function(obj, role) {
  data <- tryCatch( SeuratObject::GetAssayData(obj, assay = "RNA", layer = "counts"), error = function(e) stop("Could not access RNA counts of ", role, ": ", conditionMessage(e)) )
  if (!all(colnames(obj) %in% colnames(data))) {
    stop(role, ": count matrix does not contain all object barcodes.")
  }
  data
}
score_programs <- function(counts, cell_ids, depth, signatures, role, batch = 300L) {
  if (length(cell_ids) != length(depth) || anyNA(depth) || any(!is.finite(depth)) || any(depth <= 0)) {
    stop(role, ": depth/barcode validation failed.")
  }
  if (anyDuplicated(cell_ids) || !all(cell_ids %in% colnames(counts))) {
    stop(role, ": invalid count-matrix barcodes.")
  }
  outputs <- vector("list", length(signatures))
  names(outputs) <- names(signatures)
  for (program in names(signatures)) {
    genes <- intersect(signatures[[program]], rownames(counts))
    if (length(genes) < 10L) stop(role, ": too few genes for ", program)
    value <- numeric(length(cell_ids))
    for (first in seq.int(1L, length(cell_ids), by = batch)) {
      last <- min(first + batch - 1L, length(cell_ids))
      jj <- first:last
      block <- as.matrix(counts[genes, cell_ids[jj], drop = FALSE])
      value[jj] <- colMeans(log1p(sweep(block, 2L, depth[jj], "/") * 1e4))
    }
    outputs[[program]] <- value
    log(role, ": scored ", program, ", cells=", length(value), ", genes=", length(genes))
  }
  outputs
}
audit <- read.csv(INPUT[["frozen_gene_sets"]], check.names = FALSE)
flag_cols <- c("in_Arm1_specific", "in_Arm2_specific", "in_Shared_prebranch")
require_columns(audit, c("record_type", "resolved_external_feature", "present_in_GSE234774", flag_cols), "original Step11B frozen gene audit")
program_definitions <- c( Arm_2_specific = "in_Arm2_specific", Arm_1_specific = "in_Arm1_specific", Shared_prebranch = "in_Shared_prebranch" )
signatures <- lapply(program_definitions, function(field) {
  ok <- audit$record_type == "GENE" & !is.na(audit[[field]]) & (audit[[field]] %in% TRUE) & !is.na(audit$present_in_GSE234774) & (audit$present_in_GSE234774 %in% TRUE)
  unique(as.character(audit$resolved_external_feature[ok & !is.na(audit$resolved_external_feature) & nzchar(audit$resolved_external_feature)]))
})
if (any(lengths(signatures) < 10L)) stop("Failed to resolve original Step11B signatures.")
for (nm in names(signatures)) {
  log("ORIGINAL SIGNATURE ", nm, " n_genes=", length(signatures[[nm]]))
}
log("READ frozen Step10G full-myeloid object")
full <- readRDS(INPUT[["full_myeloid"]])
fm <- meta_from_seurat(full, "Step10G full-myeloid")
require_columns(fm, c("layer5", "biological_library", "timepoint", "myeloid_final_label", "nCount_RNA"), "Step10G full-myeloid")
if (nrow(fm) != 25205L) {
  log("WARNING: Step10G full-myeloid size differs from reported 25,205: ", nrow(fm))
}
full_meta <- make_metadata(fm, c("biological_library", "timepoint", "layer2", "layer3", "layer4", "layer5", "myeloid_final_label", "resident_lineage_status", "myhcon_res_0_6", "nCount_RNA", "nFeature_RNA"))
cross_tab <- as.data.frame( table(original_layer5 = full_meta$layer5, project_myeloid_annotation = full_meta$myeloid_final_label, useNA = "ifany"), stringsAsFactors = FALSE)
cross_tab <- cross_tab[cross_tab$Freq > 0, , drop = FALSE]
names(cross_tab)[3] <- "n_cells"
write_table(cross_tab, "02_original_author_vs_project_myeloid_annotation.csv")
full_counts <- get_counts(full, "Step10G")
for (nm in names(signatures)) {
  if (sum(signatures[[nm]] %in% rownames(full_counts)) < 10L) {
    stop("Step10G original gene availability insufficient for ", nm)
  }
}
scoring_audit <- data.frame( program = names(signatures), role = c("Target inflammatory-remodeling", "Alternative", "Shared prebranch"), source_n_genes = lengths(signatures), present_in_Step10G = vapply(signatures, function(x) sum(x %in% rownames(full_counts)), integer(1)), score_method = "mean(log1p(UMI / total_cell_UMI * 10000))", stringsAsFactors = FALSE)
write_table(scoring_audit, "03_source_signature_and_scoring_audit.csv")
full_scores <- score_programs(full_counts, full_meta$cell_barcode, suppressWarnings(as.numeric(fm$nCount_RNA)), signatures, "Step10G full myeloid")
score_names <- paste0(names(signatures), "_mean_logCP10K")
for (i in seq_along(signatures)) full_meta[[score_names[[i]]]] <- full_scores[[i]]
write_table(full_meta, "04_full_myeloid_cell_scores_original_and_project_labels.csv.gz")
write_table(summarise_scores(full_meta, c("timepoint", "biological_library", "layer5"), score_names), "05_full_myeloid_by_library_original_author_label.csv")
write_table(summarise_scores(full_meta, c("timepoint", "biological_library", "myeloid_final_label"), score_names), "06_full_myeloid_by_library_project_annotation.csv")
write_table(summarise_scores(full_meta, c("timepoint", "biological_library", "layer5", "myeloid_final_label"), score_names), "07_full_myeloid_by_library_original_and_project_labels.csv")
rm(full_counts, full, fm, full_scores, full_meta)
invisible(gc())
log("READ Step10D quality-controlled atlas to evaluate author myeloid subsets")
atlas <- readRDS(INPUT[["atlas"]])
am <- meta_from_seurat(atlas, "Step10D atlas")
require_columns(am, c("layer5", "timepoint", "biological_library", "nCount_RNA"), "Step10D atlas")
myeloid_labels <- grepl("microglia|macrophage|myeloid|monocyte|neutrophil", as.character(am$layer5), ignore.case = TRUE)
myeloid_labels[is.na(myeloid_labels)] <- FALSE
indices <- which(myeloid_labels)
if (length(indices) < 100L) {
  stop("Too few published myeloid labels remain in the Step10D QC atlas.")
}
selected_ids <- rownames(am)[indices]
ref_meta <- make_metadata(am[indices, , drop = FALSE], c("biological_library", "timepoint", "layer2", "layer3", "layer4", "layer5", "nCount_RNA", "nFeature_RNA"))
ref_counts <- get_counts(atlas, "Step10D")
ref_scores <- score_programs(ref_counts, selected_ids, suppressWarnings(as.numeric(am$nCount_RNA[indices])), signatures, "Step10D author-annotated myeloid")
for (i in seq_along(signatures)) ref_meta[[score_names[[i]]]] <- ref_scores[[i]]
write_table(ref_meta, "08_atlas_author_myeloid_cell_scores.csv.gz")
write_table(summarise_scores(ref_meta, c("timepoint", "biological_library", "layer5"), score_names), "09_atlas_author_myeloid_by_library.csv")
write_table(summarise_scores(ref_meta, c("timepoint", "layer5"), score_names), "10_atlas_author_myeloid_by_timepoint_descriptive.csv")
n_atlas <- nrow(am)
n_author_myeloid <- length(indices)
n_author_labels <- length(unique(as.character(ref_meta$layer5)))
log("Step10D QC atlas nuclei=", n_atlas, "; nuclei matching author myeloid labels=", n_author_myeloid, "; matched published layer5 classes=", n_author_labels)
rm(atlas, am, ref_counts, ref_scores, ref_meta)
invisible(gc())
log("READ original Step11B resident-only AUCell-scored object")
resident <- readRDS(INPUT[["resident_scored"]])
res_meta <- meta_from_seurat(resident, "Step11B scored resident")
auc_cols <- c("Arm_2_specific_AUCell", "Arm_1_specific_AUCell", "Shared_prebranch_AUCell")
require_columns(res_meta, c("biological_library", "timepoint", "layer5", "myeloid_final_label", auc_cols), "Step11B scored resident")
auc_meta <- make_metadata(res_meta, c("biological_library", "timepoint", "layer5", "myeloid_final_label", "nCount_RNA", "nFeature_RNA"))
for (field in auc_cols) auc_meta[[field]] <- suppressWarnings(as.numeric(res_meta[[field]]))
write_table(auc_meta, "11_existing_Step11B_resident_AUCell_with_author_labels.csv.gz")
write_table(summarise_scores(auc_meta, c("timepoint", "biological_library", "layer5"), auc_cols), "12_existing_Step11B_resident_AUCell_by_library_author_label.csv")
rm(resident, res_meta, auc_meta)
invisible(gc())
report <- c( "Step 11C | Reviewer 2 Minor Comment 1: published Tabulae Paralytica label comparisons", "", paste("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "Inputs: ORIGINAL Step10D, Step10G, Step11B artifacts in xjz/04_results/.", "No ZIP or intermediate evidence package was used as an input.", "", "Data scope:", "  Step10G: 25,205 candidate myeloid nuclei (or audited current object size),", "    with independently reviewed project labels and published layer5 labels.", paste0("  Step10D: ", n_author_myeloid, " nuclei with published myeloid labels, among ", n_atlas, " nuclei in the QC-processed atlas object."), "  Step11B: precomputed resident-like microglia AUCell scores only.", "", "Target program = Arm_2_specific (inflammatory-remodeling).", "Alternative program = Arm_1_specific.", "Shared program = Shared_prebranch.", "All Step10D/Step10G groups are scored by the same mean-logCP10K method.", "These exploratory cross-lineage expression summaries are NOT AUCell and", "should not be compared numerically with the original Step11B AUCell values.", "Score differences across author labels may reflect cell composition,", "RNA depth, and biological library effects; no lineage origin is established.", "The Step10D QC atlas is NOT the entire original 435,099-nucleus metadata.", "Per-library summaries are supplied; cells are not independent replicates.", "No differential tests or novel subtype identifications were performed.", "", "Outputs:", "  tables/01  original input audit", "  tables/02  Tabulae Paralytica original labels vs project myeloid annotations", "  tables/03  source signature definition and score-method audit", "  tables/04-07  full project myeloid comparison (cell / library levels)", "  tables/08-10  published myeloid groups in QC atlas (cell / library levels)", "  tables/11-12  original resident-only AUCell metadata and summary", "", paste("Results directory:", OUT_DIR), "Earlier objects and tables were not modified." )
writeLines(report, file.path(OUT_DIR, "10_STEP11C_FINAL_REPORT.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "00_sessionInfo.txt"))
log("SAVED 10_STEP11C_FINAL_REPORT.txt and 00_sessionInfo.txt")
files_for_zip <- list.files(OUT_DIR, recursive = TRUE, all.files = FALSE, full.names = FALSE, include.dirs = FALSE)
files_for_zip <- files_for_zip[!grepl("~$|\\.lock$", files_for_zip)]
if (file.exists(TRANSFER_ZIP)) file.remove(TRANSFER_ZIP)
zip_ok <- FALSE
if (requireNamespace("zip", quietly = TRUE)) {
  tryCatch({
    zip::zipr(zipfile = TRANSFER_ZIP, files = files_for_zip, root = OUT_DIR, include_directories = FALSE)
    zip_ok <- file.exists(TRANSFER_ZIP)
  }, error = function(e) log("TRANSFER ZIP warning: ", conditionMessage(e)))
} else {
  log("'zip' package not installed. Results were saved under xjz/04_results/.")
  log("For a transfer ZIP, install.packages('zip'), then rerun this script.")
}
if (zip_ok) {
  zip_size <- file.info(TRANSFER_ZIP)$size
  if (is.na(zip_size) || zip_size >= ZIP_LIMIT) {
    file.remove(TRANSFER_ZIP)
    log("TRANSFER ZIP EXCEEDED 512 MiB; ZIP removed; main results retained.")
  } else {
    log("TRANSFER ZIP: ", TRANSFER_ZIP, " (", round(zip_size / 1024^2, 2), " MiB)")
  }
}
log("SUCCESS: independent main results directory ", OUT_DIR)
