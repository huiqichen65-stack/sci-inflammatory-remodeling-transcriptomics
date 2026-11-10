options(stringsAsFactors = FALSE, warn = 1)
if (!requireNamespace("Matrix", quietly = TRUE)) stop("Missing package: Matrix")
if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("Missing package: SeuratObject")
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
if (!dir.exists(PROJECT_ROOT)) {
  guesses <- c(file.path(Sys.getenv("USERPROFILE"), "Desktop", "xjz"), file.path(Sys.getenv("OneDrive"), "Desktop", "xjz"), file.path(path.expand("~"), "Desktop", "xjz"))
  guesses <- guesses[nzchar(guesses) & dir.exists(guesses)]
  if (length(guesses)) PROJECT_ROOT <- guesses[1L]
}
if (!dir.exists(PROJECT_ROOT)) stop("Cannot find xjz project root. Edit PROJECT_ROOT: ", PROJECT_ROOT)
PROJECT_ROOT <- normalizePath(PROJECT_ROOT, winslash = "/", mustWork = TRUE)
STEP <- "step2D_GSE162610_target_program_lineage_specificity_audit"
OUT <- file.path(PROJECT_ROOT, "04_results", STEP)
TBL <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")
dir.create(TBL, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
rel_path <- function(...) file.path(PROJECT_ROOT, ...)
inputs <- c( atlas = rel_path("04_results", "step1F_GSE162610_broad_annotation", "GSE162610_broad_annotation_audit.rds"), myeloid = rel_path("04_results", "step2C_GSE162610_final_microglia_extraction", "GSE162610_step2_full_myeloid_with_final_lineage.rds"), cleanup = rel_path("04_results", "step3B_GSE162610_state_scoring", "GSE162610_step3B_full_with_cleanup_flags.rds"), branch = rel_path("04_results", "step5BR2_GSE162610_single_root_branch_audit", "GSE162610_step5BR2_single_root_seurat.rds"), signature = rel_path("04_results", "step11B_GSE234774_refined_main_Arm_specific_shared_AUCell_ssGSEA", "tables", "04_refined_Arm_specific_shared_signature_audit.csv") )
missing_inputs <- names(inputs)[!file.exists(inputs)]
if (length(missing_inputs)) {
  stop("Required file(s) missing in xjz: ", paste(missing_inputs, collapse = ", "), "\n", paste(inputs[missing_inputs], collapse = "\n"))
}
log_text <- c(paste("STEP2D audit launched:", Sys.time()), paste("Project:", PROJECT_ROOT), paste("Results:", OUT), "INPUTS are original xjz files; no ZIP files used or generated.", paste(names(inputs), "=", inputs, collapse = "\n"))
logmsg <- function(...) {
  x <- paste0(..., collapse = "")
  log_text <<- c(log_text, x)
}
sig_audit <- read.csv(inputs[["signature"]], check.names = FALSE)
need_cols <- c("record_type", "gene", "in_Arm2_specific")
if (!all(need_cols %in% names(sig_audit))) {
  stop("Frozen-signature audit lacks expected columns: ", paste(setdiff(need_cols, names(sig_audit)), collapse = ", "))
}
is_selected <- as.character(sig_audit$record_type) == "GENE" & tolower(as.character(sig_audit$in_Arm2_specific)) %in% c("true", "t", "1")
signature <- unique(as.character(sig_audit$gene[is_selected]))
signature <- signature[!is.na(signature) & nzchar(signature)]
if (length(signature) < 400L) stop("Unexpected frozen Arm2-specific signature length: ", length(signature))
logmsg("Frozen Arm2-specific signature genes: ", length(signature))
reviewer_markers <- c( "Tnf", "C5ar1", "C3ar1", "Myd88", "Icam1", "Fos", "Egr1", "Junb", "Atf3", "Nfkbiz", "Ccl3", "Ccl4", "Dnajb1", "Hexb", "Sall1", "Tmem119", "P2ry12", "Fcrls", "Cx3cr1", "Siglech", "Ms4a7", "Ccr2", "Ly6c2", "Fn1", "Fcgr1", "S100a8", "S100a9", "Lyve1", "Cd163", "Mrc1", "Pf4" )
requested_genes <- unique(c(signature, reviewer_markers))
logmsg("Reading original Step1F whole-atlas Seurat RDS (may use substantial RAM)...")
atlas_obj <- readRDS(inputs[["atlas"]])
if (!isS4(atlas_obj) || !"meta.data" %in% methods::slotNames(atlas_obj))
  stop("Step1F RDS is not a Seurat metadata-bearing object")
md <- atlas_obj@meta.data
required_atlas_fields <- c("sample_id", "timepoint", "nCount_RNA", "author_annotation")
if (!all(required_atlas_fields %in% names(md)))
  stop("Step1F metadata lacks: ", paste(setdiff(required_atlas_fields, names(md)), collapse = ", "))
idx_7dpi <- !is.na(md$timepoint) & as.character(md$timepoint) == "7 dpi"
if (!any(idx_7dpi)) stop("No 7 dpi cells found in Step1F atlas metadata")
atlas7 <- data.frame( cell_id = rownames(md)[idx_7dpi], sample_id = as.character(md$sample_id[idx_7dpi]), timepoint = as.character(md$timepoint[idx_7dpi]), nCount_RNA = suppressWarnings(as.numeric(md$nCount_RNA[idx_7dpi])), nFeature_RNA = if ("nFeature_RNA" %in% names(md)) suppressWarnings(as.numeric(md$nFeature_RNA[idx_7dpi])) else NA_real_, author_annotation = as.character(md$author_annotation[idx_7dpi]), stringsAsFactors = FALSE )
if (anyDuplicated(atlas7$cell_id)) stop("Step1F cell barcode duplication detected")
if (anyNA(atlas7$sample_id)) stop("Missing sample IDs in 7dpi atlas")
logmsg("Step1F 7dpi cells: ", nrow(atlas7), "; samples: ", paste(unique(atlas7$sample_id), collapse = ", "))
assay_names <- SeuratObject::Assays(atlas_obj)
assay <- if ("RNA" %in% assay_names) "RNA" else assay_names[1]
if (is.na(assay) || !nzchar(assay)) stop("No assay available in Step1F object")
raw <- tryCatch(SeuratObject::GetAssayData(atlas_obj, assay = assay, layer = "counts"), error = function(e) NULL)
if (is.null(raw)) raw <- tryCatch( SeuratObject::GetAssayData(atlas_obj, assay = assay, slot = "counts"), error = function(e) NULL)
if (is.null(raw)) stop("Could not read Step1F RNA counts. Check Seurat assay/layers.")
if (is.null(rownames(raw)) || is.null(colnames(raw))) stop("Step1F counts missing gene/cell dimnames")
present <- intersect(requested_genes, rownames(raw))
sig_present <- intersect(signature, present)
if (length(sig_present) < 400L) {
  stop("Less than 400 frozen signature genes found in Step1F counts: ", length(sig_present), "/", length(signature), ". Check original expression feature identifiers.")
}
count_ids <- intersect(atlas7$cell_id, colnames(raw))
if (length(count_ids) != nrow(atlas7)) {
  stop("7 dpi atlas cell barcodes missing from Step1F raw counts: ", nrow(atlas7) - length(count_ids))
}
raw7 <- raw[present, count_ids, drop = FALSE]
if (!inherits(raw7, "sparseMatrix")) raw7 <- Matrix::Matrix(raw7, sparse = TRUE)
raw7 <- methods::as(raw7, "dgCMatrix")
if (any(!is.finite(raw7@x))) stop("Original raw counts include NA/Inf")
rm(raw, atlas_obj, md)
invisible(gc(verbose = FALSE))
logmsg("Targeted raw counts: ", nrow(raw7), " genes x ", ncol(raw7), " cells; RNA assay: ", assay)
write.csv(data.frame(gene = signature, in_GSE162610_Step1F_RNA_counts = signature %in% rownames(raw7)), file.path(TBL, "01_frozen_Arm2_signature_gene_coverage.csv"), row.names = FALSE)
missing_markers <- setdiff(reviewer_markers, rownames(raw7))
writeLines(missing_markers, file.path(TBL, "04_marker_genes_not_in_restricted_matrix.txt"))
read_rds_meta <- function(file, cols, tag) {
  logmsg("Reading original ", tag, " RDS metadata...")
  object <- readRDS(file)
  if (!isS4(object) || !"meta.data" %in% methods::slotNames(object))
    stop(tag, " object has no @meta.data")
  d <- object@meta.data
  missing_cols <- setdiff(cols, names(d))
  if (length(missing_cols)) stop(tag, " missing metadata fields: ", paste(missing_cols, collapse = ", "))
  small <- as.data.frame(d[, cols, drop = FALSE], check.names = FALSE)
  small$cell_id <- rownames(d)
  if (is.null(small$cell_id) || anyDuplicated(small$cell_id))
    stop(tag, " has invalid/duplicated cell barcodes")
  rm(object, d)
  invisible(gc(verbose = FALSE))
  small
}
myeloid <- read_rds_meta(inputs[["myeloid"]], c("final_step2_lineage", "final_step2_status", "myeloid_cluster_r04", "AUCell_Resident_microglia", "AUCell_Infiltrating_mono_macro", "AUCell_BAM_like"), "Step2C")
cleanup <- read_rds_meta(inputs[["cleanup"]], c("step3B_cleanup_reason"), "Step3B")
branch <- read_rds_meta(inputs[["branch"]], c("step5BR2_branch_arm", "myeloid_cluster_r04"), "Step5BR2")
merge_field <- function(src, col) {
  as.character(src[[col]][match(atlas7$cell_id, src$cell_id)])
}
meta <- atlas7
meta$final_step2_lineage <- merge_field(myeloid, "final_step2_lineage")
meta$final_step2_status <- merge_field(myeloid, "final_step2_status")
meta$myeloid_cluster_r04 <- merge_field(myeloid, "myeloid_cluster_r04")
meta$step3B_cleanup_reason <- merge_field(cleanup, "step3B_cleanup_reason")
meta$frozen_branch <- merge_field(branch, "step5BR2_branch_arm")
for (cl in c("AUCell_Resident_microglia", "AUCell_Infiltrating_mono_macro", "AUCell_BAM_like")) {
  meta[[cl]] <- suppressWarnings(as.numeric(myeloid[[cl]][match(meta$cell_id, myeloid$cell_id)]))
}
rm(myeloid, cleanup, branch, atlas7)
invisible(gc(verbose = FALSE))
meta$group <- NA_character_
meta$group[!is.na(meta$final_step2_lineage) & meta$final_step2_lineage == "BAM_like_macrophage"] <- "Excluded_BAM_like"
infiltrating <- c("infiltrating_macrophage", "infiltrating_monocyte", "macrophage_biased_infiltrating_myeloid", "IFN_responsive_infiltrating_mono_macro", "MHCII_APC_like_infiltrating_macrophage")
meta$group[!is.na(meta$final_step2_lineage) & meta$final_step2_lineage %in% infiltrating] <- "Excluded_infiltrating_myeloid"
meta$group[!is.na(meta$step3B_cleanup_reason) & meta$step3B_cleanup_reason == "residual_monocyte_like_contamination"] <- "Removed_residual_monocyte_like"
meta$group[!is.na(meta$frozen_branch) & meta$frozen_branch == "Arm_1"] <- "Arm_1_alternative"
meta$group[!is.na(meta$frozen_branch) & meta$frozen_branch == "Arm_2"] <- "Arm_2_target"
selected <- meta[!is.na(meta$group), , drop = FALSE]
neu <- meta[!is.na(meta$author_annotation) & tolower(trimws(meta$author_annotation)) == "neutrophil", , drop = FALSE]
neu$group <- "Neutrophil_author_annotation"
selected <- rbind(selected, neu)
selected <- selected[selected$cell_id %in% colnames(raw7), , drop = FALSE]
if (!nrow(selected)) stop("No lineage groups selected; check group names")
levels_groups <- c("Arm_2_target", "Arm_1_alternative", "Excluded_infiltrating_myeloid", "Excluded_BAM_like", "Removed_residual_monocyte_like", "Neutrophil_author_annotation")
selected$group <- factor(selected$group, levels = levels_groups)
if (anyNA(selected$group) || anyNA(selected$sample_id)) stop("NA group/sample in selected cells")
count_tab <- as.data.frame(table(selected$sample_id, selected$group), stringsAsFactors = FALSE)
names(count_tab) <- c("sample_id", "group", "n_cells")
count_tab <- count_tab[count_tab$n_cells > 0, , drop = FALSE]
write.csv(count_tab, file.path(TBL, "02_sample_by_lineage_group_counts.csv"), row.names = FALSE)
neu_overlap <- meta[!is.na(meta$author_annotation) & tolower(trimws(meta$author_annotation)) == "neutrophil" & !is.na(meta$group), c("cell_id", "sample_id", "author_annotation", "final_step2_lineage", "frozen_branch", "group"), drop = FALSE]
write.csv(neu_overlap, file.path(TBL, "08_author_neutrophil_overlap_audit.csv"), row.names = FALSE)
logmsg("Author-Neutrophil cells overlapping a main lineage/branch group: ", nrow(neu_overlap))
logmsg("NOTE: neutrophil comparator may overlap Arm groups; never sum those rows to total.")
ids <- unique(as.character(selected$cell_id))
if (anyNA(ids) || anyDuplicated(ids)) stop("Invalid selected cell IDs")
counts_sub <- raw7[, ids, drop = FALSE]
if (!identical(colnames(counts_sub), ids)) stop("Raw cell IDs not aligned")
size_factors <- meta$nCount_RNA[match(ids, meta$cell_id)]
if (length(size_factors) != length(ids) || anyNA(size_factors) || any(!is.finite(size_factors)) || any(size_factors <= 0))
  stop("Invalid total-UMI size factors; inspect Step1F nCount_RNA")
logmat <- counts_sub %*% Matrix::Diagonal(x = 10000 / size_factors)
dimnames(logmat) <- dimnames(counts_sub)
logmat <- methods::as(logmat, "dgCMatrix")
logmat@x <- log1p(logmat@x)
if (!identical(colnames(logmat), ids)) stop("Normalized cell barcodes were lost")
sig_score <- as.numeric(Matrix::colMeans(logmat[sig_present, , drop = FALSE]))
if (length(sig_score) != length(ids) || any(!is.finite(sig_score)))
  stop("Nonfinite frozen signature score; check normalization and gene mapping")
score_match <- match(as.character(selected$cell_id), ids)
if (anyNA(score_match)) stop("Selected cells match normalized counts")
selected$target_signature_mean_logCP10k <- sig_score[score_match]
if (anyNA(selected$target_signature_mean_logCP10k))
  stop("Missing scores after barcode matching (the previous script's error)")
split_group <- split(selected, interaction(selected$sample_id, selected$group, drop = TRUE))
if (!length(split_group)) stop("Empty sample-by-group splitting")
score_stats <- do.call(rbind, lapply(split_group, function(d) {
  x <- d$target_signature_mean_logCP10k
  if (!length(x) || any(!is.finite(x))) stop("Invalid group scores")
  data.frame(sample_id = as.character(d$sample_id[1]), group = as.character(d$group[1]), n_cells = nrow(d), mean_score = mean(x), median_score = median(x), q25 = as.numeric(stats::quantile(x, 0.25, na.rm = TRUE)), q75 = as.numeric(stats::quantile(x, 0.75, na.rm = TRUE)), median_nCount_RNA = stats::median(d$nCount_RNA), median_nFeature_RNA = stats::median(d$nFeature_RNA, na.rm = TRUE))
}))
write.csv(score_stats, file.path(TBL, "03_target_signature_by_sample_lineage_group.csv"), row.names = FALSE)
markers <- intersect(reviewer_markers, rownames(counts_sub))
if (!length(markers)) stop("No reviewer-requested marker genes found in GSE162610 expression matrix")
expr <- do.call(rbind, lapply(split_group, function(d) {
  ix <- match(d$cell_id, ids)
  if (anyNA(ix)) stop("Marker expression join failure")
  p <- Matrix::rowMeans(counts_sub[markers, ix, drop = FALSE] > 0)
  av <- Matrix::rowMeans(logmat[markers, ix, drop = FALSE])
  data.frame(sample_id = as.character(d$sample_id[1]), group = as.character(d$group[1]), gene = markers, n_cells = nrow(d), pct_detected = as.numeric(p), mean_logCP10k = as.numeric(av))
}))
write.csv(expr, file.path(TBL, "05_marker_expression_by_sample_lineage_group.csv"), row.names = FALSE)
arm <- meta[!is.na(meta$frozen_branch) & meta$frozen_branch %in% c("Arm_1", "Arm_2"), , drop = FALSE]
if (!nrow(arm)) stop("No frozen Arm1/Arm2 barcodes mapped into 7dpi atlas")
cluster_ancestry <- as.data.frame(table(arm$sample_id, arm$frozen_branch, arm$myeloid_cluster_r04), stringsAsFactors = FALSE)
names(cluster_ancestry) <- c("sample_id", "frozen_branch", "Step2_myeloid_cluster", "n_cells")
cluster_ancestry <- cluster_ancestry[cluster_ancestry$n_cells > 0, , drop = FALSE]
write.csv(cluster_ancestry, file.path(TBL, "06_Arm1_Arm2_step2_cluster_provenance.csv"), row.names = FALSE)
score_columns <- c("AUCell_Resident_microglia", "AUCell_Infiltrating_mono_macro", "AUCell_BAM_like")
auc_rows <- list()
for (d in split_group) {
  for (j in seq_along(score_columns)) {
    k <- score_columns[j]
    y <- d[[k]]
    y <- y[is.finite(y)]
    auc_rows[[length(auc_rows) + 1L]] <- data.frame( sample_id = as.character(d$sample_id[1]), group = as.character(d$group[1]), lineage_signature = k, n_cells = nrow(d), n_with_score = length(y), mean_AUCell = if (length(y)) mean(y) else NA_real_, median_AUCell = if (length(y)) median(y) else NA_real_)
  }
}
auc_out <- do.call(rbind, auc_rows)
write.csv(auc_out, file.path(TBL, "07_lineage_AUCell_by_sample_group.csv"), row.names = FALSE)
plot_groups <- function(one_sample) {
  d <- selected[selected$sample_id == one_sample, , drop = FALSE]
  boxplot(target_signature_mean_logCP10k ~ group, data = d, outline = FALSE, las = 2, ylab = "Mean log1p(CP10k) across frozen Arm2 genes", main = paste("7 dpi:", one_sample), cex.axis = .76)
}
samples <- sort(unique(selected$sample_id))
pdf(file.path(FIG, "01_target_program_signature_by_sample_lineage_group.pdf"), width = max(10, 6 * length(samples)), height = 7.5, useDingbats = FALSE)
par(mfrow = c(1, length(samples)), mar = c(11, 4.5, 3.5, 1.5))
for (s in samples) plot_groups(s)
par(mfrow = c(1, 1))
dev.off()
pdf(file.path(FIG, "02_marker_expression_by_sample_lineage_group.pdf"), width = max(11, 6 * length(samples)), height = 8, useDingbats = FALSE)
par(mfrow = c(1, length(samples)), mar = c(7, 6.5, 3, 1))
for (s in samples) {
  q <- expr[expr$sample_id == s, , drop = FALSE]
  g <- levels_groups[levels_groups %in% q$group]
  gx <- match(q$group, g)
  gy <- match(q$gene, reviewer_markers)
  plot(NA, xlim = c(.5, length(g) + .5), ylim = c(.5, length(reviewer_markers) + .5), xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = paste("7 dpi:", s))
  axis(1, at = seq_along(g), labels = g, las = 2, cex.axis = .6)
  axis(2, at = seq_along(reviewer_markers), labels = reviewer_markers, las = 2, cex.axis = .65)
  color_idx <- ifelse(is.finite(q$mean_logCP10k), pmax(1, pmin(10, ceiling(q$mean_logCP10k * 4))), 1)
  points(gx, gy, pch = 21, bg = grDevices::gray.colors(10, start = 1, end = .2)[color_idx], col = "gray40", cex = .1 + 2.1 * sqrt(pmax(0, q$pct_detected)), lwd = .3)
}
par(mfrow = c(1, 1))
dev.off()
expected_tables <- c( "01_frozen_Arm2_signature_gene_coverage.csv", "02_sample_by_lineage_group_counts.csv", "03_target_signature_by_sample_lineage_group.csv", "04_marker_genes_not_in_restricted_matrix.txt", "05_marker_expression_by_sample_lineage_group.csv", "06_Arm1_Arm2_step2_cluster_provenance.csv", "07_lineage_AUCell_by_sample_group.csv", "08_author_neutrophil_overlap_audit.csv" )
checks <- data.frame( check = c("Original xjz inputs exist", "Frozen signature covers >=400 genes", "Exactly two 7dpi biological samples", "Arm1 and Arm2 both represented", "Excluded infiltrating and BAM populations represented", "No NA/NaN in scored groups", "All eight analysis tables written", "Both figures produced"), passed = c(all(file.exists(inputs)), length(sig_present) >= 400, length(samples) == 2L, all(c("Arm_1_alternative", "Arm_2_target") %in% as.character(selected$group)), all(c("Excluded_infiltrating_myeloid", "Excluded_BAM_like") %in% as.character(selected$group)), all(is.finite(selected$target_signature_mean_logCP10k)), all(file.exists(file.path(TBL, expected_tables))), all(file.exists(file.path(FIG, c("01_target_program_signature_by_sample_lineage_group.pdf", "02_marker_expression_by_sample_lineage_group.pdf"))))), stringsAsFactors = FALSE)
write.csv(checks, file.path(TBL, "09_STEP2D_audit_checks.csv"), row.names = FALSE)
writeLines(log_text, file.path(OUT, "00_run_log.txt"), useBytes = TRUE)
if (!all(checks$passed)) {
}
