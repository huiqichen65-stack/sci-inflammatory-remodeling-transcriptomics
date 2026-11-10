options(stringsAsFactors = FALSE, warn = 1)
set.seed(20261026L)
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
stopifnot(dir.exists(PROJECT_ROOT))
PROJECT_ROOT <- normalizePath(PROJECT_ROOT, winslash = "/", mustWork = TRUE)
STEP_NAME <- "step8G_GSE162610_receptor_depth_overlap_Tnf_source_sensitivity"
OUT <- file.path(PROJECT_ROOT, "04_results", STEP_NAME)
TBL <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")
RDSDIR <- file.path(OUT, "RDS")
for (d in c(OUT, TBL, FIG, RDSDIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
log_lines <- c(paste("Step8G started:", Sys.time()), paste("Root:", PROJECT_ROOT), "All inputs are original local xjz project paths. No ZIP inputs.")
logmsg <- function(...) {
  z <- paste0(..., collapse = "")
  log_lines <<- c(log_lines, paste(format(Sys.time(), "%H:%M:%S"), z))
  writeLines(log_lines, file.path(OUT, "00_run_log.txt"), useBytes = TRUE)
}
fail <- function(...) {
  logmsg("ERROR: ", paste0(..., collapse = ""))
  stop(paste0(..., collapse = ""), call. = FALSE)
}
require_pkg <- function(name) {
  if (!requireNamespace(name, quietly = TRUE)) fail("Required package not installed: ", name)
}
for (p in c("SeuratObject", "Matrix", "data.table", "ggplot2")) require_pkg(p)
src <- function(...) file.path(PROJECT_ROOT, ...)
INPUT <- c( atlas = src("04_results", "step1F_GSE162610_broad_annotation", "GSE162610_broad_annotation_audit.rds"), step2C = src("04_results", "step2C_GSE162610_final_microglia_extraction", "GSE162610_step2_full_myeloid_with_final_lineage.rds"), step8A_meta = src("04_results", "step8A_GSE162610_CellChat_input_preflight", "tables", "02_CellChat_cell_metadata.tsv.gz"), step8C_genes = src("04_results", "step8C_GSE162610_recipient_program_integration", "tables", "02_recipient_program_gene_coverage.tsv"), step8C_scores = src("04_results", "step8C_GSE162610_recipient_program_integration", "RDS", "GSE162610_step8C_recipient_program_cell_scores.rds"), step8C_effects = src("04_results", "step8C_GSE162610_recipient_program_integration", "tables", "05_7dpi_axis_program_receptor_audit_long.tsv") )
missing_input <- INPUT[!file.exists(INPUT)]
if (length(missing_input)) fail("Missing original xjz input(s):\n", paste(missing_input, collapse = "\n"))
input_bytes <- as.numeric(file.info(INPUT)$size)
md5_values <- rep(NA_character_, length(INPUT))
small_input <- which(is.finite(input_bytes) & input_bytes <= 100 * 1024^2)
if (length(small_input)) md5_values[small_input] <- unname(tools::md5sum(INPUT[small_input]))
input_audit <- data.frame(input_key = names(INPUT), path = unname(INPUT), bytes = input_bytes, md5_if_le_100MB = md5_values, stringsAsFactors = FALSE)
data.table::fwrite(input_audit, file.path(TBL, "01_input_file_and_md5_audit.csv"))
logmsg("Located all ", length(INPUT), " original project inputs.")
layer_data <- function(object, layer) {
  x <- tryCatch(SeuratObject::LayerData(object, assay = "RNA", layer = layer), error = function(e) NULL)
  if (is.null(x)) {
    x <- tryCatch(SeuratObject::GetAssayData(object, assay = "RNA", slot = layer), error = function(e) NULL)
  }
  if (is.null(x)) fail("Could not read RNA ", layer, " layer from Step1F atlas.")
  x
}
read_tsv_gz <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  utils::read.delim(con, check.names = FALSE, stringsAsFactors = FALSE)
}
write_csv <- function(x, name) data.table::fwrite(x, file.path(TBL, name), na = "NA")
write_csv_gz <- function(x, name) {
  con <- gzfile(file.path(TBL, name), open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.csv(x, con, row.names = FALSE, na = "NA")
}
positive_logical <- function(z) tolower(trimws(as.character(z))) %in% c("true", "t", "1")
meta8a <- read_tsv_gz(INPUT[["step8A_meta"]])
required_a <- c("cell", "sample_id", "timepoint_standard", "cellchat_group", "include_step8_timecourse")
if (!all(required_a %in% names(meta8a))) fail("Step8A metadata missing: ", paste(setdiff(required_a, names(meta8a)), collapse = ", "))
if (anyDuplicated(meta8a$cell)) fail("Duplicated Step8A barcodes")
meta8a$include_flag <- positive_logical(meta8a$include_step8_timecourse)
all7 <- meta8a[meta8a$timepoint_standard == "7dpi", , drop = FALSE]
if (nrow(all7) != 14505L) fail("Unexpected 7dpi cell count: ", nrow(all7), " (expected 14505).")
samples <- sort(unique(all7$sample_id))
if (!identical(samples, c("7dpi_sample1", "7dpi_sample2"))) fail("Unexpected biological sample IDs: ", paste(samples, collapse = ", "))
endo_ref <- meta8a[meta8a$cellchat_group == "Endothelial" & meta8a$include_flag & meta8a$timepoint_standard %in% c("1dpi", "3dpi", "7dpi"), , drop = FALSE]
endo_ref <- endo_ref[order(endo_ref$cellchat_group, endo_ref$timepoint_standard, endo_ref$sample_id, endo_ref$cell), , drop = FALSE]
endo7 <- endo_ref[endo_ref$timepoint_standard == "7dpi", , drop = FALSE]
if (nrow(endo7) != 1759L) fail("Expected 1759 7dpi endothelial cells, got ", nrow(endo7))
logmsg("Step8A: ", nrow(all7), " 7dpi cells, ", nrow(endo7), " 7dpi endothelial; ", nrow(endo_ref), " endothelial cells for ORIGINAL 1/3/7dpi z-score reference.")
coverage <- data.table::fread(INPUT[["step8C_genes"]], sep = "\t", data.table = FALSE)
needed_progs <- c("ENDOTHELIAL_ACTIVATION", "TNF_NFKB_RESPONSE")
if (!all(c("response_program", "n_present", "present_genes") %in% names(coverage))) fail("Step8C gene coverage format changed.")
if (!all(needed_progs %in% coverage$response_program)) fail("Missing original endothelial program(s).")
get_genes <- function(label) {
  z <- coverage[match(label, coverage$response_program), , drop = FALSE]
  g <- trimws(strsplit(as.character(z$present_genes), "\\|", fixed = FALSE)[[1]])
  g <- unique(g[nzchar(g) & !is.na(g)])
  if (length(g) != as.integer(z$n_present)) fail("Step8C gene count mismatch for ", label)
  g
}
g_act <- get_genes("ENDOTHELIAL_ACTIVATION")
g_tnf <- get_genes("TNF_NFKB_RESPONSE")
shared <- intersect(g_act, g_tnf)
act_unique <- setdiff(g_act, shared)
tnf_unique <- setdiff(g_tnf, shared)
if (!identical(c(length(g_act), length(g_tnf), length(shared)), c(14L, 20L, 8L))) {
  fail("Original program sizes/overlap are not 14/20/8. Stopping rather than silently changing gene sets.")
}
if (length(act_unique) < 4 || length(tnf_unique) < 4) fail("Overlap-free programs have too few genes.")
gene_sets <- list(ENDOTHELIAL_ACTIVATION = g_act, TNF_NFKB_RESPONSE = g_tnf, ENDOTHELIAL_ACTIVATION_NONOVERLAP = act_unique, TNF_NFKB_RESPONSE_NONOVERLAP = tnf_unique)
gene_table <- do.call(rbind, lapply(names(gene_sets), function(nm) {
  data.frame(program = nm, gene = gene_sets[[nm]], in_shared_8 = gene_sets[[nm]] %in% shared, stringsAsFactors = FALSE)
}))
write_csv(gene_table, "02_original_and_nonoverlap_program_genes.csv")
write_csv(data.frame(shared_gene = shared), "02B_eight_shared_genes.csv")
logmsg("Loading original Step1F whole-atlas RDS (may require substantial RAM)...")
atlas <- readRDS(INPUT[["atlas"]])
if (!inherits(atlas, "Seurat")) fail("Step1F input is not a Seurat object.")
if (!"RNA" %in% SeuratObject::Assays(atlas)) fail("Original Step1F atlas has no RNA assay.")
md <- atlas@meta.data
if (!all(c("nCount_RNA", "nFeature_RNA") %in% colnames(md))) fail("Step1F metadata lacks RNA nUMI/nFeature.")
if (!all(all7$cell %in% rownames(md)) || !all(endo_ref$cell %in% rownames(md))) fail("Step8A cells not found in Step1F metadata.")
all7$nUMI <- as.numeric(md[all7$cell, "nCount_RNA"])
all7$nFeature <- as.numeric(md[all7$cell, "nFeature_RNA"])
endo7$nUMI <- as.numeric(md[endo7$cell, "nCount_RNA"])
endo7$nFeature <- as.numeric(md[endo7$cell, "nFeature_RNA"])
if (any(!is.finite(endo7$nUMI) | !is.finite(endo7$nFeature) | endo7$nUMI <= 0 | endo7$nFeature <= 0)) fail("Invalid endothelial nUMI/nFeature.")
rna_counts <- layer_data(atlas, "counts")
rna_data <- layer_data(atlas, "data")
reqgenes <- unique(c(unlist(gene_sets, use.names = FALSE), "Tnf", "Tnfrsf1a"))
if (!all(reqgenes %in% rownames(rna_data))) fail("Missing normalized RNA genes: ", paste(setdiff(reqgenes, rownames(rna_data)), collapse = ", "))
if (!all(c("Tnf", "Tnfrsf1a") %in% rownames(rna_counts))) fail("RNA raw count layer lacks Tnf or Tnfrsf1a.")
if (!all(all7$cell %in% colnames(rna_counts)) || !all(endo_ref$cell %in% colnames(rna_data))) fail("Missing expected cell IDs in RNA layers.")
endo_expr <- as.matrix(rna_data[unique(unlist(gene_sets, use.names = FALSE)), endo_ref$cell, drop = FALSE])
score_genes <- function(genes) {
  a <- endo_expr[genes, , drop = FALSE]
  avg <- rowMeans(a)
  s <- apply(a, 1L, stats::sd)
  variable <- is.finite(s) & s > 0
  if (sum(variable) < 4) fail("Fewer than 4 variable genes in program.")
  b <- sweep(a[variable, , drop = FALSE], 1L, avg[variable], "-")
  b <- sweep(b, 1L, s[variable], "/")
  stats::setNames(colMeans(b), colnames(a))
}
score_lookup <- lapply(gene_sets, score_genes)
for (nm in names(score_lookup)) if (any(!is.finite(score_lookup[[nm]]))) fail("Nonfinite scores in ", nm)
all7$Tnf_counts <- as.numeric(rna_counts["Tnf", all7$cell, drop = TRUE])
all7$Tnf_logRNA <- as.numeric(rna_data["Tnf", all7$cell, drop = TRUE])
all7$Tnf_detected <- all7$Tnf_counts > 0
endo7$Tnfrsf1a_counts <- as.numeric(rna_counts["Tnfrsf1a", endo7$cell, drop = TRUE])
endo7$Tnfrsf1a_logRNA <- as.numeric(rna_data["Tnfrsf1a", endo7$cell, drop = TRUE])
endo7$receptor_positive <- as.integer(endo7$Tnfrsf1a_logRNA > 0)
if (!identical(endo7$receptor_positive == 1L, endo7$Tnfrsf1a_counts > 0)) fail("Tnfrsf1a raw/log detection mismatch.")
for (nm in names(score_lookup)) endo7[[nm]] <- unname(score_lookup[[nm]][endo7$cell])
rm(rna_counts, rna_data, endo_expr, score_lookup, atlas, md)
invisible(gc(verbose = FALSE))
logmsg("Original normalized program scores and 7dpi expression extracted from local xjz atlas.")
orig_scores <- readRDS(INPUT[["step8C_scores"]])
if (!is.data.frame(orig_scores)) fail("Step8C program score RDS is not a data frame.")
required_s <- c("cell", "sample_id", "target", "response_program", "program_score")
if (!all(required_s %in% colnames(orig_scores))) fail("Step8C original cell-score fields changed.")
orig_scores <- orig_scores[orig_scores$target == "Endothelial" & orig_scores$response_program %in% needed_progs & orig_scores$cell %in% endo7$cell, , drop = FALSE]
if (nrow(orig_scores) != 2L * nrow(endo7)) fail("Step8C original endothelial scores incomplete: ", nrow(orig_scores))
orig_scores$recomputed <- vapply(seq_len(nrow(orig_scores)), function(i) {
  endo7[[orig_scores$response_program[i]]][match(orig_scores$cell[i], endo7$cell)]
}, numeric(1))
orig_scores$abs_diff <- abs(orig_scores$program_score - orig_scores$recomputed)
score_error <- max(orig_scores$abs_diff)
if (!is.finite(score_error) || score_error > 1e-6) {
  write_csv(orig_scores, "03A_failed_cell_score_reproduction.csv")
  fail("Step8C cell-score reproduction FAILED. Max absolute difference = ", signif(score_error, 6), ". See 03A table.")
}
rm(orig_scores)
invisible(gc(verbose = FALSE))
orig_effects <- data.table::fread(INPUT[["step8C_effects"]], sep = "\t", data.table = FALSE)
orig_effects <- orig_effects[orig_effects$target == "Endothelial" & orig_effects$ligand == "Tnf" & orig_effects$receptor == "Tnfrsf1a" & orig_effects$response_program %in% needed_progs, , drop = FALSE]
if (nrow(orig_effects) != 4L) fail("Expected exactly four original Step8C comparison rows, found ", nrow(orig_effects))
repro <- do.call(rbind, lapply(needed_progs, function(nm) {
  do.call(rbind, lapply(samples, function(sid) {
    dd <- endo7[endo7$sample_id == sid, , drop = FALSE]
    pos <- dd$receptor_positive == 1L
    eff <- mean(dd[[nm]][pos]) - mean(dd[[nm]][!pos])
    old <- orig_effects[orig_effects$sample_id == sid & orig_effects$response_program == nm, , drop = FALSE]
    if (nrow(old) != 1L) fail("Missing Step8C original row: ", sid, " / ", nm)
    data.frame(sample_id = sid, program = nm, n_cells = nrow(dd), n_positive = sum(pos), n_negative = sum(!pos), original_effect = old$effect_receptor_positive_minus_negative, reproduced_effect = eff, absolute_difference = abs(eff - old$effect_receptor_positive_minus_negative), expected_positive = old$n_receptor_positive, expected_negative = old$n_receptor_negative, original_cell_score_max_abs_diff = score_error)
  }))
}))
repro$PASS <- repro$absolute_difference < 1e-6 & repro$n_positive == repro$expected_positive & repro$n_negative == repro$expected_negative
write_csv(repro, "03_original_step8C_reproduction_audit.csv")
if (!all(repro$PASS)) fail("Failed to reproduce original Step8C effects; see Table 03.")
logmsg("Original Step8C cell scores and all four endothelial effects reproduced: PASS.")
model_designs <- list(unadjusted = character(0), adjust_log_nUMI = "log_nUMI", adjust_log_nFeature = "log_nFeature", adjust_both = c("log_nUMI", "log_nFeature"))
endo7$log_nUMI <- log1p(endo7$nUMI)
endo7$log_nFeature <- log1p(endo7$nFeature)
fit_effect <- function(d, program, adjustment, covars) {
  x <- data.frame(y = d[[program]], pos = d$receptor_positive, log_nUMI = d$log_nUMI, log_nFeature = d$log_nFeature)
  x <- x[stats::complete.cases(x[, c("y", "pos", covars), drop = FALSE]), , drop = FALSE]
  if (nrow(x) < 30L || length(unique(x$pos)) != 2L) fail("Insufficient model data.")
  f <- stats::reformulate(c("pos", covars), response = "y")
  fit <- stats::lm(f, data = x)
  cc <- summary(fit)$coefficients
  if (!"pos" %in% rownames(cc)) fail("Receptor coefficient absent: ", program)
  est <- unname(cc["pos", "Estimate"])
  se <- unname(cc["pos", "Std. Error"])
  if (!is.finite(est) || !is.finite(se)) fail("Nonfinite receptor model estimate.")
  data.frame(sample_id = unique(d$sample_id), program = program, gene_set = if (grepl("NONOVERLAP", program)) "NONOVERLAP" else "ORIGINAL", covariate_model = adjustment, n_cells = nrow(x), n_pos = sum(x$pos == 1L), n_neg = sum(x$pos == 0L), beta_receptor_positive = est, se_cell_level_descriptive = se, ci95_low_cell_level_descriptive = est - 1.96 * se, ci95_high_cell_level_descriptive = est + 1.96 * se, positive_direction = est > 0, model_r_squared = summary(fit)$r.squared, stringsAsFactors = FALSE)
}
models <- do.call(rbind, lapply(samples, function(sid) {
  d <- endo7[endo7$sample_id == sid, , drop = FALSE]
  do.call(rbind, lapply(names(gene_sets), function(program) {
    do.call(rbind, lapply(names(model_designs), function(nm) {
      fit_effect(d, program, nm, model_designs[[nm]])
    }))
  }))
}))
write_csv(models, "04_receptor_depth_adjusted_within_sample_models.csv")
qc <- do.call(rbind, lapply(samples, function(sid) {
  d <- endo7[endo7$sample_id == sid, , drop = FALSE]
  do.call(rbind, lapply(c(0L, 1L), function(val) {
    x <- d[d$receptor_positive == val, , drop = FALSE]
    data.frame(sample_id = sid, receptor_positive = as.logical(val), n_cells = nrow(x), fraction_of_sample = nrow(x) / nrow(d), median_nUMI = median(x$nUMI), median_nFeature = median(x$nFeature), mean_nUMI = mean(x$nUMI), mean_nFeature = mean(x$nFeature))
  }))
}))
write_csv(qc, "05_receptor_positive_negative_depth_audit.csv")
write_csv_gz(endo7[, c("cell", "sample_id", "receptor_positive", "Tnfrsf1a_counts", "Tnfrsf1a_logRNA", "nUMI", "nFeature", names(gene_sets)), drop = FALSE], "06_endothelial_cell_level_scores_and_covariates.csv.gz")
corrs <- do.call(rbind, lapply(samples, function(sid) {
  d <- endo7[endo7$sample_id == sid, , drop = FALSE]
  data.frame(sample_id = sid, original_program_rho = suppressWarnings(stats::cor(d$ENDOTHELIAL_ACTIVATION, d$TNF_NFKB_RESPONSE, method = "spearman")), nonoverlap_program_rho = suppressWarnings(stats::cor(d$ENDOTHELIAL_ACTIVATION_NONOVERLAP, d$TNF_NFKB_RESPONSE_NONOVERLAP, method = "spearman")), n_cells = nrow(d))
}))
write_csv(corrs, "07_program_pair_correlations_before_after_deoverlap.csv")
logmsg("Loading ORIGINAL Step2C myeloid annotation RDS, metadata only...")
my_obj <- readRDS(INPUT[["step2C"]])
if (!isS4(my_obj) || !"meta.data" %in% methods::slotNames(my_obj)) fail("Step2C input has no Seurat metadata.")
my_md <- my_obj@meta.data
if (!"final_step2_lineage" %in% colnames(my_md)) fail("Step2C lacks final_step2_lineage.")
my_class <- as.character(my_md[match(all7$cell, rownames(my_md)), "final_step2_lineage"])
rm(my_obj, my_md)
invisible(gc(verbose = FALSE))
all7$step2C_lineage <- my_class
all7$reviewer_source_group <- as.character(all7$cellchat_group)
all7$reviewer_source_group[all7$cellchat_group == "Microglia_Arm_2"] <- "Target-program microglia (Arm_2)"
all7$reviewer_source_group[all7$cellchat_group == "Microglia_Arm_1"] <- "Alternative-branch microglia (Arm_1)"
all7$reviewer_source_group[all7$cellchat_group == "Microglia_other"] <- "Other microglia"
all7$reviewer_source_group[all7$cellchat_group == "Neutrophil"] <- "Neutrophil (Step8A exclusive)"
other_myeloid <- all7$cellchat_group == "Other_myeloid"
macro_cls <- c("infiltrating_macrophage", "macrophage_biased_infiltrating_myeloid", "MHCII_APC_like_infiltrating_macrophage")
all7$reviewer_source_group[other_myeloid & !is.na(my_class) & my_class %in% macro_cls] <- "Infiltrating macrophage-like"
all7$reviewer_source_group[other_myeloid & !is.na(my_class) & my_class == "infiltrating_monocyte"] <- "Infiltrating monocyte"
all7$reviewer_source_group[other_myeloid & !is.na(my_class) & my_class == "IFN_responsive_infiltrating_mono_macro"] <- "IFN-responsive infiltrating mono/mac"
all7$reviewer_source_group[other_myeloid & !is.na(my_class) & my_class == "BAM_like_macrophage"] <- "BAM-like macrophage"
all7$reviewer_source_group[other_myeloid & all7$reviewer_source_group == "Other_myeloid"] <- "Other myeloid (unresolved)"
write_csv(as.data.frame(table(step8A_group = all7$cellchat_group, step2C_lineage = ifelse(is.na(my_class), "not_in_step2C", my_class))), "08_step8A_to_step2C_lineage_crosswalk.csv")
sourcedt <- data.table::as.data.table(all7)
summary_tnf <- function(group_col) {
  d <- sourcedt[, .(n_cells = .N, n_Tnf_detected = sum(Tnf_detected), pct_Tnf_detected = 100 * mean(Tnf_detected), mean_logRNA_all_cells = mean(Tnf_logRNA), mean_logRNA_detected_only = if (sum(Tnf_detected) > 0) mean(Tnf_logRNA[Tnf_detected]) else NA_real_, mean_Tnf_UMI_per_cell = mean(Tnf_counts), total_Tnf_UMI_descriptive = sum(Tnf_counts)), by = c("sample_id", group_col)]
  data.table::setnames(d, group_col, "cell_group")
  d$group_definition <- group_col
  as.data.frame(d)
}
tnf_overview <- summary_tnf("cellchat_group")
tnf_fine <- summary_tnf("reviewer_source_group")
if (sum(tnf_overview$n_cells) != nrow(all7) || sum(tnf_fine$n_cells) != nrow(all7)) fail("Tnf source group counts not exhaustive or groups overlap.")
write_csv(tnf_overview[order(tnf_overview$sample_id, tnf_overview$cell_group), ], "09_Tnf_all_7dpi_Step8A_cell_groups_by_sample.csv")
write_csv(tnf_fine[order(tnf_fine$sample_id, tnf_fine$cell_group), ], "10_Tnf_refined_myeloid_and_other_groups_by_sample.csv")
write_csv_gz(all7[, c("cell", "sample_id", "cellchat_group", "step2C_lineage", "reviewer_source_group", "Tnf_counts", "Tnf_logRNA", "Tnf_detected", "nUMI", "nFeature"), drop = FALSE], "11_Tnf_all_7dpi_cell_level_source_audit.csv.gz")
fig_models <- models[models$covariate_model %in% c("unadjusted", "adjust_log_nUMI", "adjust_log_nFeature", "adjust_both"), , drop = FALSE]
fig_models$covariate_model <- factor(fig_models$covariate_model, levels = names(model_designs), labels = c("Unadjusted", "+ log(nUMI)", "+ log(nFeature)", "+ both"))
fig_models$program <- factor(fig_models$program, levels = names(gene_sets), labels = c("Activation (14 genes)", "TNF/NF-kB (20 genes)", "Activation (6 unique)", "TNF/NF-kB (12 unique)"))
plot1 <- ggplot2::ggplot(fig_models, ggplot2::aes(x = covariate_model, y = beta_receptor_positive, group = sample_id, color = sample_id)) + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") + ggplot2::geom_line(linewidth = 0.6) + ggplot2::geom_point(size = 2.3) + ggplot2::facet_wrap(~ program, scales = "free_y", ncol = 2) + ggplot2::labs(title = "Tnfrsf1a-detected versus undetected endothelial cells (7 dpi)", subtitle = "Within-sample effects; two biological samples; original and nonoverlapping gene sets", x = "Adjustment", y = "Receptor+ minus receptor- program-score difference", color = "Sample", caption = "Descriptive cell-level associations, not causal TNF signaling or replicate-level significance.") + ggplot2::theme_bw(base_size = 10) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1), legend.position = "bottom")
ggplot2::ggsave(file.path(FIG, "01_endothelial_depth_and_overlap_sensitivity.pdf"), plot1, width = 12, height = 8, limitsize = FALSE)
plot2 <- ggplot2::ggplot(tnf_overview, ggplot2::aes(x = sample_id, y = reorder(cell_group, mean_logRNA_all_cells), size = pct_Tnf_detected, color = mean_logRNA_all_cells)) + ggplot2::geom_point() + ggplot2::scale_size_continuous(name = "Tnf detection (%)", range = c(2, 9)) + ggplot2::labs(title = "Tnf transcription across all 7-dpi cell groups", subtitle = "Exclusive Step8A assignments, shown separately for each biological sample", x = "Biological sample", y = "Cell group", color = "Mean log RNA", caption = "RNA detection and normalized expression measure TNF protein or secretion.") + ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "right")
ggplot2::ggsave(file.path(FIG, "02_Tnf_all_7dpi_cell_groups_dotplot.pdf"), plot2, width = 11, height = 8.0, limitsize = FALSE)
plot3 <- ggplot2::ggplot(tnf_fine, ggplot2::aes(x = sample_id, y = reorder(cell_group, mean_logRNA_all_cells), size = pct_Tnf_detected, color = mean_logRNA_all_cells)) + ggplot2::geom_point() + ggplot2::scale_size_continuous(name = "Tnf detection (%)", range = c(2, 9)) + ggplot2::labs(title = "Tnf transcription: resolved myeloid sources at 7 dpi", subtitle = "Step2C lineage subdivision applies only to Step8A Other_myeloid cells", x = "Biological sample", y = "Cell group", color = "Mean log RNA", caption = "Not a protein-source or ligand secretion assay; small groups require caution.") + ggplot2::theme_bw(base_size = 10)
ggplot2::ggsave(file.path(FIG, "03_Tnf_refined_myeloid_source_dotplot.pdf"), plot3, width = 12, height = 9, limitsize = FALSE)
checks <- data.frame( check = c("All original inputs exist", "Two distinct 7dpi samples", "All 14505 7dpi cells included", "7dpi endothelial count 1759", "Original program sizes 14/20 and overlap 8", "Original Step8C effects reproduced", "Non-overlap sets at least 4 genes", "Four covariate specifications", "All exclusive Tnf source groups cover 14505 cells", "All model effects and score inputs finite"), PASS = c(all(file.exists(INPUT)), length(samples) == 2L, nrow(all7) == 14505L, nrow(endo7) == 1759L, identical(c(length(g_act), length(g_tnf), length(shared)), c(14L, 20L, 8L)), all(repro$PASS), length(act_unique) >= 4 && length(tnf_unique) >= 4, length(unique(models$covariate_model)) == 4L, sum(tnf_overview$n_cells) == 14505L && sum(tnf_fine$n_cells) == 14505L, all(is.finite(models$beta_receptor_positive)) && all(is.finite(as.matrix(endo7[, names(gene_sets), drop = FALSE])))))
write_csv(checks, "12_STEP8G_audit_checks.csv")
if (!all(checks$PASS)) fail("Step8G audit did not fully pass; inspect 12_STEP8G_audit_checks.csv")
saveRDS(list(inputs = input_audit, gene_sets = gene_sets, overlap = shared, original_reproduction = repro, models = models, depth_audit = qc, tnf_all_groups = tnf_overview, tnf_refined_groups = tnf_fine, program_correlations = corrs), file.path(RDSDIR, "GSE162610_step8G_analysis_summary.rds"), compress = TRUE)
writeLines(capture.output(utils::sessionInfo()), file.path(OUT, "00_sessionInfo.txt"))
logmsg("FINISHED. All checks PASS; output: ", OUT)
