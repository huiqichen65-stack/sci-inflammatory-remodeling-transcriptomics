options(stringsAsFactors = FALSE, width = 190)
set.seed(20261024)
ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
REF_DIR <- file.path(ROOT, "01_raw_data", "reference_gene_sets")
REF_FILE <- file.path(REF_DIR, "Marsh2022_Supplementary_Table_4.xlsx")
KEGG_LINKS_FILE <- file.path(REF_DIR, "KEGG_mmu_pathway_gene_links.tsv")
KEGG_NAMES_FILE <- file.path(REF_DIR, "KEGG_mmu_pathway_names.tsv")
SEURAT_FILE <- file.path(ROOT, "04_results", "step5BR2_GSE162610_single_root_branch_audit", "GSE162610_step5BR2_single_root_seurat.rds")
BRANCH_FILE <- file.path(ROOT, "04_results", "step5CR_GSE162610_Y29_prebranch_progressive_divergence", "GSE162610_step5CR_results.rds")
DEG_DIR <- file.path(ROOT, "04_results", "step6A_GSE162610_7dpi_branch_specific_DEG", "tables")
MASTER_FILE <- file.path(DEG_DIR, "01_all_genes_master_evidence.csv")
MODULE_FILE <- file.path(ROOT, "04_results", "step6B_GSE162610_dynamic_branch_modules", "tables", "01_dynamic_module_gene_assignments.csv")
FROZEN_FILE <- file.path(ROOT, "04_results", "step11B_GSE234774_refined_main_Arm_specific_shared_AUCell_ssGSEA", "tables", "04_refined_Arm_specific_shared_signature_audit.csv")
REGULON_FILE <- file.path(ROOT, "04_results", "step7G_GSE162610_regulon_target_program_closure", "tables", "05_priority_TF_target_network_edges.tsv")
OUTDIR <- file.path(ROOT, "04_results", "step6F_GSE162610_dissociation_activation_sensitivity")
TDIR <- file.path(OUTDIR, "tables")
FDIR <- file.path(OUTDIR, "figures")
required_packages <- c("Matrix", "SeuratObject", "readxl", "ggplot2", "irlba")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install missing R packages first: ", paste(missing_packages, collapse = ", "))
required_files <- c(REF_FILE, KEGG_LINKS_FILE, KEGG_NAMES_FILE, SEURAT_FILE, BRANCH_FILE, MASTER_FILE, MODULE_FILE, FROZEN_FILE)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Required files do not exist:\n", paste(missing_files, collapse = "\n"))
dir.create(TDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FDIR, recursive = TRUE, showWarnings = FALSE)
log_path <- file.path(OUTDIR, "00_run_log.txt")
writeLines(c(paste0("Started: ", Sys.time()), "Input sources are READ ONLY."), log_path)
msg <- function(...) {
  x <- paste(..., collapse = "")
  cat(x, "\n")
  cat(x, "\n", file = log_path, append = TRUE)
}
write_csv <- function(df, name) {
  path <- file.path(TDIR, name)
  utils::write.csv(df, path, row.names = FALSE, na = "")
  invisible(path)
}
yn <- function(x) as.character(x) %in% c("TRUE", "T", "true", "1")
score_summary <- function(x) c(median = stats::median(x), q25 = unname(stats::quantile(x, .25)), q75 = unname(stats::quantile(x, .75)), mean = mean(x))
make_report <- function(lines, file="10_STEP6F_FINAL_REPORT.txt") writeLines(lines, file.path(OUTDIR, file))
msg("Reading Marsh Supplementary Table 4: ", REF_FILE)
raw <- as.data.frame(readxl::read_excel(REF_FILE, sheet = 1, col_names = FALSE, .name_repair = "minimal"), stringsAsFactors = FALSE)
header <- trimws(as.character(unlist(raw[1, , drop = TRUE])))
expected <- c(identity = "Microglial Identity Score", all_CNS_activation = "All CNS Cells Act. Score", micro_myeloid_activation = "Micro/Myeloid Shared Act. Score")
find_column <- function(label) {
  hit <- which(header == label)
  if (length(hit) != 1L) stop("Marsh workbook header missing/ambiguous: '", label, "'. Found: ", paste(header, collapse = " | "))
  hit
}
get_genes <- function(label) {
  j <- find_column(label)
  v <- trimws(as.character(raw[-1, j, drop = TRUE]))
  v <- v[!is.na(v) & nzchar(v) & toupper(v) != "NA"]
  unique(v)
}
gene_sets <- lapply(expected, get_genes)
stopifnot(length(gene_sets$identity) == 18L, length(gene_sets$all_CNS_activation) == 19L, length(gene_sets$micro_myeloid_activation) == 25L)
activation_union <- sort(unique(c(gene_sets$all_CNS_activation, gene_sets$micro_myeloid_activation)))
if (length(activation_union) != 34L) stop("Unexpected Marsh activation union size; expected 34, found ", length(activation_union))
if (any(gene_sets$identity %in% activation_union)) {
  msg("NOTE: Some identity-list genes also appear in activation lists; only published activation memberships determine masking.")
}
msg("Marsh input: identity=", length(gene_sets$identity), ", all-CNS=", length(gene_sets$all_CNS_activation), ", micro/myeloid=", length(gene_sets$micro_myeloid_activation), ", activation union=", length(activation_union))
msg("Loading frozen trajectory Seurat and branch definition (read-only) ...")
obj <- readRDS(SEURAT_FILE)
if (!inherits(obj, "Seurat")) stop("Expected original trajectory Seurat object")
if (!"RNA" %in% names(obj@assays)) stop("RNA assay missing")
branch_obj <- readRDS(BRANCH_FILE)
if (is.null(branch_obj$arm1_cells) || is.null(branch_obj$arm2_cells)) {
  stop("Original frozen branch object must contain $arm1_cells and $arm2_cells")
}
meta <- obj@meta.data
if (!all(c("sample_id", "timepoint") %in% names(meta))) stop("Frozen Seurat metadata missing sample_id/timepoint")
cell_ids <- rownames(meta)
arm1 <- intersect(as.character(branch_obj$arm1_cells), cell_ids)
arm2 <- intersect(as.character(branch_obj$arm2_cells), cell_ids)
arm1 <- arm1[as.character(meta[arm1, "timepoint"]) == "7 dpi"]
arm2 <- arm2[as.character(meta[arm2, "timepoint"]) == "7 dpi"]
if (length(intersect(arm1, arm2))) stop("Arm definitions overlap")
if (!identical(c(length(arm1), length(arm2)), c(1399L, 1142L))) {
  stop("Branch membership mismatch with manuscript (Arm_1=1399, Arm_2=1142): observed ",length(arm1)," / ",length(arm2), ". STOP rather than silently altering frozen labels.")
}
cells <- c(arm1, arm2)
md <- meta[cells, , drop = FALSE]
md$cell_id <- cells
md$branch <- c(rep("Arm_1", length(arm1)), rep("Arm_2", length(arm2)))
md$sample_id <- as.character(md$sample_id)
if (length(unique(md$sample_id)) != 2L || any(!grepl("^7dpi", unique(md$sample_id)))) {
  stop("Expected EXACTLY two 7dpi biological sample IDs; observed: ", paste(unique(md$sample_id),collapse=","))
}
cell_counts <- as.data.frame(table(md$sample_id, md$branch), stringsAsFactors = FALSE)
names(cell_counts) <- c("sample_id", "branch", "n_cells")
write_csv(cell_counts, "02_branch_cell_counts_by_biological_sample.csv")
if (any(cell_counts$n_cells < 20L)) stop("Insufficient cells in at least one sample x branch group")
counts <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "counts")
if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")
if (!all(cells %in% colnames(counts))) stop("RNA counts missing some frozen branch cells")
counts <- counts[, cells, drop = FALSE]
rm(obj, branch_obj, meta); invisible(gc())
if (anyDuplicated(rownames(counts))) stop("Duplicate gene names in RNA count matrix")
if (any(Matrix::colSums(counts) <= 0)) stop("Zero UMI counts in branch cells")
scale_factor <- 10000 / as.numeric(Matrix::colSums(counts))
logexp <- counts %*% Matrix::Diagonal(x = scale_factor)
logexp@x <- log1p(logexp@x)
rm(counts); invisible(gc())
available <- rownames(logexp)
set_map <- lapply(gene_sets, function(v) intersect(v, available))
set_map$activation_union <- intersect(activation_union, available)
gene_audit <- do.call(rbind, lapply(names(gene_sets), function(nm) {
  data.frame(set_name=nm, gene=gene_sets[[nm]], in_GSE162610=gene_sets[[nm]] %in% available)
}))
gene_audit <- rbind(gene_audit, data.frame(set_name="activation_union", gene=activation_union, in_GSE162610=activation_union %in% available))
write_csv(gene_audit, "01_Marsh2022_published_gene_set_membership_and_coverage.csv")
coverage <- data.frame(signature=names(set_map), published_n=c(vapply(gene_sets,length,integer(1)),length(activation_union)), detected_n=vapply(set_map,length,integer(1)))
coverage$coverage_fraction <- coverage$detected_n / coverage$published_n
write_csv(coverage, "01B_Marsh2022_signature_coverage_summary.csv")
if (any(coverage$coverage_fraction < .60)) stop("Less than 60% published signature gene coverage; check species / feature names before scoring")
msg("Computing published activation module scores by sample and branch ...")
for (nm in names(set_map)) {
  md[[paste0("score_", nm)]] <- as.numeric(Matrix::colMeans(logexp[set_map[[nm]], , drop=FALSE]))
}
score_names <- paste0("score_", names(set_map))
write_csv(md[, c("cell_id","sample_id","timepoint","branch",score_names), drop=FALSE], "03_cell_level_Marsh_activation_scores.csv")
summary_records <- list()
for (sample in sort(unique(md$sample_id))) for (br in c("Arm_1", "Arm_2")) {
  take <- md$sample_id == sample & md$branch == br
  for (nm in names(set_map)) {
    s <- score_summary(md[[paste0("score_",nm)]][take])
    summary_records[[length(summary_records)+1L]] <- data.frame(sample_id=sample, branch=br, signature=nm, n_cells=sum(take), mapped_genes=length(set_map[[nm]]), median=s["median"], q25=s["q25"], q75=s["q75"], mean=s["mean"])
  }
}
score_tab <- do.call(rbind, summary_records); rownames(score_tab) <- NULL
write_csv(score_tab,"04_activation_score_sample_by_branch_summary.csv")
branch_effects <- do.call(rbind, lapply(unique(score_tab$signature), function(nm) {
  s <- score_tab[score_tab$signature==nm, ,drop=FALSE]
  a <- s[s$branch=="Arm_1",]; b <- s[s$branch=="Arm_2",]
  joined <- merge(a,b,by="sample_id",suffixes=c("_Arm1","_Arm2"))
  data.frame(signature=nm,sample_id=joined$sample_id, median_Arm2_minus_Arm1=joined$median_Arm2 - joined$median_Arm1, mean_Arm2_minus_Arm1=joined$mean_Arm2 - joined$mean_Arm1, n_Arm1=joined$n_cells_Arm1, n_Arm2=joined$n_cells_Arm2)
}))
write_csv(branch_effects,"05_paired_sample_activation_score_effects.csv")
pdat <- md[,c("sample_id","branch","score_all_CNS_activation","score_micro_myeloid_activation","score_activation_union")]
pdat_long <- do.call(rbind,lapply(c("all_CNS_activation","micro_myeloid_activation","activation_union"), function(nm) {
  x <- pdat[,c("sample_id","branch",paste0("score_",nm))]
  names(x)[3] <- "score"; x$signature <- nm; x
}))
pdat_long$branch <- factor(pdat_long$branch,levels=c("Arm_1","Arm_2"))
p <- ggplot2::ggplot(pdat_long, ggplot2::aes(x=branch,y=score,fill=branch)) + ggplot2::geom_violin(trim=TRUE,scale="width",alpha=.50,color="grey50") + ggplot2::geom_boxplot(width=.18,outlier.shape=NA,alpha=.8) + ggplot2::facet_grid(signature~sample_id,scales="free_y") + ggplot2::scale_fill_manual(values=c(Arm_1="#4B78A8",Arm_2="#C63535")) + ggplot2::labs(title="Published dissociation-associated activation signatures", subtitle="Two biological samples shown separately; individual cells are not independent replicates", x=NULL,y="Mean log1p(CP10k) across mapped signature genes") + ggplot2::theme_bw(base_size=11) + ggplot2::theme(legend.position="none")
ggplot2::ggsave(file.path(FDIR,"Fig6F_01_Marsh_activation_scores_by_sample.pdf"),p,width=10,height=10)
ggplot2::ggsave(file.path(FDIR,"Fig6F_01_Marsh_activation_scores_by_sample_review.png"),p,width=10,height=10,dpi=200)
master <- utils::read.csv(MASTER_FILE,check.names=FALSE)
if (!all(c("gene_id","high_confidence_branch_marker","strict_pseudobulk_supported", "technical_like", "cell_p", "pb_PValue") %in% names(master))) {
  stop("DEG master table columns differ from original analysis. Inspect input.")
}
master$highconf <- yn(master$high_confidence_branch_marker)
master$strict <- yn(master$strict_pseudobulk_supported)
if (sum(master$highconf) != 745L || sum(master$strict) != 343L) {
  stop("Original discovery DEG counts differ from manuscript (745 high-confidence, 343 strict). Check source files.")
}
master$in_Marsh_union <- master$gene_id %in% activation_union
master$technical_like_pattern <- grepl("^(Rpl|Rps|mt-|Gm[0-9]+)",master$gene_id)
flag_mismatch <- sum(yn(master$technical_like) != master$technical_like_pattern,na.rm=TRUE)
if (flag_mismatch != 0) stop("technical_like in master table does not match original prefix rule. CHECK SOURCE.")
modules <- utils::read.csv(MODULE_FILE,check.names=FALSE)
if (!"gene_id" %in% names(modules)) stop("Dynamic module gene table missing gene_id")
frozen <- utils::read.csv(FROZEN_FILE,check.names=FALSE)
if (!all(c("record_type","gene","in_Arm2_specific") %in% names(frozen))) stop("Frozen signature audit columns missing")
target_496 <- unique(as.character(frozen$gene[frozen$record_type=="GENE" & yn(frozen$in_Arm2_specific)]))
target_496 <- target_496[!is.na(target_496) & nzchar(target_496)]
if (length(target_496) != 496L) {
  stop("Frozen target signature must have 496 genes; identified ",length(target_496), ". Do not substitute another gene set without verification.")
}
key_groups <- list(assayed_DEG_genes=master$gene_id, high_confidence_745=master$gene_id[master$highconf], strict_pseudobulk_343=master$gene_id[master$strict], frozen_target_496=target_496, dynamic_module_genes=unique(modules$gene_id))
counts_overlap <- do.call(rbind,lapply(names(key_groups),function(name) {
  v <- unique(key_groups[[name]])
  data.frame(gene_group=name, n_original=length(v), n_in_activation_union=sum(v %in% activation_union), n_after_mask=sum(!v %in% activation_union), pct_removed=100*mean(v %in% activation_union))
}))
write_csv(counts_overlap,"06_Marsh_activation_gene_overlap_with_core_evidence.csv")
master$flag_definition <- "Gene name starts with Rpl/Rps/mt-/Gm+digits (not a dissociation marker)"
write_csv(master[master$in_Marsh_union | (master$technical_like_pattern & master$highconf), c("gene_id","in_Marsh_union","technical_like","high_confidence_branch_marker", "strict_pseudobulk_supported","final_direction","cell_avg_log2FC","pb_logFC","flag_definition")], "07_original_technical_like_flag_and_Marsh_gene_audit.csv")
write_csv(data.frame(gene=target_496, in_Marsh_activation=target_496 %in% activation_union, present_in_GSE162610=target_496 %in% available), "08_frozen_target_496_gene_mask.csv")
masked_master <- master[!master$in_Marsh_union, , drop=FALSE]
masked_master$cell_BH_recalculated_within_retained_tests <- stats::p.adjust(masked_master$cell_p,method="BH")
masked_master$pb_BH_recalculated_within_retained_tests <- stats::p.adjust(masked_master$pb_PValue,method="BH")
write_csv(masked_master[masked_master$highconf | masked_master$strict, c("gene_id","final_direction","highconf","strict","cell_avg_log2FC", "repFC_7dpi_sample1","repFC_7dpi_sample2","pb_logFC", "cell_BH_recalculated_within_retained_tests","pb_BH_recalculated_within_retained_tests")], "09_masked_retained_branch_markers_with_original_effect_estimates.csv")
orig_target <- intersect(target_496,available)
masked_target <- setdiff(orig_target,activation_union)
if (length(masked_target) < 50L) stop("Masked target signature too small; check frozen gene extraction")
md$score_frozen_target_original <- as.numeric(Matrix::colMeans(logexp[orig_target,,drop=FALSE]))
md$score_frozen_target_masked <- as.numeric(Matrix::colMeans(logexp[masked_target,,drop=FALSE]))
signature_stats <- do.call(rbind,lapply(c("score_frozen_target_original","score_frozen_target_masked"),function(sn) {
  do.call(rbind,lapply(sort(unique(md$sample_id)),function(s) {
    do.call(rbind,lapply(c("Arm_1","Arm_2"),function(b) {
      x <- md[[sn]][md$sample_id==s & md$branch==b]
      out <- score_summary(x)
      data.frame(sample_id=s,branch=b,score_definition=sn, n_cells=length(x), n_genes=if(sn=="score_frozen_target_original")length(orig_target)else length(masked_target), median=out["median"],q25=out["q25"],q75=out["q75"])
    }))
  }))
}))
write_csv(signature_stats,"10_frozen_target_signature_before_after_mask_by_sample.csv")
msg("Evaluating cross-sample branch separability before/after activation-gene masking ...")
variable_genes <- function(mat,cc,n=1500L) {
  v <- mat[,cc,drop=FALSE]
  det <- Matrix::rowSums(v>0)
  v2 <- v; v2@x <- v2@x^2
  variances <- as.numeric(Matrix::rowMeans(v2) - Matrix::rowMeans(v)^2)
  variances[!is.finite(variances)] <- 0
  variances[det<15] <- -Inf
  ix <- order(variances,decreasing=TRUE)
  rownames(v)[ix[is.finite(variances[ix])]][seq_len(min(n,sum(is.finite(variances))))]
}
centroid_predict <- function(train_s,test_s,genes) {
  train_idx <- which(md$sample_id == train_s)
  test_idx  <- which(md$sample_id == test_s)
  if (length(genes)<100L) stop("Too few training-selected variable genes")
  xtrain <- as.matrix(logexp[genes,train_idx,drop=FALSE])
  xtest <- as.matrix(logexp[genes,test_idx,drop=FALSE])
  mu <- rowMeans(xtrain)
  s <- apply(xtrain,1,stats::sd); s[!is.finite(s)|s<0.1] <- 0.1
  ztrain <- sweep(sweep(xtrain,1,mu,"-"),1,s,"/")
  ztest <- sweep(sweep(xtest,1,mu,"-"),1,s,"/")
  m1 <- rowMeans(ztrain[,md$branch[train_idx]=="Arm_1",drop=FALSE])
  m2 <- rowMeans(ztrain[,md$branch[train_idx]=="Arm_2",drop=FALSE])
  d1 <- colSums(sweep(ztest,1,m1,"-")^2)
  d2 <- colSums(sweep(ztest,1,m2,"-")^2)
  pred <- ifelse(d1<=d2,"Arm_1","Arm_2")
  truth <- md$branch[test_idx]
  data.frame(cell_id=md$cell_id[test_idx],test_sample=test_s,train_sample=train_s, actual=truth,predicted=pred,correct=truth==pred,stringsAsFactors=FALSE)
}
samples <- sort(unique(md$sample_id))
sep_records <- list(); prediction_records <- list()
for (masked in c(FALSE,TRUE)) {
  use_mat <- if(masked) logexp[!rownames(logexp)%in%activation_union,,drop=FALSE] else logexp
  for (train_s in samples) {
    test_s <- setdiff(samples,train_s)
    train_genes <- variable_genes(use_mat,md$sample_id==train_s,n=1500L)
    pred <- centroid_predict(train_s,test_s,train_genes)
    pred$analysis <- if(masked) "Marsh_activation_masked" else "All_genes"
    prediction_records[[length(prediction_records)+1L]] <- pred
    rec <- mean(pred$correct[pred$actual=="Arm_1"])
    sen <- mean(pred$correct[pred$actual=="Arm_2"])
    sep_records[[length(sep_records)+1L]] <- data.frame( analysis=unique(pred$analysis),train_sample=train_s,test_sample=test_s, n_training_genes=length(train_genes),n_test=nrow(pred), balanced_accuracy=(rec+sen)/2, sensitivity_Arm2=sen,specificity_Arm1=rec)
  }
}
pred_tab <- do.call(rbind,prediction_records)
sep_tab <- do.call(rbind,sep_records)
write_csv(pred_tab,"11_cross_sample_heldout_nearest_centroid_predictions.csv")
write_csv(sep_tab,"12_cross_sample_heldout_branch_separability.csv")
pca_meta <- list(); sep_pca <- list()
for (masked in c(FALSE,TRUE)) {
  v <- if(masked) logexp[!rownames(logexp)%in%activation_union,,drop=FALSE] else logexp
  usegenes <- variable_genes(v,rep(TRUE,ncol(v)),n=1500L)
  xx <- as.matrix(v[usegenes,,drop=FALSE]); mu <- rowMeans(xx)
  sdv <- apply(xx,1,stats::sd); sdv[sdv<0.1|!is.finite(sdv)] <- 0.1
  z <- sweep(sweep(xx,1,mu,"-"),1,sdv,"/")
  pca <- irlba::prcomp_irlba(t(z),n=10,center=FALSE,scale.=FALSE)
  pdat <- md[,c("cell_id","sample_id","branch")]
  pdat$PC1 <- pca$x[,1]; pdat$PC2 <- pca$x[,2]
  pdat$analysis <- if(masked) "Marsh activation masked" else "All genes"
  pca_meta[[length(pca_meta)+1L]] <- pdat
  for(s in samples) {
    dt <- pdat[pdat$sample_id==s,,drop=FALSE]
    m1 <- colMeans(dt[dt$branch=="Arm_1",c("PC1","PC2")]);
    m2 <- colMeans(dt[dt$branch=="Arm_2",c("PC1","PC2")])
    between <- sqrt(sum((m1-m2)^2))
    spread <- mean(c(sqrt(rowSums(sweep(as.matrix(dt[dt$branch=="Arm_1",c("PC1","PC2")]),2,m1,"-")^2)), sqrt(rowSums(sweep(as.matrix(dt[dt$branch=="Arm_2",c("PC1","PC2")]),2,m2,"-")^2))))
    sep_pca[[length(sep_pca)+1L]] <- data.frame(analysis=unique(pdat$analysis),sample_id=s, PC1_PC2_centroid_distance=between,mean_within_branch_radius=spread, separation_ratio=between/max(spread,1e-8))
  }
  rm(xx,z,v); invisible(gc())
}
pca_df <- do.call(rbind,pca_meta)
write_csv(do.call(rbind,sep_pca),"13_unsupervised_PCA_branch_centroid_separation.csv")
write_csv(pca_df,"14_unsupervised_PCA_cell_coordinates.csv")
p2 <- ggplot2::ggplot(pca_df, ggplot2::aes(x=PC1,y=PC2,color=branch)) + ggplot2::geom_point(alpha=.4,size=.50) + ggplot2::facet_grid(analysis~sample_id,scales="free") + ggplot2::scale_color_manual(values=c(Arm_1="#4B78A8",Arm_2="#C63535")) + ggplot2::labs(title="Recomputed PCA after removal of published activation genes", subtitle="Original branch labels shown; PCA is not a relearned Monocle trajectory") + ggplot2::theme_bw(base_size=11)
ggplot2::ggsave(file.path(FDIR,"Fig6F_02_masked_vs_unmasked_PCA_by_sample.pdf"),p2,width=10,height=8)
ggplot2::ggsave(file.path(FDIR,"Fig6F_02_masked_vs_unmasked_PCA_by_sample_review.png"),p2,width=10,height=8,dpi=200)
msg("Reading offline mouse KEGG annotation files ...")
read_kegg_two_columns <- function(path, label, min_rows = 100L) {
  if (!file.exists(path)) stop("Missing local KEGG ", label, " file: ", path)
  lines <- readLines(path, n = 5L, warn = FALSE, encoding = "UTF-8")
  if (!length(lines)) stop("Empty KEGG ", label, " file: ", path)
  if (any(grepl("<!DOCTYPE|<html|<head|error|not found", lines, ignore.case = TRUE))) {
    stop("KEGG ", label, " file seems to contain an HTML/error page instead of TSV: ", path)
  }
  if (!grepl("\t", lines[1], fixed = TRUE)) {
    stop("KEGG ", label, " is not TAB-delimited: ", path)
  }
  x <- tryCatch(utils::read.delim(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", fill = FALSE, check.names = FALSE, fileEncoding = "UTF-8"), error = function(e) stop("Failed to read KEGG ", label, ": ", conditionMessage(e)))
  if (ncol(x) != 2L || nrow(x) < min_rows) {
    stop("Invalid KEGG ", label, " format; expected >= ", min_rows, " rows and exactly two TAB-separated columns, got ", nrow(x), " x ", ncol(x))
  }
  x[] <- lapply(x, function(v) trimws(as.character(v)))
  if (anyNA(x) || any(!nzchar(x[[1]])) || any(!nzchar(x[[2]]))) {
    stop("Missing/empty KEGG values in ", label, ": ", path)
  }
  x
}
kegg_links_raw <- read_kegg_two_columns(KEGG_LINKS_FILE, "pathway-to-gene links", min_rows = 1000L)
kegg_names_raw <- read_kegg_two_columns(KEGG_NAMES_FILE, "pathway names", min_rows = 100L)
is_gene <- function(x) grepl("^mmu:[0-9]+$", x)
is_path <- function(x) grepl("^(path:)?mmu[0-9]{5}$", x)
links_genefirst <- mean(is_gene(kegg_links_raw[[1]]) & is_path(kegg_links_raw[[2]]))
links_pathfirst <- mean(is_path(kegg_links_raw[[1]]) & is_gene(kegg_links_raw[[2]]))
if (links_genefirst >= .99 && links_genefirst > links_pathfirst) {
  gene_col <- 1L; path_col <- 2L
  link_order <- "GENE_PATHWAY"
} else if (links_pathfirst >= .99 && links_pathfirst > links_genefirst) {
  gene_col <- 2L; path_col <- 1L
  link_order <- "PATHWAY_GENE"
} else {
  stop("KEGG links are not recognizable mouse GENE/PATHWAY pairs. ", "Check downloaded endpoints and file content; ", "gene-first fraction=", round(links_genefirst, 3), ", pathway-first fraction=", round(links_pathfirst, 3))
}
if (!all(is_path(kegg_names_raw[[1]]))) {
  stop("KEGG pathway name TSV column 1 must contain mouse pathway IDs, e.g. path:mmu04668")
}
term2gene <- unique(data.frame( term = sub("^path:", "", kegg_links_raw[[path_col]]), gene = sub("^mmu:", "", kegg_links_raw[[gene_col]]), stringsAsFactors = FALSE))
term2name <- unique(data.frame( term = sub("^path:", "", kegg_names_raw[[1]]), name = kegg_names_raw[[2]], stringsAsFactors = FALSE))
if (anyDuplicated(term2name$term)) {
  stop("The local KEGG pathway-name file has conflicting/duplicate pathway IDs")
}
if (anyDuplicated(paste(term2gene$term,term2gene$gene,sep="|"))) {
  stop("Unexpected duplicate TERM2GENE rows remained after deduplication")
}
matched_pathways <- intersect(unique(term2gene$term), term2name$term)
match_fraction <- length(matched_pathways) / length(unique(term2gene$term))
if (match_fraction < .95) {
  stop("Only ", round(100 * match_fraction, 1), "% of KEGG pathways in links have names. Check that downloaded KEGG files match.")
}
term2gene <- term2gene[term2gene$term %in% matched_pathways, , drop = FALSE]
term2name <- term2name[term2name$term %in% matched_pathways, , drop = FALSE]
kegg_genes_available <- unique(term2gene$gene)
local_ref_audit <- data.frame( metric = c("links_path", "links_bytes", "links_md5", "links_column_order", "links_raw_rows", "names_path", "names_bytes", "names_md5", "names_raw_rows", "TERM2GENE_unique_rows", "KEGG_pathways_with_names", "KEGG_unique_Entrez_genes", "pathway_name_coverage_fraction"), value = as.character(c(KEGG_LINKS_FILE, file.info(KEGG_LINKS_FILE)$size, unname(tools::md5sum(KEGG_LINKS_FILE)), link_order, nrow(kegg_links_raw), KEGG_NAMES_FILE, file.info(KEGG_NAMES_FILE)$size, unname(tools::md5sum(KEGG_NAMES_FILE)), nrow(kegg_names_raw), nrow(term2gene), length(matched_pathways), length(kegg_genes_available), match_fraction)), stringsAsFactors = FALSE)
write_csv(local_ref_audit, "15B_offline_KEGG_reference_audit.csv")
msg("Local KEGG ready: ", length(matched_pathways), " pathways, ", length(kegg_genes_available), " Entrez IDs, link order=",link_order)
for (fn in c("16_new_masked_GO_BP_KEGG_enrichment.csv", "17_significant_masked_GO_BP_KEGG_enrichment.csv", "16_new_masked_GO_BP_KEGG_enrichment_NOT_PERFORMED.csv", "15C_offline_KEGG_gene_mapping_audit.csv")) {
  target <- file.path(TDIR, fn)
  if (file.exists(target)) unlink(target)
}
msg("Running masked GO BP and OFFLINE KEGG enrichment ...")
enrich_pkgs <- c("clusterProfiler", "org.Mm.eg.db", "AnnotationDbi", "dplyr")
can_enrich <- all(vapply(enrich_pkgs, requireNamespace, logical(1), quietly = TRUE))
write_csv(data.frame(component = enrich_pkgs, installed = vapply(enrich_pkgs, requireNamespace, logical(1), quietly = TRUE)), "15_functional_enrichment_dependency_audit.csv")
functional_results <- list()
enrichment_errors <- character()
enrichment_runs <- list()
kegg_mapping_records <- list()
expected_combinations <- 2L * 2L
if (can_enrich) {
  bg <- setdiff(unique(master$gene_id), activation_union)
  bg <- bg[!is.na(bg) & nzchar(bg)]
  db <- org.Mm.eg.db::org.Mm.eg.db
  bg_map <- suppressMessages(AnnotationDbi::mapIds( db, keys = bg, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"))
  bg_entrez <- unique(unname(bg_map[!is.na(bg_map) & nzchar(bg_map)]))
  bg_kegg <- intersect(bg_entrez, kegg_genes_available)
  if (length(bg_kegg) < 100L) {
    stop("Only ",length(bg_kegg), " tested background genes map to the local KEGG annotation. ", "Check species, IDs, installed org.Mm.eg.db and the downloaded TSV files.")
  }
  for (level in c("highconf", "strict")) for (arm in c("Arm_1", "Arm_2")) {
    label <- paste(level, arm, sep = "_")
    genes <- unique(masked_master$gene_id[ masked_master[[level]] & masked_master$final_direction %in% c(arm, paste0(arm, "-biased"))])
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (length(genes) < 10L) {
      enrichment_errors <- c(enrichment_errors, paste("Too few",label,"genes for enrichment:" ,length(genes)))
      enrichment_runs[[length(enrichment_runs) + 1L]] <- data.frame(branch=arm, evidence_level=level, database="GO_BP/KEGG_OFFLINE", status="TOO_FEW_INPUT_GENES", n_output_rows=0L)
      next
    }
    go <- tryCatch(clusterProfiler::enrichGO( gene = genes, OrgDb = db, keyType = "SYMBOL", ont = "BP", universe = bg, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500), error = function(e) e)
    if (inherits(go, "error")) {
      enrichment_errors <- c(enrichment_errors, paste("GO_BP",label,conditionMessage(go)))
      enrichment_runs[[length(enrichment_runs) + 1L]] <- data.frame(branch=arm,evidence_level=level,database="GO_BP", status="ERROR",n_output_rows=0L)
    } else {
      tg <- if(is.null(go)) data.frame() else as.data.frame(go)
      if(nrow(tg)) {
        tg$evidence_level <- level; tg$branch <- arm; tg$database <- "GO_BP"
        functional_results[[length(functional_results)+1L]] <- tg
      }
      enrichment_runs[[length(enrichment_runs) + 1L]] <- data.frame(branch=arm,evidence_level=level,database="GO_BP", status=if(nrow(tg)) "RESULT" else "NO_ENRICHED_TERMS", n_output_rows=nrow(tg))
    }
    id_map <- suppressMessages(AnnotationDbi::mapIds( db, keys = genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"))
    ids <- unique(unname(id_map[!is.na(id_map) & nzchar(id_map)]))
    kegg_inputs <- intersect(ids, bg_kegg)
    kegg_mapping_records[[length(kegg_mapping_records)+1L]] <- data.frame( evidence_level = level, branch = arm, n_input_gene_symbols = length(genes), n_mapped_entrez = length(ids), n_mapped_entrez_in_local_KEGG_background = length(kegg_inputs), n_tested_gene_symbols_in_background = length(bg), n_tested_entrez_mapped = length(bg_entrez), n_tested_entrez_in_local_KEGG = length(bg_kegg), stringsAsFactors = FALSE)
    if (length(kegg_inputs) < 10L) {
      enrichment_errors <- c(enrichment_errors, paste("KEGG",label,"has fewer than 10 mapped KEGG genes; no enrichment performed"))
      enrichment_runs[[length(enrichment_runs) + 1L]] <- data.frame(branch=arm,evidence_level=level,database="KEGG_OFFLINE", status="INSUFFICIENT_MAPPED_GENES",n_output_rows=0L)
    } else {
      kk <- tryCatch(clusterProfiler::enricher( gene = kegg_inputs, universe = bg_kegg, TERM2GENE = term2gene, TERM2NAME = term2name, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500), error = function(e) e)
      if (inherits(kk, "error")) {
        enrichment_errors <- c(enrichment_errors, paste("KEGG_OFFLINE",label,conditionMessage(kk)))
        enrichment_runs[[length(enrichment_runs) + 1L]] <- data.frame(branch=arm,evidence_level=level,database="KEGG_OFFLINE", status="ERROR",n_output_rows=0L)
      } else {
        tk <- if(is.null(kk)) data.frame() else as.data.frame(kk)
        if(nrow(tk)) {
          tk$evidence_level <- level; tk$branch <- arm; tk$database <- "KEGG_OFFLINE"
          functional_results[[length(functional_results)+1L]] <- tk
        }
        enrichment_runs[[length(enrichment_runs) + 1L]] <- data.frame(branch=arm,evidence_level=level,database="KEGG_OFFLINE", status=if(nrow(tk)) "RESULT" else "NO_ENRICHED_TERMS", n_output_rows=nrow(tk))
      }
    }
  }
} else {
  enrichment_errors <- c(enrichment_errors, paste("Missing functional enrichment R packages:", paste(enrich_pkgs[!vapply(enrich_pkgs, requireNamespace, logical(1), quietly=TRUE)],collapse=", ")))
}
if (length(kegg_mapping_records)) {
  write_csv(do.call(rbind, kegg_mapping_records),"15C_offline_KEGG_gene_mapping_audit.csv")
}
if (length(enrichment_runs)) {
  enrich_run_tbl <- do.call(rbind, enrichment_runs)
} else {
  enrich_run_tbl <- data.frame(branch=character(), evidence_level=character(), database=character(),status=character(),n_output_rows=integer())
}
write_csv(enrich_run_tbl,"15D_functional_enrichment_run_status.csv")
if (length(functional_results)) {
  tbl <- dplyr::bind_rows(functional_results)
  write_csv(tbl, "16_new_masked_GO_BP_KEGG_enrichment.csv")
  write_csv(tbl[!is.na(tbl$p.adjust) & tbl$p.adjust <= 0.05, ,drop=FALSE], "17_significant_masked_GO_BP_KEGG_enrichment.csv")
} else {
  write_csv(data.frame(note = paste( "No GO BP/KEGG results produced; see 15D_functional_enrichment_run_status.csv", "and 09_functional_enrichment_status.txt")), "16_new_masked_GO_BP_KEGG_enrichment_NOT_PERFORMED.csv")
}
n_go_success <- sum(enrich_run_tbl$database=="GO_BP" & enrich_run_tbl$status %in% c("RESULT","NO_ENRICHED_TERMS"))
n_kegg_success <- sum(enrich_run_tbl$database=="KEGG_OFFLINE" & enrich_run_tbl$status %in% c("RESULT","NO_ENRICHED_TERMS"))
go_complete <- can_enrich && n_go_success == expected_combinations
kegg_complete <- can_enrich && n_kegg_success == expected_combinations
writeLines(c( "KEGG_method=clusterProfiler::enricher_OFFLINE_TERM2GENE_TERM2NAME", paste0("local_KEGG_links=",KEGG_LINKS_FILE), paste0("local_KEGG_names=",KEGG_NAMES_FILE), paste0("GO_BP_completed_comparisons=",n_go_success,"/",expected_combinations), paste0("KEGG_OFFLINE_completed_comparisons=",n_kegg_success,"/",expected_combinations), paste0("GO_BP_complete=",go_complete), paste0("KEGG_OFFLINE_complete=",kegg_complete), if(length(enrichment_errors)) enrichment_errors else "no_enrichment_runtime_errors"), file.path(OUTDIR,"09_functional_enrichment_status.txt"))
key_genes <- c("Fos","Egr1","Junb","Atf3","Nfkbiz","Ccl3","Ccl4","Dnajb1", "Tnf","Tnfrsf1a","Nfe2l2","Nfkb1","Relb","Nfkb2","C5ar1","C3ar1","Myd88","Icam1")
key <- merge(data.frame(gene_id=key_genes),master,by="gene_id",all.x=TRUE,sort=FALSE)
key$in_Marsh_union <- key$gene_id %in% activation_union
write_csv(key,"18_key_inflammatory_stress_axis_genes_mask_status.csv")
reg_status <- "NO_INPUT_FILE"
if(file.exists(REGULON_FILE)) {
  edges <- utils::read.delim(REGULON_FILE,check.names=FALSE)
  if(all(c("TF","target","highconf_marker","strict_pseudobulk_marker")%in%names(edges))) {
    edges <- edges[edges$TF %in% c("Nfkb1","Nfe2l2","Relb","Nfkb2"),,drop=FALSE]
    if(nrow(edges)) {
      rows <- do.call(rbind,lapply(split(edges,edges$TF),function(t) {
        data.frame(TF=unique(t$TF)[1],original_motif_edges=nrow(t), activation_genes_in_targets=sum(t$target %in% activation_union), original_highconf_target_edges=sum(yn(t$highconf_marker)), retained_highconf_target_edges=sum(yn(t$highconf_marker)&!t$target%in%activation_union), original_strict_target_edges=sum(yn(t$strict_pseudobulk_marker)), retained_strict_target_edges=sum(yn(t$strict_pseudobulk_marker)&!t$target%in%activation_union))
      }))
      write_csv(rows,"19_descriptive_masked_regulon_target_edge_support.csv")
      reg_status <- "PRECOMPUTED_EDGE_OVERLAP_ONLY"
    }
  } else reg_status <- "INCOMPATIBLE_EDGE_TABLE_SCHEMA"
}
checks <- data.frame( check=c("published_Marsh_table_headers","union_34_genes","two_7dpi_biological_samples", "frozen_Arm1_1399","frozen_Arm2_1142","frozen_target_signature_496", "technical_flag_matches_prefix_rule","per_sample_scoring_done", "cross_sample_prediction_done","GO_BP_recomputed","KEGG_OFFLINE_recomputed", "trajectory_topology_relearned","regulon_activity_recomputed","CellChat_recomputed"), status=c("PASS","PASS","PASS","PASS","PASS","PASS","PASS","PASS","PASS", if(go_complete) "PASS" else "NOT_COMPLETE", if(kegg_complete) "PASS" else "NOT_COMPLETE", "NOT_PERFORMED","NOT_PERFORMED","NOT_PERFORMED"), stringsAsFactors=FALSE)
write_csv(checks,"20_STEP6F_audit_checks.csv")
input_manifest <- data.frame( file=required_files,bytes=file.info(required_files)$size,md5=unname(tools::md5sum(required_files)))
write_csv(input_manifest,"00_input_file_manifest.csv")
writeLines(capture.output(sessionInfo()),file.path(OUTDIR,"00_sessionInfo.txt"))
report <- c( "GSE162610 | Reviewer 2 Major Comment 2 | Dissociation sensitivity audit", paste("Completed:",Sys.time()), "Marsh source: Supplementary Table 4; exact 3 column headings checked.", sprintf("Marsh identity/all-CNS activation/micro-myeloid activation/union = %d/%d/%d/%d", length(gene_sets$identity),length(gene_sets$all_CNS_activation),length(gene_sets$micro_myeloid_activation),length(activation_union)), paste("Activation union gene coverage:",length(set_map$activation_union),"of",length(activation_union)), paste("Frozen Arm_1 / Arm_2 cells:",length(arm1),"/",length(arm2)), "Published activation submodules scored separately in both 7dpi biological samples.", paste("High-confidence markers retained after mask:",counts_overlap$n_after_mask[counts_overlap$gene_group=="high_confidence_745"]), paste("Strict pseudobulk markers retained after mask:",counts_overlap$n_after_mask[counts_overlap$gene_group=="strict_pseudobulk_343"]), paste("Frozen target genes retained after mask:",length(masked_target)), paste("Technical_like prefix-flag mismatch:",flag_mismatch), "technical_like is a NAME-PREFIX flag, NOT a dissociation score/marker.", paste("Local KEGG file (links):",KEGG_LINKS_FILE), paste("Local KEGG file (names):",KEGG_NAMES_FILE), paste("Local KEGG pathway / Entrez coverage:",length(matched_pathways),"/",length(kegg_genes_available)), paste("GO BP comparisons completed:",n_go_success,"of",expected_combinations), paste("Offline KEGG comparisons completed:",n_kegg_success,"of",expected_combinations), "GO/KEGG audit: 09_functional_enrichment_status.txt; tables/15D_functional_enrichment_run_status.csv.", paste("Regulon assessment:",reg_status), "Cross-sample centroid analysis uses original branch labels and TRAINING-SAMPLE feature selection.", "Unsupervised PCA is rebuilt after masking; neither method re-learns Monocle trajectory topology.", "Gene-masked DEG sensitivity retains precomputed estimates; Wilcoxon/edgeR models are NOT refitted.", "No individual cells are considered independent biological replicates.", "No causal or independent validation claim is made.", "KEY OUTPUTS:", "  10_STEP6F_FINAL_REPORT.txt", "  tables/01B_Marsh2022_signature_coverage_summary.csv", "  tables/05_paired_sample_activation_score_effects.csv", "  tables/06_Marsh_activation_gene_overlap_with_core_evidence.csv", "  tables/08_frozen_target_496_gene_mask.csv", "  tables/10_frozen_target_signature_before_after_mask_by_sample.csv", "  tables/12_cross_sample_heldout_branch_separability.csv", "  tables/13_unsupervised_PCA_branch_centroid_separation.csv", "  tables/18_key_inflammatory_stress_axis_genes_mask_status.csv", "  tables/19_descriptive_masked_regulon_target_edge_support.csv (if present)", "  tables/20_STEP6F_audit_checks.csv", "  09_functional_enrichment_status.txt", "  tables/15B_offline_KEGG_reference_audit.csv", "  tables/15C_offline_KEGG_gene_mapping_audit.csv (if present)", "  tables/15D_functional_enrichment_run_status.csv", "  tables/16_new_masked_GO_BP_KEGG_enrichment.csv (if present)", "  figures/Fig6F_01_Marsh_activation_scores_by_sample_review.png", "  figures/Fig6F_02_masked_vs_unmasked_PCA_by_sample_review.png" )
make_report(report)
review_rel <- c( "10_STEP6F_FINAL_REPORT.txt", "09_functional_enrichment_status.txt", "tables/15B_offline_KEGG_reference_audit.csv", "tables/15C_offline_KEGG_gene_mapping_audit.csv", "tables/15D_functional_enrichment_run_status.csv", "tables/01B_Marsh2022_signature_coverage_summary.csv", "tables/04_activation_score_sample_by_branch_summary.csv", "tables/05_paired_sample_activation_score_effects.csv", "tables/06_Marsh_activation_gene_overlap_with_core_evidence.csv", "tables/07_original_technical_like_flag_and_Marsh_gene_audit.csv", "tables/08_frozen_target_496_gene_mask.csv", "tables/10_frozen_target_signature_before_after_mask_by_sample.csv", "tables/12_cross_sample_heldout_branch_separability.csv", "tables/13_unsupervised_PCA_branch_centroid_separation.csv", "tables/18_key_inflammatory_stress_axis_genes_mask_status.csv", "tables/19_descriptive_masked_regulon_target_edge_support.csv", "tables/20_STEP6F_audit_checks.csv", "tables/16_new_masked_GO_BP_KEGG_enrichment.csv", "tables/17_significant_masked_GO_BP_KEGG_enrichment.csv", "figures/Fig6F_01_Marsh_activation_scores_by_sample_review.png", "figures/Fig6F_02_masked_vs_unmasked_PCA_by_sample_review.png" )
review_rel <- review_rel[file.exists(file.path(OUTDIR,review_rel))]
if (requireNamespace("zip",quietly=TRUE)) {
  review_zip <- file.path(dirname(ROOT),"Reviewer2_MajorComment2_Step6F_Review_Package.zip")
  tryCatch({
    if (file.exists(review_zip)) unlink(review_zip)
    zip::zipr(zipfile=review_zip,files=review_rel,root=OUTDIR,include_directories=FALSE)
    if (!file.exists(review_zip)) stop("ZIP output does not exist after packaging")
    if (file.info(review_zip)$size > 512 * 1024^2) {
      unlink(review_zip)
      stop("Review ZIP exceeded 512 MiB; use the listed small results instead")
    }
    msg("Output archive: ",review_zip," (",round(file.info(review_zip)$size/1024^2,2)," MiB)")
  },error=function(e) msg("Review ZIP optional packaging failed: ",conditionMessage(e),
                          ". Use the files listed in the results directory."))
} else {
  msg("Optional package 'zip' is absent; results are complete. Use the files listed in the final report.")
}
msg("SUCCESS: ",file.path(OUTDIR,"10_STEP6F_FINAL_REPORT.txt"))
