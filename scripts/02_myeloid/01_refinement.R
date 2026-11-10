options(stringsAsFactors = FALSE)
set.seed(1234)
required_pkgs <- c( "Seurat", "dplyr", "ggplot2", "patchwork", "harmony", "clustree", "pheatmap" )
missing_pkgs <- required_pkgs[ !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE) ]
if (length(missing_pkgs) > 0) {
  stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))
}
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(harmony)
  library(clustree)
  library(pheatmap)
})
project_dir <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
input_rds <- file.path( project_dir, "02_objects", "step1E_GSE162610_Harmony_cluster", "GSE162610_Harmony_multires_clustered.rds" )
out_dir <- file.path( project_dir, "04_results", "step2A_GSE162610_myeloid_refinement" )
obj_dir <- file.path( project_dir, "02_objects", "step2A_GSE162610_myeloid_refinement" )
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(obj_dir, recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_rds)) {
  stop("找不到输入 RDS：\n", input_rds)
}
message("Reading: ", input_rds)
obj <- readRDS(input_rds)
message("[CHECKPOINT 1] RDS loaded successfully")
message("Object class: ", paste(class(obj), collapse = ", "))
message("Cells: ", ncol(obj), " | Features: ", nrow(obj))
assay_names <- tryCatch( as.character(SeuratObject::Assays(obj)), error = function(e) {
    message("Assays(obj) conversion failed: ", e$message)
    names(obj@assays)
  }
)
if (is.null(assay_names) || length(assay_names) == 0) {
  assay_names <- names(obj@assays)
}
assay_names <- as.character(assay_names)
message("Assays detected: ", paste(assay_names, collapse = ", "))
if (!any(assay_names == "RNA")) {
  stop("输入对象中没有 RNA assay。检测到：", paste(assay_names, collapse = ", "))
}
DefaultAssay(obj) <- "RNA"
message("[CHECKPOINT 2] DefaultAssay set to RNA")
joinlayers_available <- exists( "JoinLayers", where = asNamespace("SeuratObject"), mode = "function", inherits = FALSE )
if (isTRUE(joinlayers_available)) {
  obj <- tryCatch( {
      message("Trying JoinLayers(RNA) ...")
      SeuratObject::JoinLayers(obj, assay = "RNA")
    },
    error = function(e) {
      message("JoinLayers skipped safely: ", e$message)
      obj
    }
  )
}
message("[CHECKPOINT 3] Preparing metadata")
meta <- tryCatch( obj[[]], error = function(e) {
    message("obj[[]] failed: ", e$message)
    obj@meta.data
  }
)
if (!is.data.frame(meta)) {
  meta <- as.data.frame(meta)
}
meta_names <- as.character(colnames(meta))
if (is.null(meta_names) || length(meta_names) == 0) {
  stop("metadata 没有列名，无法继续 Step 2。")
}
message("Metadata columns: ", length(meta_names))
message("[CHECKPOINT 4] Metadata loaded successfully")
target_broad <- c("0", "17", "4", "8", "2")
candidate_cluster_cols <- c( "seurat_clusters", "RNA_snn_res.0.4", "harmony_snn_res.0.4", grep("res\\.0\\.4$", meta_names, value = TRUE), grep("cluster", meta_names, value = TRUE, ignore.case = TRUE) )
candidate_cluster_cols <- unique(as.character(candidate_cluster_cols))
candidate_cluster_cols <- candidate_cluster_cols[ vapply(candidate_cluster_cols, function(z) any(meta_names == z), logical(1)) ]
message("Candidate broad-cluster columns: ", paste(candidate_cluster_cols, collapse = ", "))
if (length(candidate_cluster_cols) == 0) {
  stop( "没有找到可识别的 cluster metadata 列。\n", "Check required" )
}
score_cluster_col <- function(x) {
  vals_raw <- meta[[x]]
  vals <- as.character(vals_raw)
  vals <- unique(vals[!is.na(vals)])
  score <- 0
  score <- score + sum(vapply(target_broad, function(z) any(vals == z), logical(1)))
  if (length(vals) == 19) score <- score + 10
  if (all(vapply(as.character(0:18), function(z) any(vals == z), logical(1)))) {
    score <- score + 20
  }
  score
}
cluster_scores <- vapply(candidate_cluster_cols, score_cluster_col, numeric(1))
broad_cluster_col_raw <- candidate_cluster_cols[which.max(cluster_scores)]
message("Best broad-cluster candidate: ", broad_cluster_col_raw, " | score = ", max(cluster_scores))
broad_vals <- as.character(meta[[broad_cluster_col_raw]])
has_all_targets <- all( vapply(target_broad, function(z) any(broad_vals == z, na.rm = TRUE), logical(1)) )
if (!has_all_targets) {
  ident_vals <- as.character(Idents(obj))
  ident_unique <- unique(ident_vals[!is.na(ident_vals)])
  idents_have_targets <- all( vapply(target_broad, function(z) any(ident_unique == z), logical(1)) )
  if (idents_have_targets) {
    obj$broad_cluster_step1 <- ident_vals
    broad_cluster_col <- "broad_cluster_step1"
    message("Using current Idents as broad_cluster_step1")
  } else {
    stop( "无法确认 Step 1 broad cluster 列。\n", "最佳候选列：", broad_cluster_col_raw, "\n", "该列未同时包含 cluster 0,17,4,8,2。\n", "候选列：", paste(candidate_cluster_cols, collapse = ", ") )
  }
} else {
  obj$broad_cluster_step1 <- broad_vals
  broad_cluster_col <- "broad_cluster_step1"
}
message("[CHECKPOINT 5] Broad cluster source confirmed: ", broad_cluster_col)
print(table(obj$broad_cluster_step1, useNA = "ifany"))
find_first_col <- function(candidates, meta_names) {
  hit <- candidates[candidates %in% meta_names]
  if (length(hit) > 0) hit[1] else NA_character_
}
sample_col <- find_first_col( c("sample_id", "sample", "orig.ident", "Sample", "sampleID"), as.character(colnames(obj[[]])) )
time_col <- find_first_col( c("timepoint", "time_point", "Timepoint", "dpi", "group", "condition"), as.character(colnames(obj[[]])) )
author_col <- if ("original_V30" %in% colnames(obj[[]])) {
  "original_V30"
} else {
  tmp <- grep("original.*V30|V30.*original|author.*annot|celltype.*author", colnames(obj[[]]), value = TRUE, ignore.case = TRUE)
  if (length(tmp) > 0) tmp[1] else NA_character_
}
if (is.na(sample_col)) {
  stop("没有找到 sample_id/orig.ident 等样本列；Step 2 Harmony 必须按 biological sample 处理。")
}
if (is.na(time_col)) {
  warning("没有自动识别到 timepoint 列；时间分布表将跳过。")
}
if (is.na(author_col)) {
  warning("没有自动识别到 original_V30；作者 annotation cross-check 将跳过。")
}
message("sample column = ", sample_col)
message("time column = ", ifelse(is.na(time_col), "NOT FOUND", time_col))
message("author column = ", ifelse(is.na(author_col), "NOT FOUND", author_col))
cells_keep <- rownames(obj[[]])[ as.character(obj$broad_cluster_step1) %in% target_broad ]
if (length(cells_keep) == 0) {
  stop("0/17/4/8/2 没有提取到任何细胞。")
}
myeloid <- subset(obj, cells = cells_keep)
myeloid$broad_origin <- factor( as.character(myeloid$broad_cluster_step1), levels = c("0", "17", "4", "8", "2") )
message("Step 2 starting pool:")
print(table(myeloid$broad_origin))
message("Total cells = ", ncol(myeloid))
write.csv( as.data.frame(table(myeloid$broad_origin)), file.path(tab_dir, "00_step2_starting_pool_broad_cluster_counts.csv"), row.names = FALSE )
DefaultAssay(myeloid) <- "RNA"
VariableFeatures(myeloid) <- character(0)
myeloid <- NormalizeData( myeloid, normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE )
myeloid <- FindVariableFeatures( myeloid, selection.method = "vst", nfeatures = 3000, verbose = TRUE )
myeloid <- ScaleData( myeloid, features = VariableFeatures(myeloid), verbose = TRUE )
myeloid <- RunPCA( myeloid, features = VariableFeatures(myeloid), npcs = 40, verbose = TRUE )
pca_stdev <- Stdev(myeloid, reduction = "pca")
pca_var <- pca_stdev^2 / sum(pca_stdev^2)
pca_cum <- cumsum(pca_var)
pca_table <- data.frame( PC = seq_along(pca_var), variance = pca_var, cumulative_variance = pca_cum )
write.csv( pca_table, file.path(tab_dir, "01_myeloid_PCA_variance.csv"), row.names = FALSE )
npcs_use <- min(20, ncol(Embeddings(myeloid, "pca")))
dims_use <- seq_len(npcs_use)
message("RunHarmony by: ", sample_col)
myeloid <- tryCatch( {
    RunHarmony( object = myeloid, group.by.vars = sample_col, reduction.use = "pca", dims.use = dims_use, theta = 1, max.iter.harmony = 20, verbose = TRUE )
  },
  error = function(e1) {
    message("Harmony signature 1 failed: ", e1$message)
    message("Trying alternate Harmony signature...")
    RunHarmony( object = myeloid, group.by.vars = sample_col, reduction = "pca", dims.use = dims_use, theta = 1, max.iter.harmony = 20, verbose = TRUE )
  }
)
if (!"harmony" %in% Reductions(myeloid)) {
  stop("Harmony reduction 未生成。")
}
myeloid <- RunUMAP( myeloid, reduction = "harmony", dims = dims_use, reduction.name = "umap.myeloid", reduction.key = "myUMAP_", verbose = TRUE )
myeloid <- FindNeighbors( myeloid, reduction = "harmony", dims = dims_use, verbose = TRUE )
resolutions <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0)
for (res in resolutions) {
  myeloid <- FindClusters( myeloid, resolution = res, algorithm = 1, random.seed = 1234, verbose = FALSE )
}
res_cols <- grep("_res\\.", colnames(myeloid[[]]), value = TRUE)
res_cols <- res_cols[grepl( paste0("(", paste(gsub("\\.", "\\\\.", resolutions), collapse = "|"), ")$"), res_cols )]
final_res_col_candidates <- grep("res\\.0\\.4$", colnames(myeloid[[]]), value = TRUE)
if (length(final_res_col_candidates) == 0) {
  stop("未找到 resolution 0.4 的聚类列。")
}
final_res_col <- tail(final_res_col_candidates, 1)
myeloid$myeloid_cluster_r04 <- as.character(myeloid[[final_res_col]][, 1])
Idents(myeloid) <- "myeloid_cluster_r04"
message("Myeloid-only clusters at resolution 0.4:")
print(table(Idents(myeloid)))
resident_genes <- c( "P2ry12", "Tmem119", "Sall1", "Fcrls", "Hexb", "Cx3cr1", "Siglech", "Gpr34", "P2ry13" )
macrophage_genes <- c( "Ccr2", "Ly6c2", "Ms4a7", "Fn1", "Fcgr1", "S100a8", "S100a9", "Chil3", "Clec4n", "F13a1", "Arg1" )
cycling_genes <- c( "Mki67", "Top2a", "Ube2c", "Birc5", "Cenpf", "Ccna2", "Aurkb" )
all_rna_genes <- rownames(myeloid[["RNA"]])
resident_use <- intersect(resident_genes, all_rna_genes)
macrophage_use <- intersect(macrophage_genes, all_rna_genes)
cycling_use <- intersect(cycling_genes, all_rna_genes)
if (length(resident_use) < 4) {
  warning("Resident marker 命中少于 4 个，Check required")
}
if (length(macrophage_use) < 4) {
  warning("Macrophage marker 命中少于 4 个，Check required")
}
myeloid <- AddModuleScore( myeloid, features = list(resident_use), name = "ResidentMG", assay = "RNA", search = FALSE )
myeloid <- AddModuleScore( myeloid, features = list(macrophage_use), name = "InfiltratingMacro", assay = "RNA", search = FALSE )
myeloid <- AddModuleScore( myeloid, features = list(cycling_use), name = "Cycling", assay = "RNA", search = FALSE )
resident_score_col <- grep("^ResidentMG[0-9]+$", colnames(myeloid[[]]), value = TRUE)[1]
macro_score_col <- grep("^InfiltratingMacro[0-9]+$", colnames(myeloid[[]]), value = TRUE)[1]
cycling_score_col <- grep("^Cycling[0-9]+$", colnames(myeloid[[]]), value = TRUE)[1]
p_cluster <- DimPlot( myeloid, reduction = "umap.myeloid", group.by = "myeloid_cluster_r04", label = TRUE, repel = TRUE ) + ggtitle("Myeloid-only reclustering | resolution 0.4")
p_origin <- DimPlot( myeloid, reduction = "umap.myeloid", group.by = "broad_origin", label = TRUE, repel = TRUE ) + ggtitle("Origin from Step 1 broad clusters")
p_sample <- DimPlot( myeloid, reduction = "umap.myeloid", group.by = sample_col ) + ggtitle("Biological sample")
if (!is.na(time_col)) {
  p_time <- DimPlot( myeloid, reduction = "umap.myeloid", group.by = time_col ) + ggtitle("SCI timepoint")
  p_umap_all <- (p_cluster | p_origin) / (p_sample | p_time)
} else {
  p_umap_all <- (p_cluster | p_origin) / p_sample
}
ggsave( file.path(fig_dir, "01_myeloid_refinement_UMAP.pdf"), p_umap_all, width = 15, height = 11 )
clustree_prefix_candidates <- unique(sub("res\\.[0-9.]+$", "res.", res_cols))
if (length(clustree_prefix_candidates) > 0) {
  clustree_prefix <- clustree_prefix_candidates[1]
  p_tree <- tryCatch( clustree(myeloid[[]], prefix = clustree_prefix), error = function(e) {
      message("clustree failed: ", e$message)
      NULL
    }
  )
  if (!is.null(p_tree)) {
    ggsave( file.path(fig_dir, "02_myeloid_clustree_resolution_sensitivity.pdf"), p_tree, width = 11, height = 8 )
  }
}
marker_order <- unique(c(resident_use, macrophage_use, cycling_use))
p_dot <- DotPlot( myeloid, features = marker_order, group.by = "myeloid_cluster_r04", assay = "RNA" ) + RotatedAxis() + ggtitle("Resident microglia vs infiltrating macrophage vs cycling markers")
ggsave( file.path(fig_dir, "03_resident_vs_macrophage_marker_DotPlot.pdf"), p_dot, width = 16, height = 7 )
score_features <- c(resident_score_col, macro_score_col, cycling_score_col)
p_scores <- FeaturePlot( myeloid, features = score_features, reduction = "umap.myeloid", ncol = 3, order = TRUE )
ggsave( file.path(fig_dir, "04_lineage_program_score_UMAP.pdf"), p_scores, width = 15, height = 5 )
DefaultAssay(myeloid) <- "RNA"
Idents(myeloid) <- "myeloid_cluster_r04"
markers <- FindAllMarkers( myeloid, assay = "RNA", only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25, test.use = "wilcox", verbose = TRUE )
write.csv( markers, file.path(tab_dir, "02_myeloid_cluster_all_markers_r04.csv"), row.names = FALSE )
fc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(markers))[1]
if (!is.na(fc_col)) {
  top20 <- markers %>% group_by(cluster) %>% arrange(desc(.data[[fc_col]]), .by_group = TRUE) %>% slice_head(n = 20) %>% ungroup()
} else {
  top20 <- markers %>% group_by(cluster) %>% slice_head(n = 20) %>% ungroup()
}
write.csv( top20, file.path(tab_dir, "03_top20_markers_per_myeloid_cluster_r04.csv"), row.names = FALSE )
score_meta <- myeloid[[]] %>% mutate(myeloid_cluster_r04 = as.character(myeloid_cluster_r04))
cluster_scores_tbl <- score_meta %>% group_by(myeloid_cluster_r04) %>% summarise( n_cells = n(), resident_score_mean = mean(.data[[resident_score_col]], na.rm = TRUE), resident_score_median = median(.data[[resident_score_col]], na.rm = TRUE), macrophage_score_mean = mean(.data[[macro_score_col]], na.rm = TRUE), macrophage_score_median = median(.data[[macro_score_col]], na.rm = TRUE), cycling_score_mean = mean(.data[[cycling_score_col]], na.rm = TRUE), cycling_score_median = median(.data[[cycling_score_col]], na.rm = TRUE), resident_minus_macro = resident_score_mean - macrophage_score_mean, .groups = "drop" )
write.csv( cluster_scores_tbl, file.path(tab_dir, "04_cluster_lineage_program_scores_r04.csv"), row.names = FALSE )
hm <- cluster_scores_tbl %>% select( myeloid_cluster_r04, resident_score_mean, macrophage_score_mean, cycling_score_mean )
hm_mat <- as.matrix(hm[, -1, drop = FALSE])
rownames(hm_mat) <- hm$myeloid_cluster_r04
hm_mat <- t(hm_mat)
pdf( file.path(fig_dir, "05_cluster_lineage_score_heatmap.pdf"), width = max(7, 0.65 * nrow(cluster_scores_tbl) + 4), height = 4.5 )
pheatmap( hm_mat, scale = "row", cluster_rows = FALSE, cluster_cols = TRUE, main = "Lineage program scores by myeloid-only cluster" )
dev.off()
cross_origin <- as.data.frame.matrix( table(myeloid$myeloid_cluster_r04, myeloid$broad_origin) )
cross_origin$myeloid_cluster_r04 <- rownames(cross_origin)
rownames(cross_origin) <- NULL
cross_sample <- as.data.frame.matrix( table(myeloid$myeloid_cluster_r04, myeloid[[sample_col]][, 1]) )
cross_sample$myeloid_cluster_r04 <- rownames(cross_sample)
rownames(cross_sample) <- NULL
cross_all <- list( broad_origin = cross_origin, sample = cross_sample )
if (!is.na(time_col)) {
  cross_time <- as.data.frame.matrix( table(myeloid$myeloid_cluster_r04, myeloid[[time_col]][, 1]) )
  cross_time$myeloid_cluster_r04 <- rownames(cross_time)
  rownames(cross_time) <- NULL
  cross_all$timepoint <- cross_time
}
if (!is.na(author_col)) {
  cross_author <- as.data.frame.matrix( table(myeloid$myeloid_cluster_r04, myeloid[[author_col]][, 1]) )
  cross_author$myeloid_cluster_r04 <- rownames(cross_author)
  rownames(cross_author) <- NULL
  cross_all$author_annotation <- cross_author
}
long_cross <- list()
long_cross[["broad_origin"]] <- as.data.frame( table( myeloid_cluster_r04 = myeloid$myeloid_cluster_r04, category = myeloid$broad_origin ) ) %>% mutate(cross_type = "broad_origin")
long_cross[["sample"]] <- as.data.frame( table( myeloid_cluster_r04 = myeloid$myeloid_cluster_r04, category = myeloid[[sample_col]][, 1] ) ) %>% mutate(cross_type = "sample")
if (!is.na(time_col)) {
  long_cross[["timepoint"]] <- as.data.frame( table( myeloid_cluster_r04 = myeloid$myeloid_cluster_r04, category = myeloid[[time_col]][, 1] ) ) %>% mutate(cross_type = "timepoint")
}
if (!is.na(author_col)) {
  long_cross[["author_annotation"]] <- as.data.frame( table( myeloid_cluster_r04 = myeloid$myeloid_cluster_r04, category = myeloid[[author_col]][, 1] ) ) %>% mutate(cross_type = "author_annotation")
}
cross_long <- bind_rows(long_cross) %>% group_by(cross_type, myeloid_cluster_r04) %>% mutate( cluster_total = sum(Freq), proportion_within_cluster = ifelse(cluster_total > 0, Freq / cluster_total, NA_real_) ) %>% ungroup()
write.csv( cross_long, file.path(tab_dir, "05_cluster_crosschecks_origin_sample_time_author.csv"), row.names = FALSE )
classification_template <- cluster_scores_tbl %>% mutate( dominant_score = case_when( resident_score_mean > macrophage_score_mean ~ "resident-biased", macrophage_score_mean > resident_score_mean ~ "macrophage-biased", TRUE ~ "balanced" ), final_lineage_class = NA_character_, reviewer_note = NA_character_ )
write.csv( classification_template, file.path(tab_dir, "06_lineage_classification_review_template.csv"), row.names = FALSE )
saveRDS( myeloid, file.path( obj_dir, "GSE162610_step2A_myeloid_refined_multires.rds" ) )
report_file <- file.path(out_dir, "07_step2A_run_report.txt")
sink(report_file)
cat("GSE162610 Step 2A | myeloid-only refinement\n")
cat("====================================================\n\n")
cat("Input RDS:\n", input_rds, "\n\n", sep = "")
cat("Broad clusters included: 0, 17, 4, 8, 2\n")
cat("Explicitly excluded from Step 2 pool: 6, 16, 9\n\n")
cat("Total starting cells: ", ncol(myeloid), "\n", sep = "")
cat("\nBroad-origin counts:\n")
print(table(myeloid$broad_origin))
cat("\nSample column: ", sample_col, "\n", sep = "")
cat("Time column: ", ifelse(is.na(time_col), "NOT FOUND", time_col), "\n", sep = "")
cat("Author annotation column: ", ifelse(is.na(author_col), "NOT FOUND", author_col), "\n", sep = "")
cat("\nPCA dims used for Harmony: 1:", npcs_use, "\n", sep = "")
cat("Harmony grouping variable: ", sample_col, "\n", sep = "")
cat("Harmony theta: 1\n")
cat("Harmony max.iter: 20\n")
cat("\nResolutions tested:\n")
print(resolutions)
cat("\nResolution 0.4 cluster counts:\n")
print(table(myeloid$myeloid_cluster_r04))
cat("\nResident markers available:\n")
print(resident_use)
cat("\nMacrophage markers available:\n")
print(macrophage_use)
cat("\nCycling markers available:\n")
print(cycling_use)
cat("\nCluster-level lineage score summary:\n")
print(cluster_scores_tbl)
cat("\nIMPORTANT:\n")
cat("- Step 2A does NOT define DAM.\n")
cat("- Step 2A does NOT define repair/scar.\n")
cat("- Step 2A only refines resident-derived microglia vs ambiguous myeloid vs infiltrating monocyte/macrophage.\n")
cat("- Borderline clusters should remain ambiguous until marker + score + origin + sample/time + author annotation evidence is reviewed together.\n")
cat("- AUCell/ssGSEA can be added in Step 2B only if the lineage boundary remains uncertain.\n")
sink()
message("====================================================")
message("STEP 2A COMPLETE")
message("Output folder: ", out_dir)
message("Object folder: ", obj_dir)
message("")
message("Key outputs:")
message("1) figures/01_myeloid_refinement_UMAP.pdf")
message("2) figures/02_myeloid_clustree_resolution_sensitivity.pdf")
message("3) figures/03_resident_vs_macrophage_marker_DotPlot.pdf")
message("4) figures/05_cluster_lineage_score_heatmap.pdf")
message("5) tables/03_top20_markers_per_myeloid_cluster_r04.csv")
message("6) tables/05_cluster_crosschecks_origin_sample_time_author.csv")
message("7) 07_step2A_run_report.txt")
message("====================================================")
