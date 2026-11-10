options(stringsAsFactors = FALSE)
options(timeout = 1200)
set.seed(20260730)
project_dir <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
in_rds <- file.path( project_dir, "02_objects", "step1C_GSE162610_scDblFinder", "GSE162610_QC_singlets.rds" )
out_dir <- file.path( project_dir, "02_objects", "step1D_GSE162610_LogNorm_PCA_preHarmony" )
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "00_run_log.txt")
writeLines( paste0("GSE162610 Step 1D started: ", Sys.time()), log_file )
if (!file.exists(in_rds)) {
  stop("找不到 Step 1C singlet 对象：\n", in_rds)
}
pkgs <- c( "Seurat", "data.table", "ggplot2", "patchwork" )
missing_pkgs <- pkgs[ !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE) ]
if (length(missing_pkgs) > 0) {
  stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))
}
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
pbmc <- readRDS(in_rds)
required_meta <- c( "sample_id", "timepoint", "nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB" )
if (!all(required_meta %in% colnames(pbmc@meta.data))) {
  stop( "对象缺少 metadata：", paste( setdiff(required_meta, colnames(pbmc@meta.data)), collapse = ", " ) )
}
DefaultAssay(pbmc) <- "RNA"
t0 <- Sys.time()
pbmc <- NormalizeData( pbmc, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE )
pbmc <- FindVariableFeatures( pbmc, selection.method = "vst", nfeatures = 3000, verbose = FALSE )
hvg <- VariableFeatures(pbmc)
if (length(hvg) != 3000) {
  warning( "实际识别到的 HVG 数为 ", length(hvg), "，不是预期的 3000。" )
}
writeLines( hvg, file.path( tab_dir, "01_HVG_3000_genes.txt" ) )
vf_info <- HVFInfo(pbmc)
vf_info$gene <- rownames(vf_info)
vf_info$is_HVG <- vf_info$gene %in% hvg
fwrite( vf_info, file.path( tab_dir, "02_variable_feature_statistics.csv" ) )
pbmc <- ScaleData( pbmc, features = hvg, verbose = FALSE )
pbmc <- RunPCA( pbmc, features = hvg, npcs = 50, verbose = FALSE )
if (!"pca" %in% Reductions(pbmc)) {
  stop("PCA reduction 未Completed生成。")
}
stdev <- Stdev(pbmc, reduction = "pca")
variance_pct <- 100 * (stdev^2) / sum(stdev^2)
cumulative_pct <- cumsum(variance_pct)
pca_variance <- data.frame( PC = seq_along(stdev), stdev = stdev, variance_pct = variance_pct, cumulative_variance_pct = cumulative_pct )
fwrite( pca_variance, file.path( tab_dir, "03_PCA_variance_explained.csv" ) )
pc80 <- which(cumulative_pct >= 80)[1]
if (is.na(pc80)) {
  pc80 <- 30
}
candidate_pc <- max(20, min(40, pc80))
writeLines( c( paste0("PC reaching 80% cumulative variance: ", pc80), paste0("Candidate PCs for Harmony/UMAP: 1:", candidate_pc), "", "Status:", "This is only an automatic reference.", "Final PC number should be decided after reviewing ElbowPlot and biological separation." ), file.path( tab_dir, "04_candidate_PC_number.txt" ) )
p_elbow <- ElbowPlot( pbmc, ndims = 50 ) + ggtitle( "GSE162610 PCA elbow plot (pre-Harmony)" )
ggsave( file.path( fig_dir, "01_PCA_ElbowPlot_50PCs.pdf" ), plot = p_elbow, width = 7, height = 5 )
p_pca_sample <- DimPlot( pbmc, reduction = "pca", group.by = "sample_id", dims = c(1, 2), pt.size = 0.15 ) + ggtitle( "PCA before Harmony - colored by sample" )
ggsave( file.path( fig_dir, "02_PCA_PC1_PC2_by_sample.pdf" ), plot = p_pca_sample, width = 10, height = 7 )
p_pca_time <- DimPlot( pbmc, reduction = "pca", group.by = "timepoint", dims = c(1, 2), pt.size = 0.15 ) + ggtitle( "PCA before Harmony - colored by timepoint" )
ggsave( file.path( fig_dir, "03_PCA_PC1_PC2_by_timepoint.pdf" ), plot = p_pca_time, width = 9, height = 7 )
dims_preHarmony <- 1:candidate_pc
pbmc <- RunUMAP( pbmc, reduction = "pca", dims = dims_preHarmony, reduction.name = "umap.preharmony", reduction.key = "UMAPPRE_", seed.use = 20260730, verbose = FALSE )
p_umap_sample <- DimPlot( pbmc, reduction = "umap.preharmony", group.by = "sample_id", pt.size = 0.10 ) + ggtitle( "Pre-Harmony UMAP - sample_id" )
ggsave( file.path( fig_dir, "04_preHarmony_UMAP_by_sample.pdf" ), plot = p_umap_sample, width = 10, height = 7 )
p_umap_time <- DimPlot( pbmc, reduction = "umap.preharmony", group.by = "timepoint", pt.size = 0.10 ) + ggtitle( "Pre-Harmony UMAP - timepoint" )
ggsave( file.path( fig_dir, "05_preHarmony_UMAP_by_timepoint.pdf" ), plot = p_umap_time, width = 9, height = 7 )
loadings <- Loadings( pbmc, reduction = "pca" )
n_top <- 30
pcs_to_report <- 1:min(20, ncol(loadings))
top_loadings_list <- lapply( pcs_to_report, function(pc) {
    vals <- loadings[, pc]
    ord <- order( abs(vals), decreasing = TRUE )[1:min(n_top, length(vals))]
    data.frame( PC = paste0("PC", pc), gene = rownames(loadings)[ord], loading = vals[ord], abs_loading = abs(vals[ord]), stringsAsFactors = FALSE )
  }
)
top_loadings <- do.call( rbind, top_loadings_list )
fwrite( top_loadings, file.path( tab_dir, "05_top30_gene_loadings_PC1_to_PC20.csv" ) )
sample_counts <- as.data.frame( table(pbmc$sample_id), stringsAsFactors = FALSE )
colnames(sample_counts) <- c( "sample_id", "n_cells" )
fwrite( sample_counts, file.path( tab_dir, "06_cell_counts_by_sample_preHarmony.csv" ) )
time_counts <- as.data.frame( table(pbmc$timepoint), stringsAsFactors = FALSE )
colnames(time_counts) <- c( "timepoint", "n_cells" )
fwrite( time_counts, file.path( tab_dir, "07_cell_counts_by_timepoint_preHarmony.csv" ) )
saveRDS( pbmc, file.path( out_dir, "GSE162610_LogNorm_PCA_preHarmony.rds" ), compress = FALSE )
all_samples_preserved <- length( unique(as.character(pbmc$sample_id)) ) == 10
all_timepoints_preserved <- length( unique(as.character(pbmc$timepoint)) ) == 4
pca_ok <- ( "pca" %in% Reductions(pbmc) && ncol(Embeddings(pbmc, "pca")) >= 30 )
umap_ok <- ( "umap.preharmony" %in% Reductions(pbmc) )
step_pass <- all( all_samples_preserved, all_timepoints_preserved, length(hvg) >= 2000, pca_ok, umap_ok )
runtime_min <- as.numeric( difftime( Sys.time(), t0, units = "mins" ) )
if (step_pass) {
} else {
}
