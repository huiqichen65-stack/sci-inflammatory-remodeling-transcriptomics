options(stringsAsFactors = FALSE)
set.seed(1234)
required_pkgs <- c( "Seurat", "SeuratObject", "dplyr", "ggplot2", "patchwork", "harmony", "clustree", "pheatmap", "Matrix" )
missing_pkgs <- required_pkgs[ !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE) ]
if (length(missing_pkgs) > 0) {
  stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))
}
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(harmony)
  library(clustree)
  library(pheatmap)
  library(Matrix)
})
project_dir <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
input_rds <- file.path( project_dir, "04_results", "step2C_GSE162610_final_microglia_extraction", "GSE162610_microglia_clean.rds" )
out_dir <- file.path( project_dir, "04_results", "step3A_GSE162610_microglia_reclustering" )
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_rds)) {
  stop( "找不到 Step 3A 正式输入对象：\n", input_rds )
}
obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) {
  stop("输入 RDS 不是 Seurat 对象。")
}
if (!"RNA" %in% names(obj@assays)) {
  stop("输入对象中没有 RNA assay。")
}
DefaultAssay(obj) <- "RNA"
expected_n <- 20747
if (ncol(obj) != expected_n) {
  warning( "Step 3A 输入细胞数与 Step 2C 预期不一致：", " expected = ", expected_n, ", observed = ", ncol(obj) )
}
meta <- obj[[]]
required_meta <- c( "sample_id", "timepoint", "step2_final_cluster" )
missing_meta <- required_meta[ !required_meta %in% colnames(meta) ]
if (length(missing_meta) > 0) {
  stop( "缺少必要 metadata：", paste(missing_meta, collapse = ", ") )
}
obj$step2_final_cluster <- as.character(obj$step2_final_cluster)
print(table(obj$sample_id))
print(table(obj$timepoint))
print(table(obj$step2_final_cluster))
joinlayers_available <- exists( "JoinLayers", where = asNamespace("SeuratObject"), mode = "function", inherits = FALSE )
if (isTRUE(joinlayers_available)) {
  obj <- tryCatch( {
      SeuratObject::JoinLayers( obj, assay = "RNA" )
    },
    error = function(e) {
      obj
    }
  )
}
DefaultAssay(obj) <- "RNA"
old_reductions <- Reductions(obj)
writeLines( old_reductions, file.path( tab_dir, "00_input_old_reductions.txt" ) )
obj <- NormalizeData( obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE )
obj <- FindVariableFeatures( obj, assay = "RNA", selection.method = "vst", nfeatures = 3000, verbose = TRUE )
hvg <- VariableFeatures(obj)
writeLines( hvg, file.path( tab_dir, "01_microglia_HVG_3000.txt" ) )
p_vf <- VariableFeaturePlot(obj)
top10_hvg <- head(hvg, 10)
p_vf_label <- LabelPoints( plot = p_vf, points = top10_hvg, repel = TRUE )
ggsave( file.path( fig_dir, "01_microglia_variable_features.pdf" ), p_vf_label, width = 8, height = 6 )
obj <- ScaleData( obj, assay = "RNA", features = hvg, verbose = TRUE )
obj <- RunPCA( obj, assay = "RNA", features = hvg, npcs = 50, reduction.name = "pca.step3A", reduction.key = "PC3A_", verbose = TRUE )
pca_stdev <- Stdev( obj, reduction = "pca.step3A" )
pca_var <- pca_stdev^2 / sum(pca_stdev^2)
pca_cum <- cumsum(pca_var)
pca_tbl <- data.frame( PC = seq_along(pca_var), variance_fraction = pca_var, cumulative_variance = pca_cum, stringsAsFactors = FALSE )
write.csv( pca_tbl, file.path( tab_dir, "02_microglia_PCA_variance.csv" ), row.names = FALSE )
pc80 <- which(pca_cum >= 0.80)[1]
pc90 <- which(pca_cum >= 0.90)[1]
pc_summary <- data.frame( criterion = c( "80_percent_cumulative_variance", "90_percent_cumulative_variance", "formal_Step3A_dims" ), PC = c( pc80, pc90, 30 ) )
write.csv( pc_summary, file.path( tab_dir, "03_candidate_PC_summary.csv" ), row.names = FALSE )
dims_use <- 1:30
pdf( file.path( fig_dir, "02_microglia_PCA_ElbowPlot_50PCs.pdf" ), width = 7, height = 6 )
print( ElbowPlot( obj, ndims = 50, reduction = "pca.step3A" ) )
dev.off()
p_pca_sample <- DimPlot( obj, reduction = "pca.step3A", group.by = "sample_id" ) + ggtitle("Microglia PCA | biological sample")
p_pca_time <- DimPlot( obj, reduction = "pca.step3A", group.by = "timepoint" ) + ggtitle("Microglia PCA | SCI timepoint")
ggsave( file.path( fig_dir, "03_microglia_PCA_sample_timepoint.pdf" ), p_pca_sample | p_pca_time, width = 14, height = 6 )
obj <- RunUMAP( obj, reduction = "pca.step3A", dims = dims_use, reduction.name = "umap.preHarmony.step3A", reduction.key = "preH3A_", seed.use = 1234, verbose = TRUE )
p_pre_sample <- DimPlot( obj, reduction = "umap.preHarmony.step3A", group.by = "sample_id" ) + ggtitle("Pre-Harmony | sample")
p_pre_time <- DimPlot( obj, reduction = "umap.preHarmony.step3A", group.by = "timepoint" ) + ggtitle("Pre-Harmony | timepoint")
ggsave( file.path( fig_dir, "04_preHarmony_microglia_UMAP.pdf" ), p_pre_sample | p_pre_time, width = 14, height = 6 )
obj <- tryCatch( {
    RunHarmony( object = obj, group.by.vars = "sample_id", reduction.use = "pca.step3A", dims.use = dims_use, theta = 1, max.iter.harmony = 20, reduction.save = "harmony.step3A", verbose = TRUE )
  },
  error = function(e1) {
    RunHarmony( object = obj, group.by.vars = "sample_id", reduction = "pca.step3A", dims.use = dims_use, theta = 1, max.iter.harmony = 20, reduction.save = "harmony.step3A", verbose = TRUE )
  }
)
if (!"harmony.step3A" %in% Reductions(obj)) {
  if ("harmony" %in% Reductions(obj)) {
    obj[["harmony.step3A"]] <- obj[["harmony"]]
  } else {
    stop("Harmony reduction 未Completed生成。")
  }
}
obj <- RunUMAP( obj, reduction = "harmony.step3A", dims = dims_use, reduction.name = "umap.microglia.step3A", reduction.key = "mg3A_", seed.use = 1234, verbose = TRUE )
p_post_sample <- DimPlot( obj, reduction = "umap.microglia.step3A", group.by = "sample_id" ) + ggtitle("Post-Harmony | biological sample")
p_post_time <- DimPlot( obj, reduction = "umap.microglia.step3A", group.by = "timepoint" ) + ggtitle("Post-Harmony | SCI timepoint")
p_post_origin <- DimPlot( obj, reduction = "umap.microglia.step3A", group.by = "step2_final_cluster", label = TRUE, repel = TRUE ) + ggtitle("Step 2 lineage provenance")
ggsave( file.path( fig_dir, "05_postHarmony_microglia_UMAP_overview.pdf" ), (p_post_sample | p_post_time) / p_post_origin, width = 15, height = 11 )
ggsave( file.path( fig_dir, "06_pre_vs_post_Harmony_microglia.pdf" ), (p_pre_sample | p_post_sample) / (p_pre_time | p_post_time), width = 15, height = 11 )
obj <- FindNeighbors( obj, reduction = "harmony.step3A", dims = dims_use, graph.name = c( "microglia_step3A_nn", "microglia_step3A_snn" ), verbose = TRUE )
resolutions <- seq( 0.1, 1.0, by = 0.1 )
cluster_number_list <- list()
for (res in resolutions) {
  cl_name <- paste0( "step3A_res.", format( res, nsmall = 1, trim = TRUE, scientific = FALSE ) )
  obj <- FindClusters( obj, graph.name = "microglia_step3A_snn", resolution = res, algorithm = 1, random.seed = 1234, cluster.name = cl_name, verbose = FALSE )
  cluster_number_list[[cl_name]] <- data.frame( resolution = res, cluster_column = cl_name, n_clusters = length( unique( as.character( obj[[cl_name]][, 1] ) ) ) )
}
cluster_number_tbl <- bind_rows( cluster_number_list )
write.csv( cluster_number_tbl, file.path( tab_dir, "04_cluster_number_by_resolution.csv" ), row.names = FALSE )
tree_meta <- obj[[]][ , paste0( "step3A_res.", format( resolutions, nsmall = 1, trim = TRUE, scientific = FALSE ) ), drop = FALSE ]
p_tree <- clustree( tree_meta, prefix = "step3A_res." )
ggsave( file.path( fig_dir, "07_microglia_clustree_resolution_0.1_to_1.0.pdf" ), p_tree, width = 12, height = 9 )
anchor_col <- "step3A_res.0.4"
if (!anchor_col %in% colnames(obj[[]])) {
  stop( "未找到 provisional resolution 0.4 metadata。" )
}
obj$step3A_cluster_r04 <- as.character( obj[[anchor_col]][, 1] )
Idents(obj) <- "step3A_cluster_r04"
print(table(obj$step3A_cluster_r04))
p_cluster <- DimPlot( obj, reduction = "umap.microglia.step3A", group.by = "step3A_cluster_r04", label = TRUE, repel = TRUE ) + ggtitle( "Clean microglia de novo reclustering | provisional resolution 0.4" )
p_cluster_time <- DimPlot( obj, reduction = "umap.microglia.step3A", group.by = "timepoint" ) + ggtitle("SCI timepoint")
p_cluster_sample <- DimPlot( obj, reduction = "umap.microglia.step3A", group.by = "sample_id" ) + ggtitle("Biological sample")
p_cluster_origin <- DimPlot( obj, reduction = "umap.microglia.step3A", group.by = "step2_final_cluster", label = TRUE, repel = TRUE ) + ggtitle("Step 2 cluster provenance")
ggsave( file.path( fig_dir, "08_provisional_r04_microglia_UMAP.pdf" ), (p_cluster | p_cluster_time) / (p_cluster_sample | p_cluster_origin), width = 15, height = 11 )
marker_programs <- list( Homeostatic = c( "P2ry12", "Tmem119", "Sall1", "Fcrls", "Hexb", "Gpr34", "P2ry13", "Siglech" ), Activated_injury = c( "Apoe", "Lpl", "Cst7", "Cd68", "Ctsb", "Ctsd", "Cd9", "Tyrobp" ), DAM_like_candidate = c( "Trem2", "Apoe", "Lpl", "Cst7", "Axl", "Clec7a", "Itgax", "Gpnmb" ), Phagocytic_lipid = c( "Lpl", "Apoe", "Cd68", "Ctsd", "Lgals3", "Fabp5", "Cst7" ), IFN_response = c( "Ifit1", "Ifit2", "Ifit3", "Isg15", "Mx1", "Rsad2", "Usp18", "Oasl1" ), Cycling_S = c( "Mcm2", "Mcm5", "Mcm6", "Pcna", "Ccne2", "Uhrf1", "Cdc6" ), Cycling_G2M = c( "Mki67", "Top2a", "Ube2c", "Birc5", "Cenpf", "Ccna2", "Aurkb" ) )
rna_genes <- rownames(obj[["RNA"]])
marker_programs_use <- lapply( marker_programs, function(gs) intersect(gs, rna_genes) )
marker_availability <- do.call( rbind, lapply( names(marker_programs), function(nm) {
      data.frame( program = nm, n_requested = length(marker_programs[[nm]]), n_available = length(marker_programs_use[[nm]]), available_genes = paste( marker_programs_use[[nm]], collapse = ";" ), stringsAsFactors = FALSE )
    }
  )
)
write.csv( marker_availability, file.path( tab_dir, "05_state_marker_availability.csv" ), row.names = FALSE )
dot_features <- unique( unlist( marker_programs_use, use.names = FALSE ) )
p_dot <- DotPlot( obj, features = dot_features, group.by = "step3A_cluster_r04", assay = "RNA" ) + RotatedAxis() + ggtitle( "Exploratory microglial state markers | provisional r=0.4" )
ggsave( file.path( fig_dir, "09_exploratory_state_marker_DotPlot_r04.pdf" ), p_dot, width = 20, height = 8 )
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "step3A_cluster_r04"
markers <- FindAllMarkers( obj, assay = "RNA", only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25, test.use = "wilcox", verbose = TRUE )
write.csv( markers, file.path( tab_dir, "06_all_cluster_markers_provisional_r04.csv" ), row.names = FALSE )
fc_col <- intersect( c( "avg_log2FC", "avg_logFC" ), colnames(markers) )[1]
if (!is.na(fc_col)) {
  top20 <- markers %>% group_by(cluster) %>% arrange( desc(.data[[fc_col]]), .by_group = TRUE ) %>% slice_head(n = 20) %>% ungroup()
} else {
  top20 <- markers %>% group_by(cluster) %>% slice_head(n = 20) %>% ungroup()
}
write.csv( top20, file.path( tab_dir, "07_top20_markers_per_cluster_provisional_r04.csv" ), row.names = FALSE )
cluster_time_counts <- as.data.frame( table( cluster = obj$step3A_cluster_r04, timepoint = obj$timepoint ) )
cluster_time_counts <- cluster_time_counts %>% group_by(cluster) %>% mutate( cluster_total = sum(Freq), proportion_within_cluster = ifelse( cluster_total > 0, Freq / cluster_total, NA_real_ ) ) %>% ungroup()
write.csv( cluster_time_counts, file.path( tab_dir, "08_cluster_by_timepoint_provisional_r04.csv" ), row.names = FALSE )
time_mat <- xtabs( Freq ~ cluster + timepoint, data = cluster_time_counts )
time_prop <- prop.table( time_mat, margin = 1 )
pdf( file.path( fig_dir, "10_cluster_timepoint_heatmap_provisional_r04.pdf" ), width = 8, height = max( 6, 0.45 * nrow(time_prop) + 3 ) )
pheatmap( time_prop, cluster_rows = FALSE, cluster_cols = FALSE, main = "Provisional microglia clusters × timepoint | row proportion" )
dev.off()
cluster_sample_counts <- as.data.frame( table( cluster = obj$step3A_cluster_r04, sample_id = obj$sample_id ) )
cluster_sample_counts <- cluster_sample_counts %>% group_by(cluster) %>% mutate( cluster_total = sum(Freq), proportion_within_cluster = ifelse( cluster_total > 0, Freq / cluster_total, NA_real_ ) ) %>% ungroup()
write.csv( cluster_sample_counts, file.path( tab_dir, "09_cluster_by_sample_provisional_r04.csv" ), row.names = FALSE )
sample_mat <- xtabs( Freq ~ cluster + sample_id, data = cluster_sample_counts )
sample_prop <- prop.table( sample_mat, margin = 1 )
pdf( file.path( fig_dir, "11_cluster_sample_heatmap_provisional_r04.pdf" ), width = 11, height = max( 6, 0.45 * nrow(sample_prop) + 3 ) )
pheatmap( sample_prop, cluster_rows = FALSE, cluster_cols = FALSE, main = "Provisional microglia clusters × biological sample | row proportion" )
dev.off()
cluster_origin_counts <- as.data.frame( table( step3A_cluster = obj$step3A_cluster_r04, step2_cluster = obj$step2_final_cluster ) )
cluster_origin_counts <- cluster_origin_counts %>% group_by(step3A_cluster) %>% mutate( cluster_total = sum(Freq), proportion_within_step3A_cluster = ifelse( cluster_total > 0, Freq / cluster_total, NA_real_ ) ) %>% ungroup()
write.csv( cluster_origin_counts, file.path( tab_dir, "10_step3A_cluster_vs_step2_provenance.csv" ), row.names = FALSE )
origin_mat <- xtabs( Freq ~ step3A_cluster + step2_cluster, data = cluster_origin_counts )
origin_prop <- prop.table( origin_mat, margin = 1 )
pdf( file.path( fig_dir, "12_step3A_cluster_vs_step2_provenance_heatmap.pdf" ), width = 8, height = max( 6, 0.45 * nrow(origin_prop) + 3 ) )
pheatmap( origin_prop, cluster_rows = FALSE, cluster_cols = FALSE, main = "New Step 3A cluster × Step 2 provenance | row proportion" )
dev.off()
cluster_summary <- data.frame( cluster = sort( unique( as.character( obj$step3A_cluster_r04 ) ) ), stringsAsFactors = FALSE )
cluster_summary$n_cells <- vapply( cluster_summary$cluster, function(cl) {
    sum( as.character( obj$step3A_cluster_r04 ) == cl )
  },
  integer(1)
)
cluster_summary$proportion <- ( cluster_summary$n_cells / ncol(obj) )
write.csv( cluster_summary, file.path( tab_dir, "11_provisional_cluster_summary_r04.csv" ), row.names = FALSE )
output_rds <- file.path( out_dir, "GSE162610_step3A_microglia_reclustered_multires.rds" )
saveRDS( obj, output_rds )
report_file <- file.path( out_dir, "12_STEP3A_microglia_reclustering_report.txt" )
sink(report_file)
print(pc_summary)
print(cluster_number_tbl)
print( table( obj$step3A_cluster_r04 ) )
print( table( obj$step2_final_cluster ) )
print( table( obj$timepoint ) )
print( table( obj$sample_id ) )
print(marker_availability)
sink()
handoff_file <- file.path( out_dir, "13_STEP3B_handoff.txt" )
sink(handoff_file)
sink()
