options(stringsAsFactors = FALSE)
options(timeout = 1200)
set.seed(20260730)
project_dir <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
in_rds <- file.path( project_dir, "02_objects", "step1D_GSE162610_LogNorm_PCA_preHarmony", "GSE162610_LogNorm_PCA_preHarmony.rds" )
out_dir <- file.path( project_dir, "02_objects", "step1E_GSE162610_Harmony_cluster" )
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "00_run_log.txt")
log_msg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  cat( paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", msg, "\n"), file = log_file, append = TRUE )
}
writeLines( paste0("GSE162610 Step 1E started: ", Sys.time()), log_file )
if (!file.exists(in_rds)) {
  stop("找不到 Step 1D 对象：\n", in_rds)
}
cran_pkgs <- c( "Seurat", "data.table", "ggplot2", "patchwork", "harmony", "clustree" )
missing_pkgs <- cran_pkgs[ !vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE) ]
if (length(missing_pkgs) > 0) {
  install.packages( missing_pkgs, repos = "https://cloud.r-project.org" )
}
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(harmony)
  library(clustree)
})
log_msg("============================================================")
log_msg("Step 1E | Harmony + formal UMAP + multi-resolution clustering")
log_msg("============================================================")
log_msg("Reading Step 1D object ...")
pbmc <- readRDS(in_rds)
DefaultAssay(pbmc) <- "RNA"
log_msg( "Input object: ", nrow(pbmc), " genes x ", ncol(pbmc), " cells" )
required_meta <- c( "sample_id", "timepoint" )
if (!all(required_meta %in% colnames(pbmc@meta.data))) {
  stop( "对象缺少 metadata：", paste( setdiff(required_meta, colnames(pbmc@meta.data)), collapse = ", " ) )
}
if (!"pca" %in% Reductions(pbmc)) {
  stop("对象中没有 PCA reduction。")
}
if (!"umap.preharmony" %in% Reductions(pbmc)) {
  warning( "对象中没有 umap.preharmony；", "将无法生成完整的 Harmony 前后 UMAP 对照。" )
}
pcs_use <- 1:20
harmony_theta <- 1
resolution_seq <- seq( 0.2, 1.2, by = 0.1 )
provisional_resolution <- 0.4
parameter_table <- data.frame( parameter = c( "Harmony group", "PCs", "Harmony theta", "Harmony max iterations", "Resolution range", "Provisional resolution" ), value = c( "sample_id", "1:20", as.character(harmony_theta), "20", "0.2 to 1.2 by 0.1", as.character(provisional_resolution) ) )
fwrite( parameter_table, file.path( tab_dir, "01_step1E_parameters.csv" ) )
log_msg("Harmony group.by.vars = sample_id")
log_msg("Harmony PCs = 1:20")
log_msg("Harmony theta = ", harmony_theta)
log_msg("timepoint will NOT be corrected.")
log_msg("Running Harmony ...")
t0 <- Sys.time()
pbmc <- RunHarmony( object = pbmc, group.by.vars = "sample_id", reduction.use = "pca", dims.use = pcs_use, theta = harmony_theta, max.iter.harmony = 20, reduction.save = "harmony", plot_convergence = FALSE, verbose = TRUE )
t_harmony <- Sys.time()
harmony_minutes <- as.numeric( difftime( t_harmony, t0, units = "mins" ) )
log_msg( "Harmony finished in ", round(harmony_minutes, 2), " min" )
if (!"harmony" %in% Reductions(pbmc)) {
  stop("Harmony 未生成 harmony reduction。")
}
harmony_dims <- ncol( Embeddings( pbmc, "harmony" ) )
log_msg( "Harmony dimensions available: ", harmony_dims )
if (harmony_dims < 20) {
  stop("Harmony 维度少于 20，无法继续使用 PCs 1:20。")
}
log_msg("Running formal Harmony UMAP ...")
pbmc <- RunUMAP( object = pbmc, reduction = "harmony", dims = pcs_use, reduction.name = "umap.harmony", reduction.key = "UMAPHARM_", seed.use = 20260730, verbose = FALSE )
if (!"umap.harmony" %in% Reductions(pbmc)) {
  stop("正式 Harmony UMAP 未Completed生成。")
}
log_msg("Harmony UMAP finished.")
log_msg("Building Harmony neighbor/SNN graph ...")
pbmc <- FindNeighbors( object = pbmc, reduction = "harmony", dims = pcs_use, graph.name = c( "harmony_nn", "harmony_snn" ), verbose = FALSE )
if (!"harmony_snn" %in% names(pbmc@graphs)) {
  stop("harmony_snn graph 未Completed生成。")
}
log_msg("Harmony neighbor graph finished.")
log_msg( "Running clustering resolutions: ", paste(resolution_seq, collapse = ", ") )
for (res in resolution_seq) {
  pbmc <- FindClusters( object = pbmc, graph.name = "harmony_snn", resolution = res, algorithm = 1, random.seed = 20260730, verbose = FALSE )
}
log_msg("Multi-resolution clustering finished.")
res_cols <- grep( "^harmony_snn_res\\.", colnames(pbmc@meta.data), value = TRUE )
if (length(res_cols) == 0) {
  stop("没有找到 harmony_snn_res.* 聚类列。")
}
cluster_number_table <- do.call( rbind, lapply( res_cols, function(cc) {
      res_value <- sub( "^harmony_snn_res\\.", "", cc )
      data.frame( resolution_column = cc, resolution = res_value, n_clusters = length( unique( pbmc@meta.data[[cc]] ) ), stringsAsFactors = FALSE )
    }
  )
)
fwrite( cluster_number_table, file.path( tab_dir, "02_cluster_number_by_resolution.csv" ) )
provisional_col <- paste0( "harmony_snn_res.", provisional_resolution )
if (!provisional_col %in% colnames(pbmc@meta.data)) {
  candidates <- res_cols[ abs( suppressWarnings( as.numeric( sub( "^harmony_snn_res\\.", "", res_cols ) ) ) - provisional_resolution ) < 1e-6 ]
  if (length(candidates) == 1) {
    provisional_col <- candidates
  } else {
    stop( "找不到 provisional resolution=", provisional_resolution, " 对应的聚类列。" )
  }
}
pbmc$cluster_provisional <- as.character( pbmc@meta.data[[provisional_col]] )
Idents(pbmc) <- "cluster_provisional"
log_msg( "Provisional display clustering = ", provisional_col, " | clusters = ", length(unique(pbmc$cluster_provisional)) )
log_msg("Generating clustree ...")
p_tree <- clustree( pbmc@meta.data, prefix = "harmony_snn_res." ) + ggtitle( "GSE162610 Harmony clustering stability" )
ggsave( filename = file.path( fig_dir, "01_clustree_resolution_0.2_to_1.2.pdf" ), plot = p_tree, width = 13, height = 10 )
p_cluster <- DimPlot( pbmc, reduction = "umap.harmony", group.by = "cluster_provisional", label = TRUE, repel = TRUE, pt.size = 0.08 ) + ggtitle( paste0( "Harmony UMAP - provisional resolution ", provisional_resolution ) )
ggsave( file.path( fig_dir, "02_Harmony_UMAP_provisional_clusters_res0.4.pdf" ), p_cluster, width = 10, height = 8 )
p_post_sample <- DimPlot( pbmc, reduction = "umap.harmony", group.by = "sample_id", pt.size = 0.08 ) + ggtitle( "Post-Harmony UMAP - sample_id" )
ggsave( file.path( fig_dir, "03_postHarmony_UMAP_by_sample.pdf" ), p_post_sample, width = 10, height = 8 )
p_post_time <- DimPlot( pbmc, reduction = "umap.harmony", group.by = "timepoint", pt.size = 0.08 ) + ggtitle( "Post-Harmony UMAP - timepoint" )
ggsave( file.path( fig_dir, "04_postHarmony_UMAP_by_timepoint.pdf" ), p_post_time, width = 9, height = 8 )
if ("umap.preharmony" %in% Reductions(pbmc)) {
  p_pre_sample <- DimPlot( pbmc, reduction = "umap.preharmony", group.by = "sample_id", pt.size = 0.06 ) + ggtitle( "Pre-Harmony - sample_id" )
  p_compare_sample <- ( p_pre_sample | p_post_sample ) + plot_annotation( title = "Sample mixing before vs after Harmony" )
  ggsave( file.path( fig_dir, "05_pre_vs_post_Harmony_by_sample.pdf" ), p_compare_sample, width = 18, height = 8 )
}
if ("umap.preharmony" %in% Reductions(pbmc)) {
  p_pre_time <- DimPlot( pbmc, reduction = "umap.preharmony", group.by = "timepoint", pt.size = 0.06 ) + ggtitle( "Pre-Harmony - timepoint" )
  p_compare_time <- ( p_pre_time | p_post_time ) + plot_annotation( title = "Timepoint structure before vs after Harmony" )
  ggsave( file.path( fig_dir, "06_pre_vs_post_Harmony_by_timepoint.pdf" ), p_compare_time, width = 17, height = 8 )
}
cluster_sample_count <- as.data.frame( table( cluster = pbmc$cluster_provisional, sample_id = pbmc$sample_id ), stringsAsFactors = FALSE )
colnames(cluster_sample_count) <- c( "cluster", "sample_id", "n_cells" )
fwrite( cluster_sample_count, file.path( tab_dir, "03_provisional_cluster_by_sample_counts.csv" ) )
cluster_sample_prop <- cluster_sample_count
cluster_totals <- aggregate( n_cells ~ cluster, data = cluster_sample_prop, FUN = sum )
colnames(cluster_totals)[2] <- "cluster_total"
cluster_sample_prop <- merge( cluster_sample_prop, cluster_totals, by = "cluster", all.x = TRUE )
cluster_sample_prop$proportion_within_cluster <- ( cluster_sample_prop$n_cells / cluster_sample_prop$cluster_total )
fwrite( cluster_sample_prop, file.path( tab_dir, "04_provisional_cluster_by_sample_proportions.csv" ) )
cluster_time_count <- as.data.frame( table( cluster = pbmc$cluster_provisional, timepoint = pbmc$timepoint ), stringsAsFactors = FALSE )
colnames(cluster_time_count) <- c( "cluster", "timepoint", "n_cells" )
fwrite( cluster_time_count, file.path( tab_dir, "05_provisional_cluster_by_timepoint_counts.csv" ) )
cluster_time_prop <- cluster_time_count
cluster_time_totals <- aggregate( n_cells ~ cluster, data = cluster_time_prop, FUN = sum )
colnames(cluster_time_totals)[2] <- "cluster_total"
cluster_time_prop <- merge( cluster_time_prop, cluster_time_totals, by = "cluster", all.x = TRUE )
cluster_time_prop$proportion_within_cluster <- ( cluster_time_prop$n_cells / cluster_time_prop$cluster_total )
fwrite( cluster_time_prop, file.path( tab_dir, "06_provisional_cluster_by_timepoint_proportions.csv" ) )
sample_cluster_prop <- cluster_sample_count
sample_totals <- aggregate( n_cells ~ sample_id, data = sample_cluster_prop, FUN = sum )
colnames(sample_totals)[2] <- "sample_total"
sample_cluster_prop <- merge( sample_cluster_prop, sample_totals, by = "sample_id", all.x = TRUE )
sample_cluster_prop$proportion_within_sample <- ( sample_cluster_prop$n_cells / sample_cluster_prop$sample_total )
p_sample_heat <- ggplot( sample_cluster_prop, aes( x = sample_id, y = cluster, fill = proportion_within_sample ) ) + geom_tile() + labs( x = NULL, y = "Provisional cluster", fill = "Fraction", title = "Cluster composition within each sample" ) + theme_classic() + theme( axis.text.x = element_text( angle = 45, hjust = 1 ) )
ggsave( file.path( fig_dir, "07_cluster_composition_by_sample_heatmap.pdf" ), p_sample_heat, width = 10, height = 8 )
time_cluster_prop <- cluster_time_count
time_totals <- aggregate( n_cells ~ timepoint, data = time_cluster_prop, FUN = sum )
colnames(time_totals)[2] <- "time_total"
time_cluster_prop <- merge( time_cluster_prop, time_totals, by = "timepoint", all.x = TRUE )
time_cluster_prop$proportion_within_timepoint <- ( time_cluster_prop$n_cells / time_cluster_prop$time_total )
p_time_heat <- ggplot( time_cluster_prop, aes( x = timepoint, y = cluster, fill = proportion_within_timepoint ) ) + geom_tile() + labs( x = NULL, y = "Provisional cluster", fill = "Fraction", title = "Cluster composition within each timepoint" ) + theme_classic()
ggsave( file.path( fig_dir, "08_cluster_composition_by_timepoint_heatmap.pdf" ), p_time_heat, width = 7, height = 8 )
cluster_size_all <- do.call( rbind, lapply( res_cols, function(cc) {
      tt <- as.data.frame( table( pbmc@meta.data[[cc]] ), stringsAsFactors = FALSE )
      colnames(tt) <- c( "cluster", "n_cells" )
      tt$resolution_column <- cc
      tt$resolution <- sub( "^harmony_snn_res\\.", "", cc )
      tt
    }
  )
)
fwrite( cluster_size_all, file.path( tab_dir, "07_cluster_sizes_all_resolutions.csv" ) )
saveRDS( pbmc, file.path( out_dir, "GSE162610_Harmony_multires_clustered.rds" ), compress = FALSE )
all_samples_preserved <- ( length(unique(as.character(pbmc$sample_id))) == 10 )
all_timepoints_preserved <- ( length(unique(as.character(pbmc$timepoint))) == 4 )
harmony_ok <- ( "harmony" %in% Reductions(pbmc) && ncol(Embeddings(pbmc, "harmony")) >= 20 )
umap_ok <- ( "umap.harmony" %in% Reductions(pbmc) )
resolution_ok <- ( length(res_cols) >= length(resolution_seq) )
step_pass <- all( all_samples_preserved, all_timepoints_preserved, harmony_ok, umap_ok, resolution_ok )
runtime_total <- as.numeric( difftime( Sys.time(), t0, units = "mins" ) )
report <- c( "============================================================", "GSE162610 STEP 1E | HARMONY + MULTI-RESOLUTION CLUSTERING REPORT", "============================================================", paste0("Date: ", Sys.time()), paste0("Total runtime (min): ", round(runtime_total, 2)), paste0("Harmony runtime (min): ", round(harmony_minutes, 2)), "", paste0("Input cells: ", ncol(pbmc)), paste0("Genes: ", nrow(pbmc)), "", "Harmony settings:", "  group.by.vars = sample_id", "  timepoint correction = NO", "  reduction.use = pca", "  dims.use = 1:20", paste0("  theta = ", harmony_theta), "  max.iter.harmony = 20", "", "Formal downstream embedding:", "  UMAP reduction = harmony", "  FindNeighbors reduction = harmony", "  dims = 1:20", "", paste0( "Resolutions calculated: ", paste( cluster_number_table$resolution, collapse = ", " ) ), "", "Cluster number by resolution:", paste0( "  res ", cluster_number_table$resolution, " -> ", cluster_number_table$n_clusters, " clusters" ), "", paste0( "Provisional display resolution: ", provisional_resolution ), paste0( "Provisional number of clusters: ", length(unique(pbmc$cluster_provisional)) ), "", paste0( "All 10 samples preserved: ", all_samples_preserved ), paste0( "All 4 timepoints preserved: ", all_timepoints_preserved ), paste0( "Harmony successful: ", harmony_ok ), paste0( "Harmony UMAP successful: ", umap_ok ), paste0( "All requested resolutions generated: ", resolution_ok ), "", "IMPORTANT:", "Resolution 0.4 is provisional only.", "Do NOT start cell-type annotation until the clustree and pre/post-Harmony plots are reviewed.", "The key acceptance criterion is better replicate mixing without erasing timepoint-associated biology.", "", paste0( "STEP 1E COMPLETE = ", step_pass ), "", "NEXT AFTER REVIEW:", "Step 1F | Choose final broad-clustering resolution and identify cluster markers for major-cell-type annotation.", "============================================================" )
writeLines( report, file.path( out_dir, "08_STEP1E_Harmony_cluster_report.txt" ) )
log_msg( "Step 1E complete = ", step_pass )
cat("\n")
cat(paste(report, collapse = "\n"))
cat("\n\n")
if (step_pass) {
  cat( "Step 1E 运行完成。分析完成。\n" )
} else {
  cat("Step 1E status: REVIEW REQUIRED\n")
}
