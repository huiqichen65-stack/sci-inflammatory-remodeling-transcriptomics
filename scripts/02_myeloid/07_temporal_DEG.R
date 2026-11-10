options(stringsAsFactors = FALSE)
set.seed(1234)
required_pkgs <- c( "Seurat", "SeuratObject", "dplyr", "tidyr", "ggplot2", "patchwork", "ggrepel", "Matrix", "pheatmap", "clusterProfiler", "org.Mm.eg.db", "enrichplot" )
missing_pkgs <- required_pkgs[ !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE) ]
if (length(missing_pkgs) > 0) {
  stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))
}
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(Matrix)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
})
if (requireNamespace("R.utils", quietly = TRUE)) {
  R.utils::setOption( "clusterProfiler.download.method", "auto" )
}
project_dir <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
input_rds <- file.path( project_dir, "04_results", "step3C_GSE162610_activated_DAM_transition_mapping", "GSE162610_step3C_transition_mapped.rds" )
out_dir <- file.path( project_dir, "04_results", "step3D_GSE162610_time_DEG_pathway_confirmation" )
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
go_dir <- file.path(out_dir, "GO_KEGG")
gsea_dir <- file.path(out_dir, "GSEA")
pb_dir <- file.path(out_dir, "pseudobulk")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pb_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_rds)) {
  stop( "Required input not found: ", input_rds )
}
obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) {
  stop("Input object is not a Seurat object.")
}
DefaultAssay(obj) <- "RNA"
required_meta <- c( "sample_id", "timepoint", "step3A_cluster_r04" )
missing_meta <- required_meta[ !required_meta %in% colnames(obj[[]]) ]
if (length(missing_meta) > 0) {
  stop( "Missing metadata: ", paste(missing_meta, collapse = ", ") )
}
obj$step3A_cluster_r04 <- as.character( obj$step3A_cluster_r04 )
time_levels <- c( "Uninjured", "1 dpi", "3 dpi", "7 dpi" )
obj$timepoint_step3D <- factor( as.character(obj$timepoint), levels = time_levels )
if (any(is.na(obj$timepoint_step3D))) {
  stop("Consistency check failed.")
}
print(table(obj$timepoint_step3D))
print(table(obj$sample_id))
comparisons <- list( C1_1dpi_vs_Uninjured = c( case = "1 dpi", ref = "Uninjured" ), C2_3dpi_vs_1dpi = c( case = "3 dpi", ref = "1 dpi" ), C3_7dpi_vs_3dpi = c( case = "7 dpi", ref = "3 dpi" ), S1_3dpi_vs_Uninjured = c( case = "3 dpi", ref = "Uninjured" ), S2_7dpi_vs_Uninjured = c( case = "7 dpi", ref = "Uninjured" ) )
primary_comparisons <- c( "C1_1dpi_vs_Uninjured", "C2_3dpi_vs_1dpi", "C3_7dpi_vs_3dpi" )
get_fc_col <- function(df) {
  hit <- intersect( c( "avg_log2FC", "avg_logFC" ), colnames(df) )
  if (length(hit) == 0) {
    stop("DEG table does not contain avg_log2FC/avg_logFC.")
  }
  hit[1]
}
safe_findmarkers <- function( sobj, case, ref ) {
  Idents(sobj) <- "timepoint_step3D"
  out <- FindMarkers( sobj, ident.1 = case, ident.2 = ref, assay = "RNA", test.use = "wilcox", logfc.threshold = 0, min.pct = 0.01, only.pos = FALSE, verbose = FALSE )
  out$gene <- rownames(out)
  fc_col <- get_fc_col(out)
  out$log2FC <- as.numeric( out[[fc_col]] )
  out <- out[ order( out$p_val_adj, -abs(out$log2FC) ), , drop = FALSE ]
  rownames(out) <- NULL
  out
}
add_deg_class <- function( df, padj_cutoff = 0.05, lfc_cutoff = 0.25 ) {
  df$DEG_class <- "NS"
  df$DEG_class[ !is.na(df$p_val_adj) & df$p_val_adj < padj_cutoff & df$log2FC >= lfc_cutoff ] <- "UP"
  df$DEG_class[ !is.na(df$p_val_adj) & df$p_val_adj < padj_cutoff & df$log2FC <= -lfc_cutoff ] <- "DOWN"
  df
}
deg_list <- list()
for (nm in names(comparisons)) {
  case <- unname( comparisons[[nm]]["case"] )
  ref <- unname( comparisons[[nm]]["ref"] )
  df <- safe_findmarkers( obj, case = case, ref = ref )
  df <- add_deg_class(df)
  df$comparison <- nm
  df$case <- case
  df$reference <- ref
  deg_list[[nm]] <- df
  write.csv( df, file.path( tab_dir, paste0( "01_DEG_", nm, ".csv" ) ), row.names = FALSE )
}
deg_all <- dplyr::bind_rows( deg_list )
write.csv( deg_all, file.path( tab_dir, "02_all_time_DEGs_combined.csv" ), row.names = FALSE )
deg_summary <- deg_all %>% dplyr::group_by( comparison, case, reference, DEG_class ) %>% dplyr::summarise( n_genes = n(), .groups = "drop" )
write.csv( deg_summary, file.path( tab_dir, "03_DEG_count_summary.csv" ), row.names = FALSE )
core_exclude <- c( "6", "8", "9" )
core_cells <- colnames(obj)[ !obj$step3A_cluster_r04 %in% core_exclude ]
obj_core <- subset( obj, cells = core_cells )
obj_core$timepoint_step3D <- factor( as.character(obj_core$timepoint), levels = time_levels )
core_deg_list <- list()
for (nm in primary_comparisons) {
  case <- unname( comparisons[[nm]]["case"] )
  ref <- unname( comparisons[[nm]]["ref"] )
  df <- safe_findmarkers( obj_core, case = case, ref = ref )
  df <- add_deg_class(df)
  df$comparison <- nm
  df$case <- case
  df$reference <- ref
  core_deg_list[[nm]] <- df
  write.csv( df, file.path( tab_dir, paste0( "04_CORE_DEG_", nm, ".csv" ) ), row.names = FALSE )
}
c1 <- deg_list[["C1_1dpi_vs_Uninjured"]]
c2 <- deg_list[["C2_3dpi_vs_1dpi"]]
c3 <- deg_list[["C3_7dpi_vs_3dpi"]]
early_sig <- c1 %>% dplyr::filter( p_val_adj < 0.05, log2FC >= 0.25 ) %>% dplyr::arrange( desc(log2FC) )
transition_base <- c2 %>% dplyr::filter( p_val_adj < 0.05, log2FC >= 0.25 ) %>% dplyr::select( gene, p_val_adj_3v1 = p_val_adj, log2FC_3v1 = log2FC )
c3_for_merge <- c3 %>% dplyr::select( gene, p_val_adj_7v3 = p_val_adj, log2FC_7v3 = log2FC )
transition_sig <- transition_base %>% dplyr::left_join( c3_for_merge, by = "gene" ) %>% dplyr::mutate( p_val_adj_7v3 = ifelse( is.na(p_val_adj_7v3), 1, p_val_adj_7v3 ), log2FC_7v3 = ifelse( is.na(log2FC_7v3), 0, log2FC_7v3 ) ) %>% dplyr::filter( log2FC_7v3 <= 0.15 | p_val_adj_7v3 >= 0.05 ) %>% dplyr::arrange( desc(log2FC_3v1) )
late_sig <- c3 %>% dplyr::filter( p_val_adj < 0.05, log2FC >= 0.25 ) %>% dplyr::arrange( desc(log2FC) )
transition_strict <- transition_base %>% dplyr::left_join( c3_for_merge, by = "gene" ) %>% dplyr::mutate( p_val_adj_7v3 = ifelse( is.na(p_val_adj_7v3), 1, p_val_adj_7v3 ), log2FC_7v3 = ifelse( is.na(log2FC_7v3), 0, log2FC_7v3 ) ) %>% dplyr::filter( p_val_adj_7v3 >= 0.05 | log2FC_7v3 <= 0 ) %>% dplyr::arrange( desc(log2FC_3v1) )
write.csv( early_sig, file.path( tab_dir, "05_Early_activation_1dpi_signature.csv" ), row.names = FALSE )
write.csv( transition_sig, file.path( tab_dir, "06_Transition_3dpi_signature.csv" ), row.names = FALSE )
write.csv( transition_strict, file.path( tab_dir, "07_Transition_3dpi_signature_STRICT.csv" ), row.names = FALSE )
write.csv( late_sig, file.path( tab_dir, "08_Late_remodeling_7dpi_signature.csv" ), row.names = FALSE )
stage_signature_summary <- data.frame( stage = c( "Early_activation_1dpi", "Transition_3dpi", "Transition_3dpi_STRICT", "Late_remodeling_7dpi" ), definition = c( "UP in 1dpi vs Uninjured", "UP in 3dpi vs 1dpi; not further strongly increased in 7dpi vs 3dpi", "UP in 3dpi vs 1dpi; 7dpi vs 3dpi non-significant or non-increasing", "UP in 7dpi vs 3dpi" ), n_genes = c( nrow(early_sig), nrow(transition_sig), nrow(transition_strict), nrow(late_sig) ), stringsAsFactors = FALSE )
write.csv( stage_signature_summary, file.path( tab_dir, "09_stage_signature_summary.csv" ), row.names = FALSE )
top_n_each <- 20
early_top <- head( early_sig$gene, top_n_each )
transition_top <- head( transition_sig$gene, top_n_each )
late_top <- head( late_sig$gene, top_n_each )
heat_genes <- unique( c( early_top, transition_top, late_top ) )
heat_genes <- intersect( heat_genes, rownames(obj) )
safe_time_map <- c( "Uninjured" = "Uninjured", "1 dpi" = "dpi1", "3 dpi" = "dpi3", "7 dpi" = "dpi7" )
obj$timepoint_step3D_safe <- unname( safe_time_map[ as.character( obj$timepoint_step3D ) ] )
obj$timepoint_step3D_safe <- factor( obj$timepoint_step3D_safe, levels = c( "Uninjured", "dpi1", "dpi3", "dpi7" ) )
if (any(is.na(obj$timepoint_step3D_safe))) {
  stop("Safe timepoint mapping failed.")
}
avg_expr <- AverageExpression( obj, assays = "RNA", group.by = "timepoint_step3D_safe", layer = "data", verbose = FALSE )[["RNA"]]
avg_expr <- as.matrix( avg_expr )
expected_safe_cols <- c( "Uninjured", "dpi1", "dpi3", "dpi7" )
if (!all(expected_safe_cols %in% colnames(avg_expr))) {
  stop( paste0( "AverageExpression safe timepoint columns missing. Observed: ", paste(colnames(avg_expr), collapse = ", ") ) )
}
avg_expr <- avg_expr[ , expected_safe_cols, drop = FALSE ]
colnames(avg_expr) <- time_levels
heat_mat <- avg_expr[ heat_genes, , drop = FALSE ]
row_z <- function(x) {
  sx <- sd( x, na.rm = TRUE )
  if (!is.finite(sx) || sx == 0) {
    return( rep( 0, length(x) ) )
  }
  ( x - mean( x, na.rm = TRUE ) ) / sx
}
heat_z <- t( apply( heat_mat, 1, row_z ) )
stage_for_gene <- rep( "Other", length(heat_genes) )
names(stage_for_gene) <- heat_genes
stage_for_gene[ heat_genes %in% early_top ] <- "Early_1dpi"
stage_for_gene[ heat_genes %in% transition_top ] <- "Transition_3dpi"
stage_for_gene[ heat_genes %in% late_top ] <- "Late_7dpi"
row_anno <- data.frame( Stage_signature = stage_for_gene[ rownames(heat_z) ], row.names = rownames(heat_z) )
pdf( file.path( fig_dir, "02_stage_signature_time_heatmap.pdf" ), width = 8, height = max( 9, 0.16 * nrow(heat_z) + 3 ) )
pheatmap( heat_z, cluster_rows = FALSE, cluster_cols = FALSE, annotation_row = row_anno, main = "Data-driven temporal signatures | row Z-score" )
dev.off()
volcano_df <- dplyr::bind_rows( lapply( primary_comparisons, function(nm) {
      x <- deg_list[[nm]]
      x$comparison_display <- switch( nm, C1_1dpi_vs_Uninjured = "1 dpi vs Uninjured", C2_3dpi_vs_1dpi = "3 dpi vs 1 dpi", C3_7dpi_vs_3dpi = "7 dpi vs 3 dpi" )
      x
    }
  )
)
volcano_df$minus_log10_padj <- -log10( pmax( volcano_df$p_val_adj, 1e-300 ) )
label_df <- volcano_df %>% dplyr::filter( DEG_class != "NS" ) %>% dplyr::group_by( comparison_display ) %>% dplyr::arrange( p_val_adj, desc( abs(log2FC) ), .by_group = TRUE ) %>% slice_head( n = 10 ) %>% dplyr::ungroup()
p_volcano <- ggplot( volcano_df, aes( x = log2FC, y = minus_log10_padj ) ) + geom_point( aes( shape = DEG_class ), alpha = 0.45, size = 0.75 ) + geom_vline( xintercept = c( -0.25, 0.25 ), linetype = "dashed" ) + geom_hline( yintercept = -log10( 0.05 ), linetype = "dashed" ) + ggrepel::geom_text_repel( data = label_df, aes( label = gene ), size = 2.6, max.overlaps = 30, min.segment.length = 0 ) + facet_wrap( ~ comparison_display, nrow = 1, scales = "free" ) + theme_bw() + labs( title = "GSE162610 microglia | real-time differential expression", x = "log2 fold change", y = "-log10 adjusted P", shape = "DEG" )
ggsave( file.path( fig_dir, "01_time_DEG_multiVolcano.pdf" ), p_volcano, width = 16, height = 6 )
symbol_to_entrez <- function( symbols ) {
  symbols <- unique( symbols[ !is.na(symbols) & symbols != "" ] )
  if (length(symbols) == 0) {
    return( data.frame() )
  }
  suppressMessages( bitr( symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db ) )
}
stage_gene_lists <- list( Early_activation_1dpi = early_sig$gene, Transition_3dpi = transition_sig$gene, Late_remodeling_7dpi = late_sig$gene )
ora_go_all <- list()
ora_kegg_all <- list()
for (stage_nm in names(stage_gene_lists)) {
  genes_symbol <- stage_gene_lists[[stage_nm]]
  id_map <- symbol_to_entrez( genes_symbol )
  if (nrow(id_map) < 10) {
    warning( "Stage ", stage_nm, ": mapped ENTREZ genes <10; ORA skipped." )
    next
  }
  entrez <- unique( id_map$ENTREZID )
  ego <- tryCatch( enrichGO( gene = entrez, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.20, readable = TRUE ), error = function(e) {
      NULL
    }
  )
  if (!is.null(ego)) {
    ego_df <- as.data.frame( ego )
    if (nrow(ego_df) > 0) {
      ego_df$stage <- stage_nm
      ora_go_all[[stage_nm]] <- ego_df
      write.csv( ego_df, file.path( go_dir, paste0( stage_nm, "_GO_BP.csv" ) ), row.names = FALSE )
      pdf( file.path( go_dir, paste0( stage_nm, "_GO_BP_dotplot.pdf" ) ), width = 10, height = 7 )
      print( dotplot( ego, showCategory = 15, label_format = 90 ) + ggtitle( paste0( stage_nm, " | GO Biological Process" ) ) )
      dev.off()
    }
  }
  ekegg <- tryCatch( enrichKEGG( gene = entrez, organism = "mmu", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.20 ), error = function(e) {
      NULL
    }
  )
  if (!is.null(ekegg)) {
    ekegg_df <- as.data.frame( ekegg )
    if (nrow(ekegg_df) > 0) {
      ekegg_df$stage <- stage_nm
      ora_kegg_all[[stage_nm]] <- ekegg_df
      write.csv( ekegg_df, file.path( go_dir, paste0( stage_nm, "_KEGG.csv" ) ), row.names = FALSE )
      pdf( file.path( go_dir, paste0( stage_nm, "_KEGG_dotplot.pdf" ) ), width = 10, height = 7 )
      print( dotplot( ekegg, showCategory = 15, label_format = 90 ) + ggtitle( paste0( stage_nm, " | KEGG" ) ) )
      dev.off()
    }
  }
}
go_combined <- if ( length( ora_go_all ) > 0 ) {
  dplyr::bind_rows( ora_go_all )
} else {
  data.frame()
}
kegg_combined <- if ( length( ora_kegg_all ) > 0 ) {
  dplyr::bind_rows( ora_kegg_all )
} else {
  data.frame()
}
write.csv( go_combined, file.path( tab_dir, "10_GO_BP_stage_combined.csv" ), row.names = FALSE )
write.csv( kegg_combined, file.path( tab_dir, "11_KEGG_stage_combined.csv" ), row.names = FALSE )
make_enrichment_overview <- function( df, title, outfile, top_n = 8 ) {
  if ( is.null(df) || nrow(df) == 0 ) {
    return( invisible( NULL ) )
  }
  df2 <- df %>% dplyr::filter( !is.na(p.adjust) ) %>% dplyr::group_by( stage ) %>% dplyr::arrange( p.adjust, .by_group = TRUE ) %>% slice_head( n = top_n ) %>% dplyr::ungroup()
  if (nrow(df2) == 0) {
    return( invisible( NULL ) )
  }
  df2$minus_log10_padj <- -log10( pmax( df2$p.adjust, 1e-300 ) )
  p <- ggplot( df2, aes( x = minus_log10_padj, y = reorder( Description, minus_log10_padj ) ) ) + geom_point( aes( size = Count, shape = stage ) ) + facet_wrap( ~ stage, scales = "free_y", ncol = 1 ) + theme_bw() + labs( title = title, x = "-log10 adjusted P", y = NULL, size = "Gene count", shape = "Stage" )
  ggsave( outfile, p, width = 11, height = max( 8, 0.35 * nrow(df2) ) )
}
make_enrichment_overview( go_combined, "GO Biological Process | data-driven SCI microglial stages", file.path( fig_dir, "03_GO_BP_stage_overview.pdf" ) )
make_enrichment_overview( kegg_combined, "KEGG | data-driven SCI microglial stages", file.path( fig_dir, "04_KEGG_stage_overview.pdf" ) )
make_ranked_entrez <- function( deg_df ) {
  df <- deg_df %>% dplyr::filter( !is.na(log2FC), !is.na(gene) ) %>% dplyr::select( gene, log2FC )
  id_map <- symbol_to_entrez( df$gene )
  if (nrow(id_map) == 0) {
    return( numeric(0) )
  }
  x <- df %>% dplyr::inner_join( id_map, by = c( "gene" = "SYMBOL" ) ) %>% dplyr::arrange( desc( abs(log2FC) ) ) %>% dplyr::distinct( ENTREZID, .keep_all = TRUE )
  gene_list <- x$log2FC
  names(gene_list) <- x$ENTREZID
  sort( gene_list, decreasing = TRUE )
}
gsea_go_all <- list()
gsea_kegg_all <- list()
for (nm in primary_comparisons) {
  rank_list <- make_ranked_entrez( deg_list[[nm]] )
  if (length(rank_list) < 100) {
    warning( "Rank list too short for ", nm, "; GSEA skipped." )
    next
  }
  ggo <- tryCatch( gseGO( geneList = rank_list, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "BH", verbose = FALSE ), error = function(e) {
      NULL
    }
  )
  if (!is.null(ggo)) {
    ggo_df <- as.data.frame( ggo )
    if (nrow(ggo_df) > 0) {
      ggo_df$comparison <- nm
      gsea_go_all[[nm]] <- ggo_df
      write.csv( ggo_df, file.path( gsea_dir, paste0( nm, "_GSEA_GO_BP.csv" ) ), row.names = FALSE )
    }
  }
  gkegg <- tryCatch( gseKEGG( geneList = rank_list, organism = "mmu", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "BH", verbose = FALSE ), error = function(e) {
      NULL
    }
  )
  if (!is.null(gkegg)) {
    gkegg_df <- as.data.frame( gkegg )
    if (nrow(gkegg_df) > 0) {
      gkegg_df$comparison <- nm
      gsea_kegg_all[[nm]] <- gkegg_df
      write.csv( gkegg_df, file.path( gsea_dir, paste0( nm, "_GSEA_KEGG.csv" ) ), row.names = FALSE )
    }
  }
}
gsea_go_combined <- if ( length( gsea_go_all ) > 0 ) {
  dplyr::bind_rows( gsea_go_all )
} else {
  data.frame()
}
gsea_kegg_combined <- if ( length( gsea_kegg_all ) > 0 ) {
  dplyr::bind_rows( gsea_kegg_all )
} else {
  data.frame()
}
write.csv( gsea_go_combined, file.path( tab_dir, "12_GSEA_GO_BP_combined.csv" ), row.names = FALSE )
write.csv( gsea_kegg_combined, file.path( tab_dir, "13_GSEA_KEGG_combined.csv" ), row.names = FALSE )
make_gsea_overview <- function( df, title, outfile, top_n = 8 ) {
  if ( is.null(df) || nrow(df) == 0 ) {
    return( invisible( NULL ) )
  }
  df2 <- df %>% dplyr::filter( !is.na(p.adjust), !is.na(NES) ) %>% dplyr::group_by( comparison ) %>% dplyr::arrange( p.adjust, desc( abs(NES) ), .by_group = TRUE ) %>% slice_head( n = top_n ) %>% dplyr::ungroup()
  if (nrow(df2) == 0) {
    return( invisible( NULL ) )
  }
  df2$comparison_display <- dplyr::recode( df2$comparison, C1_1dpi_vs_Uninjured = "1 dpi vs Uninjured", C2_3dpi_vs_1dpi = "3 dpi vs 1 dpi", C3_7dpi_vs_3dpi = "7 dpi vs 3 dpi" )
  p <- ggplot( df2, aes( x = NES, y = reorder( Description, NES ) ) ) + geom_point( aes( size = -log10( pmax( p.adjust, 1e-300 ) ), shape = comparison_display ) ) + geom_vline( xintercept = 0, linetype = "dashed" ) + facet_wrap( ~ comparison_display, scales = "free_y", ncol = 1 ) + theme_bw() + labs( title = title, x = "NES", y = NULL, size = "-log10 adj.P", shape = "Comparison" )
  ggsave( outfile, p, width = 11, height = max( 8, 0.35 * nrow(df2) ) )
}
make_gsea_overview( gsea_kegg_combined, "GSEA KEGG | real-time microglial remodeling", file.path( fig_dir, "05_GSEA_KEGG_stage_overview.pdf" ) )
make_gsea_overview( gsea_go_combined, "GSEA GO-BP | real-time microglial remodeling", file.path( fig_dir, "06_GSEA_GO_BP_stage_overview.pdf" ) )
pseudobulk_run <- FALSE
pb_deg_list <- list()
if ( requireNamespace( "edgeR", quietly = TRUE ) ) {
  get_counts <- function( sobj ) {
    x <- tryCatch( GetAssayData( sobj, assay = "RNA", layer = "counts" ), error = function(e) NULL )
    if (is.null(x)) {
      x <- tryCatch( GetAssayData( sobj, assay = "RNA", slot = "counts" ), error = function(e) NULL )
    }
    if (is.null(x)) {
      stop( "Cannot access RNA counts for pseudobulk." )
    }
    x
  }
  counts <- get_counts( obj )
  if (!inherits(counts, "dgCMatrix")) {
    counts <- as( counts, "dgCMatrix" )
  }
  sample_fac <- factor( as.character( obj$sample_id ) )
  sample_design <- Matrix::sparse.model.matrix( ~ 0 + sample_fac )
  pb_counts <- counts %*%
    sample_design
  colnames(pb_counts) <- levels( sample_fac )
  sample_meta <- obj[[]] %>% dplyr::transmute( sample_id = as.character( sample_id ), timepoint = as.character( timepoint ) ) %>% dplyr::distinct() %>% dplyr::filter( sample_id %in% colnames( pb_counts ) )
  sample_meta <- sample_meta[ match( colnames(pb_counts), sample_meta$sample_id ), , drop = FALSE ]
  if ( any( is.na( sample_meta$timepoint ) ) ) {
    stop( "Pseudobulk sample-timepoint mapping failed." )
  }
  write.csv( sample_meta, file.path( pb_dir, "00_pseudobulk_sample_metadata.csv" ), row.names = FALSE )
  run_pb_edgeR <- function( case, ref, comparison_name ) {
    keep_samples <- sample_meta$timepoint %in%
      c( ref, case )
    cc <- pb_counts[ , keep_samples, drop = FALSE ]
    mm <- sample_meta[ keep_samples, , drop = FALSE ]
    group <- factor( mm$timepoint, levels = c( ref, case ) )
    y <- edgeR::DGEList( counts = cc, group = group )
    keep_gene <- edgeR::filterByExpr( y, group = group )
    y <- y[ keep_gene, , keep.lib.sizes = FALSE ]
    y <- edgeR::calcNormFactors( y )
    design <- model.matrix( ~ group )
    y <- edgeR::estimateDisp( y, design, robust = TRUE )
    fit <- edgeR::glmQLFit( y, design, robust = TRUE )
    qlf <- edgeR::glmQLFTest( fit, coef = 2 )
    tt <- edgeR::topTags( qlf, n = Inf, sort.by = "PValue" )$table
    tt$gene <- rownames( tt )
    tt$comparison <- comparison_name
    tt$case <- case
    tt$reference <- ref
    rownames(tt) <- NULL
    tt
  }
  for (nm in primary_comparisons) {
    case <- unname( comparisons[[nm]]["case"] )
    ref <- unname( comparisons[[nm]]["ref"] )
    pb <- tryCatch( run_pb_edgeR( case = case, ref = ref, comparison_name = nm ), error = function(e) {
        NULL
      }
    )
    if (!is.null(pb)) {
      pb_deg_list[[nm]] <- pb
      write.csv( pb, file.path( pb_dir, paste0( "01_pseudobulk_edgeR_", nm, ".csv" ) ), row.names = FALSE )
    }
  }
  if ( length( pb_deg_list ) == length( primary_comparisons ) ) {
    pseudobulk_run <- TRUE
  }
} else {
}
pb_concordance <- data.frame()
if ( length( pb_deg_list ) > 0 ) {
  concordance_list <- list()
  for (nm in names(pb_deg_list)) {
    cell_df <- deg_list[[nm]] %>% dplyr::select( gene, cell_log2FC = log2FC, cell_padj = p_val_adj )
    pb_df <- pb_deg_list[[nm]] %>% dplyr::select( gene, pb_log2FC = logFC, pb_FDR = FDR )
    xx <- dplyr::inner_join( cell_df, pb_df, by = "gene" )
    xx$comparison <- nm
    xx$direction_concordant <- sign( xx$cell_log2FC ) == sign( xx$pb_log2FC )
    concordance_list[[nm]] <- xx
  }
  pb_concordance <- dplyr::bind_rows( concordance_list )
  write.csv( pb_concordance, file.path( tab_dir, "14_cell_vs_pseudobulk_DEG_concordance.csv" ), row.names = FALSE )
  p_pb <- ggplot( pb_concordance, aes( x = cell_log2FC, y = pb_log2FC ) ) + geom_point( alpha = 0.25, size = 0.65 ) + geom_hline( yintercept = 0, linetype = "dashed" ) + geom_vline( xintercept = 0, linetype = "dashed" ) + facet_wrap( ~ comparison, scales = "free", nrow = 1 ) + theme_bw() + labs( title = "Cell-level vs sample-level pseudobulk DEG direction", x = "Seurat cell-level log2FC", y = "edgeR pseudobulk log2FC" )
  ggsave( file.path( fig_dir, "07_cell_vs_pseudobulk_concordance.pdf" ), p_pb, width = 15, height = 5 )
}
signature_pb_support <- data.frame()
if ( length( pb_deg_list ) == length( primary_comparisons ) ) {
  pb1 <- pb_deg_list[["C1_1dpi_vs_Uninjured"]] %>% dplyr::select( gene, pb_logFC = logFC, pb_FDR = FDR )
  pb2 <- pb_deg_list[["C2_3dpi_vs_1dpi"]] %>% dplyr::select( gene, pb_logFC = logFC, pb_FDR = FDR )
  pb3 <- pb_deg_list[["C3_7dpi_vs_3dpi"]] %>% dplyr::select( gene, pb_logFC = logFC, pb_FDR = FDR )
  e <- early_sig %>% dplyr::select( gene, cell_log2FC = log2FC, cell_padj = p_val_adj ) %>% dplyr::left_join( pb1, by = "gene" ) %>% dplyr::mutate( stage = "Early_activation_1dpi", pseudobulk_same_direction = !is.na(pb_logFC) & pb_logFC > 0 )
  t <- transition_sig %>% dplyr::transmute( gene, cell_log2FC = log2FC_3v1, cell_padj = p_val_adj_3v1 ) %>% dplyr::left_join( pb2, by = "gene" ) %>% dplyr::mutate( stage = "Transition_3dpi", pseudobulk_same_direction = !is.na(pb_logFC) & pb_logFC > 0 )
  l <- late_sig %>% dplyr::select( gene, cell_log2FC = log2FC, cell_padj = p_val_adj ) %>% dplyr::left_join( pb3, by = "gene" ) %>% dplyr::mutate( stage = "Late_remodeling_7dpi", pseudobulk_same_direction = !is.na(pb_logFC) & pb_logFC > 0 )
  signature_pb_support <- dplyr::bind_rows( e, t, l )
  write.csv( signature_pb_support, file.path( tab_dir, "15_stage_signature_pseudobulk_support.csv" ), row.names = FALSE )
}
obj$step3D_core_state_eligible <- !obj$step3A_cluster_r04 %in%
  core_exclude
output_rds <- file.path( out_dir, "GSE162610_step3D_time_DEG_pathway_confirmed.rds" )
saveRDS( obj, output_rds )
report_file <- file.path( out_dir, "16_STEP3D_time_DEG_pathway_report.txt" )
sink(report_file)
print( table( obj$timepoint_step3D ) )
print( table( obj$sample_id ) )
print( deg_summary )
print( stage_signature_summary )
if ( nrow( pb_concordance ) > 0 ) {
  pb_summary <- pb_concordance %>% dplyr::group_by( comparison ) %>% dplyr::summarise( n_genes = n(), direction_concordance = mean( direction_concordant, na.rm = TRUE ), .groups = "drop" )
  print( pb_summary )
}
if ( nrow( signature_pb_support ) > 0 ) {
  sig_support_summary <- signature_pb_support %>% dplyr::group_by( stage ) %>% dplyr::summarise( n_genes = n(), same_direction_fraction = mean( pseudobulk_same_direction, na.rm = TRUE ), pb_FDR_lt_0.10_fraction = mean( pb_FDR < 0.10, na.rm = TRUE ), .groups = "drop" )
  print( sig_support_summary )
}
sink()
handoff_file <- file.path( out_dir, "17_STEP4_handoff.txt" )
sink(handoff_file)
sink()
