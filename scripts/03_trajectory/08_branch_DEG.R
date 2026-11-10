options(stringsAsFactors = FALSE)
set.seed(20260731)
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
STEP5CR_DIR <- file.path( PROJECT_ROOT, "04_results", "step5CR_GSE162610_Y29_prebranch_progressive_divergence" )
STEP5BR2_DIR <- file.path( PROJECT_ROOT, "04_results", "step5BR2_GSE162610_single_root_branch_audit" )
INPUT_STEP5CR <- file.path( STEP5CR_DIR, "GSE162610_step5CR_results.rds" )
INPUT_SEURAT <- file.path( STEP5BR2_DIR, "GSE162610_step5BR2_single_root_seurat.rds" )
OUTDIR <- file.path( PROJECT_ROOT, "04_results", "step6A_GSE162610_7dpi_branch_specific_DEG" )
FIGDIR <- file.path(OUTDIR, "figures")
TABDIR <- file.path(OUTDIR, "tables")
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
OUT_FIG1 <- file.path( FIGDIR, "Fig6A_01_7dpi_Arm1_Arm2_UMAP.pdf" )
OUT_FIG2 <- file.path( FIGDIR, "Fig6A_02_celllevel_volcano.pdf" )
OUT_FIG3 <- file.path( FIGDIR, "Fig6A_03_replicate_concordance.pdf" )
OUT_FIG4 <- file.path( FIGDIR, "Fig6A_04_pseudobulk_MDS.pdf" )
OUT_FIG5 <- file.path( FIGDIR, "Fig6A_05_highconfidence_pseudobulk_heatmap.pdf" )
OUT_FIG6 <- file.path( FIGDIR, "Fig6A_06_evidence_cascade.pdf" )
OUT_XLSX <- file.path( OUTDIR, "Table6A_7dpi_branch_specific_DEG.xlsx" )
OUT_RDS <- file.path( OUTDIR, "GSE162610_step6A_7dpi_branch_DEG_results.rds" )
OUT_REP <- file.path( OUTDIR, "Step6A_report.txt" )
OUT_SESS <- file.path( OUTDIR, "Step6A_sessionInfo.txt" )
TIME_COL <- "timepoint"
SAMPLE_COL <- "sample_id"
CELL_FDR_MAX <- 0.05
CELL_ABS_LOG2FC_MIN <- 0.25
REPLICATE_ABS_LOG2FC_MIN <- 0.10
PB_FDR_STRICT <- 0.10
PB_ABS_LOG2FC_STRICT <- 0.50
MIN_PCT <- 0.05
TOP_HEAT_PER_DIRECTION <- 25
TOP_LABEL_PER_DIRECTION <- 8
required_pkgs <- c( "Seurat", "SeuratObject", "Matrix", "dplyr", "tidyr", "tibble", "ggplot2", "patchwork", "edgeR", "limma", "openxlsx" )
missing_pkgs <- required_pkgs[ !vapply( required_pkgs, requireNamespace, logical(1), quietly = TRUE ) ]
if (length(missing_pkgs) > 0) {
  stop( paste0( "缺少以下R包：", paste(missing_pkgs, collapse = ", ") ), call. = FALSE )
}
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(edgeR)
})
msg <- function(...) {
  cat( sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "" )
}
capture_table <- function(x) {
  paste( capture.output( print( tibble::as_tibble(x), n = Inf, width = Inf ) ), collapse = "\n" )
}
get_norm <- function(obj) {
  assay <- DefaultAssay(obj)
  x <- tryCatch( Seurat::GetAssayData( obj, assay = assay, layer = "data" ), error = function(e) NULL )
  if (!is.null(x)) return(x)
  Seurat::GetAssayData( obj, assay = assay, slot = "data" )
}
safe_sheet <- function(wb, sheet_name, df) {
  openxlsx::addWorksheet( wb, sheetName = sheet_name )
  ddf <- tryCatch( as.data.frame(df), error = function(e) NULL )
  if ( is.null(ddf) || nrow(ddf) == 0 || ncol(ddf) == 0 ) {
    openxlsx::writeData( wb, sheet = sheet_name, x = data.frame( note = "No result under the predefined Step6A criterion.", stringsAsFactors = FALSE ) )
  } else {
    openxlsx::writeDataTable( wb, sheet = sheet_name, x = ddf, tableStyle = "TableStyleMedium2" )
  }
  openxlsx::freezePane( wb, sheet = sheet_name, firstRow = TRUE )
}
standardize_seurat_deg <- function(x) {
  x <- as.data.frame(x)
  x$gene_id <- rownames(x)
  fc_col <- intersect( c( "avg_log2FC", "avg_logFC" ), colnames(x) )[1]
  if ( length(fc_col) == 0 || is.na(fc_col) ) {
    stop( paste0( "FindMarkers结果中找不到avg_log2FC/avg_logFC。当前列：", paste(colnames(x), collapse = ", ") ), call. = FALSE )
  }
  out <- tibble::tibble( gene_id = as.character(x$gene_id), avg_log2FC = as.numeric(x[[fc_col]]), pct_Arm1 = if ("pct.1" %in% colnames(x)) as.numeric(x$pct.1) else NA_real_, pct_Arm2 = if ("pct.2" %in% colnames(x)) as.numeric(x$pct.2) else NA_real_, p_value = if ("p_val" %in% colnames(x)) as.numeric(x$p_val) else NA_real_, FDR = if ("p_val_adj" %in% colnames(x)) as.numeric(x$p_val_adj) else NA_real_ )
  out$BH_FDR <- stats::p.adjust( out$p_value, method = "BH" )
  out
}
run_cell_deg <- function( obj, ident1 = "Arm_1", ident2 = "Arm_2") {
  Seurat::Idents(obj) <- obj$step6A_branch
  res <- Seurat::FindMarkers( object = obj, ident.1 = ident1, ident.2 = ident2, assay = "RNA", test.use = "wilcox", min.pct = MIN_PCT, logfc.threshold = 0, only.pos = FALSE, verbose = FALSE )
  standardize_seurat_deg(res)
}
choose_umap <- function(obj) {
  prefs <- c( "umap.microglia.step3A", "umap.harmony", "umap.myeloid", "umap.preharmony" )
  avail <- names(obj@reductions)
  hit <- prefs[prefs %in% avail][1]
  if ( length(hit) == 0 || is.na(hit) ) {
    stop( paste0( "没有找到可复用UMAP。现有reductions：", paste(avail, collapse = ", ") ), call. = FALSE )
  }
  hit
}
msg("====================================================")
msg("GSE162610 Step 6A | 7 dpi branch-specific DEG")
msg("====================================================")
if (!file.exists(INPUT_STEP5CR)) {
  stop( paste0("找不到Step5C-R结果：", INPUT_STEP5CR), call. = FALSE )
}
if (!file.exists(INPUT_SEURAT)) {
  stop( paste0("找不到Step5B-R2 Seurat对象：", INPUT_SEURAT), call. = FALSE )
}
step5cr <- readRDS(INPUT_STEP5CR)
obj <- readRDS(INPUT_SEURAT)
DefaultAssay(obj) <- "RNA"
required_step5_fields <- c( "arm1_cells", "arm2_cells", "branchpoint" )
missing_step5 <- required_step5_fields[ !required_step5_fields %in% names(step5cr) ]
if (length(missing_step5) > 0) {
  stop( paste0( "Step5C-R结果缺少：", paste(missing_step5, collapse = ", ") ), call. = FALSE )
}
if ( !TIME_COL %in% colnames(obj[[]]) || !SAMPLE_COL %in% colnames(obj[[]]) ) {
  stop( "Seurat对象缺少timepoint或sample_id metadata。", call. = FALSE )
}
msg("Frozen Step5 branchpoint: ", step5cr$branchpoint)
timepoint <- as.character( obj[[TIME_COL, drop = TRUE]] )
names(timepoint) <- colnames(obj)
sample_id <- as.character( obj[[SAMPLE_COL, drop = TRUE]] )
names(sample_id) <- colnames(obj)
arm1_7 <- intersect( step5cr$arm1_cells, colnames(obj) )
arm1_7 <- arm1_7[ timepoint[arm1_7] == "7 dpi" ]
arm2_7 <- intersect( step5cr$arm2_cells, colnames(obj) )
arm2_7 <- arm2_7[ timepoint[arm2_7] == "7 dpi" ]
cells7 <- unique( c( arm1_7, arm2_7 ) )
obj7 <- subset( obj, cells = cells7 )
obj7$step6A_branch <- NA_character_
obj7$step6A_branch[colnames(obj7) %in% arm1_7] <- "Arm_1"
obj7$step6A_branch[colnames(obj7) %in% arm2_7] <- "Arm_2"
obj7$step6A_branch <- factor( obj7$step6A_branch, levels = c("Arm_2", "Arm_1") )
obj7$step6A_sample <- as.character( obj7[[SAMPLE_COL, drop = TRUE]] )
cell_summary <- obj7[[]] |> tibble::rownames_to_column("cell") |> dplyr::count( step6A_sample, step6A_branch, name = "n_cells" )
msg("7 dpi Arm_1 cells: ", length(arm1_7))
msg("7 dpi Arm_2 cells: ", length(arm2_7))
msg("7 dpi total branch cells: ", length(cells7))
msg("Biological samples: ", paste(sort(unique(obj7$step6A_sample)), collapse = ", "))
print(cell_summary)
sample_branch_check <- cell_summary |> dplyr::group_by(step6A_sample) |> dplyr::summarise( n_branches = dplyr::n_distinct(step6A_branch), min_cells = min(n_cells), .groups = "drop" )
if (any(sample_branch_check$n_branches < 2)) {
  stop( "至少一个7 dpi biological sample没有同时包含Arm_1和Arm_2。", call. = FALSE )
}
if (any(sample_branch_check$min_cells < 50)) {
  stop("At least one sample x branch group has <50 cells.", call. = FALSE)
}
UMAP_USE <- choose_umap(obj7)
msg("Using validated UMAP: ", UMAP_USE)
emb <- Seurat::Embeddings( obj7, reduction = UMAP_USE )
umap_df <- tibble::tibble( cell = rownames(emb), UMAP_1 = emb[, 1], UMAP_2 = emb[, 2], branch = as.character(obj7$step6A_branch), sample = as.character(obj7$step6A_sample) )
p1a <- ggplot2::ggplot( umap_df, ggplot2::aes( UMAP_1, UMAP_2, color = branch ) ) + ggplot2::geom_point( size = 0.65, alpha = 0.70 ) + ggplot2::labs( title = "7 dpi topology-defined terminal branches", subtitle = paste0( "Frozen from Step 5; branchpoint = ", step5cr$branchpoint ), color = "Branch" ) + ggplot2::theme_classic(base_size = 11)
p1b <- ggplot2::ggplot( umap_df, ggplot2::aes( UMAP_1, UMAP_2, color = branch ) ) + ggplot2::geom_point( size = 0.55, alpha = 0.70 ) + ggplot2::facet_wrap(~sample) + ggplot2::labs( title = "Branch structure reproduced by biological sample", color = "Branch" ) + ggplot2::theme_classic(base_size = 10)
grDevices::pdf( OUT_FIG1, width = 12, height = 6, useDingbats = FALSE )
print(p1a | p1b)
grDevices::dev.off()
msg("Running cell-level Wilcoxon: Arm_1 vs Arm_2...")
cell_deg <- run_cell_deg( obj7, ident1 = "Arm_1", ident2 = "Arm_2" )
cell_deg <- cell_deg |> dplyr::mutate( cell_significant = BH_FDR <= CELL_FDR_MAX & abs(avg_log2FC) >= CELL_ABS_LOG2FC_MIN, cell_direction = dplyr::case_when( avg_log2FC > 0 ~ "Arm_1", avg_log2FC < 0 ~ "Arm_2", TRUE ~ "Neutral" ) )
msg( "Cell-level significant genes: ", sum(cell_deg$cell_significant, na.rm = TRUE) )
msg( "  Arm_1-up: ", sum( cell_deg$cell_significant & cell_deg$avg_log2FC > 0, na.rm = TRUE ) )
msg( "  Arm_2-up: ", sum( cell_deg$cell_significant & cell_deg$avg_log2FC < 0, na.rm = TRUE ) )
samples7 <- sort( unique(obj7$step6A_sample) )
msg( "Running replicate-specific Wilcoxon analyses: ", paste(samples7, collapse = ", ") )
replicate_deg_list <- list()
for (smp in samples7) {
  cells_smp <- colnames(obj7)[ obj7$step6A_sample == smp ]
  obj_smp <- subset( obj7, cells = cells_smp )
  tab_smp <- table( obj_smp$step6A_branch )
  if ( !all(c("Arm_1", "Arm_2") %in% names(tab_smp)) || any(tab_smp[c("Arm_1", "Arm_2")] < 20) ) {
    stop( paste0( smp, "中至少一个branch细胞过少。" ), call. = FALSE )
  }
  dd <- run_cell_deg( obj_smp, ident1 = "Arm_1", ident2 = "Arm_2" )
  dd$sample <- smp
  replicate_deg_list[[smp]] <- dd
}
replicate_deg_long <- dplyr::bind_rows( replicate_deg_list )
replicate_fc_wide <- replicate_deg_long |> dplyr::select( gene_id, sample, avg_log2FC ) |> tidyr::pivot_wider( names_from = sample, values_from = avg_log2FC, names_prefix = "repFC_" )
rep_fc_cols <- grep( "^repFC_", colnames(replicate_fc_wide), value = TRUE )
replicate_consistency <- replicate_fc_wide |> dplyr::rowwise() |> dplyr::mutate( n_rep_available = sum( is.finite( c_across( dplyr::all_of(rep_fc_cols) ) ) ), same_direction = n_rep_available == length(rep_fc_cols) && ( all( c_across( dplyr::all_of(rep_fc_cols) ) > 0 ) || all( c_across( dplyr::all_of(rep_fc_cols) ) < 0 ) ), all_rep_absFC_ge_0.10 = n_rep_available == length(rep_fc_cols) && all( abs( c_across( dplyr::all_of(rep_fc_cols) ) ) >= REPLICATE_ABS_LOG2FC_MIN ), min_abs_repFC = ifelse( n_rep_available > 0, min( abs( c_across( dplyr::all_of(rep_fc_cols) ) ), na.rm = TRUE ), NA_real_ ) ) |> dplyr::ungroup()
msg("Building sample x branch raw-count pseudobulk...")
obj7$step6A_pb_group <- paste0( "g", obj7$step6A_sample, "__", as.character(obj7$step6A_branch) )
counts_mat <- tryCatch( Seurat::GetAssayData( obj7, assay = "RNA", layer = "counts" ), error = function(e) NULL )
if (is.null(counts_mat)) {
  counts_mat <- Seurat::GetAssayData( obj7, assay = "RNA", slot = "counts" )
}
if (ncol(counts_mat) != ncol(obj7)) {
  stop( "Raw-count matrix columns do not match the 7 dpi Seurat object.", call. = FALSE )
}
pb_group_vector <- as.character(obj7$step6A_pb_group)
names(pb_group_vector) <- colnames(obj7)
expected_groups <- sort( unique(pb_group_vector) )
pb_columns <- lapply( expected_groups, function(g) {
    cells_g <- names(pb_group_vector)[ pb_group_vector == g ]
    if (length(cells_g) == 0) {
      stop( paste0("Pseudobulk group has zero cells: ", g), call. = FALSE )
    }
    Matrix::rowSums( counts_mat[ , cells_g, drop = FALSE ] )
  }
)
pb_counts <- do.call( cbind, pb_columns )
rownames(pb_counts) <- rownames(counts_mat)
colnames(pb_counts) <- expected_groups
storage.mode(pb_counts) <- "integer"
pb_meta <- tibble::tibble( pb_group = colnames(pb_counts) ) |> tidyr::separate( pb_group, into = c("sample", "branch"), sep = "__", remove = FALSE ) |> dplyr::mutate( sample = sub("^g", "", sample), sample = factor( sample, levels = sort(unique(sample)) ), branch = factor( branch, levels = c("Arm_2", "Arm_1") ) ) |> as.data.frame()
rownames(pb_meta) <- pb_meta$pb_group
if ( anyNA(pb_meta$sample) || anyNA(pb_meta$branch) ) {
  stop( "Failed to parse exact pseudobulk sample/branch labels.", call. = FALSE )
}
msg("Pseudobulk columns: ", ncol(pb_counts))
msg( "Pseudobulk library sizes: ", paste( colnames(pb_counts), Matrix::colSums(pb_counts), sep = "=", collapse = "; " ) )
print(pb_meta)
msg("Running paired/sample-blocked edgeR: ~ sample + branch...")
y <- edgeR::DGEList( counts = pb_counts )
design <- stats::model.matrix( ~ sample + branch, data = pb_meta )
rownames(design) <- pb_meta$pb_group
keep <- edgeR::filterByExpr( y, design = design )
msg( "Pseudobulk genes retained by filterByExpr: ", sum(keep), "/", length(keep) )
y <- y[ keep, , keep.lib.sizes = FALSE ]
y <- edgeR::calcNormFactors(y)
pb_logcpm <- edgeR::cpm( y, log = TRUE, prior.count = 2 )
edgeR_method <- "QLF"
pb_fit_ok <- TRUE
pb_deg <- tryCatch( {
    y2 <- edgeR::estimateDisp( y, design, robust = TRUE )
    fit <- edgeR::glmQLFit( y2, design, robust = TRUE )
    coef_name <- grep( "^branchArm_1$", colnames(design), value = TRUE )
    if (length(coef_name) != 1) {
      stop( paste0( "无法唯一识别branchArm_1 coefficient。design列：", paste(colnames(design), collapse = ", ") ) )
    }
    qlf <- edgeR::glmQLFTest( fit, coef = coef_name )
    tt <- edgeR::topTags( qlf, n = Inf, sort.by = "PValue" )$table
    tt$gene_id <- rownames(tt)
    tibble::as_tibble(tt)
  },
  error = function(e) {
    pb_fit_ok <<- FALSE
    edgeR_method <<- "LRT_fallback"
    msg( "QLF failed; trying edgeR GLM-LRT fallback: ", conditionMessage(e) )
    y2 <- edgeR::estimateDisp( y, design, robust = FALSE )
    fit <- edgeR::glmFit( y2, design )
    coef_name <- grep( "^branchArm_1$", colnames(design), value = TRUE )
    if (length(coef_name) != 1) {
      stop( paste0( "无法唯一识别branchArm_1 coefficient。design列：", paste(colnames(design), collapse = ", ") ) )
    }
    lrt <- edgeR::glmLRT( fit, coef = coef_name )
    tt <- edgeR::topTags( lrt, n = Inf, sort.by = "PValue" )$table
    tt$gene_id <- rownames(tt)
    tibble::as_tibble(tt)
  }
)
pb_deg <- pb_deg |> dplyr::transmute( gene_id = as.character(gene_id), pb_logFC = as.numeric(logFC), pb_logCPM = as.numeric(logCPM), pb_PValue = as.numeric(PValue), pb_FDR = as.numeric(FDR), pb_direction = dplyr::case_when( pb_logFC > 0 ~ "Arm_1", pb_logFC < 0 ~ "Arm_2", TRUE ~ "Neutral" ), pb_strict = pb_FDR <= PB_FDR_STRICT & abs(pb_logFC) >= PB_ABS_LOG2FC_STRICT )
msg( "edgeR pseudobulk strict genes (FDR<=", PB_FDR_STRICT, ", |logFC|>=", PB_ABS_LOG2FC_STRICT, "): ", sum(pb_deg$pb_strict, na.rm = TRUE) )
pb_effect_list <- list()
for (smp in levels(pb_meta$sample)) {
  c1 <- rownames(pb_meta)[ pb_meta$sample == smp & pb_meta$branch == "Arm_1" ]
  c2 <- rownames(pb_meta)[ pb_meta$sample == smp & pb_meta$branch == "Arm_2" ]
  if ( length(c1) != 1 || length(c2) != 1 ) {
    next
  }
  genes_here <- rownames(pb_logcpm)
  pb_effect_list[[smp]] <- tibble::tibble( gene_id = genes_here, sample = smp, pb_within_sample_logCPM_diff = pb_logcpm[genes_here, c1] - pb_logcpm[genes_here, c2] )
}
pb_effect_long <- dplyr::bind_rows( pb_effect_list )
pb_effect_wide <- pb_effect_long |> dplyr::select( gene_id, sample, pb_within_sample_logCPM_diff ) |> tidyr::pivot_wider( names_from = sample, values_from = pb_within_sample_logCPM_diff, names_prefix = "pbDiff_" )
pb_diff_cols <- grep( "^pbDiff_", colnames(pb_effect_wide), value = TRUE )
pb_direction_check <- pb_effect_wide |> dplyr::rowwise() |> dplyr::mutate( pb_samples_same_direction = all( is.finite( c_across( dplyr::all_of(pb_diff_cols) ) ) ) && ( all( c_across( dplyr::all_of(pb_diff_cols) ) > 0 ) || all( c_across( dplyr::all_of(pb_diff_cols) ) < 0 ) ) ) |> dplyr::ungroup()
master <- cell_deg |> dplyr::rename( cell_avg_log2FC = avg_log2FC, cell_pct_Arm1 = pct_Arm1, cell_pct_Arm2 = pct_Arm2, cell_p = p_value, cell_Bonferroni = FDR, cell_BH_FDR = BH_FDR ) |> dplyr::left_join( replicate_consistency, by = "gene_id" ) |> dplyr::left_join( pb_deg, by = "gene_id" ) |> dplyr::left_join( pb_direction_check, by = "gene_id" ) |> dplyr::mutate( cell_significant = is.finite(cell_BH_FDR) & cell_BH_FDR <= CELL_FDR_MAX & abs(cell_avg_log2FC) >= CELL_ABS_LOG2FC_MIN, cell_direction = dplyr::case_when( cell_avg_log2FC > 0 ~ "Arm_1", cell_avg_log2FC < 0 ~ "Arm_2", TRUE ~ "Neutral" ), pb_same_as_cell = is.finite(pb_logFC) & sign(pb_logFC) == sign(cell_avg_log2FC), high_confidence_branch_marker = cell_significant & same_direction & all_rep_absFC_ge_0.10 & pb_same_as_cell & pb_samples_same_direction, strict_pseudobulk_supported = high_confidence_branch_marker & pb_strict, final_direction = dplyr::case_when( high_confidence_branch_marker & cell_avg_log2FC > 0 ~ "Arm_1-biased", high_confidence_branch_marker & cell_avg_log2FC < 0 ~ "Arm_2-biased", TRUE ~ "Not_high_confidence" ), technical_like = grepl( "^(Rpl|Rps|mt-|Gm[0-9]+)", gene_id, ignore.case = FALSE ) ) |> dplyr::arrange( dplyr::desc(high_confidence_branch_marker), cell_BH_FDR, dplyr::desc(abs(cell_avg_log2FC)) )
highconf <- master |> dplyr::filter( high_confidence_branch_marker )
strict_pb <- master |> dplyr::filter( strict_pseudobulk_supported )
arm1_high <- highconf |> dplyr::filter( final_direction == "Arm_1-biased" )
arm2_high <- highconf |> dplyr::filter( final_direction == "Arm_2-biased" )
msg( "High-confidence branch markers: ", nrow(highconf) )
msg( "  Arm_1-biased: ", nrow(arm1_high) )
msg( "  Arm_2-biased: ", nrow(arm2_high) )
msg( "Strict pseudobulk-supported: ", nrow(strict_pb) )
volcano_tbl <- master |> dplyr::mutate( volcano_group = dplyr::case_when( high_confidence_branch_marker & cell_avg_log2FC > 0 ~ "Arm_1-biased high-confidence", high_confidence_branch_marker & cell_avg_log2FC < 0 ~ "Arm_2-biased high-confidence", cell_significant ~ "Cell-level only", TRUE ~ "Not significant" ), neglog10FDR = -log10( pmax( cell_BH_FDR, 1e-300 ) ) ) |> dplyr::filter( is.finite(cell_avg_log2FC), is.finite(neglog10FDR) )
label_arm1 <- volcano_tbl |> dplyr::filter( volcano_group == "Arm_1-biased high-confidence", !technical_like ) |> dplyr::arrange( cell_BH_FDR, dplyr::desc(cell_avg_log2FC) ) |> dplyr::slice_head( n = TOP_LABEL_PER_DIRECTION )
label_arm2 <- volcano_tbl |> dplyr::filter( volcano_group == "Arm_2-biased high-confidence", !technical_like ) |> dplyr::arrange( cell_BH_FDR, cell_avg_log2FC ) |> dplyr::slice_head( n = TOP_LABEL_PER_DIRECTION )
label_tbl <- dplyr::bind_rows( label_arm1, label_arm2 )
p2 <- ggplot2::ggplot( volcano_tbl, ggplot2::aes( x = cell_avg_log2FC, y = neglog10FDR, color = volcano_group ) ) + ggplot2::geom_point( size = 0.75, alpha = 0.55 ) + ggplot2::geom_vline( xintercept = c( -CELL_ABS_LOG2FC_MIN, CELL_ABS_LOG2FC_MIN ), linetype = 2 ) + ggplot2::geom_hline( yintercept = -log10(CELL_FDR_MAX), linetype = 2 ) + ggplot2::geom_text( data = label_tbl, ggplot2::aes( label = gene_id ), size = 2.5, check_overlap = TRUE, show.legend = FALSE ) + ggplot2::labs( title = "7 dpi Arm_1 vs Arm_2 cell-level DEG", subtitle = "Positive log2FC = Arm_1-biased; negative = Arm_2-biased", x = "Average log2FC (Arm_1 - Arm_2)", y = "-log10(BH FDR)", color = "Evidence" ) + ggplot2::theme_classic(base_size = 11)
grDevices::pdf( OUT_FIG2, width = 8.5, height = 7, useDingbats = FALSE )
print(p2)
grDevices::dev.off()
if (length(rep_fc_cols) != 2) {
  stop( paste0( "当前Step6A预期恰好2个7 dpi biological samples，但得到：", length(rep_fc_cols) ), call. = FALSE )
}
rep_plot <- replicate_consistency |> dplyr::left_join( master |> dplyr::select( gene_id, high_confidence_branch_marker, final_direction, technical_like ), by = "gene_id" )
names(rep_plot)[ match(rep_fc_cols, names(rep_plot)) ] <- c( "sample1_log2FC", "sample2_log2FC" )
p3 <- ggplot2::ggplot( rep_plot, ggplot2::aes( x = sample1_log2FC, y = sample2_log2FC, color = final_direction ) ) + ggplot2::geom_hline( yintercept = 0, linetype = 2 ) + ggplot2::geom_vline( xintercept = 0, linetype = 2 ) + ggplot2::geom_abline( slope = 1, intercept = 0, linetype = 3 ) + ggplot2::geom_point( size = 0.8, alpha = 0.55 ) + ggplot2::labs( title = "7 dpi replicate concordance", subtitle = "Each axis is Arm_1 vs Arm_2 within one biological sample", x = paste0(rep_fc_cols[1], " avg_log2FC"), y = paste0(rep_fc_cols[2], " avg_log2FC"), color = "Final evidence" ) + ggplot2::theme_classic(base_size = 11)
grDevices::pdf( OUT_FIG3, width = 7.5, height = 7, useDingbats = FALSE )
print(p3)
grDevices::dev.off()
mds <- limma::plotMDS( y, plot = FALSE )
mds_df <- tibble::tibble( pb_group = colnames(y), Dim1 = mds$x, Dim2 = mds$y ) |> dplyr::left_join( pb_meta, by = "pb_group" )
p4 <- ggplot2::ggplot( mds_df, ggplot2::aes( x = Dim1, y = Dim2, shape = sample, label = branch ) ) + ggplot2::geom_point( size = 4 ) + ggplot2::geom_text( nudge_y = 0.08, size = 3.5 ) + ggplot2::labs( title = "7 dpi sample x branch pseudobulk MDS", subtitle = paste0( "edgeR sample-blocked design: ~ sample + branch | method = ", edgeR_method ), x = "Leading logFC dimension 1", y = "Leading logFC dimension 2", shape = "Biological sample" ) + ggplot2::theme_classic(base_size = 11)
grDevices::pdf( OUT_FIG4, width = 7, height = 6, useDingbats = FALSE )
print(p4)
grDevices::dev.off()
heat_source <- highconf |> dplyr::filter( !technical_like )
arm1_heat <- heat_source |> dplyr::filter( final_direction == "Arm_1-biased" ) |> dplyr::arrange( cell_BH_FDR, dplyr::desc(abs(cell_avg_log2FC)) ) |> dplyr::slice_head( n = TOP_HEAT_PER_DIRECTION ) |> dplyr::pull(gene_id)
arm2_heat <- heat_source |> dplyr::filter( final_direction == "Arm_2-biased" ) |> dplyr::arrange( cell_BH_FDR, dplyr::desc(abs(cell_avg_log2FC)) ) |> dplyr::slice_head( n = TOP_HEAT_PER_DIRECTION ) |> dplyr::pull(gene_id)
heat_genes <- unique( c( arm1_heat, arm2_heat ) )
if (length(heat_genes) >= 2) {
  heat_genes <- intersect( heat_genes, rownames(pb_logcpm) )
  hm <- pb_logcpm[ heat_genes, , drop = FALSE ]
  hm_z <- t( apply( hm, 1, function(v) {
        s <- stats::sd(v)
        if (!is.finite(s) || s == 0) {
          return(rep(0, length(v)))
        }
        (v - mean(v)) / s
      }
    )
  )
  heat_long <- as.data.frame(hm_z) |> tibble::rownames_to_column("gene_id") |> tidyr::pivot_longer( cols = -gene_id, names_to = "pb_group", values_to = "row_z" ) |> dplyr::left_join( pb_meta, by = "pb_group" ) |> dplyr::mutate( direction = ifelse( gene_id %in% arm1_heat, "Arm_1-biased", "Arm_2-biased" ), column_label = paste( sample, branch, sep = "\n" ) )
  gene_order <- unique( c( rev(arm1_heat), rev(arm2_heat) ) )
  heat_long$gene_id <- factor( heat_long$gene_id, levels = gene_order )
  heat_long$column_label <- factor( heat_long$column_label, levels = unique( paste( pb_meta$sample, pb_meta$branch, sep = "\n" ) ) )
  p5 <- ggplot2::ggplot( heat_long, ggplot2::aes( x = column_label, y = gene_id, fill = row_z ) ) + ggplot2::geom_tile() + ggplot2::facet_grid( direction ~ ., scales = "free_y", space = "free_y" ) + ggplot2::scale_fill_gradient2( midpoint = 0, name = "Row z-score" ) + ggplot2::labs( title = "High-confidence 7 dpi branch markers", subtitle = "Pseudobulk expression by biological sample x branch", x = NULL, y = "Gene" ) + ggplot2::theme_classic(base_size = 9) + ggplot2::theme( axis.text.y = ggplot2::element_text(size = 6) )
} else {
  p5 <- ggplot2::ggplot( tibble::tibble(x = 1, y = 1), ggplot2::aes(x, y) ) + ggplot2::annotate( "text", x = 1, y = 1, label = "Too few high-confidence non-technical genes for heatmap.", size = 5 ) + ggplot2::theme_void()
}
grDevices::pdf( OUT_FIG5, width = 8, height = 12, useDingbats = FALSE )
print(p5)
grDevices::dev.off()
evidence_counts <- tibble::tibble( category = c( "Cell-level DEG", "Replicate-direction consistent", "High-confidence branch markers", "Strict pseudobulk-supported" ), n_genes = c( sum(master$cell_significant, na.rm = TRUE), sum( master$cell_significant & master$same_direction & master$all_rep_absFC_ge_0.10, na.rm = TRUE ), nrow(highconf), nrow(strict_pb) ) )
p6 <- ggplot2::ggplot( evidence_counts, ggplot2::aes( x = category, y = n_genes ) ) + ggplot2::geom_col( width = 0.7 ) + ggplot2::geom_text( ggplot2::aes( label = n_genes ), vjust = -0.4, size = 4 ) + ggplot2::labs( title = "Step 6A branch-marker evidence cascade", subtitle = "No repair/scar annotation was used for DEG selection", x = NULL, y = "Genes" ) + ggplot2::theme_classic(base_size = 11) + ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 25, hjust = 1 ) )
grDevices::pdf( OUT_FIG6, width = 8, height = 6, useDingbats = FALSE )
print(p6)
grDevices::dev.off()
utils::write.csv( master, file.path( TABDIR, "01_all_genes_master_evidence.csv" ), row.names = FALSE )
utils::write.csv( highconf, file.path( TABDIR, "02_high_confidence_branch_markers.csv" ), row.names = FALSE )
utils::write.csv( arm1_high, file.path( TABDIR, "03_Arm1_biased_high_confidence_genes.csv" ), row.names = FALSE )
utils::write.csv( arm2_high, file.path( TABDIR, "04_Arm2_biased_high_confidence_genes.csv" ), row.names = FALSE )
utils::write.csv( strict_pb, file.path( TABDIR, "05_strict_pseudobulk_supported_genes.csv" ), row.names = FALSE )
utils::write.csv( replicate_deg_long, file.path( TABDIR, "06_replicate_specific_cell_DEG.csv" ), row.names = FALSE )
utils::write.csv( pb_deg, file.path( TABDIR, "07_edgeR_pseudobulk_DEG.csv" ), row.names = FALSE )
utils::write.csv( pb_effect_long, file.path( TABDIR, "08_within_sample_pseudobulk_effects.csv" ), row.names = FALSE )
parameter_tbl <- tibble::tibble( parameter = c( "Step5 branchpoint", "7 dpi Arm_1 cells", "7 dpi Arm_2 cells", "7 dpi biological samples", "Cell-level method", "Cell BH-FDR threshold", "Cell |avg_log2FC| threshold", "Per-sample minimum |avg_log2FC|", "Pseudobulk design", "Pseudobulk edgeR method", "Strict pseudobulk FDR threshold", "Strict pseudobulk |logFC| threshold", "Cell-level significant genes", "High-confidence branch markers", "Arm_1-biased high-confidence", "Arm_2-biased high-confidence", "Strict pseudobulk-supported genes" ), value = c( as.character(step5cr$branchpoint), length(arm1_7), length(arm2_7), paste(samples7, collapse = ", "), "Seurat Wilcoxon", CELL_FDR_MAX, CELL_ABS_LOG2FC_MIN, REPLICATE_ABS_LOG2FC_MIN, "~ sample + branch", edgeR_method, PB_FDR_STRICT, PB_ABS_LOG2FC_STRICT, sum(master$cell_significant, na.rm = TRUE), nrow(highconf), nrow(arm1_high), nrow(arm2_high), nrow(strict_pb) ) )
wb <- openxlsx::createWorkbook()
safe_sheet(wb, "Parameters", parameter_tbl)
safe_sheet(wb, "Cell_counts", cell_summary)
safe_sheet(wb, "All_genes", master)
safe_sheet(wb, "High_confidence", highconf)
safe_sheet(wb, "Arm1_biased", arm1_high)
safe_sheet(wb, "Arm2_biased", arm2_high)
safe_sheet(wb, "Strict_PB", strict_pb)
safe_sheet(wb, "Replicate_DEG", replicate_deg_long)
safe_sheet(wb, "Pseudobulk_DEG", pb_deg)
safe_sheet(wb, "PB_within_sample", pb_effect_long)
safe_sheet(wb, "Evidence_counts", evidence_counts)
openxlsx::saveWorkbook( wb, OUT_XLSX, overwrite = TRUE )
result_object <- list( step5_branchpoint = step5cr$branchpoint, arm1_7dpi_cells = arm1_7, arm2_7dpi_cells = arm2_7, cell_summary = cell_summary, cell_level_deg = cell_deg, replicate_specific_deg = replicate_deg_long, pseudobulk_counts = pb_counts, pseudobulk_metadata = pb_meta, pseudobulk_logCPM = pb_logcpm, pseudobulk_deg = pb_deg, pseudobulk_within_sample_effects = pb_effect_long, master_evidence = master, high_confidence_branch_markers = highconf, arm1_biased_high_confidence = arm1_high, arm2_biased_high_confidence = arm2_high, strict_pseudobulk_supported = strict_pb, evidence_counts = evidence_counts, thresholds = list( cell_BH_FDR = CELL_FDR_MAX, cell_abs_log2FC = CELL_ABS_LOG2FC_MIN, replicate_abs_log2FC = REPLICATE_ABS_LOG2FC_MIN, pb_FDR_strict = PB_FDR_STRICT, pb_abs_log2FC_strict = PB_ABS_LOG2FC_STRICT ) )
saveRDS( result_object, OUT_RDS, compress = TRUE )
top_arm1 <- arm1_high |> dplyr::filter(!technical_like) |> dplyr::select( gene_id, cell_avg_log2FC, cell_BH_FDR, min_abs_repFC, pb_logFC, pb_FDR, strict_pseudobulk_supported ) |> dplyr::slice_head(n = 30)
top_arm2 <- arm2_high |> dplyr::filter(!technical_like) |> dplyr::select( gene_id, cell_avg_log2FC, cell_BH_FDR, min_abs_repFC, pb_logFC, pb_FDR, strict_pseudobulk_supported ) |> dplyr::slice_head(n = 30)
report_lines <- c( "======================================================================", "GSE162610 | STEP 6A REPORT", "7 dpi Arm_1 vs Arm_2 branch-specific DEG", "======================================================================", "", "FORMAL POSITION:", "Step 6A is the first functional-comparison step after trajectory bifurcation.", "Repair/scar candidate scores were NOT used to select cells or DEGs.", "", "FROZEN TOPOLOGY:", paste0("Step5 branchpoint: ", step5cr$branchpoint), paste0("7 dpi Arm_1 cells: ", length(arm1_7)), paste0("7 dpi Arm_2 cells: ", length(arm2_7)), "", "BIOLOGICAL SAMPLE SUPPORT:", capture_table(cell_summary), "", "METHOD 1 | CELL-LEVEL:", "Seurat Wilcoxon Arm_1 vs Arm_2.", paste0( "Formal cell-level threshold: BH-FDR <= ", CELL_FDR_MAX, " and |avg_log2FC| >= ", CELL_ABS_LOG2FC_MIN ), paste0( "Cell-level significant genes: ", sum(master$cell_significant, na.rm = TRUE) ), "", "METHOD 2 | BIOLOGICAL-REPLICATE AUDIT:", "Arm_1 vs Arm_2 was repeated independently within each 7 dpi sample.", paste0( "High-confidence markers require the same direction in BOTH samples and ", "|sample-specific avg_log2FC| >= ", REPLICATE_ABS_LOG2FC_MIN, " in both." ), "", "METHOD 3 | SAMPLE-BLOCKED PSEUDOBULK:", "Raw counts were aggregated by biological sample x branch.", "edgeR design: ~ sample + branch.", paste0("edgeR method used: ", edgeR_method), "Because only two 7 dpi biological samples are available, pseudobulk p/FDR", "is treated as sensitivity support rather than the sole definition of truth.", "", "FINAL STEP6A EVIDENCE:", paste0( "High-confidence branch markers: ", nrow(highconf) ), paste0( "Arm_1-biased high-confidence genes: ", nrow(arm1_high) ), paste0( "Arm_2-biased high-confidence genes: ", nrow(arm2_high) ), paste0( "Strict pseudobulk-supported genes: ", nrow(strict_pb) ), "", "EVIDENCE CASCADE:", capture_table(evidence_counts), "", "TOP ARM_1-BIASED HIGH-CONFIDENCE GENES:", capture_table(top_arm1), "", "TOP ARM_2-BIASED HIGH-CONFIDENCE GENES:", capture_table(top_arm2), "", "INTERPRETATION RULE:", "At this step retain neutral names Arm_1-biased and Arm_2-biased.", "Do NOT rename them repair-biased/scar-biased until Step6 pathway/module", "analysis shows coherent functional separation.", "", "STATUS:", "After Step6A PASS, proceed to Step6B dynamic branch modules, then GO/KEGG/GSEA.", "", "FILES TO REVIEW:", OUT_FIG1, OUT_FIG2, OUT_FIG3, OUT_FIG4, OUT_FIG5, OUT_FIG6, OUT_XLSX, OUT_RDS, OUT_REP )
writeLines( report_lines, OUT_REP, useBytes = TRUE )
writeLines( capture.output(sessionInfo()), OUT_SESS, useBytes = TRUE )
msg("====================================================")
msg("STEP 6A COMPLETE")
msg("====================================================")
msg("7 dpi Arm_1 cells: ", length(arm1_7))
msg("7 dpi Arm_2 cells: ", length(arm2_7))
msg( "Cell-level significant genes: ", sum(master$cell_significant, na.rm = TRUE) )
msg( "High-confidence branch markers: ", nrow(highconf) )
msg( "Arm_1-biased high-confidence: ", nrow(arm1_high) )
msg( "Arm_2-biased high-confidence: ", nrow(arm2_high) )
msg( "Strict pseudobulk-supported genes: ", nrow(strict_pb) )
msg("")
msg("Key outputs:")
msg("1) figures/Fig6A_01_7dpi_Arm1_Arm2_UMAP.pdf")
msg("2) figures/Fig6A_02_celllevel_volcano.pdf")
msg("3) figures/Fig6A_03_replicate_concordance.pdf")
msg("4) figures/Fig6A_04_pseudobulk_MDS.pdf")
msg("5) figures/Fig6A_05_highconfidence_pseudobulk_heatmap.pdf")
msg("6) figures/Fig6A_06_evidence_cascade.pdf")
msg("7) Step6A_report.txt")
msg("")
msg("精确基因列表见 Table6A_7dpi_branch_specific_DEG.xlsx")
