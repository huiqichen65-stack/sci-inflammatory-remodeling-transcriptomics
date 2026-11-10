options(stringsAsFactors = FALSE)
set.seed(20260731)
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
STEP5BR_DIR <- file.path( PROJECT_ROOT, "04_results", "step5BR_GSE162610_injury_response_trajectory_refinement" )
STEP5BR2_DIR <- file.path( PROJECT_ROOT, "04_results", "step5BR2_GSE162610_single_root_branch_audit" )
STEP5CV2_DIR <- file.path( PROJECT_ROOT, "04_results", "step5C_GSE162610_Y29_branch_dependent_genes_v2" )
OUTDIR <- file.path( PROJECT_ROOT, "04_results", "step5CR_GSE162610_Y29_prebranch_progressive_divergence" )
dir.create(STEP5CV2_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
INPUT_CDS <- file.path( STEP5BR_DIR, "GSE162610_step5BR_injury_monocle3_cds.rds" )
INPUT_SEURAT <- file.path( STEP5BR2_DIR, "GSE162610_step5BR2_single_root_seurat.rds" )
CACHE_FILE <- file.path( STEP5CV2_DIR, "Step5C_v2_expensive_models_and_LRT_cache.RData" )
OUT_FIG1 <- file.path(OUTDIR, "Fig5CR_01_shared_prebranch_definition.pdf")
OUT_FIG2 <- file.path(OUTDIR, "Fig5CR_02_progressive_divergence_statistics.pdf")
OUT_FIG3 <- file.path(OUTDIR, "Fig5CR_03_progressive_divergence_heatmap.pdf")
OUT_FIG4 <- file.path(OUTDIR, "Fig5CR_04_top_gene_branchpoint_trajectories.pdf")
OUT_XLSX <- file.path(OUTDIR, "Table5CR_Y29_prebranch_progressive_divergence.xlsx")
OUT_RDS <- file.path(OUTDIR, "GSE162610_step5CR_results.rds")
OUT_REP <- file.path(OUTDIR, "Step5CR_report.txt")
OUT_SESS <- file.path(OUTDIR, "Step5CR_sessionInfo.txt")
TIME_COL <- "timepoint"
SAMPLE_COL <- "sample_id"
BRANCH_COL <- "step5BR2_branch_arm"
PT_COL <- "step5BR_pseudotime"
TRANSITION_COL <- "shared_transition_3dpi_highconf"
PRIMARY_BRANCHPOINT <- "Y_29"
PREBRANCH_GRAPH_RADIUS <- 4
Q_MAX <- 0.05
MIN_DIVERGENCE_GAIN <- 0.05
MIN_ABS_7DPI_DIFF <- 0.10
N_BINS <- 12
TOP_HEAT_PER_DIRECTION <- 25
TOP_SMOOTH_PER_DIRECTION <- 6
CORES <- if (.Platform$OS.type == "windows") 1L else 4L
required_pkgs <- c( "Seurat", "SeuratObject", "Matrix", "monocle3", "SingleCellExperiment", "SummarizedExperiment", "igraph", "splines", "dplyr", "tidyr", "tibble", "ggplot2", "patchwork", "openxlsx" )
missing_pkgs <- required_pkgs[ !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE) ]
if (length(missing_pkgs) > 0) {
  stop( paste0( "缺少R包：", paste(missing_pkgs, collapse = ", ") ), call. = FALSE )
}
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(monocle3)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(igraph)
  library(splines)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})
msg <- function(...) {
  cat( sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "" )
}
capture_table <- function(x) {
  paste( capture.output( print(tibble::as_tibble(x), n = Inf, width = Inf) ), collapse = "\n" )
}
get_norm <- function(obj) {
  assay <- DefaultAssay(obj)
  x <- tryCatch( Seurat::GetAssayData(obj, assay = assay, layer = "data"), error = function(e) NULL )
  if (!is.null(x)) return(x)
  Seurat::GetAssayData(obj, assay = assay, slot = "data")
}
get_closest_vertex <- function(cds, pg) {
  aux <- cds@principal_graph_aux[["UMAP"]]
  raw <- aux$pr_graph_cell_proj_closest_vertex
  if (is.null(raw)) {
    stop("Monocle3 closest-principal-vertex mapping not found.", call. = FALSE)
  }
  if (is.matrix(raw) || is.data.frame(raw)) {
    vv <- raw[, 1]
    names(vv) <- rownames(raw)
  } else {
    vv <- raw
  }
  if (is.null(names(vv))) {
    names(vv) <- colnames(cds)
  }
  vv_chr <- as.character(vv)
  pg_names <- igraph::V(pg)$name
  if (all(vv_chr %in% pg_names)) {
    return(stats::setNames(vv_chr, names(vv)))
  }
  y <- paste0("Y_", vv_chr)
  if (all(y %in% pg_names)) {
    return(stats::setNames(y, names(vv)))
  }
  idx <- suppressWarnings(as.integer(vv_chr))
  if ( all(!is.na(idx)) && all(idx >= 1) && all(idx <= igraph::vcount(pg)) ) {
    return( stats::setNames( pg_names[idx], names(vv) ) )
  }
  stop("Cannot translate Monocle3 closest-vertex IDs.", call. = FALSE)
}
standardize_compare <- function(x, prefix) {
  x <- tibble::as_tibble(x, rownames = "gene_id_from_rownames")
  if (!"gene_id" %in% colnames(x)) {
    if ("id" %in% colnames(x)) {
      x$gene_id <- as.character(x$id)
    } else {
      x$gene_id <- as.character(x$gene_id_from_rownames)
    }
  }
  p_candidates <- c("p_value", "pval", "p.value", "lrtest_p_value")
  p_col <- p_candidates[p_candidates %in% colnames(x)][1]
  if (length(p_col) == 0 || is.na(p_col)) {
    stop( paste0( "compare_models没有识别到p值列。当前列：", paste(colnames(x), collapse = ", ") ), call. = FALSE )
  }
  pp <- as.numeric(x[[p_col]])
  q_candidates <- c("q_value", "qval", "q.value")
  q_col <- q_candidates[q_candidates %in% colnames(x)][1]
  if (length(q_col) > 0 && !is.na(q_col)) {
    qq <- as.numeric(x[[q_col]])
  } else {
    qq <- stats::p.adjust(pp, method = "BH")
  }
  out <- tibble::tibble( gene_id = as.character(x$gene_id) )
  out[[paste0(prefix, "_p_value")]] <- pp
  out[[paste0(prefix, "_q_value")]] <- qq
  out
}
safe_sheet <- function(wb, name, df) {
  openxlsx::addWorksheet(wb, name)
  ddf <- tryCatch(as.data.frame(df), error = function(e) NULL)
  if ( is.null(ddf) || nrow(ddf) == 0 || ncol(ddf) == 0 ) {
    openxlsx::writeData( wb, name, data.frame( note = "No result under the predefined Step5C-R criterion." ) )
  } else {
    openxlsx::writeDataTable( wb, name, ddf, tableStyle = "TableStyleMedium2" )
  }
  openxlsx::freezePane(wb, name, firstRow = TRUE)
}
row_z <- function(mat) {
  mu <- rowMeans(mat, na.rm = TRUE)
  ss <- apply(mat, 1, stats::sd, na.rm = TRUE)
  ss[!is.finite(ss) | ss == 0] <- 1
  z <- sweep(mat, 1, mu, "-")
  sweep(z, 1, ss, "/")
}
null_plot <- function(txt) {
  ggplot2::ggplot( tibble::tibble(x = 1, y = 1), ggplot2::aes(x, y) ) + ggplot2::annotate("text", x = 1, y = 1, label = txt, size = 5) + ggplot2::theme_void()
}
msg("====================================================")
msg("GSE162610 Step 5C-R | Y_29 prebranch divergence audit")
msg("====================================================")
core_objects <- c( "full_fit", "null_fit", "additive_fit", "branch_dep_tbl", "shape_tbl" )
have_memory <- all( vapply( core_objects, exists, logical(1), envir = .GlobalEnv ) )
if (have_memory) {
  msg("检测到当前R session中的Step5C-v2耗时对象。")
  msg("立即保存缓存，避免以后重跑3291-gene negative-binomial models...")
  save( full_fit, null_fit, additive_fit, branch_dep_tbl, shape_tbl, file = CACHE_FILE, compress = TRUE )
  msg("缓存已保存：", CACHE_FILE)
} else if (file.exists(CACHE_FILE)) {
  msg("当前session对象不完整，正在读取已保存缓存...")
  load( CACHE_FILE, envir = .GlobalEnv )
} else {
  stop(paste0("Missing Step5C model objects: ", paste(core_objects, collapse = ", ")), call. = FALSE)
}
formal_genes <- branch_dep_tbl |> dplyr::filter( is.finite(branch_dependence_q_value), branch_dependence_q_value <= 0.05 ) |> dplyr::pull(gene_id) |> unique()
msg( "Step5C-v2正式branch-dependent genes带入本步：", length(formal_genes) )
if (length(formal_genes) < 10) {
  stop("正式branch-dependent genes少于10个。", call. = FALSE)
}
if (!file.exists(INPUT_CDS)) {
  stop(paste0("找不到：", INPUT_CDS), call. = FALSE)
}
if (!file.exists(INPUT_SEURAT)) {
  stop(paste0("找不到：", INPUT_SEURAT), call. = FALSE)
}
cds <- readRDS(INPUT_CDS)
obj <- readRDS(INPUT_SEURAT)
DefaultAssay(obj) <- "RNA"
req <- c( TIME_COL, SAMPLE_COL, BRANCH_COL, PT_COL, TRANSITION_COL )
miss <- req[!req %in% colnames(obj[[]])]
if (length(miss) > 0) {
  stop( paste0( "缺少metadata：", paste(miss, collapse = ", ") ), call. = FALSE )
}
pg <- monocle3::principal_graph(cds)[["UMAP"]]
if ( is.null(pg) || !PRIMARY_BRANCHPOINT %in% igraph::V(pg)$name ) {
  stop("Step5B-R graph中找不到Y_29。", call. = FALSE)
}
cds_cells <- colnames(cds)
obj_cells <- colnames(obj)
if (!all(cds_cells %in% obj_cells)) {
  stop("部分CDS细胞不在Step5B-R2 Seurat对象中。", call. = FALSE)
}
branch_arm <- as.character(obj[[BRANCH_COL, drop = TRUE]])
names(branch_arm) <- obj_cells
timepoint <- as.character(obj[[TIME_COL, drop = TRUE]])
names(timepoint) <- obj_cells
sample_id <- as.character(obj[[SAMPLE_COL, drop = TRUE]])
names(sample_id) <- obj_cells
pt <- as.numeric(obj[[PT_COL, drop = TRUE]])
names(pt) <- obj_cells
shared_transition <- as.logical(obj[[TRANSITION_COL, drop = TRUE]])
names(shared_transition) <- obj_cells
closest_vertex <- get_closest_vertex(cds, pg)
closest_vertex <- closest_vertex[cds_cells]
dmat <- igraph::distances( pg, v = PRIMARY_BRANCHPOINT, to = igraph::V(pg), mode = "all" )
vertex_distance <- as.numeric( dmat[1, igraph::V(pg)$name] )
names(vertex_distance) <- igraph::V(pg)$name
cell_distance <- vertex_distance[ closest_vertex[cds_cells] ]
names(cell_distance) <- cds_cells
prebranch_cells <- cds_cells[ timepoint[cds_cells] == "3 dpi" & shared_transition[cds_cells] & branch_arm[cds_cells] %in% c("Pre_branch_root_side", "Branchpoint_vertex") & cell_distance[cds_cells] <= PREBRANCH_GRAPH_RADIUS ]
arm1_cells <- cds_cells[ branch_arm[cds_cells] == "Arm_1" & timepoint[cds_cells] %in% c("3 dpi", "7 dpi") ]
arm2_cells <- cds_cells[ branch_arm[cds_cells] == "Arm_2" & timepoint[cds_cells] %in% c("3 dpi", "7 dpi") ]
msg("Shared prebranch cells: ", length(prebranch_cells))
msg("Arm_1 downstream cells: ", length(arm1_cells))
msg("Arm_2 downstream cells: ", length(arm2_cells))
if (length(prebranch_cells) < 50) {
  stop(paste0("High-confidence 3dpi prebranch cells near Y_29: ", length(prebranch_cells), " (<50)."), call. = FALSE)
}
analysis_cells <- unique( c( prebranch_cells, arm1_cells, arm2_cells ) )
state <- rep(NA_character_, length(analysis_cells))
names(state) <- analysis_cells
state[prebranch_cells] <- "Shared_prebranch"
state[arm1_cells] <- "Arm_1"
state[arm2_cells] <- "Arm_2"
signed_progress <- rep(NA_real_, length(analysis_cells))
names(signed_progress) <- analysis_cells
signed_progress[prebranch_cells] <- -cell_distance[prebranch_cells]
signed_progress[arm1_cells] <- cell_distance[arm1_cells]
signed_progress[arm2_cells] <- cell_distance[arm2_cells]
signed_progress <- signed_progress / max(abs(signed_progress), na.rm = TRUE)
post_progress <- pmax(signed_progress, 0)
branch_contrast <- rep(0, length(analysis_cells))
names(branch_contrast) <- analysis_cells
branch_contrast[arm1_cells] <- 0.5
branch_contrast[arm2_cells] <- -0.5
analysis_meta <- tibble::tibble( cell = analysis_cells, state = state[analysis_cells], timepoint = timepoint[analysis_cells], sample = sample_id[analysis_cells], retained_step5BR_pseudotime = pt[analysis_cells], graph_distance_from_Y29 = cell_distance[analysis_cells], signed_progress = signed_progress[analysis_cells], post_progress = post_progress[analysis_cells], branch_contrast = branch_contrast[analysis_cells] )
analysis_summary <- analysis_meta |> dplyr::count( state, timepoint, sample, name = "n_cells" )
umap <- SingleCellExperiment::reducedDim(cds, "UMAP")
umap_df <- tibble::tibble( cell = rownames(umap), UMAP_1 = umap[, 1], UMAP_2 = umap[, 2] ) |> dplyr::mutate( display_state = dplyr::case_when( cell %in% prebranch_cells ~ "Shared_prebranch", cell %in% arm1_cells ~ "Arm_1", cell %in% arm2_cells ~ "Arm_2", TRUE ~ "Other" ) )
p1a <- ggplot2::ggplot( umap_df, ggplot2::aes( UMAP_1, UMAP_2, color = display_state ) ) + ggplot2::geom_point(size = 0.45, alpha = 0.55) + ggplot2::scale_color_manual( values = c( "Other" = "#D9D9D9", "Shared_prebranch" = "#6A3D9A", "Arm_1" = "#1F78B4", "Arm_2" = "#E31A1C" ) ) + ggplot2::labs( title = "Y_29 shared prebranch and corrected downstream arms", subtitle = paste0( "Prebranch: 3 dpi high-confidence shared transition, graph radius <= ", PREBRANCH_GRAPH_RADIUS ) ) + ggplot2::theme_classic(base_size = 11)
p1b <- ggplot2::ggplot( analysis_meta, ggplot2::aes( x = signed_progress, fill = state ) ) + ggplot2::geom_histogram( bins = 50, position = "identity", alpha = 0.55 ) + ggplot2::geom_vline(xintercept = 0, linetype = 2) + ggplot2::labs( title = "Branchpoint-centered graph progression", subtitle = "Negative = prebranch; positive = downstream from Y_29", x = "Scaled signed graph progress", y = "Cells" ) + ggplot2::theme_classic(base_size = 11)
p1c <- ggplot2::ggplot( analysis_summary, ggplot2::aes( x = state, y = n_cells, fill = timepoint ) ) + ggplot2::geom_col() + ggplot2::facet_wrap(~sample, scales = "free_y") + ggplot2::labs( title = "Biological-sample support", x = NULL, y = "Cells" ) + ggplot2::theme_classic(base_size = 9) + ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 25, hjust = 1 ) )
p1d <- ggplot2::ggplot( analysis_meta |> dplyr::filter(state == "Shared_prebranch"), ggplot2::aes( x = graph_distance_from_Y29, y = retained_step5BR_pseudotime ) ) + ggplot2::geom_point(size = 0.7, alpha = 0.5) + ggplot2::labs( title = "Local 3 dpi shared-transition cells before Y_29", x = "Graph distance to Y_29", y = "Retained Step5B-R pseudotime" ) + ggplot2::theme_classic(base_size = 10)
grDevices::pdf( OUT_FIG1, width = 14, height = 10, useDingbats = FALSE )
print( (p1a | p1b) / (p1c | p1d) + patchwork::plot_annotation( title = "Step 5C-R | Shared activated/DAM-like prebranch definition" ) )
grDevices::dev.off()
formal_genes <- intersect( formal_genes, rownames(cds) )
cds_local <- cds[ formal_genes, analysis_cells, drop = FALSE ]
meta_order <- analysis_meta[ match( colnames(cds_local), analysis_meta$cell ), , drop = FALSE ]
if ( !identical( as.character(meta_order$cell), colnames(cds_local) ) ) {
  stop("Local metadata order mismatch.", call. = FALSE)
}
SummarizedExperiment::colData(cds_local)$step5CR_sample <- factor(meta_order$sample)
SummarizedExperiment::colData(cds_local)$step5CR_signed_progress <- meta_order$signed_progress
SummarizedExperiment::colData(cds_local)$step5CR_post_progress <- meta_order$post_progress
SummarizedExperiment::colData(cds_local)$step5CR_branch_contrast <- meta_order$branch_contrast
FULL_FORMULA <- paste0( "~ step5CR_sample + ", "splines::ns(step5CR_signed_progress, df = 3) + ", "step5CR_branch_contrast + ", "step5CR_branch_contrast:splines::ns(step5CR_post_progress, df = 2)" )
ADDITIVE_FORMULA <- paste0( "~ step5CR_sample + ", "splines::ns(step5CR_signed_progress, df = 3) + ", "step5CR_branch_contrast" )
NULL_FORMULA <- paste0( "~ step5CR_sample + ", "splines::ns(step5CR_signed_progress, df = 3)" )
msg("Genes modeled: ", nrow(cds_local))
msg("Cells modeled: ", ncol(cds_local))
msg("Windows-safe cores: ", CORES)
msg("FULL: ", FULL_FORMULA)
msg("ADDITIVE: ", ADDITIVE_FORMULA)
msg("NULL: ", NULL_FORMULA)
msg("Fitting FULL prebranch-anchored model...")
fit_full_CR <- monocle3::fit_models( cds_local, model_formula_str = FULL_FORMULA, expression_family = "negbinomial", clean_model = TRUE, cores = CORES )
msg("Fitting ADDITIVE model...")
fit_add_CR <- monocle3::fit_models( cds_local, model_formula_str = ADDITIVE_FORMULA, expression_family = "negbinomial", clean_model = TRUE, cores = CORES )
msg("Fitting NULL common-trajectory model...")
fit_null_CR <- monocle3::fit_models( cds_local, model_formula_str = NULL_FORMULA, expression_family = "negbinomial", clean_model = TRUE, cores = CORES )
fit_status <- tibble::tibble( gene_id = formal_genes, full_status = as.character( fit_full_CR$status[ match(formal_genes, fit_full_CR$gene_id) ] ), additive_status = as.character( fit_add_CR$status[ match(formal_genes, fit_add_CR$gene_id) ] ), null_status = as.character( fit_null_CR$status[ match(formal_genes, fit_null_CR$gene_id) ] ) )
msg( "FULL model OK: ", sum(fit_status$full_status == "OK", na.rm = TRUE), "/", nrow(fit_status) )
msg( "ADDITIVE model OK: ", sum(fit_status$additive_status == "OK", na.rm = TRUE), "/", nrow(fit_status) )
msg( "NULL model OK: ", sum(fit_status$null_status == "OK", na.rm = TRUE), "/", nrow(fit_status) )
msg("PRIMARY progressive-shape test: FULL vs ADDITIVE...")
shape_raw_CR <- monocle3::compare_models( fit_full_CR, fit_add_CR )
shape_test <- standardize_compare( shape_raw_CR, "progressive_shape" )
msg("SECONDARY anchored-branch test: FULL vs NULL...")
branch_raw_CR <- monocle3::compare_models( fit_full_CR, fit_null_CR )
branch_test <- standardize_compare( branch_raw_CR, "anchored_branch" )
norm_mat <- get_norm(obj)
formal_genes <- intersect( formal_genes, rownames(norm_mat) )
post_meta <- analysis_meta |> dplyr::filter( state %in% c("Arm_1", "Arm_2") )
q25 <- stats::quantile( post_meta$post_progress, 0.25, na.rm = TRUE )
q75 <- stats::quantile( post_meta$post_progress, 0.75, na.rm = TRUE )
early_arm1 <- post_meta$cell[ post_meta$state == "Arm_1" & post_meta$post_progress <= q25 ]
early_arm2 <- post_meta$cell[ post_meta$state == "Arm_2" & post_meta$post_progress <= q25 ]
late_arm1 <- post_meta$cell[ post_meta$state == "Arm_1" & post_meta$post_progress >= q75 ]
late_arm2 <- post_meta$cell[ post_meta$state == "Arm_2" & post_meta$post_progress >= q75 ]
mean_expr <- function(genes, cells) {
  if (length(cells) == 0) {
    return(rep(NA_real_, length(genes)))
  }
  as.numeric( Matrix::rowMeans( norm_mat[ genes, cells, drop = FALSE ] ) )
}
early_diff <- mean_expr(formal_genes, early_arm1) - mean_expr(formal_genes, early_arm2)
late_diff <- mean_expr(formal_genes, late_arm1) - mean_expr(formal_genes, late_arm2)
seven_arm1 <- arm1_cells[ timepoint[arm1_cells] == "7 dpi" ]
seven_arm2 <- arm2_cells[ timepoint[arm2_cells] == "7 dpi" ]
seven_diff <- mean_expr(formal_genes, seven_arm1) - mean_expr(formal_genes, seven_arm2)
divergence_tbl <- tibble::tibble( gene_id = formal_genes, early_arm1_minus_arm2 = early_diff, late_arm1_minus_arm2 = late_diff, abs_early_difference = abs(early_diff), abs_late_difference = abs(late_diff), divergence_gain = abs(late_diff) - abs(early_diff), mean_7dpi_arm1_minus_arm2 = seven_diff )
seven_samples <- sort( unique( sample_id[ c(seven_arm1, seven_arm2) ] ) )
replicate_effects <- dplyr::bind_rows( lapply( seven_samples, function(smp) {
      c1 <- seven_arm1[ sample_id[seven_arm1] == smp ]
      c2 <- seven_arm2[ sample_id[seven_arm2] == smp ]
      if (length(c1) == 0 || length(c2) == 0) {
        return(NULL)
      }
      dd <- mean_expr(formal_genes, c1) - mean_expr(formal_genes, c2)
      tibble::tibble( gene_id = formal_genes, sample = smp, arm1_minus_arm2 = dd )
    }
  )
)
replicate_consistency <- replicate_effects |> dplyr::group_by(gene_id) |> dplyr::summarise( n_samples = dplyr::n(), n_positive = sum(arm1_minus_arm2 > 0), n_negative = sum(arm1_minus_arm2 < 0), direction_consistent = n_positive == n_samples | n_negative == n_samples, .groups = "drop" )
master <- tibble::tibble( gene_id = formal_genes ) |> dplyr::left_join( branch_dep_tbl |> dplyr::select( gene_id, step5C_v2_p = branch_dependence_p_value, step5C_v2_q = branch_dependence_q_value ), by = "gene_id" ) |> dplyr::left_join( shape_test, by = "gene_id" ) |> dplyr::left_join( branch_test, by = "gene_id" ) |> dplyr::left_join( divergence_tbl, by = "gene_id" ) |> dplyr::left_join( replicate_consistency, by = "gene_id" ) |> dplyr::mutate( progressive_shape_significant = is.finite(progressive_shape_q_value) & progressive_shape_q_value <= Q_MAX, anchored_branch_significant = is.finite(anchored_branch_q_value) & anchored_branch_q_value <= Q_MAX, terminal_direction = dplyr::case_when( mean_7dpi_arm1_minus_arm2 > 0 ~ "Arm_1", mean_7dpi_arm1_minus_arm2 < 0 ~ "Arm_2", TRUE ~ "Neutral" ), strong_progressive = progressive_shape_significant & divergence_gain >= MIN_DIVERGENCE_GAIN & abs(mean_7dpi_arm1_minus_arm2) >= MIN_ABS_7DPI_DIFF & direction_consistent ) |> dplyr::arrange( progressive_shape_q_value, dplyr::desc(divergence_gain) )
progressive_tbl <- master |> dplyr::filter(progressive_shape_significant)
strong_tbl <- master |> dplyr::filter(strong_progressive)
msg( "Progressive-shape genes q<=0.05: ", nrow(progressive_tbl) )
msg( "Anchored branch-deviation genes q<=0.05: ", sum(master$anchored_branch_significant, na.rm = TRUE) )
msg( "Strong progressive genes: ", nrow(strong_tbl) )
msg( "Minimum progressive-shape q: ", signif( min(master$progressive_shape_q_value, na.rm = TRUE), 6 ) )
plot_source <- if (nrow(strong_tbl) >= 10) {
  strong_tbl
} else {
  progressive_tbl
}
arm1_top <- plot_source |> dplyr::filter(terminal_direction == "Arm_1") |> dplyr::arrange( progressive_shape_q_value, dplyr::desc(divergence_gain) ) |> dplyr::slice_head(n = TOP_HEAT_PER_DIRECTION) |> dplyr::pull(gene_id)
arm2_top <- plot_source |> dplyr::filter(terminal_direction == "Arm_2") |> dplyr::arrange( progressive_shape_q_value, dplyr::desc(divergence_gain) ) |> dplyr::slice_head(n = TOP_HEAT_PER_DIRECTION) |> dplyr::pull(gene_id)
heat_genes <- unique( c(arm1_top, arm2_top) )
smooth_genes <- unique( c( head(arm1_top, TOP_SMOOTH_PER_DIRECTION), head(arm2_top, TOP_SMOOTH_PER_DIRECTION) ) )
plot_bin_tbl <- tibble::tibble()
if (length(heat_genes) >= 2) {
  pre_mean <- mean_expr( heat_genes, prebranch_cells )
  pre_tbl <- tibble::tibble( gene_id = heat_genes, state = "Shared_prebranch", progress_bin = 0L, progress_mid = 0, mean_expression = pre_mean )
  post_meta2 <- post_meta |> dplyr::mutate( progress_bin = cut( post_progress, breaks = seq( 0, max(post_progress, na.rm = TRUE), length.out = N_BINS + 1 ), include.lowest = TRUE, labels = FALSE ) )
  post_list <- list()
  for (aa in c("Arm_1", "Arm_2")) {
    bins_aa <- sort( unique( post_meta2$progress_bin[ post_meta2$state == aa ] ) )
    for (bb in bins_aa) {
      cc <- post_meta2$cell[ post_meta2$state == aa & post_meta2$progress_bin == bb ]
      if (length(cc) == 0) next
      vv <- mean_expr(heat_genes, cc)
      pmid <- median( post_meta2$post_progress[ post_meta2$state == aa & post_meta2$progress_bin == bb ], na.rm = TRUE )
      post_list[[paste0(aa, "_", bb)]] <- tibble::tibble( gene_id = heat_genes, state = aa, progress_bin = as.integer(bb), progress_mid = pmid, mean_expression = vv )
    }
  }
  plot_bin_tbl <- dplyr::bind_rows( pre_tbl, dplyr::bind_rows(post_list) )
}
p2a <- ggplot2::ggplot( master, ggplot2::aes( x = divergence_gain, y = -log10( pmax(progressive_shape_q_value, 1e-300) ) ) ) + ggplot2::geom_point(size = 0.75, alpha = 0.45) + ggplot2::geom_hline( yintercept = -log10(Q_MAX), linetype = 2 ) + ggplot2::geom_vline( xintercept = MIN_DIVERGENCE_GAIN, linetype = 2 ) + ggplot2::labs( title = "Prebranch-anchored progressive divergence", subtitle = "Primary LRT: FULL vs ADDITIVE", x = "Late-minus-early absolute Arm1/Arm2 difference", y = "-log10(progressive-shape q)" ) + ggplot2::theme_classic(base_size = 11)
count_tbl <- tibble::tibble( category = c( "Step5C-v2 branch-dependent", "Anchored branch-deviation", "Progressive-shape", "Strong progressive" ), n_genes = c( length(formal_genes), sum(master$anchored_branch_significant, na.rm = TRUE), nrow(progressive_tbl), nrow(strong_tbl) ) )
p2b <- ggplot2::ggplot( count_tbl, ggplot2::aes( x = category, y = n_genes ) ) + ggplot2::geom_col(width = 0.7) + ggplot2::geom_text( ggplot2::aes(label = n_genes), vjust = -0.35, size = 3.5 ) + ggplot2::labs( title = "Evidence cascade", x = NULL, y = "Genes" ) + ggplot2::theme_classic(base_size = 10) + ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 25, hjust = 1 ) )
p2c <- ggplot2::ggplot( master, ggplot2::aes( x = abs_early_difference, y = abs_late_difference ) ) + ggplot2::geom_point(size = 0.7, alpha = 0.4) + ggplot2::geom_abline( slope = 1, intercept = 0, linetype = 2 ) + ggplot2::labs( title = "Does Arm1/Arm2 separation increase after Y_29?", x = "Early downstream absolute difference", y = "Late downstream absolute difference" ) + ggplot2::theme_classic(base_size = 11)
rep_counts <- master |> dplyr::filter(progressive_shape_significant) |> dplyr::count( terminal_direction, direction_consistent, name = "n_genes" )
if (nrow(rep_counts) > 0) {
  p2d <- ggplot2::ggplot( rep_counts, ggplot2::aes( x = terminal_direction, y = n_genes, fill = direction_consistent ) ) + ggplot2::geom_col() + ggplot2::labs( title = "7 dpi replicate direction", x = "Terminal direction", y = "Progressive genes", fill = "Consistent\nin both samples" ) + ggplot2::theme_classic(base_size = 10)
} else {
  p2d <- null_plot( "No progressive-shape genes passed q <= 0.05." )
}
grDevices::pdf( OUT_FIG2, width = 14, height = 10, useDingbats = FALSE )
print( (p2a | p2b) / (p2c | p2d) + patchwork::plot_annotation( title = "Step 5C-R | Y_29 progressive-divergence statistics" ) )
grDevices::dev.off()
if ( length(heat_genes) >= 2 && nrow(plot_bin_tbl) > 0 ) {
  heat_wide <- plot_bin_tbl |> dplyr::mutate( state_bin = paste0( state, "__", sprintf("%02d", progress_bin) ) ) |> dplyr::select( gene_id, state_bin, mean_expression ) |> tidyr::pivot_wider( names_from = state_bin, values_from = mean_expression )
  heat_mat <- as.matrix( heat_wide[ , setdiff( colnames(heat_wide), "gene_id" ), drop = FALSE ] )
  rownames(heat_mat) <- heat_wide$gene_id
  heat_mat[is.na(heat_mat)] <- 0
  heat_z <- row_z(heat_mat)
  heat_long <- as.data.frame(heat_z) |> tibble::rownames_to_column("gene_id") |> tidyr::pivot_longer( cols = -gene_id, names_to = "state_bin", values_to = "z_expression" ) |> tidyr::separate( state_bin, into = c("state", "bin"), sep = "__" ) |> dplyr::mutate( bin = as.integer(bin), direction = ifelse( gene_id %in% arm1_top, "Arm_1-biased", "Arm_2-biased" ) )
  gene_order <- unique( c( rev(arm1_top), rev(arm2_top) ) )
  heat_long$gene_id <- factor( heat_long$gene_id, levels = gene_order )
  heat_long$state <- factor( heat_long$state, levels = c( "Shared_prebranch", "Arm_1", "Arm_2" ) )
  p3 <- ggplot2::ggplot( heat_long, ggplot2::aes( x = bin, y = gene_id, fill = z_expression ) ) + ggplot2::geom_tile() + ggplot2::facet_grid( direction ~ state, scales = "free", space = "free" ) + ggplot2::scale_fill_gradient2( midpoint = 0, name = "Row z-score" ) + ggplot2::labs( title = "Shared prebranch -> Arm_1 / Arm_2 gene-program separation", x = "Branchpoint-centered progression bin", y = "Gene" ) + ggplot2::theme_classic(base_size = 8) + ggplot2::theme( axis.text.y = ggplot2::element_text(size = 5.5) )
} else {
  p3 <- null_plot( "No significant prebranch-anchored progressive genes.\nDo not force progressive molecular fate divergence." )
}
grDevices::pdf( OUT_FIG3, width = 13, height = 14, useDingbats = FALSE )
print(p3)
grDevices::dev.off()
if ( length(smooth_genes) >= 2 && nrow(plot_bin_tbl) > 0 ) {
  smooth_plot <- plot_bin_tbl |> dplyr::filter( gene_id %in% smooth_genes ) |> dplyr::mutate( x = ifelse( state == "Shared_prebranch", -0.05, progress_mid ), direction = ifelse( gene_id %in% arm1_top, "Arm_1-biased", "Arm_2-biased" ) )
  shared <- smooth_plot |> dplyr::filter( state == "Shared_prebranch" )
  line_tbl <- dplyr::bind_rows( shared |> dplyr::mutate(lineage = "Arm_1"), shared |> dplyr::mutate(lineage = "Arm_2"), smooth_plot |> dplyr::filter( state %in% c("Arm_1", "Arm_2") ) |> dplyr::mutate(lineage = state) )
  p4 <- ggplot2::ggplot( line_tbl, ggplot2::aes( x = x, y = mean_expression, group = lineage, linetype = lineage ) ) + ggplot2::geom_line(linewidth = 0.8) + ggplot2::geom_point(size = 1.3) + ggplot2::geom_vline( xintercept = 0, linetype = 2, linewidth = 0.4 ) + ggplot2::facet_wrap( direction ~ gene_id, scales = "free_y", ncol = 3 ) + ggplot2::labs( title = "Top genes from one shared prebranch state into two Y_29 arms", subtitle = "Shared prebranch point duplicated for visualization only", x = "Progress away from Y_29", y = "Mean normalized expression", linetype = "Lineage" ) + ggplot2::theme_classic(base_size = 9)
} else {
  p4 <- null_plot( "No progressive-divergence genes available for trajectory plotting." )
}
grDevices::pdf( OUT_FIG4, width = 13, height = 11, useDingbats = FALSE )
print(p4)
grDevices::dev.off()
result_object <- list( branchpoint = PRIMARY_BRANCHPOINT, prebranch_graph_radius = PREBRANCH_GRAPH_RADIUS, prebranch_cells = prebranch_cells, arm1_cells = arm1_cells, arm2_cells = arm2_cells, formal_step5Cv2_genes = formal_genes, progressive_shape_genes = progressive_tbl$gene_id, strong_progressive_genes = strong_tbl$gene_id, master_table = master, analysis_metadata = analysis_meta, formulas = list( full = FULL_FORMULA, additive = ADDITIVE_FORMULA, null = NULL_FORMULA ) )
saveRDS( result_object, OUT_RDS, compress = TRUE )
save( fit_full_CR, fit_add_CR, fit_null_CR, shape_test, branch_test, file = file.path( OUTDIR, "Step5CR_model_cache.RData" ), compress = TRUE )
parameter_tbl <- tibble::tibble( parameter = c( "Primary branchpoint", "Prebranch graph radius", "Shared prebranch cells", "Arm_1 downstream cells", "Arm_2 downstream cells", "Step5C-v2 formal genes tested", "Model family", "Progressive-shape q threshold", "Progressive-shape genes", "Anchored branch-deviation genes", "Strong progressive genes", "Minimum progressive-shape p", "Minimum progressive-shape q", "R cores" ), value = c( PRIMARY_BRANCHPOINT, PREBRANCH_GRAPH_RADIUS, length(prebranch_cells), length(arm1_cells), length(arm2_cells), length(formal_genes), "negbinomial", Q_MAX, nrow(progressive_tbl), sum(master$anchored_branch_significant, na.rm = TRUE), nrow(strong_tbl), min(master$progressive_shape_p_value, na.rm = TRUE), min(master$progressive_shape_q_value, na.rm = TRUE), CORES ) )
wb <- openxlsx::createWorkbook()
safe_sheet(wb, "Parameters", parameter_tbl)
safe_sheet(wb, "Analysis_cells", analysis_meta)
safe_sheet(wb, "Cell_summary", analysis_summary)
safe_sheet(wb, "Model_fit_status", fit_status)
safe_sheet(wb, "Progressive_shape_LRT", shape_test)
safe_sheet(wb, "Anchored_branch_LRT", branch_test)
safe_sheet(wb, "Master_568_genes", master)
safe_sheet(wb, "Progressive_genes", progressive_tbl)
safe_sheet(wb, "Strong_progressive", strong_tbl)
safe_sheet(wb, "7dpi_replicate_effects", replicate_effects)
safe_sheet(wb, "Binned_expression", plot_bin_tbl)
openxlsx::saveWorkbook( wb, OUT_XLSX, overwrite = TRUE )
top_progressive <- master |> dplyr::filter(progressive_shape_significant) |> dplyr::select( gene_id, progressive_shape_q_value, anchored_branch_q_value, early_arm1_minus_arm2, late_arm1_minus_arm2, divergence_gain, mean_7dpi_arm1_minus_arm2, direction_consistent, strong_progressive ) |> dplyr::slice_head(n = 30)
report_lines <- c( "======================================================================", "GSE162610 | STEP 5C-R REPORT", "Y_29 shared prebranch -> Arm_1/Arm_2 progressive divergence audit", "======================================================================", "", "FORMAL POSITION:", "This is still Step 5. GO/KEGG/GSEA remain Step 6.", "", "PRIOR RESULT:", paste0( "Step5C-v2 formal branch-dependent genes carried forward: ", length(formal_genes) ), "", "SHARED PREBRANCH:", paste0( "Definition = 3 dpi + shared_transition_3dpi_highconf + ", "Pre_branch_root_side/Branchpoint_vertex + graph distance <= ", PREBRANCH_GRAPH_RADIUS, " from Y_29" ), paste0( "Shared prebranch cells: ", length(prebranch_cells) ), paste0( "Arm_1 downstream cells: ", length(arm1_cells) ), paste0( "Arm_2 downstream cells: ", length(arm2_cells) ), "", "CELL SUPPORT:", capture_table(analysis_summary), "", "MODEL:", paste0("FULL: ", FULL_FORMULA), paste0("ADDITIVE: ", ADDITIVE_FORMULA), paste0("NULL: ", NULL_FORMULA), "", "PRIMARY TEST:", "FULL vs ADDITIVE asks whether Arm1-vs-Arm2 separation changes as cells", "progress away from the shared Y_29 prebranch state.", "", "SECONDARY TEST:", "FULL vs NULL asks whether any branch-specific deviation exists after Y_29.", "", "RESULT:", paste0( "Progressive-shape genes q<=0.05: ", nrow(progressive_tbl) ), paste0( "Anchored branch-deviation genes q<=0.05: ", sum(master$anchored_branch_significant, na.rm = TRUE) ), paste0( "Strong progressive genes: ", nrow(strong_tbl) ), paste0( "Minimum progressive-shape p: ", signif( min(master$progressive_shape_p_value, na.rm = TRUE), 6 ) ), paste0( "Minimum progressive-shape q: ", signif( min(master$progressive_shape_q_value, na.rm = TRUE), 6 ) ), "", "TOP PROGRESSIVE GENES:", capture_table(top_progressive), "", "DECISION RULE:", "Strong progressive fate-divergence support requires a non-trivial set of", "significant progressive-shape genes, increasing Arm1/Arm2 separation from", "early to late downstream progression, and consistent late direction in both", "7 dpi biological samples.", "", "If this fails, retain the valid conclusion of Y_29-associated branch-biased", "transcriptional states, but do not claim progressive molecular fate divergence.", "", "STATUS:", "If Step5C-R passes, proceed to Step5D Monocle2 DDRTree sensitivity.", "", "OUTPUTS:", OUT_FIG1, OUT_FIG2, OUT_FIG3, OUT_FIG4, OUT_XLSX, OUT_RDS, OUT_REP, OUT_SESS )
writeLines( report_lines, OUT_REP, useBytes = TRUE )
writeLines( capture.output(sessionInfo()), OUT_SESS, useBytes = TRUE )
msg("====================================================")
msg("STEP 5C-R COMPLETE")
msg("====================================================")
msg("Shared prebranch cells: ", length(prebranch_cells))
msg("Step5C-v2 genes tested: ", length(formal_genes))
msg("Progressive-shape genes q<=0.05: ", nrow(progressive_tbl))
msg( "Anchored branch-deviation genes q<=0.05: ", sum(master$anchored_branch_significant, na.rm = TRUE) )
msg("Strong progressive genes: ", nrow(strong_tbl))
msg( "Minimum progressive-shape q: ", signif( min(master$progressive_shape_q_value, na.rm = TRUE), 6 ) )
msg("")
msg("Key outputs:")
msg("1) Fig5CR_01_shared_prebranch_definition.pdf")
msg("2) Fig5CR_02_progressive_divergence_statistics.pdf")
msg("3) Fig5CR_03_progressive_divergence_heatmap.pdf")
msg("4) Fig5CR_04_top_gene_branchpoint_trajectories.pdf")
msg("5) Table5CR_Y29_prebranch_progressive_divergence.xlsx")
msg("6) Step5CR_report.txt")
