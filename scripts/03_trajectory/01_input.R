options(stringsAsFactors = FALSE)
set.seed(20260731)
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
INPUT_RDS <- file.path( PROJECT_ROOT, "04_results", "step4A_GSE162610_repair_scar_candidate_program_discovery", "GSE162610_step4A_candidate_program_scored.rds" )
OUTDIR <- file.path( PROJECT_ROOT, "04_results", "step5A_GSE162610_balanced_root_trajectory_input" )
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
OUT_FIG1 <- file.path(OUTDIR, "Fig5A_01_uninjured_root_candidate_audit.pdf")
OUT_FIG2 <- file.path(OUTDIR, "Fig5A_02_balanced_root_selection.pdf")
OUT_FIG3 <- file.path(OUTDIR, "Fig5A_03_trajectory_input_overview.pdf")
OUT_XLSX <- file.path(OUTDIR, "Table5A_root_and_trajectory_input.xlsx")
OUT_RDS  <- file.path(OUTDIR, "GSE162610_step5A_balanced_root_trajectory_input.rds")
OUT_REP  <- file.path(OUTDIR, "Step5A_report.txt")
OUT_SESS <- file.path(OUTDIR, "Step5A_sessionInfo.txt")
TIME_COL       <- "timepoint"
SAMPLE_COL     <- "sample_id"
CLUSTER_COL    <- "step3A_cluster_r04"
UMAP_REDUCTION <- "umap.microglia.step3A"
ROOT_CANDIDATE_CLUSTERS <- c("0", "3")
TRAJECTORY_CLUSTERS <- c("0", "1", "2", "3", "4", "5", "10")
CYCLING_CLUSTERS <- c("6", "8")
IFN_CLUSTER <- "9"
ROOT_TOP_PERCENTILE <- 0.70
ROOT_TARGET_PER_SAMPLE <- 300
ROOT_MIN_PER_SAMPLE <- 80
HOMEOSTATIC_GENES <- c( "P2ry12", "P2ry13", "Tmem119", "Sall1", "Fcrls", "Hexb", "Cx3cr1", "Gpr34", "Siglech", "Olfml3", "Slc2a5", "Adgrg1" )
ACTIVATED_DAM_GENES <- c( "Apoe", "Lpl", "Spp1", "Gpnmb", "Lgals3", "Fabp5", "Cst7", "Itgax", "Axl", "Clec7a", "Ccl2", "Ccl3", "Ccl4", "Il1b", "Tnf", "Nfkbia" )
required_pkgs <- c( "Seurat", "SeuratObject", "Matrix", "dplyr", "tidyr", "tibble", "ggplot2", "patchwork", "openxlsx" )
missing_pkgs <- required_pkgs[ !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE) ]
if (length(missing_pkgs) > 0) {
  stop(paste0("Missing R packages: ", paste(missing_pkgs, collapse = ", ")), call. = FALSE)
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
})
TIME_COLORS <- c( "Uninjured" = "#66C2A5", "1 dpi"     = "#FC8D62", "3 dpi"     = "#8DA0CB", "7 dpi"     = "#E78AC3" )
CLUSTER_COLORS <- c( "0"  = "#4E79A7", "1"  = "#E15759", "2"  = "#59A14F", "3"  = "#F28E2B", "4"  = "#B07AA1", "5"  = "#76B7B2", "10" = "#2F4B7C" )
ROOT_STATUS_COLORS <- c( "Not_selected" = "#D9D9D9", "Balanced_root" = "#D73027" )
ROOT_SCORE_COLORS <- c( "#313695", "#74ADD1", "#FFFFBF", "#FDAE61", "#A50026" )
msg <- function(...) {
  cat( sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "" )
}
capture_table <- function(x) {
  x <- tibble::as_tibble(x)
  paste( capture.output( print( x, n = Inf, width = Inf ) ), collapse = "\n" )
}
get_data_matrix <- function(obj) {
  assay <- DefaultAssay(obj)
  x <- tryCatch( Seurat::GetAssayData( obj, assay = assay, layer = "data" ), error = function(e) NULL )
  if (!is.null(x)) return(x)
  Seurat::GetAssayData( obj, assay = assay, slot = "data" )
}
calc_gene_z_program <- function(mat, genes) {
  genes <- intersect( genes, rownames(mat) )
  if (length(genes) == 0) {
    stop("No genes available for requested program.", call. = FALSE)
  }
  small <- as.matrix( mat[genes, , drop = FALSE] )
  gene_mean <- rowMeans( small, na.rm = TRUE )
  gene_sd <- apply( small, 1, stats::sd, na.rm = TRUE )
  gene_sd[ !is.finite(gene_sd) | gene_sd == 0 ] <- 1
  zmat <- sweep( small, 1, gene_mean, "-" )
  zmat <- sweep( zmat, 1, gene_sd, "/" )
  score <- colMeans( zmat, na.rm = TRUE )
  list( score = as.numeric(score), genes = genes )
}
msg("====================================================")
msg("GSE162610 Step 5A | balanced root + trajectory input")
msg("====================================================")
if (!file.exists(INPUT_RDS)) {
  stop( paste0( "找不到Step 4A scored RDS：\n", INPUT_RDS, "\n\nCheck required" ), call. = FALSE )
}
msg("Reading: ", INPUT_RDS)
obj <- readRDS(INPUT_RDS)
if (!inherits(obj, "Seurat")) {
  stop("Input RDS is not a Seurat object.", call. = FALSE)
}
DefaultAssay(obj) <- "RNA"
msg( "Loaded object: ", ncol(obj), " cells; ", nrow(obj), " features." )
required_meta <- c( TIME_COL, SAMPLE_COL, CLUSTER_COL, "transition_status", "shared_transition_3dpi_highconf" )
missing_meta <- required_meta[ !required_meta %in% colnames(obj[[]]) ]
if (length(missing_meta) > 0) {
  stop( paste0( "缺少正式metadata：", paste(missing_meta, collapse = ", ") ), call. = FALSE )
}
if (!UMAP_REDUCTION %in% Seurat::Reductions(obj)) {
  stop( paste0( "缺少正式Step3A UMAP：", UMAP_REDUCTION, "\n本脚本不会重新计算或替换UMAP。" ), call. = FALSE )
}
obj$step3A_cluster_r04 <- as.character( obj$step3A_cluster_r04 )
time_vec <- as.character( obj[[TIME_COL, drop = TRUE]] )
sample_vec <- as.character( obj[[SAMPLE_COL, drop = TRUE]] )
cluster_vec <- as.character( obj[[CLUSTER_COL, drop = TRUE]] )
trajectory_keep <- ( cluster_vec %in% TRAJECTORY_CLUSTERS )
trajectory_cells <- colnames(obj)[ trajectory_keep ]
if (length(trajectory_cells) < 1000) {
  stop( "Trajectory input contains <1000 cells; check Step3 cluster metadata.", call. = FALSE )
}
msg( "Trajectory input cells (clusters 0/1/2/3/4/5/10): ", length(trajectory_cells) )
trajectory_time_summary <- tibble::tibble( cell = trajectory_cells, timepoint = time_vec[ match(trajectory_cells, colnames(obj)) ], sample = sample_vec[ match(trajectory_cells, colnames(obj)) ], cluster = cluster_vec[ match(trajectory_cells, colnames(obj)) ] ) |> dplyr::count( timepoint, sample, cluster, name = "n_cells" )
root_candidate_keep <- ( time_vec == "Uninjured" & cluster_vec %in% ROOT_CANDIDATE_CLUSTERS & trajectory_keep )
root_candidate_cells <- colnames(obj)[ root_candidate_keep ]
if (length(root_candidate_cells) == 0) {
  stop( "No Uninjured cluster 0/3 root candidates detected.", call. = FALSE )
}
root_candidate_samples <- sort( unique( sample_vec[ match( root_candidate_cells, colnames(obj) ) ] ) )
msg( "Uninjured root candidate samples: ", paste(root_candidate_samples, collapse = ", ") )
if (length(root_candidate_samples) < 3) {
  warning( paste0( "Detected only ", length(root_candidate_samples), " uninjured samples. Expected 3." ) )
}
root_candidate_count <- tibble::tibble( cell = root_candidate_cells, sample = sample_vec[ match(root_candidate_cells, colnames(obj)) ], cluster = cluster_vec[ match(root_candidate_cells, colnames(obj)) ] ) |> dplyr::count( sample, cluster, name = "n_cells" )
msg("Root candidate cells: ", length(root_candidate_cells))
msg("Calculating biological root score...")
expr <- get_data_matrix(obj)
if (!identical( colnames(expr), colnames(obj) )) {
  expr <- expr[ , colnames(obj), drop = FALSE ]
}
expr_root <- expr[ , root_candidate_cells, drop = FALSE ]
homeo_avail <- intersect( HOMEOSTATIC_GENES, rownames(expr_root) )
dam_avail <- intersect( ACTIVATED_DAM_GENES, rownames(expr_root) )
msg( "Homeostatic genes available: ", length(homeo_avail), "/", length(HOMEOSTATIC_GENES) )
msg( "Activated/DAM genes available: ", length(dam_avail), "/", length(ACTIVATED_DAM_GENES) )
if (length(homeo_avail) < 6) {
  stop( paste0( "可用homeostatic genes <6：", paste(homeo_avail, collapse = ", ") ), call. = FALSE )
}
if (length(dam_avail) < 6) {
  stop( paste0( "可用activated/DAM genes <6：", paste(dam_avail, collapse = ", ") ), call. = FALSE )
}
homeo_score_obj <- calc_gene_z_program( expr_root, homeo_avail )
dam_score_obj <- calc_gene_z_program( expr_root, dam_avail )
root_score_tbl <- tibble::tibble( cell = root_candidate_cells, sample = sample_vec[ match(root_candidate_cells, colnames(obj)) ], cluster = cluster_vec[ match(root_candidate_cells, colnames(obj)) ], homeostatic_program_score = homeo_score_obj$score, activated_DAM_penalty_score = dam_score_obj$score ) |> dplyr::mutate( root_score = homeostatic_program_score - activated_DAM_penalty_score ) |> dplyr::group_by(sample) |> dplyr::mutate( root_score_percentile = dplyr::percent_rank(root_score) ) |> dplyr::ungroup()
highconf_root_candidates <- root_score_tbl |> dplyr::filter( root_score_percentile >= ROOT_TOP_PERCENTILE )
highconf_counts <- highconf_root_candidates |> dplyr::count( sample, name = "n_highconf" )
if (nrow(highconf_counts) != length(root_candidate_samples)) {
  stop( "At least one uninjured sample has no high-confidence root candidates.", call. = FALSE )
}
root_n_per_sample <- min( ROOT_TARGET_PER_SAMPLE, min(highconf_counts$n_highconf) )
msg( "Balanced root cells selected per uninjured sample: ", root_n_per_sample )
if (root_n_per_sample < ROOT_MIN_PER_SAMPLE) {
  stop(paste0("Balanced root has ", root_n_per_sample, " cells/sample (<", ROOT_MIN_PER_SAMPLE, ")."), call. = FALSE)
}
selected_root_tbl <- highconf_root_candidates |> dplyr::group_by(sample) |> dplyr::arrange( dplyr::desc(root_score), .by_group = TRUE ) |> dplyr::mutate( within_sample_rank = dplyr::row_number() ) |> dplyr::filter( within_sample_rank <= root_n_per_sample ) |> dplyr::ungroup()
selected_root_cells <- selected_root_tbl$cell
selected_root_summary <- selected_root_tbl |> dplyr::count( sample, cluster, name = "n_root" )
sample_balance_check <- selected_root_tbl |> dplyr::count( sample, name = "n_root" )
if (length(unique(sample_balance_check$n_root)) != 1) {
  stop( "Internal error: selected root is not exactly sample-balanced.", call. = FALSE )
}
msg( "Final balanced root cells: ", length(selected_root_cells) )
obj$step5A_use_for_trajectory <- trajectory_keep
obj$step5A_root_candidate <- FALSE
obj$step5A_root_candidate[ match(root_candidate_cells, colnames(obj)) ] <- TRUE
obj$step5A_balanced_root <- FALSE
obj$step5A_balanced_root[ match(selected_root_cells, colnames(obj)) ] <- TRUE
obj$step5A_homeostatic_program_score <- NA_real_
obj$step5A_activated_DAM_penalty_score <- NA_real_
obj$step5A_root_score <- NA_real_
obj$step5A_root_score_percentile <- NA_real_
candidate_idx <- match( root_score_tbl$cell, colnames(obj) )
obj$step5A_homeostatic_program_score[candidate_idx] <- root_score_tbl$homeostatic_program_score
obj$step5A_activated_DAM_penalty_score[candidate_idx] <- root_score_tbl$activated_DAM_penalty_score
obj$step5A_root_score[candidate_idx] <- root_score_tbl$root_score
obj$step5A_root_score_percentile[candidate_idx] <- root_score_tbl$root_score_percentile
prior_state <- rep( "Other_trajectory_cell", ncol(obj) )
prior_state[ time_vec == "1 dpi" & cluster_vec %in% c("2", "5", "10") ] <- "Early_activated_reference"
prior_state[ as.logical(obj$shared_transition_3dpi_highconf) ] <- "Shared_transition_reference"
prior_state[ time_vec == "7 dpi" & cluster_vec == "1" & as.character(obj$transition_status) == "Mature_DAM_like" ] <- "Established_DAM_reference"
prior_state[ obj$step5A_balanced_root ] <- "Balanced_homeostatic_root"
obj$step5A_prior_reference_state <- prior_state
obj_traj <- subset( obj, cells = trajectory_cells )
root_score_summary_by_sample <- root_score_tbl |> dplyr::group_by(sample) |> dplyr::summarise( n_candidates = dplyr::n(), root_score_mean = mean(root_score, na.rm = TRUE), root_score_median = median(root_score, na.rm = TRUE), root_score_q75 = stats::quantile( root_score, 0.75, na.rm = TRUE ), .groups = "drop" )
selected_root_score_summary <- selected_root_tbl |> dplyr::group_by(sample) |> dplyr::summarise( n_selected = dplyr::n(), homeostatic_mean = mean( homeostatic_program_score, na.rm = TRUE ), activated_DAM_mean = mean( activated_DAM_penalty_score, na.rm = TRUE ), root_score_mean = mean(root_score, na.rm = TRUE), root_score_median = median(root_score, na.rm = TRUE), .groups = "drop" )
trajectory_summary_time <- obj_traj[[]] |> dplyr::count( .data[[TIME_COL]], name = "n_cells" )
colnames(trajectory_summary_time)[ colnames(trajectory_summary_time) == TIME_COL ] <- "timepoint"
trajectory_summary_sample <- obj_traj[[]] |> dplyr::count( .data[[TIME_COL]], .data[[SAMPLE_COL]], name = "n_cells" )
colnames(trajectory_summary_sample)[ colnames(trajectory_summary_sample) == TIME_COL ] <- "timepoint"
colnames(trajectory_summary_sample)[ colnames(trajectory_summary_sample) == SAMPLE_COL ] <- "sample"
trajectory_summary_cluster <- obj_traj[[]] |> dplyr::count( .data[[TIME_COL]], .data[[CLUSTER_COL]], name = "n_cells" )
colnames(trajectory_summary_cluster)[ colnames(trajectory_summary_cluster) == TIME_COL ] <- "timepoint"
colnames(trajectory_summary_cluster)[ colnames(trajectory_summary_cluster) == CLUSTER_COL ] <- "cluster"
umap <- Seurat::Embeddings( obj, reduction = UMAP_REDUCTION )
umap_df <- tibble::tibble( cell = rownames(umap), UMAP_1 = umap[, 1], UMAP_2 = umap[, 2], timepoint = time_vec[ match(rownames(umap), colnames(obj)) ], sample = sample_vec[ match(rownames(umap), colnames(obj)) ], cluster = cluster_vec[ match(rownames(umap), colnames(obj)) ], trajectory_input = rownames(umap) %in% trajectory_cells, root_candidate = rownames(umap) %in% root_candidate_cells, balanced_root = rownames(umap) %in% selected_root_cells )
root_umap_df <- umap_df |> dplyr::filter( root_candidate ) |> dplyr::left_join( root_score_tbl |> dplyr::select( cell, homeostatic_program_score, activated_DAM_penalty_score, root_score, root_score_percentile ), by = "cell" )
traj_umap_df <- umap_df |> dplyr::filter( trajectory_input )
all_samples <- sort( unique( traj_umap_df$sample ) )
sample_palette <- grDevices::hcl.colors( length(all_samples), palette = "Dark 3" )
names(sample_palette) <- all_samples
uninj_samples <- sort( unique( root_umap_df$sample ) )
uninj_palette <- grDevices::hcl.colors( length(uninj_samples), palette = "Dark 3" )
names(uninj_palette) <- uninj_samples
msg("Creating Fig5A_01_uninjured_root_candidate_audit.pdf...")
theme5 <- ggplot2::theme_classic(base_size = 11) + ggplot2::theme( axis.text = ggplot2::element_text( color = "#222222" ), axis.title = ggplot2::element_text( color = "#222222" ), plot.title = ggplot2::element_text( face = "bold", color = "#111111" ), plot.subtitle = ggplot2::element_text( color = "#333333" ) )
p1a <- ggplot2::ggplot( root_umap_df, ggplot2::aes( UMAP_1, UMAP_2, color = factor(cluster) ) ) + ggplot2::geom_point( size = 0.55, alpha = 0.75 ) + ggplot2::scale_color_manual( values = CLUSTER_COLORS, name = "Step3A cluster" ) + ggplot2::labs( title = "Uninjured homeostatic candidate pool", subtitle = "Only formal clusters 0 and 3" ) + theme5
p1b <- ggplot2::ggplot( root_umap_df, ggplot2::aes( UMAP_1, UMAP_2, color = root_score ) ) + ggplot2::geom_point( size = 0.55, alpha = 0.80 ) + ggplot2::scale_color_gradientn( colors = ROOT_SCORE_COLORS, name = "Root score" ) + ggplot2::labs( title = "Biological root score", subtitle = "Homeostatic program minus activated/DAM penalty" ) + theme5
p1c <- ggplot2::ggplot( root_umap_df, ggplot2::aes( UMAP_1, UMAP_2, color = factor(sample) ) ) + ggplot2::geom_point( size = 0.55, alpha = 0.75 ) + ggplot2::scale_color_manual( values = uninj_palette, name = "Uninjured sample" ) + ggplot2::labs( title = "Candidate pool by biological sample" ) + theme5
score_violin_df <- root_score_tbl |> dplyr::mutate( selected = ifelse( cell %in% selected_root_cells, "Balanced_root", "Not_selected" ) )
p1d <- ggplot2::ggplot( score_violin_df, ggplot2::aes( x = factor(sample), y = root_score, fill = factor(sample) ) ) + ggplot2::geom_violin( trim = TRUE, scale = "width", alpha = 0.8, linewidth = 0.25 ) + ggplot2::geom_boxplot( width = 0.12, outlier.shape = NA, fill = "#FFFFFF", linewidth = 0.25 ) + ggplot2::scale_fill_manual( values = uninj_palette, guide = "none" ) + ggplot2::labs( title = "Root-score distribution by uninjured sample", x = "Biological sample", y = "Root score" ) + theme5 + ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 25, hjust = 1 ) )
grDevices::pdf( OUT_FIG1, width = 13, height = 10, useDingbats = FALSE )
print( (p1a | p1b) / (p1c | p1d) + patchwork::plot_annotation( title = "Step 5A | Homeostatic-root candidate audit" ) )
grDevices::dev.off()
msg("Creating Fig5A_02_balanced_root_selection.pdf...")
root_overlay_df <- traj_umap_df |> dplyr::mutate( root_status = ifelse( balanced_root, "Balanced_root", "Not_selected" ) )
p2a <- ggplot2::ggplot( root_overlay_df, ggplot2::aes( UMAP_1, UMAP_2 ) ) + ggplot2::geom_point( data = root_overlay_df |> dplyr::filter( root_status == "Not_selected" ), color = ROOT_STATUS_COLORS[["Not_selected"]], size = 0.35, alpha = 0.30 ) + ggplot2::geom_point( data = root_overlay_df |> dplyr::filter( root_status == "Balanced_root" ), color = ROOT_STATUS_COLORS[["Balanced_root"]], size = 0.80, alpha = 0.90 ) + ggplot2::labs( title = "Balanced homeostatic root on formal Step3A UMAP", subtitle = "Red cells will orient pseudotime in Step 5B", x = "UMAP_1", y = "UMAP_2" ) + theme5
selected_root_plot_df <- root_umap_df |> dplyr::filter( cell %in% selected_root_cells )
p2b <- ggplot2::ggplot( selected_root_plot_df, ggplot2::aes( UMAP_1, UMAP_2, color = factor(sample) ) ) + ggplot2::geom_point( size = 0.85, alpha = 0.80 ) + ggplot2::scale_color_manual( values = uninj_palette, name = "Root sample" ) + ggplot2::labs( title = "Selected root cells by biological sample", subtitle = paste0( root_n_per_sample, " cells selected per sample" ) ) + theme5
selected_score_long <- selected_root_tbl |> dplyr::select( cell, sample, homeostatic_program_score, activated_DAM_penalty_score, root_score ) |> tidyr::pivot_longer( cols = c( homeostatic_program_score, activated_DAM_penalty_score, root_score ), names_to = "metric", values_to = "score" ) |> dplyr::mutate( metric = dplyr::recode( metric, homeostatic_program_score = "Homeostatic program", activated_DAM_penalty_score = "Activated/DAM penalty", root_score = "Final root score" ) )
p2c <- ggplot2::ggplot( selected_score_long, ggplot2::aes( x = factor(sample), y = score, fill = factor(sample) ) ) + ggplot2::geom_boxplot( width = 0.55, outlier.shape = NA, alpha = 0.80 ) + ggplot2::facet_wrap( ~metric, nrow = 1, scales = "free_y" ) + ggplot2::scale_fill_manual( values = uninj_palette, guide = "none" ) + ggplot2::labs( title = "Selected-root biological profile", x = "Uninjured sample", y = "Program score" ) + theme5 + ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 25, hjust = 1 ) )
root_composition_plot <- selected_root_tbl |> dplyr::count( sample, cluster, name = "n_cells" ) |> dplyr::group_by(sample) |> dplyr::mutate( fraction = n_cells / sum(n_cells) ) |> dplyr::ungroup()
p2d <- ggplot2::ggplot( root_composition_plot, ggplot2::aes( x = factor(sample), y = fraction, fill = factor(cluster) ) ) + ggplot2::geom_col( width = 0.70 ) + ggplot2::scale_fill_manual( values = CLUSTER_COLORS, name = "Step3A cluster" ) + ggplot2::labs( title = "Balanced-root cluster composition", x = "Uninjured sample", y = "Fraction of selected root" ) + theme5 + ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 25, hjust = 1 ) )
grDevices::pdf( OUT_FIG2, width = 14, height = 10, useDingbats = FALSE )
print( (p2a | p2b) / (p2c | p2d) + patchwork::plot_annotation( title = "Step 5A | Sample-balanced homeostatic-root selection" ) )
grDevices::dev.off()
msg("Creating Fig5A_03_trajectory_input_overview.pdf...")
p3a <- ggplot2::ggplot( traj_umap_df, ggplot2::aes( UMAP_1, UMAP_2, color = factor(timepoint) ) ) + ggplot2::geom_point( size = 0.40, alpha = 0.60 ) + ggplot2::scale_color_manual( values = TIME_COLORS, name = "Timepoint" ) + ggplot2::labs( title = "Formal trajectory input by true time" ) + theme5
p3b <- ggplot2::ggplot( traj_umap_df, ggplot2::aes( UMAP_1, UMAP_2, color = factor(cluster) ) ) + ggplot2::geom_point( size = 0.40, alpha = 0.60 ) + ggplot2::scale_color_manual( values = CLUSTER_COLORS, name = "Step3A cluster" ) + ggplot2::labs( title = "Formal trajectory input by Step3A state" ) + theme5
p3c <- ggplot2::ggplot( root_overlay_df, ggplot2::aes( UMAP_1, UMAP_2 ) ) + ggplot2::geom_point( color = "#D9D9D9", size = 0.35, alpha = 0.30 ) + ggplot2::geom_point( data = root_overlay_df |> dplyr::filter( balanced_root ), color = "#D73027", size = 0.85, alpha = 0.95 ) + ggplot2::labs( title = "Root orientation anchor", subtitle = "Root is balanced; graph will be learned independently in Step 5B" ) + theme5
composition_df <- obj_traj[[]] |> dplyr::count( .data[[TIME_COL]], .data[[SAMPLE_COL]], name = "n_cells" )
colnames(composition_df)[ colnames(composition_df) == TIME_COL ] <- "timepoint"
colnames(composition_df)[ colnames(composition_df) == SAMPLE_COL ] <- "sample"
p3d <- ggplot2::ggplot( composition_df, ggplot2::aes( x = factor(timepoint), y = n_cells, fill = factor(sample) ) ) + ggplot2::geom_col( position = "stack", width = 0.70 ) + ggplot2::scale_fill_manual( values = sample_palette, name = "Biological sample" ) + ggplot2::labs( title = "Trajectory-input composition by true time and sample", x = "True timepoint", y = "Cells" ) + theme5
grDevices::pdf( OUT_FIG3, width = 14, height = 10, useDingbats = FALSE )
print( (p3a | p3b) / (p3c | p3d) + patchwork::plot_annotation( title = "Step 5A | Formal trajectory-input construction" ) )
grDevices::dev.off()
msg("Writing Table5A_root_and_trajectory_input.xlsx...")
program_gene_audit <- tibble::tibble( program = c( rep( "Homeostatic", length(HOMEOSTATIC_GENES) ), rep( "Activated_DAM_penalty", length(ACTIVATED_DAM_GENES) ) ), gene = c( HOMEOSTATIC_GENES, ACTIVATED_DAM_GENES ) ) |> dplyr::mutate( available_in_RNA = gene %in% rownames(expr) )
root_parameters <- tibble::tibble( parameter = c( "Root candidate timepoint", "Root candidate clusters", "Root top percentile threshold", "Target root cells per sample", "Minimum allowed root cells per sample", "Actual balanced root cells per sample", "Trajectory clusters", "Excluded cycling clusters", "Excluded IFN cluster", "Formal UMAP" ), value = c( "Uninjured", paste(ROOT_CANDIDATE_CLUSTERS, collapse = ";"), as.character(ROOT_TOP_PERCENTILE), as.character(ROOT_TARGET_PER_SAMPLE), as.character(ROOT_MIN_PER_SAMPLE), as.character(root_n_per_sample), paste(TRAJECTORY_CLUSTERS, collapse = ";"), paste(CYCLING_CLUSTERS, collapse = ";"), IFN_CLUSTER, UMAP_REDUCTION ) )
selected_root_export <- selected_root_tbl |> dplyr::arrange( sample, within_sample_rank )
wb <- openxlsx::createWorkbook()
add_sheet <- function(sheet_name, df) {
  openxlsx::addWorksheet( wb, sheetName = sheet_name )
  openxlsx::writeDataTable( wb, sheet = sheet_name, x = df, tableStyle = "TableStyleMedium2" )
  openxlsx::freezePane( wb, sheet = sheet_name, firstRow = TRUE )
}
add_sheet( "Root_parameters", root_parameters )
add_sheet( "Program_gene_audit", program_gene_audit )
add_sheet( "Root_candidate_counts", root_candidate_count )
add_sheet( "Highconf_counts", highconf_counts )
add_sheet( "Root_score_by_sample", root_score_summary_by_sample )
add_sheet( "Selected_root_summary", selected_root_summary )
add_sheet( "Selected_root_scores", selected_root_score_summary )
add_sheet( "Selected_root_cells", selected_root_export )
add_sheet( "Trajectory_by_time", trajectory_summary_time )
add_sheet( "Trajectory_by_sample", trajectory_summary_sample )
add_sheet( "Trajectory_by_cluster", trajectory_summary_cluster )
add_sheet( "Trajectory_full_audit", trajectory_time_summary )
openxlsx::saveWorkbook( wb, OUT_XLSX, overwrite = TRUE )
saveRDS( obj_traj, OUT_RDS, compress = TRUE )
msg("Writing Step5A_report.txt...")
old_cluster3_root <- tibble::tibble( sample = sample_vec[ time_vec == "Uninjured" & cluster_vec == "3" ] ) |> dplyr::count( sample, name = "old_cluster3_root_cells" )
report_lines <- c( "======================================================================", "GSE162610 | STEP 5A REPORT", "Sample-balanced homeostatic root + formal trajectory input", "======================================================================", "", "SCOPE:", "This step does NOT run Monocle.", "This step only constructs a biologically justified balanced root and the", "formal non-cycling trajectory input for Step 5B.", "", "------------------------------- INPUT ---------------------------------", paste0("Input RDS: ", INPUT_RDS), paste0("Formal cluster field: ", CLUSTER_COL), paste0("Formal UMAP: ", UMAP_REDUCTION), paste0("Input cells: ", ncol(obj)), "", "----------------------- OLD ROOT IMBALANCE AUDIT ----------------------", "Uninjured + cluster 3 counts in the current object:", capture_table(old_cluster3_root), "", "The old cluster-3-only definition is NOT used as the final root.", "", "------------------------ NEW ROOT DEFINITION --------------------------", paste0( "Candidate pool: Uninjured clusters ", paste(ROOT_CANDIDATE_CLUSTERS, collapse = "/") ), paste0( "Candidate cells: ", length(root_candidate_cells) ), paste0( "Homeostatic genes available: ", length(homeo_avail), "/", length(HOMEOSTATIC_GENES) ), paste0( "Activated/DAM penalty genes available: ", length(dam_avail), "/", length(ACTIVATED_DAM_GENES) ), paste0( "Within-sample high-confidence percentile threshold: ", ROOT_TOP_PERCENTILE ), paste0( "Balanced selected root cells/sample: ", root_n_per_sample ), paste0( "Total balanced root cells: ", length(selected_root_cells) ), "", "Selected root counts by sample/cluster:", capture_table(selected_root_summary), "", "Selected root biological-score summary:", capture_table(selected_root_score_summary), "", "----------------------- TRAJECTORY INPUT ------------------------------", paste0( "Included Step3A clusters: ", paste(TRAJECTORY_CLUSTERS, collapse = ", ") ), paste0( "Excluded cycling clusters: ", paste(CYCLING_CLUSTERS, collapse = ", ") ), paste0( "Excluded IFN side-state cluster: ", IFN_CLUSTER ), paste0( "Formal trajectory input cells: ", ncol(obj_traj) ), "", "Trajectory cells by true time:", capture_table(trajectory_summary_time), "", "---------------------- ANTI-CIRCULARITY RULE --------------------------", "Repair/scar scores from Step 4 are NOT used to:", "  - select root cells;", "  - select trajectory cells;", "  - construct a graph;", "  - assign any future branch.", "", "They are carried forward only so that AFTER an independent Step 5B graph", "is learned, the Step 4 functional programs can be projected onto branches.", "", "-------------------------- STEP 5A DECISION ---------------------------", "PASS criteria:", "  1. all expected uninjured biological samples contribute root candidates;", "  2. exactly the same number of high-confidence root cells is selected", "     from every uninjured sample;", "  3. selected root cells show high homeostatic and low activated/DAM scores;", "  4. trajectory input spans Uninjured / 1 dpi / 3 dpi / 7 dpi;", "  5. cycling and IFN side states are excluded from the primary trajectory.", "", "Only after Step 5A PASS should Step 5B Monocle3 be run.", "", "Formal outputs:", paste0("  ", OUT_FIG1), paste0("  ", OUT_FIG2), paste0("  ", OUT_FIG3), paste0("  ", OUT_XLSX), paste0("  ", OUT_RDS), paste0("  ", OUT_REP), paste0("  ", OUT_SESS), "" )
writeLines( report_lines, OUT_REP, useBytes = TRUE )
writeLines( capture.output( sessionInfo() ), OUT_SESS, useBytes = TRUE )
msg("====================================================")
msg("STEP 5A COMPLETE")
msg("====================================================")
msg("Output directory: ", OUTDIR)
msg("")
msg("Balanced root cells/sample: ", root_n_per_sample)
msg("Total balanced root cells: ", length(selected_root_cells))
msg("Trajectory input cells: ", ncol(obj_traj))
msg("")
msg("Key outputs:")
msg("  1) Fig5A_01_uninjured_root_candidate_audit.pdf")
msg("  2) Fig5A_02_balanced_root_selection.pdf")
msg("  3) Fig5A_03_trajectory_input_overview.pdf")
msg("  4) Table5A_root_and_trajectory_input.xlsx")
msg("  5) Step5A_report.txt")
msg("")
msg("Trajectory input status: complete")
