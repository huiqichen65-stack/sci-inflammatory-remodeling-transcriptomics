options(stringsAsFactors = FALSE)
set.seed(20260801)
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
STEP6A_DIR <- file.path( PROJECT_ROOT, "04_results", "step6A_GSE162610_7dpi_branch_specific_DEG" )
INPUT_STEP6A <- file.path( STEP6A_DIR, "GSE162610_step6A_7dpi_branch_DEG_results.rds" )
OUTDIR <- file.path( PROJECT_ROOT, "04_results", "step6E_GSE162610_pseudobulk_ssGSEA" )
FIGDIR <- file.path(OUTDIR, "figures")
TABDIR <- file.path(OUTDIR, "tables")
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
OUT_FIG1 <- file.path( FIGDIR, "Fig6E_01_pseudobulk_ssGSEA_heatmap.pdf" )
OUT_FIG2 <- file.path( FIGDIR, "Fig6E_02_paired_branch_scores.pdf" )
OUT_FIG3 <- file.path( FIGDIR, "Fig6E_03_within_sample_branch_effects.pdf" )
OUT_FIG4 <- file.path( FIGDIR, "Fig6E_04_program_direction_summary.pdf" )
OUT_XLSX <- file.path( OUTDIR, "Table6E_pseudobulk_ssGSEA.xlsx" )
OUT_RDS <- file.path( OUTDIR, "GSE162610_step6E_pseudobulk_ssGSEA_results.rds" )
OUT_REP <- file.path( OUTDIR, "Step6E_report.txt" )
OUT_SESS <- file.path( OUTDIR, "Step6E_sessionInfo.txt" )
MIN_DETECTED_GENES <- 5L
MAX_GENESET_SIZE <- 500L
ARM_COLORS <- c( "Arm_1" = "#2F6B9A", "Arm_2" = "#B64B45" )
SAMPLE_SHAPES <- c( "7dpi_sample1" = 16, "7dpi_sample2" = 17 )
HEAT_COLORS <- grDevices::colorRampPalette( c( "#2166AC", "#F7F7F7", "#B2182B" ) )(101)
required_pkgs <- c( "GSVA", "BiocParallel", "org.Mm.eg.db", "AnnotationDbi", "dplyr", "tidyr", "tibble", "ggplot2", "patchwork", "pheatmap", "openxlsx" )
missing_pkgs <- required_pkgs[ !vapply( required_pkgs, requireNamespace, logical(1), quietly = TRUE ) ]
if (length(missing_pkgs) > 0) {
  stop(paste0("Missing R packages: ", paste(missing_pkgs, collapse = ", ")), call. = FALSE)
}
suppressPackageStartupMessages({
  library(GSVA)
  library(BiocParallel)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})
capture_table <- function(x) {
  paste( capture.output( print( tibble::as_tibble(x), n = Inf, width = Inf ) ), collapse = "\n" )
}
safe_sheet <- function(wb, sheet_name, df) {
  openxlsx::addWorksheet( wb, sheetName = sheet_name )
  ddf <- tryCatch( as.data.frame(df), error = function(e) NULL )
  if ( is.null(ddf) || nrow(ddf) == 0 || ncol(ddf) == 0 ) {
    openxlsx::writeData( wb, sheet = sheet_name, x = data.frame( note = "No result under the predefined Step6E criterion.", stringsAsFactors = FALSE ) )
  } else {
    for (cc in colnames(ddf)) {
      if (is.list(ddf[[cc]])) {
        ddf[[cc]] <- vapply( ddf[[cc]], function(z) {
            paste( as.character(z), collapse = "; " )
          },
          character(1)
        )
      }
    }
    openxlsx::writeDataTable( wb, sheet = sheet_name, x = ddf, tableStyle = "TableStyleMedium2" )
  }
  openxlsx::freezePane( wb, sheet = sheet_name, firstRow = TRUE )
}
row_zscore <- function(mat) {
  mat <- as.matrix(mat)
  row_mean <- rowMeans( mat, na.rm = TRUE )
  row_sd <- apply( mat, 1, stats::sd, na.rm = TRUE )
  row_sd[ !is.finite(row_sd) | row_sd == 0 ] <- 1
  out <- sweep( mat, 1, row_mean, "-" )
  out <- sweep( out, 1, row_sd, "/" )
  out
}
run_ssgsea_compatible <- function( expression_matrix, gene_sets) {
  exports <- getNamespaceExports( "GSVA" )
  if ("ssgseaParam" %in% exports) {
    param <- GSVA::ssgseaParam( exprData = as.matrix(expression_matrix), geneSets = gene_sets, minSize = MIN_DETECTED_GENES, maxSize = MAX_GENESET_SIZE, alpha = 0.25, normalize = TRUE )
    serial_param <- BiocParallel::SerialParam( progressbar = FALSE )
    ans <- tryCatch( GSVA::gsva( param, verbose = FALSE, BPPARAM = serial_param ), error = function(e) {
        GSVA::gsva( param, verbose = FALSE )
      }
    )
  } else {
    ans <- GSVA::gsva( expr = as.matrix(expression_matrix), gset.idx.list = gene_sets, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE, min.sz = MIN_DETECTED_GENES, max.sz = MAX_GENESET_SIZE, ssgsea.norm = TRUE, parallel.sz = 1, verbose = FALSE )
  }
  as.matrix(ans)
}
get_go_gene_sets <- function( specifications, expressed_genes) {
  go_ids <- unique( specifications$GO_ID )
  mapping <- AnnotationDbi::select( org.Mm.eg.db, keys = go_ids, columns = "SYMBOL", keytype = "GOALL" )
  mapping <- as.data.frame( mapping )
  if ( !"GOALL" %in% colnames(mapping) || !"SYMBOL" %in% colnames(mapping) ) {
    stop( paste0( "GOALL mapping returned unexpected columns: ", paste( colnames(mapping), collapse = ", " ) ), call. = FALSE )
  }
  mapping <- mapping |> dplyr::filter( !is.na(GOALL), !is.na(SYMBOL), nzchar(SYMBOL) ) |> dplyr::distinct( GOALL, SYMBOL ) |> dplyr::inner_join( specifications, by = c( "GOALL" = "GO_ID" ) )
  gene_sets <- split( mapping$SYMBOL, mapping$program )
  gene_sets <- lapply( gene_sets, function(x) {
      intersect( unique( as.character(x) ), expressed_genes )
    }
  )
  list( gene_sets = gene_sets, mapping = mapping )
}
if (!file.exists(INPUT_STEP6A)) {
  stop( paste0( "Required input not found: ", INPUT_STEP6A ), call. = FALSE )
}
step6a <- readRDS( INPUT_STEP6A )
required_fields <- c( "pseudobulk_logCPM", "pseudobulk_metadata" )
missing_fields <- required_fields[ !required_fields %in% names(step6a) ]
if (length(missing_fields) > 0) {
  stop( paste0( "Required input not found: ", paste( missing_fields, collapse = ", " ) ), call. = FALSE )
}
pb_logcpm <- as.matrix( step6a$pseudobulk_logCPM )
pb_meta <- as.data.frame( step6a$pseudobulk_metadata )
storage.mode(pb_logcpm) <- "double"
if ( is.null(rownames(pb_logcpm)) || is.null(colnames(pb_logcpm)) ) {
  stop( "Required input not found: ", call. = FALSE )
}
required_meta_cols <- c( "pb_group", "sample", "branch" )
missing_meta_cols <- required_meta_cols[ !required_meta_cols %in% colnames(pb_meta) ]
if (length(missing_meta_cols) > 0) {
  stop( paste0( "Missing metadata: ", paste( missing_meta_cols, collapse = ", " ) ), call. = FALSE )
}
pb_meta$pb_group <- as.character( pb_meta$pb_group )
pb_meta$sample <- as.character( pb_meta$sample )
pb_meta$branch <- as.character( pb_meta$branch )
if ( !all( colnames(pb_logcpm) %in% pb_meta$pb_group ) ) {
  stop( "Validation check.", call. = FALSE )
}
pb_meta <- pb_meta[ match( colnames(pb_logcpm), pb_meta$pb_group ), , drop = FALSE ]
if (anyNA(pb_meta$pb_group)) {
  stop( "Validation check.", call. = FALSE )
}
pb_meta$sample <- factor( pb_meta$sample, levels = sort( unique( pb_meta$sample ) ) )
pb_meta$branch <- factor( pb_meta$branch, levels = c( "Arm_2", "Arm_1" ) )
pair_audit <- pb_meta |> dplyr::count( sample, branch, name = "n_pseudobulk" )
if ( any( pair_audit$n_pseudobulk != 1 ) ) {
  stop( "Requirement not met.", call. = FALSE )
}
sample_branch_counts <- pair_audit |> dplyr::group_by( sample ) |> dplyr::summarise( n_branches = dplyr::n_distinct( branch ), .groups = "drop" )
if ( any( sample_branch_counts$n_branches != 2 ) ) {
  stop( "Required result not found.", call. = FALSE )
}
pb_meta <- pb_meta |> dplyr::arrange( sample, branch )
pb_logcpm <- pb_logcpm[ , pb_meta$pb_group, drop = FALSE ]
print( pb_meta[ , required_meta_cols, drop = FALSE ] )
go_specs <- tibble::tribble( ~GO_ID,      ~program,                         ~category,                 ~expected_direction, "GO:0006119", "Oxidative_phosphorylation",      "Energy_metabolism",       "Arm_1", "GO:0042775", "Mitochondrial_ATP_synthesis",    "Energy_metabolism",       "Arm_1", "GO:0006909", "Phagocytosis",                   "Clearance_lysosome",       "Arm_1", "GO:0006911", "Phagocytosis_engulfment",        "Clearance_lysosome",       "Arm_1", "GO:0007040", "Lysosome_organization",          "Clearance_lysosome",       "Arm_1", "GO:0016042", "Lipid_catabolic_process",        "Lipid_handling",           "Arm_1", "GO:0033344", "Cholesterol_efflux",              "Lipid_handling",           "Arm_1", "GO:0042060", "Wound_healing",                  "Tissue_support",           "Arm_1", "GO:0001525", "Angiogenesis",                   "Tissue_support",           "Arm_1", "GO:0006954", "Inflammatory_response",          "Inflammation",             "Arm_2", "GO:0001816", "Cytokine_production",            "Inflammation",             "Arm_2", "GO:0033209", "TNF_mediated_signaling",         "Inflammation",             "Arm_2", "GO:0051092", "NFkB_activation",                "Inflammation",             "Arm_2", "GO:0032496", "Response_to_LPS",                "Innate_immune_response",    "Arm_2", "GO:0006956", "Complement_activation",          "Innate_immune_response",    "Arm_2", "GO:0097190", "Apoptotic_signaling",            "Stress_cell_death",        "Arm_2", "GO:0006979", "Response_to_oxidative_stress",   "Stress_cell_death",        "Arm_2", "GO:0030198", "Extracellular_matrix_organization", "Matrix_TGFb",           "Arm_2", "GO:0007179", "TGFb_receptor_signaling",        "Matrix_TGFb",              "Arm_2" )
go_data <- get_go_gene_sets( specifications = go_specs, expressed_genes = rownames(pb_logcpm) )
go_gene_sets <- go_data$gene_sets
state_gene_sets <- list( Homeostatic_microglia = c( "P2ry12", "Tmem119", "Sall1", "Fcrls", "Hexb", "Cx3cr1", "Gpr34", "Olfml3", "Slc2a5", "Siglech" ), Activated_DAM_like = c( "Apoe", "Trem2", "Tyrobp", "Lpl", "Cst7", "Clec7a", "Axl", "Itgax", "Lgals3", "Ctsd", "Ctsb", "Fabp5", "Spp1", "Cd9" ) )
state_gene_sets <- lapply( state_gene_sets, function(x) {
    intersect( unique(x), rownames(pb_logcpm) )
  }
)
gene_sets_all <- c( state_gene_sets, go_gene_sets )
program_metadata <- dplyr::bind_rows( tibble::tibble( program = names(state_gene_sets), category = "State_context", expected_direction = c( "Arm_1", "Neutral" ), source = "Pre-specified curated state markers" ), go_specs |> dplyr::select( program, category, expected_direction ) |> dplyr::mutate( source = "Mouse GOALL annotation" ) ) |> dplyr::distinct( program, .keep_all = TRUE )
gene_set_audit <- tibble::tibble( program = names(gene_sets_all), detected_genes = vapply( gene_sets_all, length, integer(1) ), detected_gene_list = vapply( gene_sets_all, function(x) {
      paste( sort( unique(x) ), collapse = "; " )
    },
    character(1)
  )
) |> dplyr::left_join(
    program_metadata, by = "program"
  ) |> dplyr::mutate(
    retained_for_ssGSEA = detected_genes >= MIN_DETECTED_GENES & detected_genes <= MAX_GENESET_SIZE
  ) |> dplyr::arrange(
    category, program
  )
print( gene_set_audit |> dplyr::select( program, category, detected_genes, retained_for_ssGSEA ) )
gene_sets_use <- gene_sets_all[ gene_set_audit$program[ gene_set_audit$retained_for_ssGSEA ] ]
if (length(gene_sets_use) < 8) {
  stop( paste0( "Only ", length(gene_sets_use), " gene sets passed the detected-size criterion; continue." ), call. = FALSE )
}
ssgsea_scores <- run_ssgsea_compatible( expression_matrix = pb_logcpm, gene_sets = gene_sets_use )
if ( !all( colnames(pb_logcpm) %in% colnames(ssgsea_scores) ) ) {
  stop( "ssGSEA score columns match the input pseudobulk profiles.", call. = FALSE )
}
ssgsea_scores <- ssgsea_scores[ , colnames(pb_logcpm), drop = FALSE ]
score_long <- as.data.frame( t(ssgsea_scores) ) |> tibble::rownames_to_column( "pb_group" ) |> tidyr::pivot_longer( cols = -pb_group, names_to = "program", values_to = "ssGSEA_score" ) |> dplyr::left_join( pb_meta |> dplyr::select( pb_group, sample, branch ), by = "pb_group" ) |> dplyr::left_join( program_metadata, by = "program" ) |> dplyr::group_by( program ) |> dplyr::mutate( centered_score = ssGSEA_score - mean( ssGSEA_score, na.rm = TRUE ) ) |> dplyr::ungroup()
score_wide <- score_long |> dplyr::select( program, category, expected_direction, sample, branch, ssGSEA_score ) |> tidyr::pivot_wider( names_from = branch, values_from = ssGSEA_score )
required_score_cols <- c( "Arm_1", "Arm_2" )
if ( !all( required_score_cols %in% colnames(score_wide) ) ) {
  stop( "Cannot construct paired Arm_1 and Arm_2 ssGSEA scores.", call. = FALSE )
}
branch_effects <- score_wide |> dplyr::mutate( Arm_1_minus_Arm_2 = Arm_1 - Arm_2, observed_direction = dplyr::case_when( Arm_1_minus_Arm_2 > 0 ~ "Arm_1", Arm_1_minus_Arm_2 < 0 ~ "Arm_2", TRUE ~ "Neutral" ) )
program_summary <- branch_effects |> dplyr::group_by( program, category, expected_direction ) |> dplyr::summarise( n_samples = dplyr::n(), mean_Arm_1_minus_Arm_2 = mean( Arm_1_minus_Arm_2, na.rm = TRUE ), min_effect = min( Arm_1_minus_Arm_2, na.rm = TRUE ), max_effect = max( Arm_1_minus_Arm_2, na.rm = TRUE ), min_abs_effect = min( abs( Arm_1_minus_Arm_2 ), na.rm = TRUE ), same_direction_all_samples = all( Arm_1_minus_Arm_2 > 0 ) || all( Arm_1_minus_Arm_2 < 0 ), observed_direction = dplyr::case_when( mean_Arm_1_minus_Arm_2 > 0 ~ "Arm_1", mean_Arm_1_minus_Arm_2 < 0 ~ "Arm_2", TRUE ~ "Neutral" ), .groups = "drop" ) |> dplyr::mutate( matches_expected_direction = expected_direction == "Neutral" | ( same_direction_all_samples & observed_direction == expected_direction ) ) |> dplyr::arrange( category, dplyr::desc( mean_Arm_1_minus_Arm_2 ) )
sample_effect_wide <- branch_effects |> dplyr::select( program, sample, Arm_1_minus_Arm_2 ) |> tidyr::pivot_wider( names_from = sample, values_from = Arm_1_minus_Arm_2, names_prefix = "effect_" )
program_summary <- program_summary |> dplyr::left_join( sample_effect_wide, by = "program" )
print( program_summary |> dplyr::select( program, category, expected_direction, mean_Arm_1_minus_Arm_2, same_direction_all_samples, observed_direction, matches_expected_direction ) )
score_z <- row_zscore( ssgsea_scores )
heat_order <- program_summary |> dplyr::arrange( category, dplyr::desc( mean_Arm_1_minus_Arm_2 ) ) |> dplyr::pull( program )
heat_order <- heat_order[ heat_order %in% rownames(score_z) ]
score_z <- score_z[ heat_order, pb_meta$pb_group, drop = FALSE ]
annotation_col <- pb_meta |> dplyr::select( sample, branch ) |> as.data.frame()
rownames(annotation_col) <- pb_meta$pb_group
annotation_row <- program_summary |> dplyr::select( program, category, expected_direction ) |> dplyr::distinct() |> as.data.frame()
rownames(annotation_row) <- annotation_row$program
annotation_row$program <- NULL
annotation_row <- annotation_row[ rownames(score_z), , drop = FALSE ]
annotation_colors <- list( branch = ARM_COLORS, sample = c( "7dpi_sample1" = "#6A3D9A", "7dpi_sample2" = "#1B9E77" ), expected_direction = c( "Arm_1" = "#2F6B9A", "Arm_2" = "#B64B45", "Neutral" = "#7F7F7F" ) )
pheatmap::pheatmap( mat = score_z, color = HEAT_COLORS, breaks = seq( -2.5, 2.5, length.out = length(HEAT_COLORS) + 1 ), cluster_rows = FALSE, cluster_cols = FALSE, annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors, show_colnames = TRUE, show_rownames = TRUE, fontsize = 9, fontsize_row = 8, fontsize_col = 8, angle_col = 45, border_color = NA, main = paste0( "7 dpi pseudobulk ssGSEA | row z-score\n", "Gene sets are independent of Arm_1/Arm_2 DEG selection" ), filename = OUT_FIG1, width = 12, height = 11 )
score_long$branch <- factor( score_long$branch, levels = c( "Arm_2", "Arm_1" ) )
p2 <- ggplot2::ggplot( score_long, ggplot2::aes( x = branch, y = centered_score, group = sample ) ) + ggplot2::geom_hline( yintercept = 0, linewidth = 0.3, color = "grey75" ) + ggplot2::geom_line( color = "grey55", linewidth = 0.6 ) + ggplot2::geom_point( ggplot2::aes( color = branch, shape = sample ), size = 2.5 ) + ggplot2::facet_wrap( ~program, scales = "free_y", ncol = 4 ) + ggplot2::scale_color_manual( values = ARM_COLORS ) + ggplot2::scale_shape_manual( values = SAMPLE_SHAPES ) + ggplot2::labs( title = "Paired 7 dpi pseudobulk ssGSEA scores", subtitle = "Each line connects Arm_2 and Arm_1 within one biological sample", x = NULL, y = "Program-centered ssGSEA score", color = NULL, shape = "Biological sample" ) + ggplot2::theme_classic( base_size = 9 ) + ggplot2::theme( legend.position = "top", strip.background = ggplot2::element_rect( fill = "grey95", color = NA ), strip.text = ggplot2::element_text( face = "bold", size = 7.5 ), axis.text.x = ggplot2::element_text( angle = 25, hjust = 1 ) )
ggplot2::ggsave( filename = OUT_FIG2, plot = p2, width = 14, height = 12, units = "in" )
effect_order <- program_summary |> dplyr::arrange( mean_Arm_1_minus_Arm_2 ) |> dplyr::pull( program )
branch_effects$program <- factor( branch_effects$program, levels = effect_order )
p3 <- ggplot2::ggplot( branch_effects, ggplot2::aes( x = Arm_1_minus_Arm_2, y = program, shape = sample, color = observed_direction ) ) + ggplot2::geom_vline( xintercept = 0, linetype = 2, linewidth = 0.45 ) + ggplot2::geom_point( size = 2.8, alpha = 0.9 ) + ggplot2::facet_grid( category ~ ., scales = "free_y", space = "free_y" ) + ggplot2::scale_color_manual( values = c( "Arm_1" = ARM_COLORS[["Arm_1"]], "Arm_2" = ARM_COLORS[["Arm_2"]], "Neutral" = "#7F7F7F" ) ) + ggplot2::scale_shape_manual( values = SAMPLE_SHAPES ) + ggplot2::labs( title = "Within-sample branch effects", subtitle = "Positive = Arm_1 higher; negative = Arm_2 higher", x = "Arm_1 - Arm_2 ssGSEA score", y = NULL, color = "Observed direction", shape = "Biological sample" ) + ggplot2::theme_classic( base_size = 9 ) + ggplot2::theme( legend.position = "top", strip.background = ggplot2::element_rect( fill = "grey95", color = NA ), strip.text.y = ggplot2::element_text( angle = 0, face = "bold" ) )
ggplot2::ggsave( filename = OUT_FIG3, plot = p3, width = 12, height = 11, units = "in" )
summary_plot <- program_summary |> dplyr::mutate( program = factor( program, levels = effect_order ), direction_color = dplyr::case_when( observed_direction == "Arm_1" ~ "Arm_1", observed_direction == "Arm_2" ~ "Arm_2", TRUE ~ "Neutral" ) )
p4 <- ggplot2::ggplot( summary_plot, ggplot2::aes( x = mean_Arm_1_minus_Arm_2, y = program, color = direction_color ) ) + ggplot2::geom_vline( xintercept = 0, linetype = 2, linewidth = 0.45 ) + ggplot2::geom_errorbarh( ggplot2::aes( xmin = min_effect, xmax = max_effect ), height = 0, linewidth = 0.7 ) + ggplot2::geom_point( ggplot2::aes( shape = same_direction_all_samples, size = min_abs_effect ), alpha = 0.95 ) + ggplot2::facet_grid( category ~ ., scales = "free_y", space = "free_y" ) + ggplot2::scale_color_manual( values = c( "Arm_1" = ARM_COLORS[["Arm_1"]], "Arm_2" = ARM_COLORS[["Arm_2"]], "Neutral" = "#7F7F7F" ) ) + ggplot2::scale_shape_manual( values = c( "TRUE" = 16, "FALSE" = 1 ) ) + ggplot2::labs( title = "Pseudobulk program direction summary", subtitle = "Point = mean effect; horizontal range = two biological samples", x = "Mean Arm_1 - Arm_2 ssGSEA score", y = NULL, color = "Observed direction", shape = "Same direction\nin both samples", size = "Minimum absolute\nsample effect" ) + ggplot2::theme_classic( base_size = 9 ) + ggplot2::theme( legend.position = "top", strip.background = ggplot2::element_rect( fill = "grey95", color = NA ), strip.text.y = ggplot2::element_text( angle = 0, face = "bold" ) )
ggplot2::ggsave( filename = OUT_FIG4, plot = p4, width = 12, height = 11, units = "in" )
score_matrix_export <- as.data.frame( ssgsea_scores ) |> tibble::rownames_to_column( "program" )
utils::write.csv( gene_set_audit, file.path( TABDIR, "01_gene_set_audit.csv" ), row.names = FALSE )
utils::write.csv( score_long, file.path( TABDIR, "02_pseudobulk_ssGSEA_scores_long.csv" ), row.names = FALSE )
utils::write.csv( branch_effects, file.path( TABDIR, "03_within_sample_branch_effects.csv" ), row.names = FALSE )
utils::write.csv( program_summary, file.path( TABDIR, "04_program_direction_summary.csv" ), row.names = FALSE )
parameter_tbl <- tibble::tibble( parameter = c( "Species", "Input expression", "Profiles", "Biological samples", "Branches", "ssGSEA implementation", "Minimum detected genes per set", "Maximum gene-set size", "Positive branch effect", "Formal p-values", "Branch-DEG-derived gene sets", "GO annotation key type" ), value = c( "Mus musculus", "Step6A edgeR TMM-normalized logCPM", ncol(pb_logcpm), paste( levels(pb_meta$sample), collapse = ", " ), "Arm_1 and Arm_2", if ( "ssgseaParam" %in% getNamespaceExports("GSVA") ) {
      "GSVA::ssgseaParam"
    } else {
      "GSVA legacy direct interface"
    },
    MIN_DETECTED_GENES, MAX_GENESET_SIZE, "Arm_1 higher than Arm_2", "Not calculated: only two biological samples", "None", "GOALL"
  )
)
wb <- openxlsx::createWorkbook()
safe_sheet( wb, "Parameters", parameter_tbl )
safe_sheet( wb, "Pseudobulk_metadata", pb_meta )
safe_sheet( wb, "Gene_set_audit", gene_set_audit )
safe_sheet( wb, "GO_specs", go_specs )
safe_sheet( wb, "Score_matrix", score_matrix_export )
safe_sheet( wb, "Scores_long", score_long )
safe_sheet( wb, "Branch_effects", branch_effects )
safe_sheet( wb, "Program_summary", program_summary )
openxlsx::saveWorkbook( wb, OUT_XLSX, overwrite = TRUE )
result_object <- list( pseudobulk_logCPM = pb_logcpm, pseudobulk_metadata = pb_meta, GO_specifications = go_specs, gene_sets = gene_sets_use, gene_set_audit = gene_set_audit, ssGSEA_scores = ssgsea_scores, score_long = score_long, branch_effects = branch_effects, program_summary = program_summary, thresholds = list( minimum_detected_genes = MIN_DETECTED_GENES, maximum_gene_set_size = MAX_GENESET_SIZE ) )
saveRDS( result_object, OUT_RDS, compress = TRUE )
robust_programs <- program_summary |> dplyr::filter( same_direction_all_samples ) |> dplyr::arrange( dplyr::desc( abs( mean_Arm_1_minus_Arm_2 ) ) )
expected_supported <- program_summary |> dplyr::filter( expected_direction != "Neutral", same_direction_all_samples, observed_direction == expected_direction ) |> dplyr::arrange( category, dplyr::desc( abs( mean_Arm_1_minus_Arm_2 ) ) )
discordant_programs <- program_summary |> dplyr::filter( !same_direction_all_samples ) |> dplyr::arrange( category, program )
writeLines( capture.output( sessionInfo() ), OUT_SESS, useBytes = TRUE )
