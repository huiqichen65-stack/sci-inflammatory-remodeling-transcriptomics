options( stringsAsFactors = FALSE, warn = 1, Seurat.object.assay.version = "v5", future.globals.maxSize = 24 * 1024^3 )
set.seed( 20260803L )
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
MAIN_TIMEPOINTS <- c( "Uninjured", "7 dpi", "14 dpi", "1 month", "2 months" )
TIMEPOINT_LEVELS <- MAIN_TIMEPOINTS
FROZEN_CLUSTER_COLUMN <- "myhcon_res_0_6"
FROZEN_UMAP <- "umap_myeloid_harmony"
PRIMARY_CLUSTERS <- c( "0", "1", "2", "3", "7" )
EXPECTED_PRIMARY_INPUT_FEATURES <- 26918L
EXPECTED_PRIMARY_INPUT_NUCLEI <- 20179L
MIN_LIBRARIES_PER_MAIN_TIMEPOINT <- 2L
MIN_SIGNATURE_GENES_PRESENT <- 20L
MIN_SIGNATURE_MAPPING_FRACTION <- 0.70
AUC_MAX_RANK_FRACTION <- 0.05
SSGSEA_ALPHA <- 0.25
PSEUDOBULK_PRIOR_COUNT <- 2
ARM1_STEP6B_MODULES <- c( "M3_Shared_prebranch" )
ARM2_STEP6B_MODULES <- c( "M1_Arm2_terminal", "M2_Arm2_terminal" )
TECHNICAL_GENE_PATTERN <- "^(mt-|Rpl|Rps|Hb[ab])"
PRIMARY_OBJECT_PATH <- file.path( PROJECT_ROOT, "04_results", "step10G_GSE234774_frozen_myeloid_annotations_resident_objects", "RDS", "GSE234774_step10G_primary_noncycling_resident_microglia_BPCells.rds" )
STEP6A_DIR <- file.path( PROJECT_ROOT, "04_results", "step6A_GSE162610_7dpi_branch_specific_DEG" )
STEP6A_ARM1_FILE <- file.path( STEP6A_DIR, "tables", "03_Arm1_biased_high_confidence_genes.csv" )
STEP6A_ARM2_FILE <- file.path( STEP6A_DIR, "tables", "04_Arm2_biased_high_confidence_genes.csv" )
STEP6B_ASSIGNMENT_FILE <- file.path( PROJECT_ROOT, "04_results", "step6B_GSE162610_dynamic_branch_modules", "tables", "01_dynamic_module_gene_assignments.csv" )
OUTPUT_DIR <- file.path( PROJECT_ROOT, "04_results", "step11A_GSE234774_main_timepoint_Arm_AUCell_ssGSEA" )
FIGURES_DIR <- file.path( OUTPUT_DIR, "figures" )
TABLES_DIR <- file.path( OUTPUT_DIR, "tables" )
RDS_DIR <- file.path( OUTPUT_DIR, "RDS" )
for (directory in c( OUTPUT_DIR, FIGURES_DIR, TABLES_DIR, RDS_DIR )) {
  dir.create( directory, recursive = TRUE, showWarnings = FALSE )
}
RUN_LOG <- file.path( OUTPUT_DIR, "00_run_log.txt" )
if (file.exists( RUN_LOG )) {
  file.remove( RUN_LOG )
}
SCORED_OBJECT_PATH <- file.path( RDS_DIR, "GSE234774_step11A_main_timepoint_primary_resident_Arm_scored_BPCells.rds" )
RESULTS_RDS_PATH <- file.path( RDS_DIR, "GSE234774_step11A_Arm_AUCell_ssGSEA_results.rds" )
write_csv_safe <- function( x, filename ) {
  path <- file.path( TABLES_DIR, filename )
  utils::write.csv( x, path, row.names = FALSE, na = "" )
  invisible( path )
}
stop_with_report <- function(reason) {
  report_path <- file.path( OUTPUT_DIR, "10_STEP11A_FINAL_REPORT.txt" )
  stop( reason, call. = FALSE )
}
safe_package_version <- function( package_name ) {
  if (!requireNamespace( package_name, quietly = TRUE )) {
    return( NA_character_ )
  }
  as.character( utils::packageVersion( package_name ) )
}
require_columns <- function( data, columns, label ) {
  missing_columns <- setdiff( columns, colnames( data ) )
  if ( length( missing_columns ) > 0L ) {
    stop_with_report( paste0( label, " is missing required column(s): ", paste( missing_columns, collapse = ", " ), "." ) )
  }
}
detect_gene_column <- function( data, label ) {
  candidates <- c( "gene_id", "gene", "Gene", "symbol", "SYMBOL" )
  detected <- intersect( candidates, colnames( data ) )
  if ( length( detected ) == 0L ) {
    stop_with_report( paste0( "Could not identify a gene-symbol column in ", label, "." ) )
  }
  detected[[1L]]
}
clean_gene_vector <- function(x) {
  output <- unique( trimws( as.character( x ) ) )
  output[ !is.na( output ) & nzchar( output ) ]
}
resolve_features <- function( requested_features, available_features ) {
  available_upper <- toupper( available_features )
  first_index <- !duplicated( available_upper )
  lookup <- setNames( available_features[ first_index ], available_upper[ first_index ] )
  resolved <- unname( lookup[ toupper( requested_features ) ] )
  names( resolved ) <- requested_features
  resolved
}
safe_zscore <- function(x) {
  value_sd <- stats::sd( x, na.rm = TRUE )
  if ( !is.finite( value_sd ) || value_sd == 0 ) {
    return( rep( 0, length( x ) ) )
  }
  output <- ( x - mean( x, na.rm = TRUE ) ) / value_sd
  output[ !is.finite( output ) ] <- 0
  output
}
to_dgCMatrix <- function( matrix_object, label ) {
  if (inherits( matrix_object, "dgCMatrix" )) {
    return( matrix_object )
  }
  converted <- tryCatch( methods::as( matrix_object, "dgCMatrix" ), error = function(error_object) {
      error_object
    }
  )
  if (inherits( converted, "error" )) {
    intermediate <- tryCatch( methods::as( matrix_object, "CsparseMatrix" ), error = function(error_object) {
        error_object
      }
    )
    if (!inherits( intermediate, "error" )) {
      converted <- tryCatch( methods::as( intermediate, "dgCMatrix" ), error = function(error_object) {
          error_object
        }
      )
    }
  }
  if (inherits( converted, "error" )) {
    stop_with_report( paste0( "Could not convert ", label, " to dgCMatrix without dense conversion: ", conditionMessage( converted ) ) )
  }
  converted
}
run_ssgsea_version_safe <- function( expression_matrix, gene_sets ) {
  gsva_exports <- getNamespaceExports( "GSVA" )
  if ( "ssgseaParam" %in% gsva_exports ) {
    parameter_arguments <- list( exprData = expression_matrix, geneSets = gene_sets, minSize = MIN_SIGNATURE_GENES_PRESENT, maxSize = Inf, alpha = SSGSEA_ALPHA, normalize = TRUE, checkNA = "yes", use = "all.obs", ondisk = "no", verbose = TRUE )
    parameter_formals <- names( formals( GSVA::ssgseaParam ) )
    if ( !"..." %in% parameter_formals ) {
      parameter_arguments <- parameter_arguments[ names( parameter_arguments ) %in% parameter_formals ]
    }
    parameter_object <- do.call( GSVA::ssgseaParam, parameter_arguments )
    result <- GSVA::gsva( parameter_object, verbose = TRUE, BPPARAM = BiocParallel::SerialParam() )
  } else {
    legacy_arguments <- list( expr = expression_matrix, gset.idx.list = gene_sets, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE, ssgsea.norm = TRUE, min.sz = MIN_SIGNATURE_GENES_PRESENT, max.sz = Inf, parallel.sz = 1L, verbose = TRUE )
    legacy_formals <- names( formals( GSVA::gsva ) )
    if ( !"..." %in% legacy_formals ) {
      legacy_arguments <- legacy_arguments[ names( legacy_arguments ) %in% legacy_formals ]
    }
    result <- do.call( GSVA::gsva, legacy_arguments )
  }
  as.matrix( result )
}
plot_umap_score <- function( plot_data, score_column, title_text, legend_title ) {
  current_data <- plot_data
  current_data$score <- current_data[[score_column]]
  ggplot2::ggplot( current_data, ggplot2::aes( x = UMAP_1, y = UMAP_2 ) ) + ggplot2::geom_point( ggplot2::aes( colour = score ), size = 0.28, alpha = 0.72 ) + ggplot2::scale_colour_viridis_c() + ggplot2::coord_equal() + ggplot2::labs( title = title_text, x = "Myeloid Harmony UMAP 1", y = "Myeloid Harmony UMAP 2", colour = legend_title ) + ggplot2::theme_bw()
}
safe_cor_test <- function( x, y ) {
  complete <- stats::complete.cases( x, y )
  x_complete <- x[ complete ]
  y_complete <- y[ complete ]
  if ( length( x_complete ) < 4L || stats::sd( x_complete ) == 0 || stats::sd( y_complete ) == 0 ) {
    return( list( rho = NA_real_, p_value = NA_real_, n = length( x_complete ) ) )
  }
  test <- suppressWarnings( stats::cor.test( x_complete, y_complete, method = "spearman", exact = FALSE ) )
  list( rho = unname( test$estimate ), p_value = test$p.value, n = length( x_complete ) )
}
required_packages <- c( "Matrix", "data.table", "ggplot2", "patchwork", "scales", "Seurat", "SeuratObject", "AUCell", "GSVA", "edgeR", "BiocParallel" )
missing_packages <- required_packages[ !vapply( required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1) ) ]
if ( length( missing_packages ) > 0L ) {
  stop_with_report( paste0( "Missing required package(s): ", paste( missing_packages, collapse = ", " ), "." ) )
}
suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(Seurat)
  library(SeuratObject)
  library(AUCell)
  library(GSVA)
  library(edgeR)
  library(BiocParallel)
})
package_versions <- data.frame( package = c( "R", required_packages ), version = c( paste0( R.version$major, ".", R.version$minor ), vapply( required_packages, safe_package_version, FUN.VALUE = character(1) ) ), stringsAsFactors = FALSE )
write_csv_safe( package_versions, "00_package_versions.csv" )
input_files <- c( primary_resident_object = PRIMARY_OBJECT_PATH, Step6A_Arm1_highconf = STEP6A_ARM1_FILE, Step6A_Arm2_highconf = STEP6A_ARM2_FILE, Step6B_module_assignments = STEP6B_ASSIGNMENT_FILE )
input_file_audit <- data.frame( input_name = names( input_files ), path = unname( input_files ), exists = file.exists( input_files ), size_bytes = vapply( input_files, function(path) {
      if ( file.exists( path ) ) {
        file.info( path )$size
      } else {
        NA_real_
      }
    },
    numeric(1)
  ), stringsAsFactors = FALSE
)
write_csv_safe( input_file_audit, "01_input_file_audit.csv" )
if ( !all( input_file_audit$exists ) ) {
  stop_with_report( paste0( "Missing required input file(s): ", paste( input_file_audit$input_name[ !input_file_audit$exists ], collapse = ", " ), "." ) )
}
primary_object <- readRDS( PRIMARY_OBJECT_PATH )
if (!inherits( primary_object, "Seurat" )) {
  stop_with_report( "The Step 10G primary input is not a Seurat object." )
}
if ( nrow( primary_object ) != EXPECTED_PRIMARY_INPUT_FEATURES || ncol( primary_object ) != EXPECTED_PRIMARY_INPUT_NUCLEI ) {
  stop_with_report( paste0( "Unexpected Step 10G primary-object dimensions. Expected ", EXPECTED_PRIMARY_INPUT_FEATURES, " x ", EXPECTED_PRIMARY_INPUT_NUCLEI, "; detected ", nrow( primary_object ), " x ", ncol( primary_object ), "." ) )
}
required_metadata <- c( "biological_library", "timepoint", "myeloid_final_label", "arm_analysis_set", FROZEN_CLUSTER_COLUMN )
missing_metadata <- setdiff( required_metadata, colnames( primary_object@meta.data ) )
if ( length( missing_metadata ) > 0L ) {
  stop_with_report( paste0( "Primary object is missing metadata column(s): ", paste( missing_metadata, collapse = ", " ), "." ) )
}
if (!FROZEN_UMAP %in% names( primary_object@reductions ) ) {
  stop_with_report( paste0( "Frozen UMAP reduction is missing: ", FROZEN_UMAP, "." ) )
}
input_cluster_values <- sort( unique( as.character( primary_object@meta.data[ , FROZEN_CLUSTER_COLUMN ] ) ), method = "radix" )
if (!setequal( input_cluster_values, PRIMARY_CLUSTERS )) {
  stop_with_report( paste0( "Primary object clusters match 0,1,2,3,7. Detected: ", paste( input_cluster_values, collapse = ", " ), "." ) )
}
if ( any( as.character( primary_object$arm_analysis_set ) != "PRIMARY_NONCYCLING_RESIDENT" ) ) {
  stop_with_report( "The primary object contains cells outside PRIMARY_NONCYCLING_RESIDENT." )
}
primary_timepoint_values <- as.character( primary_object$timepoint )
main_cells <- colnames( primary_object )[ primary_timepoint_values %in% MAIN_TIMEPOINTS ]
if ( length( main_cells ) == 0L ) {
  stop_with_report( "No nuclei were retained for the five main validation timepoints." )
}
main_object <- subset( x = primary_object, cells = main_cells )
main_object$timepoint <- factor( as.character( main_object$timepoint ), levels = TIMEPOINT_LEVELS, ordered = TRUE )
if ( any( as.character( main_object$timepoint ) %in% c( "1 dpi", "4 dpi" ) ) ) {
  stop_with_report( "The main analysis still contains 1 dpi or 4 dpi cells." )
}
if ( "Cycling_microglia" %in% as.character( main_object$myeloid_final_label ) ) {
  stop_with_report( "Cycling microglia is present in the main analysis object." )
}
if ( !"IFN_responsive_microglia" %in% as.character( main_object$myeloid_final_label ) ) {
  stop_with_report( "IFN-responsive microglia cluster 7 is absent from the main analysis." )
}
main_metadata <- main_object@meta.data
library_manifest <- unique( data.table::data.table( biological_library = as.character( main_metadata$biological_library ), timepoint = as.character( main_metadata$timepoint ) ) )
library_mapping_audit <- library_manifest[ , .( n_timepoint_labels = uniqueN( timepoint ), timepoint_labels = paste( sort( unique( timepoint ) ), collapse = ";" ) ), by = biological_library ]
if ( any( library_mapping_audit$n_timepoint_labels != 1L ) ) {
  stop_with_report( "At least one biological_library maps to multiple main timepoints." )
}
library_cell_counts <- data.table::data.table( biological_library = as.character( main_metadata$biological_library ), timepoint = as.character( main_metadata$timepoint ) )[ , .( n_cells = .N ), by = .( biological_library, timepoint ) ]
library_cell_counts$timepoint <- factor( library_cell_counts$timepoint, levels = TIMEPOINT_LEVELS, ordered = TRUE )
data.table::setorder( library_cell_counts, timepoint, biological_library )
timepoint_replication <- library_cell_counts[ , .( n_libraries = uniqueN( biological_library ), n_cells = sum( n_cells ), minimum_library_cells = min( n_cells ), maximum_library_cells = max( n_cells ) ), by = timepoint ]
missing_main_timepoints <- setdiff( MAIN_TIMEPOINTS, as.character( timepoint_replication$timepoint ) )
if ( length( missing_main_timepoints ) > 0L ) {
  stop_with_report( paste0( "Missing main validation timepoint(s): ", paste( missing_main_timepoints, collapse = ", " ), "." ) )
}
low_replication_timepoints <- as.character( timepoint_replication$timepoint[ timepoint_replication$n_libraries < MIN_LIBRARIES_PER_MAIN_TIMEPOINT ] )
if ( length( low_replication_timepoints ) > 0L ) {
  stop_with_report( paste0( "Main validation timepoint(s) have fewer than ", MIN_LIBRARIES_PER_MAIN_TIMEPOINT, " biological libraries: ", paste( low_replication_timepoints, collapse = ", " ), "." ) )
}
main_object_audit <- rbind( data.frame( record_type = "MAIN_OBJECT", timepoint = "ALL", biological_library = "ALL", n_libraries = length( unique( library_cell_counts$biological_library ) ), n_cells = ncol( main_object ), minimum_library_cells = min( library_cell_counts$n_cells ), maximum_library_cells = max( library_cell_counts$n_cells ), stringsAsFactors = FALSE ), data.frame( record_type = "TIMEPOINT", timepoint = as.character( timepoint_replication$timepoint ), biological_library = "ALL", n_libraries = timepoint_replication$n_libraries, n_cells = timepoint_replication$n_cells, minimum_library_cells = timepoint_replication$minimum_library_cells, maximum_library_cells = timepoint_replication$maximum_library_cells, stringsAsFactors = FALSE ), data.frame( record_type = "LIBRARY", timepoint = as.character( library_cell_counts$timepoint ), biological_library = library_cell_counts$biological_library, n_libraries = 1L, n_cells = library_cell_counts$n_cells, minimum_library_cells = library_cell_counts$n_cells, maximum_library_cells = library_cell_counts$n_cells, stringsAsFactors = FALSE ) )
write_csv_safe( main_object_audit, "02_main_timepoint_object_and_library_audit.csv" )
step6a_arm1 <- data.table::fread( STEP6A_ARM1_FILE, data.table = FALSE, check.names = FALSE, showProgress = FALSE )
step6a_arm2 <- data.table::fread( STEP6A_ARM2_FILE, data.table = FALSE, check.names = FALSE, showProgress = FALSE )
step6b_assignments <- data.table::fread( STEP6B_ASSIGNMENT_FILE, data.table = FALSE, check.names = FALSE, showProgress = FALSE )
arm1_gene_column <- detect_gene_column( step6a_arm1, "Step6A Arm_1 high-confidence marker table" )
arm2_gene_column <- detect_gene_column( step6a_arm2, "Step6A Arm_2 high-confidence marker table" )
require_columns( step6b_assignments, c( "gene_id", "module_name" ), "Step6B dynamic module assignment table" )
step6a_arm1_genes <- clean_gene_vector( step6a_arm1[ , arm1_gene_column ] )
step6a_arm2_genes <- clean_gene_vector( step6a_arm2[ , arm2_gene_column ] )
step6b_assignments$gene_id <- trimws( as.character( step6b_assignments$gene_id ) )
step6b_assignments$module_name <- trimws( as.character( step6b_assignments$module_name ) )
step6b_assignments <- step6b_assignments[ !is.na( step6b_assignments$gene_id ) & nzchar( step6b_assignments$gene_id ) & !is.na( step6b_assignments$module_name ) & nzchar( step6b_assignments$module_name ), , drop = FALSE ]
observed_modules <- unique( step6b_assignments$module_name )
required_modules <- c( ARM1_STEP6B_MODULES, ARM2_STEP6B_MODULES )
missing_modules <- setdiff( required_modules, observed_modules )
if ( length( missing_modules ) > 0L ) {
  stop_with_report( paste0( "Step6B assignment table is missing frozen module(s): ", paste( missing_modules, collapse = ", " ), "." ) )
}
step6b_arm1_genes <- clean_gene_vector( step6b_assignments$gene_id[ step6b_assignments$module_name %in% ARM1_STEP6B_MODULES ] )
step6b_arm2_genes <- clean_gene_vector( step6b_assignments$gene_id[ step6b_assignments$module_name %in% ARM2_STEP6B_MODULES ] )
arm1_frozen_requested <- unique( c( step6a_arm1_genes, step6b_arm1_genes ) )
arm2_frozen_requested <- unique( c( step6a_arm2_genes, step6b_arm2_genes ) )
if ( length( arm1_frozen_requested ) < MIN_SIGNATURE_GENES_PRESENT || length( arm2_frozen_requested ) < MIN_SIGNATURE_GENES_PRESENT ) {
  stop_with_report( "A frozen Arm signature contains too few genes before external mapping." )
}
available_features <- rownames( main_object )
arm1_resolution <- resolve_features( requested_features = arm1_frozen_requested, available_features = available_features )
arm2_resolution <- resolve_features( requested_features = arm2_frozen_requested, available_features = available_features )
arm1_present <- unique( unname( arm1_resolution[ !is.na( arm1_resolution ) ] ) )
arm2_present <- unique( unname( arm2_resolution[ !is.na( arm2_resolution ) ] ) )
arm1_mapping_fraction <- length( arm1_present ) / length( arm1_frozen_requested )
arm2_mapping_fraction <- length( arm2_present ) / length( arm2_frozen_requested )
if ( length( arm1_present ) < MIN_SIGNATURE_GENES_PRESENT || length( arm2_present ) < MIN_SIGNATURE_GENES_PRESENT || arm1_mapping_fraction < MIN_SIGNATURE_MAPPING_FRACTION || arm2_mapping_fraction < MIN_SIGNATURE_MAPPING_FRACTION ) {
  stop_with_report( paste0( "Frozen-signature mapping failed. Arm_1 present=", length( arm1_present ), "/", length( arm1_frozen_requested ), " (", round( arm1_mapping_fraction, 3 ), "); Arm_2 present=", length( arm2_present ), "/", length( arm2_frozen_requested ), " (", round( arm2_mapping_fraction, 3 ), ")." ) )
}
overlap_requested <- intersect( arm1_frozen_requested, arm2_frozen_requested )
overlap_present <- intersect( arm1_present, arm2_present )
arm_gene_sets <- list( Arm_1_frozen = arm1_present, Arm_2_frozen = arm2_present )
all_signature_genes <- sort( unique( c( arm1_frozen_requested, arm2_frozen_requested ) ) )
signature_gene_audit <- data.frame( record_type = "GENE", gene = all_signature_genes, component_Step6A_Arm1 = all_signature_genes %in% step6a_arm1_genes, component_Step6B_M3 = all_signature_genes %in% step6b_arm1_genes, component_Step6A_Arm2 = all_signature_genes %in% step6a_arm2_genes, component_Step6B_M1_M2 = all_signature_genes %in% step6b_arm2_genes, in_Arm1_frozen = all_signature_genes %in% arm1_frozen_requested, in_Arm2_frozen = all_signature_genes %in% arm2_frozen_requested, overlap_between_frozen_arms = all_signature_genes %in% overlap_requested, resolved_external_feature = unname( resolve_features( requested_features = all_signature_genes, available_features = available_features ) ), present_in_GSE234774 = !is.na( resolve_features( requested_features = all_signature_genes, available_features = available_features ) ), technical_gene_pattern = grepl( TECHNICAL_GENE_PATTERN, all_signature_genes, ignore.case = TRUE ), input_gene_count = NA_integer_, present_gene_count = NA_integer_, mapping_fraction = NA_real_, stringsAsFactors = FALSE )
signature_summary_rows <- data.frame( record_type = "SUMMARY", gene = c( "Step6A_Arm1_highconf", "Step6B_M3_Shared_prebranch", "Arm_1_frozen_union", "Step6A_Arm2_highconf", "Step6B_M1_M2_Arm2_terminal", "Arm_2_frozen_union", "Arm1_Arm2_overlap" ), component_Step6A_Arm1 = NA, component_Step6B_M3 = NA, component_Step6A_Arm2 = NA, component_Step6B_M1_M2 = NA, in_Arm1_frozen = NA, in_Arm2_frozen = NA, overlap_between_frozen_arms = NA, resolved_external_feature = "", present_in_GSE234774 = NA, technical_gene_pattern = NA, input_gene_count = c( length( step6a_arm1_genes ), length( step6b_arm1_genes ), length( arm1_frozen_requested ), length( step6a_arm2_genes ), length( step6b_arm2_genes ), length( arm2_frozen_requested ), length( overlap_requested ) ), present_gene_count = c( sum( step6a_arm1_genes %in% available_features ), sum( step6b_arm1_genes %in% available_features ), length( arm1_present ), sum( step6a_arm2_genes %in% available_features ), sum( step6b_arm2_genes %in% available_features ), length( arm2_present ), length( overlap_present ) ), mapping_fraction = c( mean( step6a_arm1_genes %in% available_features ), mean( step6b_arm1_genes %in% available_features ), arm1_mapping_fraction, mean( step6a_arm2_genes %in% available_features ), mean( step6b_arm2_genes %in% available_features ), arm2_mapping_fraction, ifelse( length( overlap_requested ) == 0L, NA_real_, length( overlap_present ) / length( overlap_requested ) ) ), stringsAsFactors = FALSE )
signature_gene_audit <- rbind( signature_summary_rows, signature_gene_audit )
write_csv_safe( signature_gene_audit, "04_frozen_Arm_signature_gene_audit.csv" )
Seurat::DefaultAssay( main_object ) <- "RNA"
rna_layers <- SeuratObject::Layers( object = main_object, assay = "RNA", search = NA )
counts_layers <- rna_layers[ grepl( "^counts($|\\.)", rna_layers ) ]
data_layers <- rna_layers[ grepl( "^data($|\\.)", rna_layers ) ]
if ( length( counts_layers ) != 1L || length( data_layers ) != 1L ) {
  stop_with_report( paste0( "Expected one RNA counts layer and one RNA data layer; detected ", length( counts_layers ), " and ", length( data_layers ), "." ) )
}
counts_layer_name <- counts_layers[[1L]]
data_layer_name <- data_layers[[1L]]
counts_matrix_raw <- SeuratObject::LayerData( object = main_object, assay = "RNA", layer = counts_layer_name, fast = FALSE )
counts_matrix <- to_dgCMatrix( counts_matrix_raw, "main-timepoint RNA counts" )
rm( counts_matrix_raw )
invisible( gc() )
if ( nrow( counts_matrix ) != nrow( main_object ) || ncol( counts_matrix ) != ncol( main_object ) || !identical( rownames( counts_matrix ), rownames( main_object ) ) || !identical( colnames( counts_matrix ), colnames( main_object ) ) ) {
  stop_with_report( "Sparse counts matrix does not match the main Seurat object." )
}
auc_max_rank <- max( 1L, ceiling( AUC_MAX_RANK_FRACTION * nrow( counts_matrix ) ) )
set.seed( 20260803L )
cell_rankings <- AUCell::AUCell_buildRankings( exprMat = counts_matrix, featureType = "genes", plotStats = FALSE, splitByBlocks = TRUE, BPPARAM = BiocParallel::SerialParam(), keepZeroesAsNA = FALSE, verbose = TRUE )
cell_auc_object <- AUCell::AUCell_calcAUC( geneSets = arm_gene_sets, rankings = cell_rankings, nCores = 1L, normAUC = TRUE, aucMaxRank = auc_max_rank, verbose = TRUE )
cell_auc_matrix <- as.matrix( AUCell::getAUC( cell_auc_object ) )
required_auc_rows <- c( "Arm_1_frozen", "Arm_2_frozen" )
if ( !all( required_auc_rows %in% rownames( cell_auc_matrix ) ) || !identical( colnames( cell_auc_matrix ), colnames( main_object ) ) || any( !is.finite( cell_auc_matrix[ required_auc_rows, , drop = FALSE ] ) ) ) {
  stop_with_report( "AUCell output dimensions, cell order, or finite-value checks failed." )
}
main_object$Arm_1_AUCell <- as.numeric( cell_auc_matrix[ "Arm_1_frozen", colnames( main_object ) ] )
main_object$Arm_2_AUCell <- as.numeric( cell_auc_matrix[ "Arm_2_frozen", colnames( main_object ) ] )
main_object$Arm_balance_AUCell <- main_object$Arm_1_AUCell - main_object$Arm_2_AUCell
rm( cell_rankings, cell_auc_object )
invisible( gc() )
cell_score_table <- data.table::data.table( cell = colnames( main_object ), biological_library = as.character( main_object$biological_library ), timepoint = as.character( main_object$timepoint ), myeloid_final_label = as.character( main_object$myeloid_final_label ), Arm_1_AUCell = main_object$Arm_1_AUCell, Arm_2_AUCell = main_object$Arm_2_AUCell, Arm_balance_AUCell = main_object$Arm_balance_AUCell )
library_auc_summary <- cell_score_table[ , .( n_cells = .N, Arm_1_AUCell_median = median( Arm_1_AUCell ), Arm_1_AUCell_mean = mean( Arm_1_AUCell ), Arm_2_AUCell_median = median( Arm_2_AUCell ), Arm_2_AUCell_mean = mean( Arm_2_AUCell ), Arm_balance_AUCell_median = median( Arm_balance_AUCell ), Arm_balance_AUCell_mean = mean( Arm_balance_AUCell ) ), by = .( biological_library, timepoint ) ]
library_order_table <- unique( library_auc_summary[ , c( "biological_library", "timepoint" ), with = FALSE ] )
library_order_table$timepoint <- factor( library_order_table$timepoint, levels = TIMEPOINT_LEVELS, ordered = TRUE )
data.table::setorder( library_order_table, timepoint, biological_library )
library_order <- library_order_table$biological_library
cell_library_factor <- factor( as.character( main_object$biological_library ), levels = library_order )
if ( anyNA( cell_library_factor ) ) {
  stop_with_report( "At least one main-analysis cell could not be mapped to a pseudobulk library." )
}
library_indicator <- Matrix::sparseMatrix( i = seq_len( ncol( main_object ) ), j = as.integer( cell_library_factor ), x = 1, dims = c( ncol( main_object ), length( library_order ) ), dimnames = list( colnames( main_object ), library_order ) )
pseudobulk_counts_sparse <- counts_matrix %*%
  library_indicator
if ( ncol( pseudobulk_counts_sparse ) != length( library_order ) || !identical( colnames( pseudobulk_counts_sparse ), library_order ) ) {
  stop_with_report( "Pseudobulk count aggregation failed the library-order check." )
}
pseudobulk_library_totals <- Matrix::colSums( pseudobulk_counts_sparse )
if ( any( pseudobulk_library_totals <= 0 ) ) {
  stop_with_report( "At least one pseudobulk library has zero total counts." )
}
pseudobulk_keep_genes <- Matrix::rowSums( pseudobulk_counts_sparse ) > 0
pseudobulk_counts <- as.matrix( pseudobulk_counts_sparse[ pseudobulk_keep_genes, , drop = FALSE ] )
rm( pseudobulk_counts_sparse, library_indicator )
invisible( gc() )
sample_manifest <- data.frame( biological_library = library_order, timepoint = as.character( library_order_table$timepoint[ match( library_order, library_order_table$biological_library ) ] ), stringsAsFactors = FALSE )
rownames( sample_manifest ) <- sample_manifest$biological_library
dge <- edgeR::DGEList( counts = pseudobulk_counts, samples = sample_manifest )
dge <- edgeR::calcNormFactors( dge, method = "TMM" )
pseudobulk_logcpm <- edgeR::cpm( dge, log = TRUE, prior.count = PSEUDOBULK_PRIOR_COUNT )
if ( any( !is.finite( pseudobulk_logcpm ) ) || !identical( colnames( pseudobulk_logcpm ), library_order ) ) {
  stop_with_report( "TMM-normalized pseudobulk logCPM failed finite-value or library-order checks." )
}
ssgsea_gene_sets <- lapply( arm_gene_sets, function(genes) {
    intersect( genes, rownames( pseudobulk_logcpm ) )
  }
)
if ( any( vapply( ssgsea_gene_sets, length, integer(1) ) < MIN_SIGNATURE_GENES_PRESENT ) ) {
  stop_with_report( "At least one frozen signature has too few genes in the pseudobulk logCPM matrix." )
}
ssgsea_scores <- run_ssgsea_version_safe( expression_matrix = pseudobulk_logcpm, gene_sets = ssgsea_gene_sets )
if ( !all( required_auc_rows %in% rownames( ssgsea_scores ) ) || !all( library_order %in% colnames( ssgsea_scores ) ) ) {
  stop_with_report( "ssGSEA output does not contain both frozen signatures and all libraries." )
}
ssgsea_scores <- ssgsea_scores[ required_auc_rows, library_order, drop = FALSE ]
if ( any( !is.finite( ssgsea_scores ) ) ) {
  stop_with_report( "ssGSEA returned non-finite values." )
}
library_scores <- merge( as.data.frame( library_auc_summary ), data.frame( biological_library = library_order, Arm_1_ssGSEA = as.numeric( ssgsea_scores[ "Arm_1_frozen", library_order ] ), Arm_2_ssGSEA = as.numeric( ssgsea_scores[ "Arm_2_frozen", library_order ] ), stringsAsFactors = FALSE ), by = "biological_library", all.x = TRUE, sort = FALSE )
library_scores <- library_scores[ match( library_order, library_scores$biological_library ), , drop = FALSE ]
library_scores$Arm_balance_ssGSEA <- library_scores$Arm_1_ssGSEA - library_scores$Arm_2_ssGSEA
library_scores$timepoint <- factor( library_scores$timepoint, levels = TIMEPOINT_LEVELS, ordered = TRUE )
library_scores$Arm_1_AUCell_median_z <- safe_zscore( library_scores$Arm_1_AUCell_median )
library_scores$Arm_2_AUCell_median_z <- safe_zscore( library_scores$Arm_2_AUCell_median )
library_scores$Arm_balance_AUCell_median_z <- safe_zscore( library_scores$Arm_balance_AUCell_median )
library_scores$Arm_1_ssGSEA_z <- safe_zscore( library_scores$Arm_1_ssGSEA )
library_scores$Arm_2_ssGSEA_z <- safe_zscore( library_scores$Arm_2_ssGSEA )
library_scores$Arm_balance_ssGSEA_z <- safe_zscore( library_scores$Arm_balance_ssGSEA )
arm1_concordance <- safe_cor_test( library_scores$Arm_1_AUCell_median, library_scores$Arm_1_ssGSEA )
arm2_concordance <- safe_cor_test( library_scores$Arm_2_AUCell_median, library_scores$Arm_2_ssGSEA )
balance_concordance <- safe_cor_test( library_scores$Arm_balance_AUCell_median, library_scores$Arm_balance_ssGSEA )
concordance_table <- data.frame( signature = c( "Arm_1", "Arm_2", "Arm_balance" ), AUCell_library_statistic = "median_per_cell_AUC", ssGSEA_library_statistic = "TMM_logCPM_pseudobulk_ssGSEA", spearman_rho = c( arm1_concordance$rho, arm2_concordance$rho, balance_concordance$rho ), p_value = c( arm1_concordance$p_value, arm2_concordance$p_value, balance_concordance$p_value ), n_libraries = c( arm1_concordance$n, arm2_concordance$n, balance_concordance$n ), stringsAsFactors = FALSE )
library_long_parts <- list( data.frame( record_type = "LIBRARY_SCORE", biological_library = library_scores$biological_library, timepoint = as.character( library_scores$timepoint ), n_cells = library_scores$n_cells, n_libraries = 1L, method = "AUCell_library_median", signature = "Arm_1", score = library_scores$Arm_1_AUCell_median, score_z = library_scores$Arm_1_AUCell_median_z, stringsAsFactors = FALSE ), data.frame( record_type = "LIBRARY_SCORE", biological_library = library_scores$biological_library, timepoint = as.character( library_scores$timepoint ), n_cells = library_scores$n_cells, n_libraries = 1L, method = "AUCell_library_median", signature = "Arm_2", score = library_scores$Arm_2_AUCell_median, score_z = library_scores$Arm_2_AUCell_median_z, stringsAsFactors = FALSE ), data.frame( record_type = "LIBRARY_SCORE", biological_library = library_scores$biological_library, timepoint = as.character( library_scores$timepoint ), n_cells = library_scores$n_cells, n_libraries = 1L, method = "AUCell_library_median", signature = "Arm_balance_Arm1_minus_Arm2", score = library_scores$Arm_balance_AUCell_median, score_z = library_scores$Arm_balance_AUCell_median_z, stringsAsFactors = FALSE ), data.frame( record_type = "LIBRARY_SCORE", biological_library = library_scores$biological_library, timepoint = as.character( library_scores$timepoint ), n_cells = library_scores$n_cells, n_libraries = 1L, method = "ssGSEA_pseudobulk", signature = "Arm_1", score = library_scores$Arm_1_ssGSEA, score_z = library_scores$Arm_1_ssGSEA_z, stringsAsFactors = FALSE ), data.frame( record_type = "LIBRARY_SCORE", biological_library = library_scores$biological_library, timepoint = as.character( library_scores$timepoint ), n_cells = library_scores$n_cells, n_libraries = 1L, method = "ssGSEA_pseudobulk", signature = "Arm_2", score = library_scores$Arm_2_ssGSEA, score_z = library_scores$Arm_2_ssGSEA_z, stringsAsFactors = FALSE ), data.frame( record_type = "LIBRARY_SCORE", biological_library = library_scores$biological_library, timepoint = as.character( library_scores$timepoint ), n_cells = library_scores$n_cells, n_libraries = 1L, method = "ssGSEA_pseudobulk", signature = "Arm_balance_Arm1_minus_Arm2", score = library_scores$Arm_balance_ssGSEA, score_z = library_scores$Arm_balance_ssGSEA_z, stringsAsFactors = FALSE ) )
library_score_long <- do.call( rbind, library_long_parts )
library_score_long$timepoint <- factor( library_score_long$timepoint, levels = TIMEPOINT_LEVELS, ordered = TRUE )
timepoint_score_summary <- data.table::as.data.table( library_score_long )[ , .( n_libraries = uniqueN( biological_library ), n_cells = sum( n_cells ), mean_score = mean( score ), median_score = median( score ), minimum_score = min( score ), maximum_score = max( score ), sd_score = stats::sd( score ), mean_score_z = mean( score_z ), median_score_z = median( score_z ) ), by = .( timepoint, method, signature ) ]
timepoint_score_rows <- data.frame( record_type = "TIMEPOINT_SUMMARY", biological_library = "ALL", timepoint = as.character( timepoint_score_summary$timepoint ), n_cells = timepoint_score_summary$n_cells, n_libraries = timepoint_score_summary$n_libraries, method = timepoint_score_summary$method, signature = timepoint_score_summary$signature, score = timepoint_score_summary$median_score, score_z = timepoint_score_summary$median_score_z, mean_score = timepoint_score_summary$mean_score, median_score = timepoint_score_summary$median_score, minimum_score = timepoint_score_summary$minimum_score, maximum_score = timepoint_score_summary$maximum_score, sd_score = timepoint_score_summary$sd_score, mean_score_z = timepoint_score_summary$mean_score_z, median_score_z = timepoint_score_summary$median_score_z, stringsAsFactors = FALSE )
library_score_rows <- data.frame( library_score_long, mean_score = NA_real_, median_score = NA_real_, minimum_score = NA_real_, maximum_score = NA_real_, sd_score = NA_real_, mean_score_z = NA_real_, median_score_z = NA_real_, stringsAsFactors = FALSE )
combined_score_table <- rbind( library_score_rows, timepoint_score_rows )
combined_score_table$timepoint <- factor( combined_score_table$timepoint, levels = TIMEPOINT_LEVELS, ordered = TRUE )
combined_score_table <- combined_score_table[ order( combined_score_table$record_type, combined_score_table$method, combined_score_table$signature, combined_score_table$timepoint, combined_score_table$biological_library ), , drop = FALSE ]
combined_score_table$timepoint <- as.character( combined_score_table$timepoint )
write_csv_safe( combined_score_table, "05_main_validation_library_and_timepoint_scores.csv" )
write_csv_safe( library_scores, "06_library_scores_wide.csv" )
write_csv_safe( concordance_table, "07_AUCell_ssGSEA_concordance.csv" )
umap_embeddings <- Seurat::Embeddings( object = main_object, reduction = FROZEN_UMAP )
umap_plot_data <- data.frame( UMAP_1 = umap_embeddings[ , 1L ], UMAP_2 = umap_embeddings[ , 2L ], Arm_1_AUCell = main_object$Arm_1_AUCell, Arm_2_AUCell = main_object$Arm_2_AUCell, Arm_balance_AUCell = main_object$Arm_balance_AUCell, timepoint = main_object$timepoint, stringsAsFactors = FALSE )
p_arm1_umap <- plot_umap_score( plot_data = umap_plot_data, score_column = "Arm_1_AUCell", title_text = "Frozen Arm_1 AUCell", legend_title = "Normalized AUC" )
p_arm2_umap <- plot_umap_score( plot_data = umap_plot_data, score_column = "Arm_2_AUCell", title_text = "Frozen Arm_2 AUCell", legend_title = "Normalized AUC" )
p_balance_umap <- plot_umap_score( plot_data = umap_plot_data, score_column = "Arm_balance_AUCell", title_text = "AUCell Arm balance", legend_title = "Arm_1 - Arm_2" )
p_timepoint_umap <- ggplot2::ggplot( umap_plot_data, ggplot2::aes( x = UMAP_1, y = UMAP_2, colour = timepoint ) ) + ggplot2::geom_point( size = 0.28, alpha = 0.72 ) + ggplot2::coord_equal() + ggplot2::labs( title = "Main validation timepoints", x = "Myeloid Harmony UMAP 1", y = "Myeloid Harmony UMAP 2", colour = "Timepoint" ) + ggplot2::theme_bw()
ggplot2::ggsave( filename = file.path( FIGURES_DIR, "01_AUCell_Arm_scores_UMAP.pdf" ), plot = ( p_arm1_umap | p_arm2_umap ) / ( p_balance_umap | p_timepoint_umap ), width = 19, height = 14, limitsize = FALSE )
plot_library_long <- library_score_long
plot_library_long$timepoint <- factor( as.character( plot_library_long$timepoint ), levels = TIMEPOINT_LEVELS, ordered = TRUE )
p_library_timecourse <- ggplot2::ggplot( plot_library_long, ggplot2::aes( x = timepoint, y = score, group = biological_library ) ) + ggplot2::geom_point( ggplot2::aes( shape = biological_library ), size = 2.4, position = ggplot2::position_jitter( width = 0.08, height = 0 ) ) + ggplot2::stat_summary( ggplot2::aes( group = 1 ), fun = median, geom = "line", linewidth = 0.8 ) + ggplot2::stat_summary( ggplot2::aes( group = 1 ), fun = median, geom = "point", size = 3.0 ) + ggplot2::facet_grid( rows = ggplot2::vars( method ), cols = ggplot2::vars( signature ), scales = "free_y" ) + ggplot2::labs( title = "Frozen Arm scores by biological library", subtitle = "Points are biological libraries; line is the timepoint median", x = NULL, y = "Method-specific score", shape = "Biological library" ) + ggplot2::theme_bw() + ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 ), legend.position = "none" )
ggplot2::ggsave( filename = file.path( FIGURES_DIR, "02_library_AUCell_ssGSEA_timecourse.pdf" ), plot = p_library_timecourse, width = 18, height = 10, limitsize = FALSE )
concordance_plot_data <- rbind( data.frame( biological_library = library_scores$biological_library, timepoint = library_scores$timepoint, signature = "Arm_1", AUCell_z = library_scores$Arm_1_AUCell_median_z, ssGSEA_z = library_scores$Arm_1_ssGSEA_z, stringsAsFactors = FALSE ), data.frame( biological_library = library_scores$biological_library, timepoint = library_scores$timepoint, signature = "Arm_2", AUCell_z = library_scores$Arm_2_AUCell_median_z, ssGSEA_z = library_scores$Arm_2_ssGSEA_z, stringsAsFactors = FALSE ), data.frame( biological_library = library_scores$biological_library, timepoint = library_scores$timepoint, signature = "Arm_balance", AUCell_z = library_scores$Arm_balance_AUCell_median_z, ssGSEA_z = library_scores$Arm_balance_ssGSEA_z, stringsAsFactors = FALSE ) )
concordance_labels <- data.frame( signature = concordance_table$signature, label = paste0( "Spearman rho=", ifelse( is.na( concordance_table$spearman_rho ), "NA", format( round( concordance_table$spearman_rho, 3 ), nsmall = 3 ) ), "\nn=", concordance_table$n_libraries ), x = -Inf, y = Inf, stringsAsFactors = FALSE )
p_concordance <- ggplot2::ggplot( concordance_plot_data, ggplot2::aes( x = AUCell_z, y = ssGSEA_z, shape = timepoint ) ) + ggplot2::geom_hline( yintercept = 0, linetype = 3 ) + ggplot2::geom_vline( xintercept = 0, linetype = 3 ) + ggplot2::geom_point( size = 2.8 ) + ggplot2::geom_smooth( method = "lm", se = FALSE, linewidth = 0.7 ) + ggplot2::geom_text( data = concordance_labels, ggplot2::aes( x = x, y = y, label = label ), inherit.aes = FALSE, hjust = -0.05, vjust = 1.1, size = 3.2 ) + ggplot2::facet_wrap( ggplot2::vars( signature ), scales = "free" ) + ggplot2::labs( title = "Library-level AUCell and pseudobulk ssGSEA concordance", subtitle = "Scores are standardized within each method before plotting", x = "AUCell library-median Z score", y = "Pseudobulk ssGSEA Z score", shape = "Timepoint" ) + ggplot2::theme_bw()
ggplot2::ggsave( filename = file.path( FIGURES_DIR, "03_AUCell_ssGSEA_library_concordance.pdf" ), plot = p_concordance, width = 15, height = 6.5, limitsize = FALSE )
main_object@misc$GSE234774_step11A <- list( main_timepoints = MAIN_TIMEPOINTS, statistical_unit = "biological_library", frozen_Arm1_definition = "Step6A Arm1 high-confidence union Step6B M3_Shared_prebranch", frozen_Arm2_definition = "Step6A Arm2 high-confidence union Step6B M1/M2 Arm2_terminal", Arm1_requested_genes = arm1_frozen_requested, Arm2_requested_genes = arm2_frozen_requested, Arm1_present_genes = arm1_present, Arm2_present_genes = arm2_present, overlapping_present_genes = overlap_present, AUCell_aucMaxRank = auc_max_rank, AUCell_normAUC = TRUE, pseudobulk_normalization = "edgeR TMM logCPM", pseudobulk_prior_count = PSEUDOBULK_PRIOR_COUNT, ssGSEA_alpha = SSGSEA_ALPHA, ssGSEA_normalize = TRUE, cell_level_p_values_run = FALSE, sensitivity_analysis_run = FALSE, cycling_analysis_run = FALSE, Arm_signatures_redefined = FALSE, completed_at = as.character( Sys.time() ) )
saveRDS( main_object, SCORED_OBJECT_PATH, compress = FALSE )
results_object <- list( main_timepoints = MAIN_TIMEPOINTS, library_manifest = as.data.frame( library_cell_counts ), signature_gene_sets_requested = list( Arm_1_frozen = arm1_frozen_requested, Arm_2_frozen = arm2_frozen_requested ), signature_gene_sets_present = arm_gene_sets, signature_overlap_requested = overlap_requested, signature_overlap_present = overlap_present, AUCell_aucMaxRank = auc_max_rank, cell_AUCell = cell_auc_matrix, library_scores = library_scores, pseudobulk_counts = pseudobulk_counts, pseudobulk_logCPM = pseudobulk_logcpm, ssGSEA_scores = ssgsea_scores, concordance = concordance_table, timepoint_summary = as.data.frame( timepoint_score_summary ), completed_at = as.character( Sys.time() ) )
saveRDS( results_object, RESULTS_RDS_PATH, compress = TRUE )
reloaded_scored_object <- tryCatch( readRDS( SCORED_OBJECT_PATH ), error = function(error_object) {
    error_object
  }
)
scored_object_reload_pass <- !inherits( reloaded_scored_object, "error" ) && inherits( reloaded_scored_object, "Seurat" ) && ncol( reloaded_scored_object ) == ncol( main_object ) && identical( colnames( reloaded_scored_object ), colnames( main_object ) ) && all( c( "Arm_1_AUCell", "Arm_2_AUCell", "Arm_balance_AUCell" ) %in% colnames( reloaded_scored_object@meta.data ) )
rm( reloaded_scored_object )
invisible( gc() )
reloaded_results <- tryCatch( readRDS( RESULTS_RDS_PATH ), error = function(error_object) {
    error_object
  }
)
results_reload_pass <- !inherits( reloaded_results, "error" ) && is.list( reloaded_results ) && all( c( "library_scores", "ssGSEA_scores", "cell_AUCell", "signature_gene_sets_present" ) %in% names( reloaded_results ) )
rm( reloaded_results )
invisible( gc() )
required_figure_paths <- file.path( FIGURES_DIR, c( "01_AUCell_Arm_scores_UMAP.pdf", "02_library_AUCell_ssGSEA_timecourse.pdf", "03_AUCell_ssGSEA_library_concordance.pdf" ) )
required_table_paths <- file.path( TABLES_DIR, c( "04_frozen_Arm_signature_gene_audit.csv", "05_main_validation_library_and_timepoint_scores.csv" ) )
audit_checks <- data.frame( check = c( "Frozen Step 10G primary object dimensions and cluster membership are correct", "Only the five prespecified main validation timepoints are included", "Every main timepoint has at least two biological libraries", "IFN-responsive cluster 7 is retained and cycling cluster 6 is absent", "Frozen GSE162610 Step6A/Step6B source files were loaded without redefining Arm labels", "Both frozen Arm signatures meet external mapping thresholds", "AUCell returned finite Arm_1 and Arm_2 scores for every main-analysis cell", "Pseudobulk aggregation produced one nonzero profile per biological library", "ssGSEA returned finite Arm_1 and Arm_2 scores for every biological library", "AUCell and ssGSEA library-level concordance statistics were produced", "All three required figures and two core tables were created", "Scored Seurat object and complete results RDS reload successfully", "Cell-level P values were calculated", "Non-IFN sensitivity, 1 dpi, 4 dpi, or cycling analyses were run", "Frozen Arm_1/Arm_2 signatures were redefined" ), required_for_pass = c( TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE ), status = c( "PASS", ifelse( setequal( unique( as.character( main_object$timepoint ) ), MAIN_TIMEPOINTS ), "PASS", "STOP" ), ifelse( all( timepoint_replication$n_libraries >= MIN_LIBRARIES_PER_MAIN_TIMEPOINT ), "PASS", "STOP" ), ifelse( "IFN_responsive_microglia" %in% as.character( main_object$myeloid_final_label ) && !"Cycling_microglia" %in% as.character( main_object$myeloid_final_label ), "PASS", "STOP" ), "PASS", ifelse( length( arm1_present ) >= MIN_SIGNATURE_GENES_PRESENT && length( arm2_present ) >= MIN_SIGNATURE_GENES_PRESENT && arm1_mapping_fraction >= MIN_SIGNATURE_MAPPING_FRACTION && arm2_mapping_fraction >= MIN_SIGNATURE_MAPPING_FRACTION, "PASS", "STOP" ), ifelse( ncol( cell_auc_matrix ) == ncol( main_object ) && all( is.finite( cell_auc_matrix[ required_auc_rows, , drop = FALSE ] ) ), "PASS", "STOP" ), ifelse( ncol( pseudobulk_counts ) == length( library_order ) && all( colSums( pseudobulk_counts ) > 0 ), "PASS", "STOP" ), ifelse( ncol( ssgsea_scores ) == length( library_order ) && all( is.finite( ssgsea_scores ) ), "PASS", "STOP" ), ifelse( nrow( concordance_table ) == 3L && all( concordance_table$n_libraries == length( library_order ) ), "PASS", "STOP" ), ifelse( all( file.exists( required_figure_paths ) ) && all( file.exists( required_table_paths ) ), "PASS", "STOP" ), ifelse( scored_object_reload_pass && results_reload_pass, "PASS", "STOP" ), "CONFIRMED_NOT_RUN", "CONFIRMED_NOT_RUN", "CONFIRMED_NOT_REDEFINED" ), detail = c( paste0( nrow( primary_object ), " features x ", ncol( primary_object ), " nuclei; clusters=", paste( input_cluster_values, collapse = "," ) ), paste( MAIN_TIMEPOINTS, collapse = ", " ), paste( paste0( as.character( timepoint_replication$timepoint ), "=", timepoint_replication$n_libraries, " libraries" ), collapse = "; " ), paste0( "labels=", paste( sort( unique( as.character( main_object$myeloid_final_label ) ) ), collapse = ";" ) ), paste0( "Arm1=Step6A Arm1 highconf + ", paste( ARM1_STEP6B_MODULES, collapse = "+" ), "; Arm2=Step6A Arm2 highconf + ", paste( ARM2_STEP6B_MODULES, collapse = "+" ) ), paste0( "Arm1 ", length( arm1_present ), "/", length( arm1_frozen_requested ), " (", round( arm1_mapping_fraction, 3 ), "); Arm2 ", length( arm2_present ), "/", length( arm2_frozen_requested ), " (", round( arm2_mapping_fraction, 3 ), "); retained overlap=", length( overlap_present ) ), paste0( ncol( cell_auc_matrix ), " cells; aucMaxRank=", auc_max_rank ), paste0( ncol( pseudobulk_counts ), " libraries; total-count range=", min( colSums( pseudobulk_counts ) ), "-", max( colSums( pseudobulk_counts ) ) ), paste0( ncol( ssgsea_scores ), " libraries; GSVA version=", safe_package_version( "GSVA" ) ), paste( paste0( concordance_table$signature, " rho=", round( concordance_table$spearman_rho, 3 ) ), collapse = "; " ), paste( c( basename( required_figure_paths ), basename( required_table_paths ) ), collapse = "; " ), paste0( "scored object reload=", scored_object_reload_pass, "; results RDS reload=", results_reload_pass ), "No cell-level inferential test was run.", "Main timepoints only; no sensitivity or cycling analysis.", "Frozen source genes were mapped and scored without data-driven reselection." ), stringsAsFactors = FALSE )
write_csv_safe( audit_checks, "09_STEP11A_audit_checks.csv" )
required_failures <- audit_checks$required_for_pass & audit_checks$status != "PASS"
final_decision <- if ( any( required_failures ) ) {
  "STOP"
} else {
  "PASS_WITH_MAIN_ARM_VALIDATION_REVIEW"
}
timepoint_lines <- apply( as.data.frame( timepoint_replication ), 1L, function(row) {
    paste0( row[["timepoint"]], ": ", row[["n_libraries"]], " biological libraries; ", row[["n_cells"]], " nuclei; library-cell range ", row[["minimum_library_cells"]], "-", row[["maximum_library_cells"]] )
  }
)
concordance_lines <- apply( concordance_table, 1L, function(row) {
    paste0( row[["signature"]], ": Spearman rho=", ifelse( is.na( as.numeric( row[["spearman_rho"]] ) ), "NA", round( as.numeric( row[["spearman_rho"]] ), 3 ) ), "; p=", ifelse( is.na( as.numeric( row[["p_value"]] ) ), "NA", signif( as.numeric( row[["p_value"]] ), 4 ) ), "; n=", row[["n_libraries"]] )
  }
)
if ( final_decision == "STOP" ) {
  failed <- audit_checks[ audit_checks$required_for_pass & audit_checks$status != "PASS", , drop = FALSE ]
  for ( index in seq_len( nrow( failed ) ) ) {
  }
}
report_path <- file.path( OUTPUT_DIR, "10_STEP11A_FINAL_REPORT.txt" )
if ( final_decision == "STOP" ) {
  stop( "GSE234774 Step 11A decision is STOP. Review the final report.", call. = FALSE )
}
