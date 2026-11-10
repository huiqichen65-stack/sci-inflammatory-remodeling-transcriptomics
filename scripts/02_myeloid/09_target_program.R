options(stringsAsFactors = FALSE)
set.seed(20260731)
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
STEP3D_RDS <- file.path( PROJECT_ROOT, "04_results", "step3D_GSE162610_time_DEG_pathway_confirmation", "GSE162610_step3D_time_DEG_pathway_confirmed.rds" )
STEP3E_DIR <- file.path( PROJECT_ROOT, "04_results", "step3E_GSE162610_noncycling_transition_refinement" )
OUTDIR <- file.path( PROJECT_ROOT, "04_results", "step4A_GSE162610_repair_scar_candidate_program_discovery" )
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
OUT_FIG1 <- file.path(OUTDIR, "Fig4A_01_7dpi_program_UMAP.pdf")
OUT_FIG2 <- file.path(OUTDIR, "Fig4A_02_program_distribution.pdf")
OUT_FIG3 <- file.path(OUTDIR, "Fig4A_03_repair_scar_continuous_map.pdf")
OUT_XLSX <- file.path(OUTDIR, "Table4A_program_construction_and_scores.xlsx")
OUT_RDS  <- file.path(OUTDIR, "GSE162610_step4A_candidate_program_scored.rds")
OUT_REP  <- file.path(OUTDIR, "Step4A_report.txt")
OUT_SESS <- file.path(OUTDIR, "Step4A_sessionInfo.txt")
TIME_COL       <- "timepoint"
SAMPLE_COL     <- "sample_id"
CLUSTER_COL    <- "step3A_cluster_r04"
UMAP_REDUCTION <- "umap.microglia.step3A"
CYCLING_CLUSTERS <- c("6", "8")
IFN_CLUSTER      <- "9"
CORE_ACTIVATED_DAM_CLUSTERS <- c("0", "1", "4")
EXPECTED_RETAINED_STEP3_CLUSTERS <- c("0", "1", "2", "3", "4", "5", "6", "8", "9", "10")
MIN_GENE_DETECTION_PROP <- 0.02
MIN_GENES_PER_COMPONENT <- 8
MAX_GENES_PER_COMPONENT <- 250
MAX_SCATTER_CELLS <- 6000
required_pkgs <- c( "Seurat", "SeuratObject", "Matrix", "dplyr", "tidyr", "tibble", "ggplot2", "patchwork", "openxlsx", "AnnotationDbi", "org.Mm.eg.db", "GO.db" )
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
COL_CLUSTER <- c( "0"  = "#4E79A7", "1"  = "#E15759", "2"  = "#59A14F", "3"  = "#F28E2B", "4"  = "#B07AA1", "5"  = "#76B7B2", "6"  = "#EDC948", "8"  = "#FF9DA7", "9"  = "#BAB0AC", "10" = "#2F4B7C" )
COL_SAMPLE <- c( "7dpi_sample1" = "#3B82F6", "7dpi_sample2" = "#F97316" )
COL_REPAIR <- c("#F7FBFF", "#6BAED6", "#08306B")
COL_SCAR   <- c("#FFF5F0", "#FB6A4A", "#99000D")
COL_DIFF   <- c("#B2182B", "#F7F7F7", "#2166AC")
msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
}
safe_z <- function(x) {
  x <- as.numeric(x)
  sx <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sx) || sx == 0) return(rep(0, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / sx)
}
get_data_matrix <- function(obj) {
  assay <- DefaultAssay(obj)
  x <- tryCatch( Seurat::GetAssayData(obj, assay = assay, layer = "data"), error = function(e) NULL )
  if (!is.null(x)) return(x)
  Seurat::GetAssayData(obj, assay = assay, slot = "data")
}
gene_detection_prop <- function(mat) {
  Matrix::rowMeans(mat > 0)
}
gene_variance_sparse <- function(mat) {
  m1 <- Matrix::rowMeans(mat)
  m2 <- Matrix::rowMeans(mat ^ 2)
  vv <- m2 - m1 ^ 2
  vv[vv < 0] <- 0
  vv
}
mean_gene_score <- function(mat, genes) {
  genes <- intersect(genes, rownames(mat))
  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(mat)))
  }
  as.numeric(Matrix::colMeans(mat[genes, , drop = FALSE]))
}
build_go_bp_term_table <- function() {
  term_vec <- AnnotationDbi::Term(GO.db::GOTERM)
  ont_vec  <- AnnotationDbi::Ontology(GO.db::GOTERM)
  tibble::tibble( GO_ID = names(term_vec), TERM = unname(term_vec), ONTOLOGY = unname(ont_vec[names(term_vec)]) ) |> dplyr::filter( !is.na(TERM), ONTOLOGY == "BP" )
}
go_component_genes <- function( regex, go_term_table, expressed_genes, detection, variance_vec, min_detection = 0.02, max_genes = 250) {
  matched_terms <- go_term_table |> dplyr::filter( grepl(regex, TERM, ignore.case = TRUE, perl = TRUE) )
  if (nrow(matched_terms) == 0) {
    return(list( genes = character(), terms = character(), n_terms = 0 ))
  }
  go_ids <- unique(matched_terms$GO_ID)
  map <- suppressMessages( tryCatch( AnnotationDbi::select( org.Mm.eg.db::org.Mm.eg.db, keys = go_ids, keytype = "GOALL", columns = c("SYMBOL", "GOALL") ), error = function(e) NULL ) )
  if (is.null(map) || nrow(map) == 0) {
    map <- suppressMessages( tryCatch( AnnotationDbi::select( org.Mm.eg.db::org.Mm.eg.db, keys = go_ids, keytype = "GO", columns = c("SYMBOL", "GO") ), error = function(e) NULL ) )
  }
  if (is.null(map) || nrow(map) == 0 || !"SYMBOL" %in% colnames(map)) {
    return(list( genes = character(), terms = unique(matched_terms$TERM), n_terms = nrow(matched_terms) ))
  }
  genes <- unique(na.omit(as.character(map$SYMBOL)))
  genes <- intersect(genes, expressed_genes)
  if (length(genes) > 0) {
    genes <- genes[detection[genes] >= min_detection]
  }
  if (length(genes) > max_genes) {
    ord <- order( variance_vec[genes], decreasing = TRUE, na.last = NA )
    genes <- genes[ord][seq_len(max_genes)]
  }
  list( genes = unique(genes), terms = unique(matched_terms$TERM), n_terms = nrow(matched_terms) )
}
score_component_set <- function(expr_mat, component_list, prefix) {
  out <- lapply(names(component_list), function(component_name) {
    genes <- component_list[[component_name]]
    raw_score <- mean_gene_score( expr_mat, genes )
    tibble::tibble( cell = colnames(expr_mat), !!paste0(prefix, "_", component_name, "_raw") := raw_score, !!paste0(prefix, "_", component_name, "_z") := safe_z(raw_score) )
  })
  Reduce( function(x, y) dplyr::left_join(x, y, by = "cell"), out )
}
msg("====================================================")
msg("GSE162610 Step 4A | repair/scar candidate programs")
msg("====================================================")
msg("Reading formal Step 3D object...")
if (!file.exists(STEP3D_RDS)) {
  stop( paste0( "找不到Step 3D正式输入RDS：\n", STEP3D_RDS ), call. = FALSE )
}
obj <- readRDS(STEP3D_RDS)
if (!inherits(obj, "Seurat")) {
  stop("Step 3D输入不是Seurat object。", call. = FALSE)
}
DefaultAssay(obj) <- "RNA"
msg("Loaded Step 3D object: ", ncol(obj), " cells; ", nrow(obj), " features.")
required_meta <- c( TIME_COL, SAMPLE_COL, CLUSTER_COL, "transition_status", "shared_transition_3dpi_highconf" )
missing_meta <- required_meta[ !required_meta %in% colnames(obj[[]]) ]
if (length(missing_meta) > 0) {
  stop( paste0( "Step 3D对象缺少Step 3正式metadata：\n ", paste(missing_meta, collapse = ", "), "\n\n本脚本不会自动猜测其他字段。" ), call. = FALSE )
}
available_reductions <- Seurat::Reductions(obj)
if (!UMAP_REDUCTION %in% available_reductions) {
  stop( paste0( "Step 3D对象缺少正式microglia UMAP：", UMAP_REDUCTION, "\n现有reductions：", paste(available_reductions, collapse = ", "), "\n\n本脚本不会重新计算UMAP，也不会使用其他降维代替。" ), call. = FALSE )
}
obj$step3A_cluster_r04 <- as.character(obj$step3A_cluster_r04)
observed_clusters_all <- sort( unique(obj$step3A_cluster_r04) )
unexpected_clusters <- setdiff( observed_clusters_all, EXPECTED_RETAINED_STEP3_CLUSTERS )
if (length(unexpected_clusters) > 0) {
  warning(paste0("Unexpected Step 3B clusters: ", paste(unexpected_clusters, collapse = ", ")))
}
msg("Authoritative cluster field: ", CLUSTER_COL)
msg("Authoritative UMAP reduction: ", UMAP_REDUCTION)
msg("Observed Step 3 clusters: ", paste(observed_clusters_all, collapse = ", "))
time_vec    <- as.character(obj[[TIME_COL, drop = TRUE]])
sample_vec  <- as.character(obj[[SAMPLE_COL, drop = TRUE]])
cluster_vec <- as.character(obj[[CLUSTER_COL, drop = TRUE]])
all_7dpi_cells <- colnames(obj)[ time_vec == "7 dpi" ]
if (length(all_7dpi_cells) == 0) {
  stop("Step 3D对象中未找到timepoint == '7 dpi'的细胞。", call. = FALSE)
}
primary_keep <- ( time_vec == "7 dpi" & !(cluster_vec %in% CYCLING_CLUSTERS) & cluster_vec != IFN_CLUSTER )
primary_cells <- colnames(obj)[primary_keep]
sensitivity_noncycling_keep <- ( time_vec == "7 dpi" & !(cluster_vec %in% CYCLING_CLUSTERS) )
sensitivity_noncycling_cells <- colnames(obj)[ sensitivity_noncycling_keep ]
core_adam_keep <- ( primary_keep & cluster_vec %in% CORE_ACTIVATED_DAM_CLUSTERS )
core_adam_cells <- colnames(obj)[core_adam_keep]
seven_samples <- sort( unique(sample_vec[time_vec == "7 dpi"]) )
observed_primary_clusters <- sort( unique(cluster_vec[primary_keep]) )
msg("7 dpi all resident microglia: ", length(all_7dpi_cells))
msg( "Primary 7 dpi landscape after excluding cycling 6/8 and IFN 9: ", length(primary_cells) )
msg( "Core activated/DAM reference cells (clusters 0/1/4): ", length(core_adam_cells) )
msg( "Primary Step 3 clusters: ", paste(observed_primary_clusters, collapse = ", ") )
msg( "7 dpi biological samples: ", paste(seven_samples, collapse = ", ") )
if (length(primary_cells) < 100) {
  stop("Primary 7 dpi landscape has <100 cells.", call. = FALSE)
}
if (length(seven_samples) < 2) {
  warning("7 dpi biological samples <2，不能形成replicate-level结论。")
}
obj$step4A_analysis_status <- "not_primary"
obj$step4A_analysis_status[ time_vec == "7 dpi" ] <- "7dpi_all"
obj$step4A_analysis_status[ sensitivity_noncycling_keep ] <- "7dpi_noncycling"
obj$step4A_analysis_status[ primary_keep ] <- "7dpi_primary_noncycling_nonIFN"
obj$step4A_core_activated_DAM_reference <- core_adam_keep
msg("Preparing normalized RNA expression matrix...")
expr_all <- get_data_matrix(obj)
if (!identical(colnames(expr_all), colnames(obj))) {
  expr_all <- expr_all[, colnames(obj), drop = FALSE]
}
expr_primary <- expr_all[, primary_cells, drop = FALSE]
det <- gene_detection_prop(expr_primary)
varv <- gene_variance_sparse(expr_primary)
names(det)  <- rownames(expr_primary)
names(varv) <- rownames(expr_primary)
expressed_genes <- names(det)[ det >= MIN_GENE_DETECTION_PROP ]
msg( "Genes detected in >= ", sprintf("%.1f%%", MIN_GENE_DETECTION_PROP * 100), " of primary 7 dpi cells: ", length(expressed_genes) )
msg("Building GO-BP functional candidate components...")
go_bp <- build_go_bp_term_table()
repair_component_regex <- list( Phagocytosis_Efferocytosis = paste0( "phagocyt|", "efferocyt|", "engulfment of apoptotic|", "clearance of apoptotic" ), Lysosome_Endocytosis = paste0( "lysosom|", "endocytosis|", "endocytic" ), Lipid_Handling_Catabolism = paste0( "lipid catabolic|", "lipid homeostasis|", "lipid transport|", "lipoprotein|", "cholesterol efflux|", "cholesterol transport|", "fatty acid catabolic|", "fatty acid oxidation|", "fatty acid transport" ), Inflammation_Resolution = paste0( "negative regulation of inflammatory response|", "negative regulation of cytokine production|", "resolution of inflammation" ), Neuro_Myelin_Trophic_Support = paste0( "axon regeneration|", "neuron projection regeneration|", "myelination|", "axon ensheathment|", "oligodendrocyte differentiation|", "neurotrophin" ), Vascular_Endothelial_Support = paste0( "angiogenesis|", "blood vessel development|", "blood vessel morphogenesis|", "endothelial cell migration|", "vascular endothelial" ) )
scar_component_regex <- list( Cytokine_Chemokine_Inflammation = paste0( "positive regulation of inflammatory response|", "inflammatory response|", "chemokine production|", "chemokine-mediated signaling|", "cytokine-mediated signaling|", "positive regulation of cytokine production|", "tumor necrosis factor|", "interleukin-1" ), Complement_Coagulation = paste0( "complement activation|", "regulation of complement|", "blood coagulation|", "coagulation" ), ECM_Matrix_Remodeling = paste0( "extracellular matrix organization|", "extracellular matrix assembly|", "extracellular matrix disassembly|", "collagen fibril organization|", "collagen metabolic|", "matrix organization" ), Fibrosis_TGFbeta = paste0( "fibroblast proliferation|", "fibroblast migration|", "fibrosis|", "transforming growth factor beta|", "TGF-beta" ), Chronic_Innate_Signaling = paste0( "toll-like receptor signaling|", "NOD.*signaling|", "innate immune response|", "NF-kappaB.*signaling" ) )
build_component_list <- function(regex_list, program_name) {
  out <- lapply(names(regex_list), function(component_name) {
    tmp <- go_component_genes( regex = regex_list[[component_name]], go_term_table = go_bp, expressed_genes = expressed_genes, detection = det, variance_vec = varv, min_detection = MIN_GENE_DETECTION_PROP, max_genes = MAX_GENES_PER_COMPONENT )
    if (length(tmp$genes) < MIN_GENES_PER_COMPONENT) {
      warning( paste0( program_name, " / ", component_name, " retained only ", length(tmp$genes), " genes after expression filtering." ) )
    }
    list( program = program_name, component = component_name, regex = regex_list[[component_name]], genes = tmp$genes, terms = tmp$terms, n_terms = tmp$n_terms )
  })
  names(out) <- names(regex_list)
  out
}
repair_components <- build_component_list( repair_component_regex, "Repair_related_candidate" )
scar_components <- build_component_list( scar_component_regex, "Scar_related_candidate" )
repair_union_raw <- unique( unlist( lapply(repair_components, `[[`, "genes"), use.names = FALSE ) )
scar_union_raw <- unique( unlist( lapply(scar_components, `[[`, "genes"), use.names = FALSE ) )
cross_program_overlap <- intersect( repair_union_raw, scar_union_raw )
repair_union_unique <- setdiff( repair_union_raw, cross_program_overlap )
scar_union_unique <- setdiff( scar_union_raw, cross_program_overlap )
repair_components_score <- lapply( repair_components, function(x) setdiff(x$genes, cross_program_overlap) )
scar_components_score <- lapply( scar_components, function(x) setdiff(x$genes, cross_program_overlap) )
msg("Raw repair program genes: ", length(repair_union_raw))
msg("Raw scar program genes: ", length(scar_union_raw))
msg( "Repair/scar overlapping genes removed from both composites: ", length(cross_program_overlap) )
msg("Unique repair composite genes: ", length(repair_union_unique))
msg("Unique scar composite genes: ", length(scar_union_unique))
if (length(repair_union_unique) < 20) {
  warning("Final repair composite contains <20 unique genes.")
}
if (length(scar_union_unique) < 20) {
  warning("Final scar composite contains <20 unique genes.")
}
msg("Scoring primary 7 dpi cells...")
repair_scores <- score_component_set( expr_primary, repair_components_score, "repair" )
scar_scores <- score_component_set( expr_primary, scar_components_score, "scar" )
score_tbl <- dplyr::left_join( repair_scores, scar_scores, by = "cell" )
repair_z_cols <- grep( "^repair_.*_z$", colnames(score_tbl), value = TRUE )
scar_z_cols <- grep( "^scar_.*_z$", colnames(score_tbl), value = TRUE )
score_tbl$repair_candidate_score <- rowMeans( as.matrix( score_tbl[, repair_z_cols, drop = FALSE] ), na.rm = TRUE )
score_tbl$scar_candidate_score <- rowMeans( as.matrix( score_tbl[, scar_z_cols, drop = FALSE] ), na.rm = TRUE )
score_tbl$repair_minus_scar <- score_tbl$repair_candidate_score - score_tbl$scar_candidate_score
score_tbl <- score_tbl |> dplyr::mutate( cluster = cluster_vec[match(cell, colnames(obj))], sample = sample_vec[match(cell, colnames(obj))], timepoint = time_vec[match(cell, colnames(obj))], transition_status = as.character(obj$transition_status[match(cell, colnames(obj))]), core_activated_DAM_reference = cluster %in% CORE_ACTIVATED_DAM_CLUSTERS )
obj$repair_candidate_score <- NA_real_
obj$scar_candidate_score   <- NA_real_
obj$repair_minus_scar      <- NA_real_
idx <- match(score_tbl$cell, colnames(obj))
obj$repair_candidate_score[idx] <- score_tbl$repair_candidate_score
obj$scar_candidate_score[idx]   <- score_tbl$scar_candidate_score
obj$repair_minus_scar[idx]      <- score_tbl$repair_minus_scar
component_cols_to_store <- c( repair_z_cols, scar_z_cols )
for (nm in component_cols_to_store) {
  full_vec <- rep(NA_real_, ncol(obj))
  full_vec[idx] <- score_tbl[[nm]]
  obj[[nm]] <- full_vec
}
msg("Calculating sample- and cluster-level summaries...")
sample_summary <- score_tbl |> dplyr::group_by(sample) |> dplyr::summarise( n_cells = dplyr::n(), repair_mean = mean(repair_candidate_score, na.rm = TRUE), repair_median = median(repair_candidate_score, na.rm = TRUE), scar_mean = mean(scar_candidate_score, na.rm = TRUE), scar_median = median(scar_candidate_score, na.rm = TRUE), repair_minus_scar_mean = mean(repair_minus_scar, na.rm = TRUE), .groups = "drop" )
cluster_summary <- score_tbl |> dplyr::group_by(cluster) |> dplyr::summarise( n_cells = dplyr::n(), repair_mean = mean(repair_candidate_score, na.rm = TRUE), repair_median = median(repair_candidate_score, na.rm = TRUE), scar_mean = mean(scar_candidate_score, na.rm = TRUE), scar_median = median(scar_candidate_score, na.rm = TRUE), repair_minus_scar_mean = mean(repair_minus_scar, na.rm = TRUE), .groups = "drop" )
cluster_sample_summary <- score_tbl |> dplyr::group_by(sample, cluster) |> dplyr::summarise( n_cells = dplyr::n(), repair_mean = mean(repair_candidate_score, na.rm = TRUE), repair_median = median(repair_candidate_score, na.rm = TRUE), scar_mean = mean(scar_candidate_score, na.rm = TRUE), scar_median = median(scar_candidate_score, na.rm = TRUE), repair_minus_scar_mean = mean(repair_minus_scar, na.rm = TRUE), .groups = "drop" )
component_cluster_summary <- score_tbl |> dplyr::select( cell, cluster, sample, dplyr::all_of(component_cols_to_store) ) |> tidyr::pivot_longer( cols = dplyr::all_of(component_cols_to_store), names_to = "component", values_to = "z_score" ) |> dplyr::group_by(cluster, component) |> dplyr::summarise( mean_z = mean(z_score, na.rm = TRUE), median_z = median(z_score, na.rm = TRUE), .groups = "drop" )
replicate_consistency <- tibble::tibble()
if (length(unique(cluster_sample_summary$sample)) >= 2) {
  samps <- sort( unique(cluster_sample_summary$sample) )
  sample_pairs <- combn( samps, 2, simplify = FALSE )
  replicate_consistency <- lapply( sample_pairs, function(pair) {
      a <- cluster_sample_summary |> dplyr::filter(sample == pair[1]) |> dplyr::select( cluster, repair_a = repair_median, scar_a = scar_median )
      b <- cluster_sample_summary |> dplyr::filter(sample == pair[2]) |> dplyr::select( cluster, repair_b = repair_median, scar_b = scar_median )
      ab <- dplyr::inner_join( a, b, by = "cluster" )
      tibble::tibble( sample_A = pair[1], sample_B = pair[2], shared_clusters = nrow(ab), repair_cluster_spearman = ifelse( nrow(ab) >= 3, suppressWarnings( stats::cor( ab$repair_a, ab$repair_b, method = "spearman" ) ), NA_real_ ), scar_cluster_spearman = ifelse( nrow(ab) >= 3, suppressWarnings( stats::cor( ab$scar_a, ab$scar_b, method = "spearman" ) ), NA_real_ ) )
    }
  ) |> dplyr::bind_rows()
}
score_cor_all <- suppressWarnings( stats::cor( score_tbl$repair_candidate_score, score_tbl$scar_candidate_score, method = "spearman", use = "complete.obs" ) )
score_cor_by_sample <- score_tbl |> dplyr::group_by(sample) |> dplyr::summarise( n_cells = dplyr::n(), repair_scar_spearman = suppressWarnings( stats::cor( repair_candidate_score, scar_candidate_score, method = "spearman", use = "complete.obs" ) ), .groups = "drop" )
transition_summary <- score_tbl |> dplyr::group_by(transition_status) |> dplyr::summarise( n_cells = dplyr::n(), repair_median = median(repair_candidate_score, na.rm = TRUE), scar_median = median(scar_candidate_score, na.rm = TRUE), repair_minus_scar_median = median(repair_minus_scar, na.rm = TRUE), .groups = "drop" ) |> dplyr::arrange(dplyr::desc(n_cells))
msg("Using formal Step 3A UMAP: ", UMAP_REDUCTION)
umap <- Seurat::Embeddings( obj, reduction = UMAP_REDUCTION )
if (ncol(umap) < 2) {
  stop( paste0( UMAP_REDUCTION, "少于2个维度。" ), call. = FALSE )
}
umap_df <- tibble::tibble( cell = rownames(umap), UMAP_1 = umap[, 1], UMAP_2 = umap[, 2] ) |> dplyr::filter( cell %in% primary_cells ) |> dplyr::left_join( score_tbl |> dplyr::select( cell, cluster, sample, transition_status, repair_candidate_score, scar_candidate_score, repair_minus_scar, core_activated_DAM_reference ), by = "cell" )
missing_cluster_colors <- setdiff( sort(unique(umap_df$cluster)), names(COL_CLUSTER) )
if (length(missing_cluster_colors) > 0) {
  extra_cols <- grDevices::hcl.colors( length(missing_cluster_colors), palette = "Dark 3" )
  names(extra_cols) <- missing_cluster_colors
  COL_CLUSTER <- c(COL_CLUSTER, extra_cols)
}
missing_sample_colors <- setdiff( sort(unique(umap_df$sample)), names(COL_SAMPLE) )
if (length(missing_sample_colors) > 0) {
  extra_sample_cols <- grDevices::hcl.colors( length(missing_sample_colors), palette = "Dark 3" )
  names(extra_sample_cols) <- missing_sample_colors
  COL_SAMPLE <- c(COL_SAMPLE, extra_sample_cols)
}
msg("Creating Fig4A_01_7dpi_program_UMAP.pdf...")
theme_step4 <- ggplot2::theme_classic(base_size = 11) + ggplot2::theme( axis.text = ggplot2::element_text(color = "#222222"), axis.title = ggplot2::element_text(color = "#222222"), plot.title = ggplot2::element_text( face = "bold", color = "#111111" ), plot.subtitle = ggplot2::element_text( color = "#333333" ), legend.title = ggplot2::element_text( face = "bold" ) )
p_cluster <- ggplot2::ggplot( umap_df, ggplot2::aes( x = UMAP_1, y = UMAP_2, color = factor(cluster) ) ) + ggplot2::geom_point( size = 0.45, alpha = 0.75 ) + ggplot2::scale_color_manual( values = COL_CLUSTER, name = "Step 3A cluster" ) + ggplot2::labs( title = "7 dpi primary resident-microglia landscape", subtitle = "Formal Step3A clusters; cycling 6/8 and IFN 9 excluded" ) + theme_step4
p_repair <- ggplot2::ggplot( umap_df, ggplot2::aes( UMAP_1, UMAP_2, color = repair_candidate_score ) ) + ggplot2::geom_point( size = 0.45, alpha = 0.80 ) + ggplot2::scale_color_gradientn( colors = COL_REPAIR, name = "Repair-related\ncandidate score" ) + ggplot2::labs( title = "Repair-related candidate program", subtitle = "Continuous functional score; not a branch assignment" ) + theme_step4
p_scar <- ggplot2::ggplot( umap_df, ggplot2::aes( UMAP_1, UMAP_2, color = scar_candidate_score ) ) + ggplot2::geom_point( size = 0.45, alpha = 0.80 ) + ggplot2::scale_color_gradientn( colors = COL_SCAR, name = "Scar-related\ncandidate score" ) + ggplot2::labs( title = "Scar-related candidate program", subtitle = "Continuous functional score; not a branch assignment" ) + theme_step4
p_diff <- ggplot2::ggplot( umap_df, ggplot2::aes( UMAP_1, UMAP_2, color = repair_minus_scar ) ) + ggplot2::geom_point( size = 0.45, alpha = 0.80 ) + ggplot2::scale_color_gradient2( low = COL_DIFF[1], mid = COL_DIFF[2], high = COL_DIFF[3], midpoint = 0, name = "Repair - scar" ) + ggplot2::labs( title = "Relative functional direction", subtitle = "Visualization only; no cell classification" ) + theme_step4
grDevices::pdf( OUT_FIG1, width = 13, height = 10, useDingbats = FALSE )
print( (p_cluster | p_repair) / (p_scar | p_diff) + patchwork::plot_annotation( title = "Step 4A | Candidate functional programs in 7 dpi microglia" ) )
grDevices::dev.off()
msg("Creating Fig4A_02_program_distribution.pdf...")
score_long <- score_tbl |> dplyr::select( cell, cluster, sample, repair_candidate_score, scar_candidate_score ) |> tidyr::pivot_longer( cols = c( repair_candidate_score, scar_candidate_score ), names_to = "program", values_to = "score" ) |> dplyr::mutate( program = dplyr::recode( program, repair_candidate_score = "Repair-related candidate", scar_candidate_score = "Scar-related candidate" ) )
p_violin_cluster <- ggplot2::ggplot( score_long, ggplot2::aes( x = factor(cluster), y = score, fill = factor(cluster) ) ) + ggplot2::geom_violin( scale = "width", trim = TRUE, linewidth = 0.25, alpha = 0.80 ) + ggplot2::geom_boxplot( width = 0.12, outlier.shape = NA, linewidth = 0.25, fill = "#FFFFFF" ) + ggplot2::facet_wrap( ~program, ncol = 1, scales = "free_y" ) + ggplot2::scale_fill_manual( values = COL_CLUSTER, guide = "none" ) + ggplot2::labs( x = "Formal Step 3A microglia cluster", y = "Continuous candidate-program score", title = "Program-score distribution across 7 dpi Step 3A clusters" ) + theme_step4
sample_score_long <- score_tbl |> dplyr::select( cell, sample, repair_candidate_score, scar_candidate_score ) |> tidyr::pivot_longer( cols = c( repair_candidate_score, scar_candidate_score ), names_to = "program", values_to = "score" ) |> dplyr::mutate( program = dplyr::recode( program, repair_candidate_score = "Repair-related candidate", scar_candidate_score = "Scar-related candidate" ) )
p_violin_sample <- ggplot2::ggplot( sample_score_long, ggplot2::aes( x = factor(sample), y = score, fill = factor(sample) ) ) + ggplot2::geom_violin( scale = "width", trim = TRUE, linewidth = 0.25, alpha = 0.80 ) + ggplot2::geom_boxplot( width = 0.12, outlier.shape = NA, linewidth = 0.25, fill = "#FFFFFF" ) + ggplot2::facet_wrap( ~program, ncol = 1, scales = "free_y" ) + ggplot2::scale_fill_manual( values = COL_SAMPLE, guide = "none" ) + ggplot2::labs( x = "7 dpi biological sample", y = "Continuous candidate-program score", title = "Program-score reproducibility across 7 dpi replicates" ) + theme_step4 + ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 25, hjust = 1 ) )
component_plot_df <- component_cluster_summary |> dplyr::mutate( component = gsub( "^repair_", "Repair | ", component ), component = gsub( "^scar_", "Scar | ", component ), component = gsub( "_z$", "", component ), component = gsub( "_", " ", component ) )
p_component_heat <- ggplot2::ggplot( component_plot_df, ggplot2::aes( x = factor(cluster), y = component, fill = mean_z ) ) + ggplot2::geom_tile( color = "#FFFFFF", linewidth = 0.30 ) + ggplot2::scale_fill_gradient2( low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0, name = "Mean z" ) + ggplot2::labs( x = "Formal Step 3A microglia cluster", y = NULL, title = "Functional-component profile" ) + ggplot2::theme_minimal(base_size = 10) + ggplot2::theme( axis.text = ggplot2::element_text( color = "#222222" ), axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 ), plot.title = ggplot2::element_text( face = "bold", color = "#111111" ), panel.grid = ggplot2::element_blank() )
grDevices::pdf( OUT_FIG2, width = 15, height = 11, useDingbats = FALSE )
print( (p_violin_cluster | p_violin_sample) / p_component_heat + patchwork::plot_layout( heights = c(1.35, 1) ) + patchwork::plot_annotation( title = "Step 4A | Candidate-program distribution and replicate support" ) )
grDevices::dev.off()
msg("Creating Fig4A_03_repair_scar_continuous_map.pdf...")
scatter_df <- score_tbl
if (nrow(scatter_df) > MAX_SCATTER_CELLS) {
  set.seed(20260731)
  n_strata <- nrow( dplyr::distinct( scatter_df, sample, cluster ) )
  n_each <- max( 20, ceiling( MAX_SCATTER_CELLS / max(1, n_strata) ) )
  scatter_df <- scatter_df |> dplyr::mutate(.step4A_random = stats::runif(dplyr::n())) |> dplyr::group_by( sample, cluster ) |> dplyr::arrange( .step4A_random, .by_group = TRUE ) |> dplyr::slice_head( n = n_each ) |> dplyr::ungroup() |> dplyr::select( -.step4A_random )
  if (nrow(scatter_df) > MAX_SCATTER_CELLS) {
    scatter_df <- dplyr::slice_sample( scatter_df, n = MAX_SCATTER_CELLS )
  }
}
p_xy_cluster <- ggplot2::ggplot( scatter_df, ggplot2::aes( x = repair_candidate_score, y = scar_candidate_score, color = factor(cluster) ) ) + ggplot2::geom_hline( yintercept = 0, linewidth = 0.35, linetype = 2, color = "#666666" ) + ggplot2::geom_vline( xintercept = 0, linewidth = 0.35, linetype = 2, color = "#666666" ) + ggplot2::geom_point( size = 0.75, alpha = 0.55 ) + ggplot2::scale_color_manual( values = COL_CLUSTER, name = "Step 3A cluster" ) + ggplot2::labs( x = "Repair-related candidate score", y = "Scar-related candidate score", title = "Two continuous functional coordinates", subtitle = "No quadrant is interpreted as a terminal cell fate" ) + theme_step4
p_xy_sample <- ggplot2::ggplot( scatter_df, ggplot2::aes( x = repair_candidate_score, y = scar_candidate_score, color = factor(sample) ) ) + ggplot2::geom_hline( yintercept = 0, linewidth = 0.35, linetype = 2, color = "#666666" ) + ggplot2::geom_vline( xintercept = 0, linewidth = 0.35, linetype = 2, color = "#666666" ) + ggplot2::geom_point( size = 0.75, alpha = 0.55 ) + ggplot2::scale_color_manual( values = COL_SAMPLE, name = "7 dpi sample" ) + ggplot2::labs( x = "Repair-related candidate score", y = "Scar-related candidate score", title = "Replicate overlay", subtitle = "Replicates should reproduce the same relative functional structure" ) + theme_step4
p_xy_facet <- ggplot2::ggplot( scatter_df, ggplot2::aes( x = repair_candidate_score, y = scar_candidate_score, color = factor(cluster) ) ) + ggplot2::geom_hline( yintercept = 0, linewidth = 0.30, linetype = 2, color = "#777777" ) + ggplot2::geom_vline( xintercept = 0, linewidth = 0.30, linetype = 2, color = "#777777" ) + ggplot2::geom_point( size = 0.60, alpha = 0.50 ) + ggplot2::facet_wrap( ~sample ) + ggplot2::scale_color_manual( values = COL_CLUSTER, guide = "none" ) + ggplot2::labs( x = "Repair-related candidate score", y = "Scar-related candidate score", title = "Within-replicate continuous state maps" ) + theme_step4
grDevices::pdf( OUT_FIG3, width = 14, height = 10, useDingbats = FALSE )
print( (p_xy_cluster | p_xy_sample) / p_xy_facet + patchwork::plot_annotation( title = "Step 4A | Repair-related versus scar-related functional coordinates" ) )
grDevices::dev.off()
msg("Building formal Excel result workbook...")
component_gene_table <- dplyr::bind_rows( lapply( names(repair_components), function(nm) {
      genes_raw <- repair_components[[nm]]$genes
      tibble::tibble( program = "Repair_related_candidate", component = nm, gene = genes_raw, used_in_composite = genes_raw %in% repair_components_score[[nm]], excluded_cross_program_overlap = genes_raw %in% cross_program_overlap, detection_prop_7dpi_primary = det[genes_raw], expression_variance_7dpi_primary = varv[genes_raw] )
    }
  ), lapply(
    names(scar_components), function(nm) {
      genes_raw <- scar_components[[nm]]$genes
      tibble::tibble( program = "Scar_related_candidate", component = nm, gene = genes_raw, used_in_composite = genes_raw %in% scar_components_score[[nm]], excluded_cross_program_overlap = genes_raw %in% cross_program_overlap, detection_prop_7dpi_primary = det[genes_raw], expression_variance_7dpi_primary = varv[genes_raw] )
    }
  )
)
go_term_audit <- dplyr::bind_rows( lapply( names(repair_components), function(nm) {
      tibble::tibble( program = "Repair_related_candidate", component = nm, GO_term = repair_components[[nm]]$terms )
    }
  ), lapply(
    names(scar_components), function(nm) {
      tibble::tibble( program = "Scar_related_candidate", component = nm, GO_term = scar_components[[nm]]$terms )
    }
  )
)
program_summary <- tibble::tibble( metric = c( "Step3D input cells", "All 7 dpi resident microglia", "Primary 7 dpi cells after excluding 6/8/9", "Core activated/DAM reference cells (0/1/4)", "Primary Step3A clusters", "Detected 7 dpi biological samples", "Raw repair union genes", "Raw scar union genes", "Cross-program overlap removed from BOTH scores", "Unique repair composite genes", "Unique scar composite genes", "Overall repair-scar Spearman" ), value = c( as.character(ncol(obj)), as.character(length(all_7dpi_cells)), as.character(length(primary_cells)), as.character(length(core_adam_cells)), paste(observed_primary_clusters, collapse = ";"), paste(seven_samples, collapse = ";"), as.character(length(repair_union_raw)), as.character(length(scar_union_raw)), as.character(length(cross_program_overlap)), as.character(length(repair_union_unique)), as.character(length(scar_union_unique)), format(round(score_cor_all, 4), nsmall = 4) ) )
overlap_table <- tibble::tibble( gene = sort(cross_program_overlap), note = paste0( "Excluded from BOTH composite scores because ", "the gene maps to both functional directions." ) )
cluster_audit <- tibble::tibble( cluster = observed_clusters_all, expected_after_step3B_cleanup = observed_clusters_all %in% EXPECTED_RETAINED_STEP3_CLUSTERS, excluded_from_primary_step4A = observed_clusters_all %in% c(CYCLING_CLUSTERS, IFN_CLUSTER), biological_note = dplyr::case_when( observed_clusters_all == "0" ~ "homeostatic-like/adapted", observed_clusters_all == "1" ~ "DAM-like/phagocytic", observed_clusters_all == "2" ~ "acute inflammatory", observed_clusters_all == "3" ~ "homeostatic", observed_clusters_all == "4" ~ "activated/partial-DAM", observed_clusters_all == "5" ~ "Ccl12/Lcn2 early injury-responsive", observed_clusters_all == "6" ~ "G2/M cycling; primary Step4A excluded", observed_clusters_all == "8" ~ "S phase cycling; primary Step4A excluded", observed_clusters_all == "9" ~ "IFN side-state; primary Step4A excluded", observed_clusters_all == "10" ~ "acute inflammatory/stress", TRUE ~ "unexpected/review" ) )
wb <- openxlsx::createWorkbook()
add_sheet <- function(sheet_name, df) {
  openxlsx::addWorksheet( wb, sheetName = sheet_name )
  openxlsx::writeDataTable( wb, sheet = sheet_name, x = df, tableStyle = "TableStyleMedium2" )
  openxlsx::freezePane( wb, sheet = sheet_name, firstRow = TRUE )
}
add_sheet("Program_summary", program_summary)
add_sheet("Step3_cluster_audit", cluster_audit)
add_sheet("Component_genes", component_gene_table)
add_sheet("Cross_program_overlap", overlap_table)
add_sheet("GO_term_audit", go_term_audit)
add_sheet("Cell_scores", score_tbl)
add_sheet("Sample_summary", sample_summary)
add_sheet("Cluster_summary", cluster_summary)
add_sheet("Cluster_sample_summary", cluster_sample_summary)
add_sheet("Replicate_consistency", replicate_consistency)
add_sheet("Score_correlation_sample", score_cor_by_sample)
add_sheet("Transition_status_summary", transition_summary)
add_sheet("Component_cluster_mean", component_cluster_summary)
openxlsx::saveWorkbook( wb, OUT_XLSX, overwrite = TRUE )
saveRDS( obj, OUT_RDS, compress = TRUE )
msg("Writing Step4A_report.txt...")
report_lines <- c( "======================================================================", "GSE162610 | STEP 4A REPORT", "Repair-related vs scar-related candidate PROGRAM construction", "======================================================================", "", "STEP 4A SCOPE:", "Continuous functional programs only.", "NO repair/scar cell classification.", "NO trajectory bifurcation claim.", "NO Monocle.", "", "----------------------- AUTHORITATIVE STEP 3 INPUT ----------------------", paste0("Input RDS: ", STEP3D_RDS), paste0("Timepoint field: ", TIME_COL), paste0("Biological sample field: ", SAMPLE_COL), paste0("Formal Step3 cluster field: ", CLUSTER_COL), paste0("Formal Step3 UMAP: ", UMAP_REDUCTION), "", "------------------------------ CELL AUDIT ------------------------------", paste0("Step3D object cells: ", ncol(obj)), paste0( "Observed Step3A clusters: ", paste(observed_clusters_all, collapse = ", ") ), paste0( "Unexpected clusters relative to formal Step3B retained set: ", ifelse( length(unexpected_clusters) == 0, "None", paste(unexpected_clusters, collapse = ", ") ) ), paste0("All 7 dpi resident microglia: ", length(all_7dpi_cells)), paste0( "Primary 7 dpi cells after excluding cycling 6/8 and IFN 9: ", length(primary_cells) ), paste0( "Primary Step3A clusters: ", paste(observed_primary_clusters, collapse = ", ") ), paste0( "Core activated/DAM reference cells (clusters 0/1/4): ", length(core_adam_cells) ), paste0( "7 dpi biological samples: ", paste(seven_samples, collapse = ", ") ), "", "--------------------------- PROGRAM AUDIT ------------------------------", paste0("Raw repair gene union: ", length(repair_union_raw)), paste0("Raw scar gene union: ", length(scar_union_raw)), paste0( "Cross-program overlap removed from BOTH composite scores: ", length(cross_program_overlap) ), paste0( "Final unique repair composite genes: ", length(repair_union_unique) ), paste0( "Final unique scar composite genes: ", length(scar_union_unique) ), "", "Composite-score design:", " 1. Programs are built from functional GO-BP directions specified in the", " approved Step 4 project plan.", " 2. Genes must be expressed in the correctly defined primary 7 dpi", " resident-microglia landscape.", " 3. Broad GO components are capped by expression variance.", " 4. Genes assigned to BOTH repair/scar directions are excluded from BOTH", " composite scores.", " 5. Each functional component is independently scored and z-standardized;", " the composite is the mean of component z-scores.", "", "---------------------- CONTINUOUS RELATIONSHIP -------------------------", paste0( "Overall repair/scar Spearman correlation: ", format(round(score_cor_all, 4), nsmall = 4) ), "", "Per-sample repair/scar correlations:", paste( capture.output( print(score_cor_by_sample, n = Inf) ), collapse = "\n" ), "", "Important interpretation:", " - Positive correlation may indicate shared activated/DAM biology.", " - Weak/negative correlation may indicate more independent functional axes.", " - Correlation alone NEVER proves or disproves a fate bifurcation.", " - No quadrant or score threshold is used to label cells in Step 4A.", "", "------------------------- REPLICATE CHECK ------------------------------", if (nrow(replicate_consistency) > 0) {
    paste( capture.output( print(replicate_consistency, n = Inf) ), collapse = "\n" )
  } else {
    "Replicate-consistency correlation could not be calculated."
  },
  "", "------------------------- STEP 4A DECISION -----------------------------", "Judge Step 4A from:", "  1. both candidate programs are measurable;", "  2. their UMAP/cluster distributions show structured gradients;", "  3. broad relative patterns are reproducible in both 7 dpi replicates;", "  4. program behavior remains biologically interpretable without assigning", "     repair/scar cells;", "  5. if strong shared activation remains, refine it in Step 4B before", "     trajectory interpretation.", "", "DO NOT RUN MONOCLE UNTIL STEP 5.", "", "Formal outputs:", paste0("  ", OUT_FIG1), paste0("  ", OUT_FIG2), paste0("  ", OUT_FIG3), paste0("  ", OUT_XLSX), paste0("  ", OUT_RDS), paste0("  ", OUT_REP), paste0("  ", OUT_SESS), ""
)
writeLines( report_lines, OUT_REP, useBytes = TRUE )
writeLines( capture.output( sessionInfo() ), OUT_SESS, useBytes = TRUE )
msg("====================================================")
msg("STEP 4A COMPLETED.")
msg("Output directory: ", OUTDIR)
msg("====================================================")
msg("")
msg("Key outputs:")
msg("  1) Fig4A_01_7dpi_program_UMAP.pdf")
msg("  2) Fig4A_02_program_distribution.pdf")
msg("  3) Fig4A_03_repair_scar_continuous_map.pdf")
msg("  4) Table4A_program_construction_and_scores.xlsx")
msg("  5) Step4A_report.txt")
msg("")
msg("Target-program discovery status: complete")
