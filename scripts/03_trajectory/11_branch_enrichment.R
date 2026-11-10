options(stringsAsFactors = FALSE)
options(timeout = max(600, getOption("timeout", 60)))
set.seed(20260801)
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
STEP6A_DIR <- file.path( PROJECT_ROOT, "04_results", "step6A_GSE162610_7dpi_branch_specific_DEG" )
INPUT_STEP6A <- file.path( STEP6A_DIR, "GSE162610_step6A_7dpi_branch_DEG_results.rds" )
OUTDIR <- file.path( PROJECT_ROOT, "04_results", "step6CR_GSE162610_branch_DEG_GO_KEGG" )
FIGDIR <- file.path(OUTDIR, "figures")
TABDIR <- file.path(OUTDIR, "tables")
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
OUT_FIG1 <- file.path( FIGDIR, "Fig6CR_01_branch_GO_BP_dotplot.pdf" )
OUT_FIG2 <- file.path( FIGDIR, "Fig6CR_02_branch_KEGG_dotplot.pdf" )
OUT_FIG3 <- file.path( FIGDIR, "Fig6CR_03_enrichment_quality_and_term_counts.pdf" )
OUT_FIG4 <- file.path( FIGDIR, "Fig6CR_04_strict_pseudobulk_sensitivity.pdf" )
OUT_XLSX <- file.path( OUTDIR, "Table6CR_branch_DEG_GO_KEGG.xlsx" )
OUT_RDS <- file.path( OUTDIR, "GSE162610_step6CR_branch_DEG_GO_KEGG_results.rds" )
OUT_REP <- file.path( OUTDIR, "Step6CR_report.txt" )
OUT_SESS <- file.path( OUTDIR, "Step6CR_sessionInfo.txt" )
PADJ_MAX <- 0.05
EXPLORATORY_PADJ_MAX <- 0.10
MIN_TERM_COUNT <- 3L
MIN_GENESET_SIZE <- 10L
MAX_GENESET_SIZE <- 500L
TOP_GO_BP_PER_BRANCH <- 12L
TOP_KEGG_PER_BRANCH <- 10L
TOP_SENSITIVITY_PER_BRANCH <- 8L
ARM_COLORS <- c( "Arm_1-biased" = "#2F6B9A", "Arm_2-biased" = "#B64B45" )
required_pkgs <- c( "clusterProfiler", "org.Mm.eg.db", "AnnotationDbi", "dplyr", "tidyr", "tibble", "ggplot2", "patchwork", "stringr", "openxlsx" )
missing_pkgs <- required_pkgs[ !vapply( required_pkgs, requireNamespace, logical(1), quietly = TRUE ) ]
if (length(missing_pkgs) > 0) {
  stop(paste0("Missing R packages: ", paste(missing_pkgs, collapse = ", ")), call. = FALSE)
}
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(stringr)
})
msg <- function(...) {
  cat( sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "" )
}
capture_table <- function(x) {
  paste( capture.output( print( tibble::as_tibble(x), n = Inf, width = Inf ) ), collapse = "\n" )
}
safe_sheet <- function(wb, sheet_name, df) {
  openxlsx::addWorksheet( wb, sheetName = sheet_name )
  ddf <- tryCatch( as.data.frame(df), error = function(e) NULL )
  if ( is.null(ddf) || nrow(ddf) == 0 || ncol(ddf) == 0 ) {
    openxlsx::writeData( wb, sheet = sheet_name, x = data.frame( note = "No result under the predefined Step6C-R criterion.", stringsAsFactors = FALSE ) )
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
ratio_to_numeric <- function(x) {
  vapply( strsplit( as.character(x), "/", fixed = TRUE ), function(z) {
      if (length(z) != 2) {
        return(NA_real_)
      }
      as.numeric(z[1]) / as.numeric(z[2])
    },
    numeric(1)
  )
}
wrap_term <- function(x, width = 52) {
  stringr::str_wrap( x, width = width )
}
null_plot <- function(text) {
  ggplot2::ggplot( tibble::tibble(x = 1, y = 1), ggplot2::aes(x, y) ) + ggplot2::annotate( "text", x = 1, y = 1, label = text, size = 5 ) + ggplot2::theme_void()
}
convert_symbols <- function(symbols) {
  symbols <- unique( as.character(symbols) )
  symbols <- symbols[ !is.na(symbols) & nzchar(symbols) ]
  mapping <- clusterProfiler::bitr( symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db ) |> dplyr::distinct( SYMBOL, ENTREZID )
  list( input_symbols = symbols, mapping = mapping, mapped_symbols = unique(mapping$SYMBOL), mapped_entrez = unique(mapping$ENTREZID), mapping_rate = length(unique(mapping$SYMBOL)) / max(length(symbols), 1) )
}
standardize_result <- function( enrichment_object, branch, evidence_level, database, ontology = NA_character_) {
  if ( is.null(enrichment_object) || nrow(as.data.frame(enrichment_object)) == 0 ) {
    return(tibble::tibble())
  }
  as.data.frame( enrichment_object ) |> tibble::as_tibble() |> dplyr::mutate( branch = branch, evidence_level = evidence_level, database = database, ontology = ontology, GeneRatio_numeric = ratio_to_numeric( GeneRatio ), BgRatio_numeric = ratio_to_numeric( BgRatio ), neglog10_FDR = -log10( pmax( p.adjust, 1e-300 ) ), significant = p.adjust <= PADJ_MAX, exploratory = p.adjust <= EXPLORATORY_PADJ_MAX ) |> dplyr::relocate( branch, evidence_level, database, ontology )
}
run_go_one <- function( gene_entrez, universe_entrez, branch, evidence_level, ontology) {
  if (length(gene_entrez) < MIN_GENESET_SIZE) {
    return(tibble::tibble())
  }
  msg( "Running GO ", ontology, ": ", branch, " | ", evidence_level, " | genes=", length(unique(gene_entrez)) )
  obj <- tryCatch( clusterProfiler::enrichGO( gene = unique(gene_entrez), universe = unique(universe_entrez), OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = ontology, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = MIN_GENESET_SIZE, maxGSSize = MAX_GENESET_SIZE, readable = TRUE ), error = function(e) {
      warning( paste0( "enrichGO failed for ", branch, " / ", evidence_level, " / ", ontology, ": ", conditionMessage(e) ) )
      NULL
    }
  )
  if (is.null(obj)) {
    msg( "GO ", ontology, " returned NULL: ", branch, " | ", evidence_level )
    return(tibble::tibble())
  }
  result <- standardize_result( enrichment_object = obj, branch = branch, evidence_level = evidence_level, database = "GO", ontology = ontology ) |> dplyr::filter( is.finite(p.adjust), Count >= MIN_TERM_COUNT )
  msg( "GO ", ontology, " complete: ", branch, " | ", evidence_level, " | tested terms=", nrow(result), " | FDR<=0.05=", sum(result$significant, na.rm = TRUE) )
  result
}
run_kegg_once <- function( gene_entrez, universe_entrez) {
  clusterProfiler::enrichKEGG( gene = unique(gene_entrez), universe = unique(universe_entrez), organism = "mmu", keyType = "kegg", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = MIN_GENESET_SIZE, maxGSSize = MAX_GENESET_SIZE, use_internal_data = FALSE )
}
run_kegg_one <- function( gene_entrez, universe_entrez, branch, evidence_level, max_attempts = 3L) {
  if (length(gene_entrez) < MIN_GENESET_SIZE) {
    return(tibble::tibble())
  }
  waits <- c(0, 10, 30)
  obj <- NULL
  for (attempt in seq_len(max_attempts)) {
    if (waits[attempt] > 0) {
      msg( "Waiting ", waits[attempt], " seconds before KEGG retry..." )
      Sys.sleep( waits[attempt] )
    }
    msg( "Running KEGG ORA: ", branch, " | ", evidence_level, " | attempt ", attempt, "/", max_attempts )
    ans <- tryCatch( run_kegg_once( gene_entrez = gene_entrez, universe_entrez = universe_entrez ), error = function(e) e )
    if (!inherits(ans, "error")) {
      obj <- ans
      break
    }
    msg( "KEGG attempt failed: ", conditionMessage(ans) )
  }
  if (is.null(obj)) {
    warning( paste0( "KEGG ORA failed after ", max_attempts, " attempts for ", branch, " / ", evidence_level, ". GO results remain valid." ) )
    return(tibble::tibble())
  }
  obj_readable <- tryCatch( clusterProfiler::setReadable( obj, OrgDb = org.Mm.eg.db, keyType = "ENTREZID" ), error = function(e) obj )
  standardize_result( enrichment_object = obj_readable, branch = branch, evidence_level = evidence_level, database = "KEGG", ontology = NA_character_ ) |> dplyr::filter( is.finite(p.adjust), Count >= MIN_TERM_COUNT )
}
top_terms <- function( df, n_per_branch) {
  if (nrow(df) == 0) {
    return(df)
  }
  df |> dplyr::group_by( branch ) |> dplyr::arrange( p.adjust, dplyr::desc(Count), .by_group = TRUE ) |> dplyr::slice_head( n = n_per_branch ) |> dplyr::ungroup()
}
make_branch_dotplot <- function( df, title_text, subtitle_text, n_per_branch) {
  sig <- df |> dplyr::filter( significant )
  display_level <- "FDR <= 0.05"
  if (nrow(sig) == 0) {
    sig <- df |> dplyr::filter( exploratory )
    display_level <- "No FDR <= 0.05 terms; displaying exploratory FDR <= 0.10"
  }
  if (nrow(sig) == 0) {
    return( null_plot( paste0( "No terms passed exploratory FDR <= ", EXPLORATORY_PADJ_MAX, "." ) ) )
  }
  plot_df <- top_terms( sig, n_per_branch ) |> dplyr::mutate( Description_wrapped = wrap_term( Description, 55 ) )
  ggplot2::ggplot( plot_df, ggplot2::aes( x = GeneRatio_numeric, y = stats::reorder( Description_wrapped, GeneRatio_numeric ), size = Count, color = branch ) ) + ggplot2::geom_point( alpha = 0.90 ) + ggplot2::facet_wrap( ~branch, scales = "free_y", ncol = 2 ) + ggplot2::scale_color_manual( values = ARM_COLORS ) + ggplot2::labs( title = title_text, subtitle = paste0( subtitle_text, "\n", display_level ), x = "Gene ratio", y = NULL, size = "Genes", color = NULL ) + ggplot2::theme_classic( base_size = 10 ) + ggplot2::theme( legend.position = "top", strip.background = ggplot2::element_rect( fill = "grey95", color = NA ), strip.text = ggplot2::element_text( face = "bold" ) )
}
msg("====================================================")
msg("GSE162610 Step 6C-R v3 FAST | branch DEG GO / KEGG")
msg("====================================================")
if (!file.exists(INPUT_STEP6A)) {
  stop( paste0( "找不到Step6A结果：", INPUT_STEP6A ), call. = FALSE )
}
step6a <- readRDS( INPUT_STEP6A )
required_fields <- c( "master_evidence", "pseudobulk_deg" )
missing_fields <- required_fields[ !required_fields %in% names(step6a) ]
if (length(missing_fields) > 0) {
  stop( paste0( "Step6A结果缺少：", paste( missing_fields, collapse = ", " ) ), call. = FALSE )
}
master <- as.data.frame( step6a$master_evidence )
pb_deg <- as.data.frame( step6a$pseudobulk_deg )
required_master_cols <- c( "gene_id", "final_direction", "strict_pseudobulk_supported" )
missing_master_cols <- required_master_cols[ !required_master_cols %in% colnames(master) ]
if (length(missing_master_cols) > 0) {
  stop( paste0( "Step6A master_evidence缺少：", paste( missing_master_cols, collapse = ", " ) ), call. = FALSE )
}
if (!"gene_id" %in% colnames(pb_deg)) {
  stop( "Step6A pseudobulk_deg缺少gene_id。", call. = FALSE )
}
direction_levels <- c( "Arm_1-biased", "Arm_2-biased" )
primary_master <- master |> dplyr::filter( final_direction %in% direction_levels ) |> dplyr::distinct( gene_id, .keep_all = TRUE )
strict_master <- primary_master |> dplyr::filter( strict_pseudobulk_supported )
primary_gene_lists <- split( primary_master$gene_id, primary_master$final_direction )
strict_gene_lists <- split( strict_master$gene_id, strict_master$final_direction )
for (branch in direction_levels) {
  if (!branch %in% names(primary_gene_lists)) {
    primary_gene_lists[[branch]] <- character()
  }
  if (!branch %in% names(strict_gene_lists)) {
    strict_gene_lists[[branch]] <- character()
  }
}
background_symbols <- unique( as.character( pb_deg$gene_id ) )
background_symbols <- background_symbols[ !is.na(background_symbols) & nzchar(background_symbols) ]
msg( "Primary high-confidence branch markers: ", nrow(primary_master) )
msg( "Arm_1-biased primary genes: ", length( primary_gene_lists[["Arm_1-biased"]] ) )
msg( "Arm_2-biased primary genes: ", length( primary_gene_lists[["Arm_2-biased"]] ) )
msg( "Strict pseudobulk-supported branch markers: ", nrow(strict_master) )
msg( "Background genes retained by edgeR: ", length(background_symbols) )
background_map <- convert_symbols( background_symbols )
background_entrez <- background_map$mapped_entrez
if (length(background_entrez) < 1000) {
  stop( "Mapped enrichment background contains fewer than 1,000 Entrez IDs.", call. = FALSE )
}
mapping_records <- list( tibble::tibble( branch = "Background", evidence_level = "edgeR_background", input_symbols = length( background_map$input_symbols ), mapped_symbols = length( background_map$mapped_symbols ), mapped_entrez = length( background_map$mapped_entrez ), mapping_rate = background_map$mapping_rate ) )
mapped_gene_sets <- list( Primary = list(), Strict_pseudobulk = list() )
for (evidence_level in names(mapped_gene_sets)) {
  source_lists <- if ( evidence_level == "Primary" ) {
    primary_gene_lists
  } else {
    strict_gene_lists
  }
  for (branch in direction_levels) {
    mp <- convert_symbols( source_lists[[branch]] )
    mapped_gene_sets[[evidence_level]][[branch]] <- mp
    mapping_records[[ paste( evidence_level, branch, sep = "__" ) ]] <- tibble::tibble( branch = branch, evidence_level = evidence_level, input_symbols = length( mp$input_symbols ), mapped_symbols = length( mp$mapped_symbols ), mapped_entrez = length( mp$mapped_entrez ), mapping_rate = mp$mapping_rate )
  }
}
mapping_audit <- dplyr::bind_rows( mapping_records )
print(mapping_audit)
go_results <- list()
kegg_results <- list()
for (evidence_level in names(mapped_gene_sets)) {
  for (branch in direction_levels) {
    entrez <- mapped_gene_sets[[evidence_level]][[branch]]$mapped_entrez
    for (ont in c("BP", "CC", "MF")) {
      key <- paste( evidence_level, branch, ont, sep = "__" )
      go_results[[key]] <- run_go_one( gene_entrez = entrez, universe_entrez = background_entrez, branch = branch, evidence_level = evidence_level, ontology = ont )
    }
    key <- paste( evidence_level, branch, "KEGG", sep = "__" )
    kegg_results[[key]] <- run_kegg_one( gene_entrez = entrez, universe_entrez = background_entrez, branch = branch, evidence_level = evidence_level, max_attempts = 3L )
  }
}
all_go <- dplyr::bind_rows( go_results )
all_kegg <- dplyr::bind_rows( kegg_results )
primary_go <- all_go |> dplyr::filter( evidence_level == "Primary" )
strict_go <- all_go |> dplyr::filter( evidence_level == "Strict_pseudobulk" )
primary_kegg <- all_kegg |> dplyr::filter( evidence_level == "Primary" )
strict_kegg <- all_kegg |> dplyr::filter( evidence_level == "Strict_pseudobulk" )
term_counts <- dplyr::bind_rows( all_go |> dplyr::filter( significant ) |> dplyr::count( evidence_level, branch, ontology, name = "n_significant_terms" ), all_kegg |> dplyr::filter( significant ) |> dplyr::count( evidence_level, branch, name = "n_significant_terms" ) |> dplyr::mutate( ontology = "KEGG" ) )
primary_strict_overlap <- dplyr::bind_rows( all_go |> dplyr::filter( significant ) |> dplyr::select( evidence_level, branch, database, ontology, ID, Description, p.adjust, Count ), all_kegg |> dplyr::filter( significant ) |> dplyr::select( evidence_level, branch, database, ontology, ID, Description, p.adjust, Count ) ) |> dplyr::group_by( branch, database, ontology, ID, Description ) |> dplyr::summarise( found_in_primary = any( evidence_level == "Primary" ), found_in_strict = any( evidence_level == "Strict_pseudobulk" ), best_FDR = min( p.adjust, na.rm = TRUE ), maximum_count = max( Count, na.rm = TRUE ), .groups = "drop" ) |> dplyr::mutate( reproduced_in_strict = found_in_primary & found_in_strict ) |> dplyr::arrange( branch, best_FDR )
p1 <- make_branch_dotplot( df = primary_go |> dplyr::filter( ontology == "BP" ), title_text = "Branch-level GO Biological Process enrichment", subtitle_text = paste0( "Primary high-confidence DEGs; background = ", length(background_symbols), " edgeR-tested genes" ), n_per_branch = TOP_GO_BP_PER_BRANCH )
grDevices::pdf( OUT_FIG1, width = 13, height = 10, useDingbats = FALSE )
print(p1)
grDevices::dev.off()
p2 <- make_branch_dotplot( df = primary_kegg, title_text = "Branch-level KEGG enrichment", subtitle_text = "Primary high-confidence DEGs; mouse KEGG (mmu)", n_per_branch = TOP_KEGG_PER_BRANCH )
grDevices::pdf( OUT_FIG2, width = 13, height = 9, useDingbats = FALSE )
print(p2)
grDevices::dev.off()
mapping_plot_df <- mapping_audit |> dplyr::filter( branch != "Background" ) |> dplyr::mutate( group_label = paste( branch, evidence_level, sep = "\n" ) )
p3a <- ggplot2::ggplot( mapping_plot_df, ggplot2::aes( x = group_label, y = 100 * mapping_rate, fill = branch ) ) + ggplot2::geom_col( width = 0.68 ) + ggplot2::geom_text( ggplot2::aes( label = paste0( round( 100 * mapping_rate, 1 ), "%" ) ), vjust = -0.35, size = 3.2 ) + ggplot2::scale_fill_manual( values = ARM_COLORS ) + ggplot2::coord_cartesian( ylim = c(0, 105) ) + ggplot2::labs( title = "Mouse SYMBOL-to-Entrez mapping", x = NULL, y = "Mapped genes (%)", fill = NULL ) + ggplot2::theme_classic( base_size = 10 ) + ggplot2::theme( legend.position = "top", axis.text.x = ggplot2::element_text( angle = 20, hjust = 1 ) )
if (nrow(term_counts) > 0) {
  p3b <- ggplot2::ggplot( term_counts, ggplot2::aes( x = ontology, y = n_significant_terms, fill = branch ) ) + ggplot2::geom_col( position = ggplot2::position_dodge( width = 0.8 ), width = 0.7 ) + ggplot2::facet_wrap( ~evidence_level, ncol = 1, scales = "free_y" ) + ggplot2::scale_fill_manual( values = ARM_COLORS ) + ggplot2::labs( title = "Significant terms", subtitle = paste0( "FDR <= ", PADJ_MAX, "; term Count >= ", MIN_TERM_COUNT ), x = NULL, y = "Terms", fill = NULL ) + ggplot2::theme_classic( base_size = 10 ) + ggplot2::theme( legend.position = "top" )
} else {
  p3b <- null_plot( "No GO/KEGG terms passed FDR <= 0.05." )
}
grDevices::pdf( OUT_FIG3, width = 13, height = 7, useDingbats = FALSE )
print( p3a | p3b )
grDevices::dev.off()
strict_display <- dplyr::bind_rows( strict_go |> dplyr::filter( ontology == "BP" ), strict_kegg )
if (nrow(strict_display) > 0) {
  strict_sig <- strict_display |> dplyr::filter( significant )
  sensitivity_label <- "FDR <= 0.05"
  if (nrow(strict_sig) == 0) {
    strict_sig <- strict_display |> dplyr::filter( exploratory )
    sensitivity_label <- "No FDR <= 0.05 terms; displaying exploratory FDR <= 0.10"
  }
  if (nrow(strict_sig) > 0) {
    strict_plot_df <- strict_sig |> dplyr::group_by( branch, database ) |> dplyr::arrange( p.adjust, dplyr::desc(Count), .by_group = TRUE ) |> dplyr::slice_head( n = TOP_SENSITIVITY_PER_BRANCH ) |> dplyr::ungroup() |> dplyr::mutate( Description_wrapped = wrap_term( Description, 55 ) )
    p4 <- ggplot2::ggplot( strict_plot_df, ggplot2::aes( x = GeneRatio_numeric, y = stats::reorder( Description_wrapped, GeneRatio_numeric ), size = Count, color = branch ) ) + ggplot2::geom_point( alpha = 0.90 ) + ggplot2::facet_grid( database ~ branch, scales = "free_y", space = "free_y" ) + ggplot2::scale_color_manual( values = ARM_COLORS ) + ggplot2::labs( title = "Strict pseudobulk-supported enrichment sensitivity", subtitle = sensitivity_label, x = "Gene ratio", y = NULL, size = "Genes", color = NULL ) + ggplot2::theme_classic( base_size = 9 ) + ggplot2::theme( legend.position = "top", strip.background = ggplot2::element_rect( fill = "grey95", color = NA ), strip.text = ggplot2::element_text( face = "bold" ) )
  } else {
    p4 <- null_plot( "No strict pseudobulk-supported GO BP/KEGG terms passed FDR <= 0.10." )
  }
} else {
  p4 <- null_plot( "No strict pseudobulk-supported GO BP/KEGG results were returned." )
}
grDevices::pdf( OUT_FIG4, width = 13, height = 11, useDingbats = FALSE )
print(p4)
grDevices::dev.off()
utils::write.csv( all_go, file.path( TABDIR, "01_all_branch_GO_results.csv" ), row.names = FALSE )
utils::write.csv( all_kegg, file.path( TABDIR, "02_all_branch_KEGG_results.csv" ), row.names = FALSE )
utils::write.csv( primary_strict_overlap, file.path( TABDIR, "03_primary_strict_term_overlap.csv" ), row.names = FALSE )
utils::write.csv( mapping_audit, file.path( TABDIR, "04_symbol_entrez_mapping_audit.csv" ), row.names = FALSE )
parameter_tbl <- tibble::tibble( parameter = c( "Species", "GO OrgDb", "KEGG organism", "Primary gene definition", "Strict sensitivity definition", "Background definition", "Background symbols", "Mapped background symbols", "Background mapping rate", "Primary Arm_1 genes", "Primary Arm_2 genes", "Strict Arm_1 genes", "Strict Arm_2 genes", "FDR threshold", "Exploratory FDR threshold", "Minimum term count" ), value = c( "Mus musculus", "org.Mm.eg.db", "mmu", "Step6A final_direction Arm_1-biased / Arm_2-biased", "Primary markers additionally strict_pseudobulk_supported", "All genes retained in Step6A edgeR pseudobulk analysis", length(background_symbols), length(background_map$mapped_symbols), background_map$mapping_rate, length( primary_gene_lists[["Arm_1-biased"]] ), length( primary_gene_lists[["Arm_2-biased"]] ), length( strict_gene_lists[["Arm_1-biased"]] ), length( strict_gene_lists[["Arm_2-biased"]] ), PADJ_MAX, EXPLORATORY_PADJ_MAX, MIN_TERM_COUNT ) )
primary_go_top <- primary_go |> dplyr::filter( significant ) |> dplyr::group_by( branch, ontology ) |> dplyr::arrange( p.adjust, .by_group = TRUE ) |> dplyr::slice_head( n = 20 ) |> dplyr::ungroup()
primary_kegg_top <- primary_kegg |> dplyr::filter( significant ) |> dplyr::group_by( branch ) |> dplyr::arrange( p.adjust, .by_group = TRUE ) |> dplyr::slice_head( n = 20 ) |> dplyr::ungroup()
strict_go_top <- strict_go |> dplyr::filter( significant ) |> dplyr::group_by( branch, ontology ) |> dplyr::arrange( p.adjust, .by_group = TRUE ) |> dplyr::slice_head( n = 20 ) |> dplyr::ungroup()
strict_kegg_top <- strict_kegg |> dplyr::filter( significant ) |> dplyr::group_by( branch ) |> dplyr::arrange( p.adjust, .by_group = TRUE ) |> dplyr::slice_head( n = 20 ) |> dplyr::ungroup()
wb <- openxlsx::createWorkbook()
safe_sheet(wb, "Parameters", parameter_tbl)
safe_sheet(wb, "Mapping_audit", mapping_audit)
safe_sheet(wb, "Term_counts", term_counts)
safe_sheet(wb, "Primary_GO_all", primary_go)
safe_sheet(wb, "Primary_GO_top", primary_go_top)
safe_sheet(wb, "Primary_KEGG_all", primary_kegg)
safe_sheet(wb, "Primary_KEGG_top", primary_kegg_top)
safe_sheet(wb, "Strict_GO_all", strict_go)
safe_sheet(wb, "Strict_GO_top", strict_go_top)
safe_sheet(wb, "Strict_KEGG_all", strict_kegg)
safe_sheet(wb, "Strict_KEGG_top", strict_kegg_top)
safe_sheet(wb, "Primary_strict_overlap", primary_strict_overlap)
safe_sheet(wb, "Primary_markers", primary_master)
safe_sheet(wb, "Strict_markers", strict_master)
openxlsx::saveWorkbook( wb, OUT_XLSX, overwrite = TRUE )
result_object <- list( primary_markers = primary_master, strict_markers = strict_master, primary_gene_lists = primary_gene_lists, strict_gene_lists = strict_gene_lists, background_symbols = background_symbols, background_mapping = background_map$mapping, mapping_audit = mapping_audit, primary_GO = primary_go, primary_KEGG = primary_kegg, strict_GO = strict_go, strict_KEGG = strict_kegg, term_counts = term_counts, primary_strict_overlap = primary_strict_overlap, thresholds = list( FDR = PADJ_MAX, exploratory_FDR = EXPLORATORY_PADJ_MAX, minimum_term_count = MIN_TERM_COUNT, minimum_gene_set_size = MIN_GENESET_SIZE, maximum_gene_set_size = MAX_GENESET_SIZE ) )
saveRDS( result_object, OUT_RDS, compress = TRUE )
report_primary_go <- primary_go |> dplyr::filter( significant, ontology == "BP" ) |> dplyr::group_by( branch ) |> dplyr::arrange( p.adjust, dplyr::desc(Count), .by_group = TRUE ) |> dplyr::slice_head( n = 20 ) |> dplyr::ungroup() |> dplyr::select( branch, ID, Description, GeneRatio, Count, p.adjust, geneID )
report_primary_kegg <- primary_kegg |> dplyr::filter( significant ) |> dplyr::group_by( branch ) |> dplyr::arrange( p.adjust, dplyr::desc(Count), .by_group = TRUE ) |> dplyr::slice_head( n = 20 ) |> dplyr::ungroup() |> dplyr::select( branch, ID, Description, GeneRatio, Count, p.adjust, geneID )
report_strict <- primary_strict_overlap |> dplyr::filter( found_in_strict ) |> dplyr::group_by( branch, database ) |> dplyr::arrange( best_FDR, .by_group = TRUE ) |> dplyr::slice_head( n = 20 ) |> dplyr::ungroup()
report_lines <- c( "======================================================================", "GSE162610 | STEP 6C-R REPORT", "Branch-level GO / KEGG enrichment", "======================================================================", "", "FORMAL POSITION:", "Step6C-R directly compares the functional enrichment of Step6A", "GO semantic-similarity simplify() was intentionally not run because it", "caused the previous Windows R session to stall. Raw GO enrichment is", "retained; redundant terms can be curated after result inspection.", "Arm_1-biased and Arm_2-biased high-confidence branch markers.", "Repair/scar labels were NOT used for gene selection or pathway naming.", "", "PRIMARY GENE SETS:", paste0( "Arm_1-biased genes: ", length( primary_gene_lists[["Arm_1-biased"]] ) ), paste0( "Arm_2-biased genes: ", length( primary_gene_lists[["Arm_2-biased"]] ) ), "", "STRICT PSEUDOBULK SENSITIVITY:", paste0( "Arm_1-biased strict genes: ", length( strict_gene_lists[["Arm_1-biased"]] ) ), paste0( "Arm_2-biased strict genes: ", length( strict_gene_lists[["Arm_2-biased"]] ) ), "", "BACKGROUND:", paste0( "Step6A edgeR-tested genes: ", length(background_symbols) ), paste0( "Mapped background symbols: ", length( background_map$mapped_symbols ), " (", round( 100 * background_map$mapping_rate, 2 ), "%)" ), "The whole genome was not used as the enrichment universe.", "", "MAPPING AUDIT:", capture_table(mapping_audit), "", "SIGNIFICANT TERM COUNTS:", capture_table(term_counts), "", "TOP PRIMARY GO BP TERMS:", capture_table(report_primary_go), "", "TOP PRIMARY KEGG TERMS:", capture_table(report_primary_kegg), "", "TERMS REPRODUCED IN STRICT PSEUDOBULK SENSITIVITY:", capture_table(report_strict), "", "INTERPRETATION RULE:", "Primary interpretation should prioritize:", "  1. coherent pathway families rather than one isolated term;", "  2. terms supported by several genes;", "  3. functional contrast between Arm_1-biased and Arm_2-biased genes;", "  4. reproduction in the strict pseudobulk sensitivity analysis;", "  5. consistency with Step5C-R dynamics and Step6A/6B branch evidence.", "", "STEP6C-R PASS CRITERIA:", "  1. adequate mouse SYMBOL-to-Entrez mapping;", "  2. interpretable GO BP or KEGG terms for at least one branch;", "  3. a coherent directional functional contrast between branches;", "  4. key conclusions not dependent on a single isolated term;", "  5. preferably, major themes reproduced by strict pseudobulk genes.", "", "STATUS:", "After Step6C-R interpretation, run sample/pseudobulk ssGSEA using", "predefined external functional gene sets, then proceed to Step7 pySCENIC.", "", "OUTPUTS:", OUT_FIG1, OUT_FIG2, OUT_FIG3, OUT_FIG4, OUT_XLSX, OUT_RDS, OUT_REP, OUT_SESS )
writeLines( report_lines, OUT_REP, useBytes = TRUE )
writeLines( capture.output( sessionInfo() ), OUT_SESS, useBytes = TRUE )
msg("====================================================")
msg("STEP 6C-R COMPLETE")
msg("====================================================")
msg( "Primary GO significant terms: ", sum( primary_go$significant, na.rm = TRUE ) )
msg( "Primary KEGG significant terms: ", sum( primary_kegg$significant, na.rm = TRUE ) )
msg( "Strict GO significant terms: ", sum( strict_go$significant, na.rm = TRUE ) )
msg( "Strict KEGG significant terms: ", sum( strict_kegg$significant, na.rm = TRUE ) )
msg("")
msg("Key outputs:")
msg("1) figures/Fig6CR_01_branch_GO_BP_dotplot.pdf")
msg("2) figures/Fig6CR_02_branch_KEGG_dotplot.pdf")
msg("3) figures/Fig6CR_03_enrichment_quality_and_term_counts.pdf")
msg("4) figures/Fig6CR_04_strict_pseudobulk_sensitivity.pdf")
msg("5) Step6CR_report.txt")
msg("")
msg("如需精确通路基因，再补充 Table6CR_branch_DEG_GO_KEGG.xlsx。")
