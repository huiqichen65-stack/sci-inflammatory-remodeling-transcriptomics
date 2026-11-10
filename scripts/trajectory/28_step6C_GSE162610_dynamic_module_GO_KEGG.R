# ==============================================================================
# 28_step6C_GSE162610_dynamic_module_GO_KEGG.R
# ==============================================================================
# STEP 6C
# GSE162610 | Dynamic-module GO / KEGG enrichment
#
# INPUT
#   Step6B dynamic modules:
#     - M1/M2/M4/M5: Arm_2-associated terminal modules
#     - M3: shared-prebranch program relatively retained in Arm_1
#
# PRIMARY PRINCIPLES
#   1. Mouse annotation:
#        org.Mm.eg.db
#        KEGG organism = "mmu"
#
#   2. Tested-gene background:
#        The enrichment universe is restricted to the genes that ACTUALLY
#        entered Step6B dynamic clustering, not the whole genome.
#
#   3. Two analysis levels:
#        A) Per dynamic module
#        B) Neutral program-level summaries derived from Step6B:
#             Arm2_core_terminal
#             Arm1_retained_shared
#             Arm2_exploratory
#
#   4. Repair/scar labels are NOT imposed in this script.
#      Functional naming is considered only after inspecting enrichment.
#
#   5. GSEA is NOT performed here; it remains Step6D.
#
# OUTPUT UPLOAD LIMIT
#   The final console asks for only five files.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(20260801)

# ==============================================================================
# 0. PATHS / PARAMETERS
# ==============================================================================

PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
STEP6B_DIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step6B_GSE162610_dynamic_branch_modules"
)

INPUT_STEP6B <- file.path(
  STEP6B_DIR,
  "GSE162610_step6B_dynamic_modules_results.rds"
)

OUTDIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step6C_GSE162610_dynamic_module_GO_KEGG"
)

FIGDIR <- file.path(OUTDIR, "figures")
TABDIR <- file.path(OUTDIR, "tables")

dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)

OUT_FIG1 <- file.path(
  FIGDIR,
  "Fig6C_01_module_GO_BP_dotplot.pdf"
)

OUT_FIG2 <- file.path(
  FIGDIR,
  "Fig6C_02_module_KEGG_dotplot.pdf"
)

OUT_FIG3 <- file.path(
  FIGDIR,
  "Fig6C_03_enrichment_quality_and_term_counts.pdf"
)

OUT_FIG4 <- file.path(
  FIGDIR,
  "Fig6C_04_program_level_functional_summary.pdf"
)

OUT_XLSX <- file.path(
  OUTDIR,
  "Table6C_dynamic_module_GO_KEGG.xlsx"
)

OUT_RDS <- file.path(
  OUTDIR,
  "GSE162610_step6C_GO_KEGG_results.rds"
)

OUT_REP <- file.path(
  OUTDIR,
  "Step6C_report.txt"
)

OUT_SESS <- file.path(
  OUTDIR,
  "Step6C_sessionInfo.txt"
)

PADJ_MAX <- 0.05
MIN_TERM_COUNT <- 3L
MIN_GENESET_SIZE <- 10L
MAX_GENESET_SIZE <- 500L

GO_SIMPLIFY_CUTOFF <- 0.70

TOP_GO_BP_PER_MODULE <- 8L
TOP_KEGG_PER_MODULE <- 8L

TOP_PROGRAM_GO <- 10L
TOP_PROGRAM_KEGG <- 8L

# ==============================================================================
# 1. PACKAGES
# ==============================================================================

required_pkgs <- c(
  "clusterProfiler",
  "org.Mm.eg.db",
  "AnnotationDbi",
  "enrichplot",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "patchwork",
  "stringr",
  "openxlsx"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "缺少以下R包：",
      paste(missing_pkgs, collapse = ", "),
      "\n可使用BiocManager::install()安装Bioconductor包。"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(stringr)
})

# ==============================================================================
# 2. HELPERS
# ==============================================================================

msg <- function(...) {
  cat(
    sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
    ...,
    "\n",
    sep = ""
  )
}

capture_table <- function(x) {
  paste(
    capture.output(
      print(
        tibble::as_tibble(x),
        n = Inf,
        width = Inf
      )
    ),
    collapse = "\n"
  )
}

safe_sheet <- function(wb, sheet_name, df) {

  openxlsx::addWorksheet(
    wb,
    sheetName = sheet_name
  )

  ddf <- tryCatch(
    as.data.frame(df),
    error = function(e) NULL
  )

  if (
    is.null(ddf) ||
    nrow(ddf) == 0 ||
    ncol(ddf) == 0
  ) {

    openxlsx::writeData(
      wb,
      sheet = sheet_name,
      x = data.frame(
        note = "No result under the predefined Step6C criterion.",
        stringsAsFactors = FALSE
      )
    )

  } else {

    openxlsx::writeDataTable(
      wb,
      sheet = sheet_name,
      x = ddf,
      tableStyle = "TableStyleMedium2"
    )
  }

  openxlsx::freezePane(
    wb,
    sheet = sheet_name,
    firstRow = TRUE
  )
}

ratio_to_numeric <- function(x) {
  vapply(
    strsplit(as.character(x), "/", fixed = TRUE),
    function(z) {
      if (length(z) != 2) return(NA_real_)
      as.numeric(z[1]) / as.numeric(z[2])
    },
    numeric(1)
  )
}

wrap_term <- function(x, width = 48) {
  stringr::str_wrap(
    x,
    width = width
  )
}

convert_symbols <- function(symbols) {

  symbols <- unique(
    as.character(symbols)
  )

  symbols <- symbols[
    !is.na(symbols) &
      nzchar(symbols)
  ]

  mapped <- clusterProfiler::bitr(
    symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )

  mapped <- mapped |>
    dplyr::distinct(
      SYMBOL,
      ENTREZID
    )

  list(
    input_symbols = symbols,
    mapping = mapped,
    mapped_symbols = unique(mapped$SYMBOL),
    mapped_entrez = unique(mapped$ENTREZID),
    mapping_rate =
      length(unique(mapped$SYMBOL)) /
      max(length(symbols), 1)
  )
}

standardize_enrich_df <- function(
    enrich_object,
    group_name,
    analysis_level,
    database,
    ontology = NA_character_) {

  if (
    is.null(enrich_object) ||
    nrow(as.data.frame(enrich_object)) == 0
  ) {
    return(tibble::tibble())
  }

  x <- as.data.frame(enrich_object)

  if (nrow(x) == 0) {
    return(tibble::tibble())
  }

  x |>
    tibble::as_tibble() |>
    dplyr::mutate(
      group_name = group_name,
      analysis_level = analysis_level,
      database = database,
      ontology = ontology,
      GeneRatio_numeric =
        ratio_to_numeric(GeneRatio),
      BgRatio_numeric =
        ratio_to_numeric(BgRatio),
      neglog10_FDR =
        -log10(
          pmax(
            p.adjust,
            1e-300
          )
        )
    ) |>
    dplyr::relocate(
      group_name,
      analysis_level,
      database,
      ontology
    )
}

run_go_one <- function(
    gene_entrez,
    universe_entrez,
    group_name,
    analysis_level,
    ontology) {

  if (length(gene_entrez) < MIN_GENESET_SIZE) {
    return(tibble::tibble())
  }

  obj <- tryCatch(
    clusterProfiler::enrichGO(
      gene = unique(gene_entrez),
      universe = unique(universe_entrez),
      OrgDb = org.Mm.eg.db,
      keyType = "ENTREZID",
      ont = ontology,
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      minGSSize = MIN_GENESET_SIZE,
      maxGSSize = MAX_GENESET_SIZE,
      readable = TRUE
    ),
    error = function(e) {
      warning(
        paste0(
          "enrichGO failed for ",
          group_name,
          " / ",
          ontology,
          ": ",
          conditionMessage(e)
        )
      )
      NULL
    }
  )

  if (is.null(obj)) {
    return(tibble::tibble())
  }

  # Reduce highly redundant GO terms when possible.
  obj_simple <- tryCatch(
    clusterProfiler::simplify(
      obj,
      cutoff = GO_SIMPLIFY_CUTOFF,
      by = "p.adjust",
      select_fun = min,
      measure = "Wang"
    ),
    error = function(e) obj
  )

  standardize_enrich_df(
    enrich_object = obj_simple,
    group_name = group_name,
    analysis_level = analysis_level,
    database = "GO",
    ontology = ontology
  ) |>
    dplyr::filter(
      is.finite(p.adjust),
      p.adjust <= PADJ_MAX,
      Count >= MIN_TERM_COUNT
    )
}

run_kegg_one <- function(
    gene_entrez,
    universe_entrez,
    group_name,
    analysis_level) {

  if (length(gene_entrez) < MIN_GENESET_SIZE) {
    return(tibble::tibble())
  }

  obj <- tryCatch(
    clusterProfiler::enrichKEGG(
      gene = unique(gene_entrez),
      universe = unique(universe_entrez),
      organism = "mmu",
      keyType = "kegg",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      minGSSize = MIN_GENESET_SIZE,
      maxGSSize = MAX_GENESET_SIZE,
      use_internal_data = FALSE
    ),
    error = function(e) {
      warning(
        paste0(
          "enrichKEGG failed for ",
          group_name,
          ": ",
          conditionMessage(e)
        )
      )
      NULL
    }
  )

  if (is.null(obj)) {
    return(tibble::tibble())
  }

  # Convert KEGG Entrez IDs back to mouse symbols when possible.
  obj_readable <- tryCatch(
    clusterProfiler::setReadable(
      obj,
      OrgDb = org.Mm.eg.db,
      keyType = "ENTREZID"
    ),
    error = function(e) obj
  )

  standardize_enrich_df(
    enrich_object = obj_readable,
    group_name = group_name,
    analysis_level = analysis_level,
    database = "KEGG",
    ontology = NA_character_
  ) |>
    dplyr::filter(
      is.finite(p.adjust),
      p.adjust <= PADJ_MAX,
      Count >= MIN_TERM_COUNT
    )
}

top_terms <- function(
    df,
    n_per_group) {

  if (nrow(df) == 0) {
    return(df)
  }

  df |>
    dplyr::group_by(
      group_name
    ) |>
    dplyr::arrange(
      p.adjust,
      dplyr::desc(Count),
      .by_group = TRUE
    ) |>
    dplyr::slice_head(
      n = n_per_group
    ) |>
    dplyr::ungroup()
}

null_plot <- function(txt) {
  ggplot2::ggplot(
    tibble::tibble(x = 1, y = 1),
    ggplot2::aes(x, y)
  ) +
    ggplot2::annotate(
      "text",
      x = 1,
      y = 1,
      label = txt,
      size = 5
    ) +
    ggplot2::theme_void()
}

# ==============================================================================
# 3. LOAD STEP6B
# ==============================================================================

msg("====================================================")
msg("GSE162610 Step 6C | dynamic-module GO / KEGG")
msg("====================================================")

if (!file.exists(INPUT_STEP6B)) {
  stop(
    paste0(
      "找不到Step6B结果：",
      INPUT_STEP6B
    ),
    call. = FALSE
  )
}

step6b <- readRDS(
  INPUT_STEP6B
)

required_fields <- c(
  "genes_entering_clustering",
  "gene_modules",
  "module_annotation",
  "selected_k"
)

missing_fields <- required_fields[
  !required_fields %in% names(step6b)
]

if (length(missing_fields) > 0) {
  stop(
    paste0(
      "Step6B结果缺少：",
      paste(
        missing_fields,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

gene_modules <- as.data.frame(
  step6b$gene_modules
)

module_annotation <- as.data.frame(
  step6b$module_annotation
)

required_gene_cols <- c(
  "gene_id",
  "module_name"
)

missing_gene_cols <- required_gene_cols[
  !required_gene_cols %in% colnames(gene_modules)
]

if (length(missing_gene_cols) > 0) {
  stop(
    paste0(
      "Step6B gene_modules缺少：",
      paste(
        missing_gene_cols,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

msg("Selected dynamic-module K: ", step6b$selected_k)
msg(
  "Genes entering Step6B clustering: ",
  length(step6b$genes_entering_clustering)
)

# ==============================================================================
# 4. MOUSE SYMBOL -> ENTREZ MAPPING / TESTED BACKGROUND
# ==============================================================================

background_symbols <- unique(
  as.character(
    step6b$genes_entering_clustering
  )
)

background_map <- convert_symbols(
  background_symbols
)

background_entrez <- background_map$mapped_entrez

msg(
  "Background mouse symbol mapping: ",
  length(background_map$mapped_symbols),
  "/",
  length(background_symbols),
  " (",
  round(
    100 * background_map$mapping_rate,
    1
  ),
  "%)"
)

if (length(background_entrez) < 100) {
  stop(
    "Mapped enrichment background contains fewer than 100 Entrez IDs.",
    call. = FALSE
  )
}

mapping_table <- background_map$mapping

# ==============================================================================
# 5. PER-MODULE GENE LISTS / MAPPING AUDIT
# ==============================================================================

module_gene_lists <- split(
  gene_modules$gene_id,
  gene_modules$module_name
)

module_gene_lists <- lapply(
  module_gene_lists,
  unique
)

module_mapping_list <- lapply(
  names(module_gene_lists),
  function(mod) {

    mp <- convert_symbols(
      module_gene_lists[[mod]]
    )

    tibble::tibble(
      group_name = mod,
      analysis_level = "Module",
      input_symbols =
        length(mp$input_symbols),
      mapped_symbols =
        length(mp$mapped_symbols),
      mapped_entrez =
        length(mp$mapped_entrez),
      mapping_rate =
        mp$mapping_rate
    )
  }
)

module_mapping_audit <- dplyr::bind_rows(
  module_mapping_list
)

print(module_mapping_audit)

# ==============================================================================
# 6. DERIVE NEUTRAL PROGRAM-LEVEL GENE SETS
# ==============================================================================

required_annotation_cols <- c(
  "module_name",
  "temporal_label",
  "mean_Arm1_minus_Arm2",
  "same_direction_all_samples",
  "robust_terminal_branch_effect"
)

missing_annotation_cols <- required_annotation_cols[
  !required_annotation_cols %in%
    colnames(module_annotation)
]

if (length(missing_annotation_cols) > 0) {
  stop(
    paste0(
      "Step6B module_annotation缺少：",
      paste(
        missing_annotation_cols,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

robust_arm2_modules <- module_annotation |>
  dplyr::filter(
    temporal_label == "Arm2_terminal",
    robust_terminal_branch_effect
  ) |>
  dplyr::pull(
    module_name
  )

arm1_retained_modules <- module_annotation |>
  dplyr::filter(
    temporal_label == "Shared_prebranch",
    same_direction_all_samples,
    mean_Arm1_minus_Arm2 > 0
  ) |>
  dplyr::pull(
    module_name
  )

exploratory_arm2_modules <- module_annotation |>
  dplyr::filter(
    temporal_label == "Arm2_terminal",
    !robust_terminal_branch_effect
  ) |>
  dplyr::pull(
    module_name
  )

program_gene_lists <- list()

if (length(robust_arm2_modules) > 0) {
  program_gene_lists$Arm2_core_terminal <- unique(
    gene_modules$gene_id[
      gene_modules$module_name %in%
        robust_arm2_modules
    ]
  )
}

if (length(arm1_retained_modules) > 0) {
  program_gene_lists$Arm1_retained_shared <- unique(
    gene_modules$gene_id[
      gene_modules$module_name %in%
        arm1_retained_modules
    ]
  )
}

if (length(exploratory_arm2_modules) > 0) {
  program_gene_lists$Arm2_exploratory <- unique(
    gene_modules$gene_id[
      gene_modules$module_name %in%
        exploratory_arm2_modules
    ]
  )
}

program_definition <- tibble::tibble(
  program_name = c(
    "Arm2_core_terminal",
    "Arm1_retained_shared",
    "Arm2_exploratory"
  ),
  source_modules = c(
    paste(
      robust_arm2_modules,
      collapse = "; "
    ),
    paste(
      arm1_retained_modules,
      collapse = "; "
    ),
    paste(
      exploratory_arm2_modules,
      collapse = "; "
    )
  )
) |>
  dplyr::filter(
    nzchar(source_modules)
  )

msg("Program-level definitions:")
print(program_definition)

program_mapping_audit <- dplyr::bind_rows(
  lapply(
    names(program_gene_lists),
    function(prg) {

      mp <- convert_symbols(
        program_gene_lists[[prg]]
      )

      tibble::tibble(
        group_name = prg,
        analysis_level = "Program",
        input_symbols =
          length(mp$input_symbols),
        mapped_symbols =
          length(mp$mapped_symbols),
        mapped_entrez =
          length(mp$mapped_entrez),
        mapping_rate =
          mp$mapping_rate
      )
    }
  )
)

mapping_audit <- dplyr::bind_rows(
  module_mapping_audit,
  program_mapping_audit
)

# ==============================================================================
# 7. RUN PER-MODULE GO / KEGG
# ==============================================================================

msg("Running per-module GO BP/CC/MF and KEGG...")

module_go_results <- list()
module_kegg_results <- list()

for (mod in names(module_gene_lists)) {

  mp <- convert_symbols(
    module_gene_lists[[mod]]
  )

  for (ont in c("BP", "CC", "MF")) {

    key <- paste(
      mod,
      ont,
      sep = "__"
    )

    module_go_results[[key]] <- run_go_one(
      gene_entrez = mp$mapped_entrez,
      universe_entrez = background_entrez,
      group_name = mod,
      analysis_level = "Module",
      ontology = ont
    )
  }

  module_kegg_results[[mod]] <- run_kegg_one(
    gene_entrez = mp$mapped_entrez,
    universe_entrez = background_entrez,
    group_name = mod,
    analysis_level = "Module"
  )
}

module_go <- dplyr::bind_rows(
  module_go_results
)

module_kegg <- dplyr::bind_rows(
  module_kegg_results
)

# ==============================================================================
# 8. RUN PROGRAM-LEVEL GO / KEGG
# ==============================================================================

msg("Running neutral program-level GO BP and KEGG...")

program_go_results <- list()
program_kegg_results <- list()

for (prg in names(program_gene_lists)) {

  mp <- convert_symbols(
    program_gene_lists[[prg]]
  )

  program_go_results[[prg]] <- run_go_one(
    gene_entrez = mp$mapped_entrez,
    universe_entrez = background_entrez,
    group_name = prg,
    analysis_level = "Program",
    ontology = "BP"
  )

  program_kegg_results[[prg]] <- run_kegg_one(
    gene_entrez = mp$mapped_entrez,
    universe_entrez = background_entrez,
    group_name = prg,
    analysis_level = "Program"
  )
}

program_go <- dplyr::bind_rows(
  program_go_results
)

program_kegg <- dplyr::bind_rows(
  program_kegg_results
)

# ==============================================================================
# 9. TERM COUNTS / QUALITY TABLES
# ==============================================================================

module_go_counts <- if (nrow(module_go) > 0) {
  module_go |>
    dplyr::count(
      group_name,
      ontology,
      name = "n_significant_terms"
    )
} else {
  tibble::tibble()
}

module_kegg_counts <- if (nrow(module_kegg) > 0) {
  module_kegg |>
    dplyr::count(
      group_name,
      name = "n_significant_terms"
    ) |>
    dplyr::mutate(
      ontology = "KEGG"
    )
} else {
  tibble::tibble()
}

term_counts <- dplyr::bind_rows(
  module_go_counts,
  module_kegg_counts
)

# ==============================================================================
# 10. FIGURE 1 | MODULE GO BP DOTPLOT
# ==============================================================================

module_go_bp <- module_go |>
  dplyr::filter(
    ontology == "BP"
  )

top_module_go_bp <- top_terms(
  module_go_bp,
  TOP_GO_BP_PER_MODULE
)

if (nrow(top_module_go_bp) > 0) {

  top_module_go_bp <- top_module_go_bp |>
    dplyr::mutate(
      Description_wrapped =
        wrap_term(
          Description,
          52
        )
    )

  p1 <- ggplot2::ggplot(
    top_module_go_bp,
    ggplot2::aes(
      x = GeneRatio_numeric,
      y = stats::reorder(
        Description_wrapped,
        GeneRatio_numeric
      ),
      size = Count,
      color = neglog10_FDR
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.85
    ) +
    ggplot2::facet_wrap(
      ~group_name,
      scales = "free_y",
      ncol = 2
    ) +
    ggplot2::labs(
      title = "GO Biological Process enrichment by dynamic module",
      subtitle = "Background = 803 genes entering Step6B clustering; GO terms simplified for redundancy",
      x = "Gene ratio",
      y = NULL,
      size = "Genes",
      color = "-log10(FDR)"
    ) +
    ggplot2::theme_classic(
      base_size = 9
    )

} else {

  p1 <- null_plot(
    "No module-level GO BP terms passed FDR <= 0.05 and Count >= 3."
  )
}

grDevices::pdf(
  OUT_FIG1,
  width = 13,
  height = 12,
  useDingbats = FALSE
)

print(p1)
grDevices::dev.off()

# ==============================================================================
# 11. FIGURE 2 | MODULE KEGG DOTPLOT
# ==============================================================================

top_module_kegg <- top_terms(
  module_kegg,
  TOP_KEGG_PER_MODULE
)

if (nrow(top_module_kegg) > 0) {

  top_module_kegg <- top_module_kegg |>
    dplyr::mutate(
      Description_wrapped =
        wrap_term(
          Description,
          52
        )
    )

  p2 <- ggplot2::ggplot(
    top_module_kegg,
    ggplot2::aes(
      x = GeneRatio_numeric,
      y = stats::reorder(
        Description_wrapped,
        GeneRatio_numeric
      ),
      size = Count,
      color = neglog10_FDR
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.85
    ) +
    ggplot2::facet_wrap(
      ~group_name,
      scales = "free_y",
      ncol = 2
    ) +
    ggplot2::labs(
      title = "KEGG enrichment by dynamic module",
      subtitle = "Mouse KEGG (mmu); tested-gene background from Step6B",
      x = "Gene ratio",
      y = NULL,
      size = "Genes",
      color = "-log10(FDR)"
    ) +
    ggplot2::theme_classic(
      base_size = 9
    )

} else {

  p2 <- null_plot(
    paste0(
      "No module-level KEGG terms passed the predefined threshold.\n",
      "If warnings indicate KEGG connection failure, GO results remain valid."
    )
  )
}

grDevices::pdf(
  OUT_FIG2,
  width = 13,
  height = 11,
  useDingbats = FALSE
)

print(p2)
grDevices::dev.off()

# ==============================================================================
# 12. FIGURE 3 | QUALITY / TERM COUNTS
# ==============================================================================

p3a <- ggplot2::ggplot(
  mapping_audit,
  ggplot2::aes(
    x = group_name,
    y = 100 * mapping_rate,
    fill = analysis_level
  )
) +
  ggplot2::geom_col(
    width = 0.7
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = paste0(
        round(
          100 * mapping_rate,
          1
        ),
        "%"
      )
    ),
    vjust = -0.35,
    size = 3
  ) +
  ggplot2::coord_cartesian(
    ylim = c(
      0,
      105
    )
  ) +
  ggplot2::labs(
    title = "Mouse SYMBOL-to-Entrez mapping quality",
    x = NULL,
    y = "Mapped symbols (%)",
    fill = "Level"
  ) +
  ggplot2::theme_classic(
    base_size = 10
  ) +
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 30,
        hjust = 1
      )
  )

if (nrow(term_counts) > 0) {

  p3b <- ggplot2::ggplot(
    term_counts,
    ggplot2::aes(
      x = group_name,
      y = n_significant_terms,
      fill = ontology
    )
  ) +
    ggplot2::geom_col(
      position = "dodge"
    ) +
    ggplot2::labs(
      title = "Significant enrichment terms by module",
      subtitle = paste0(
        "FDR <= ",
        PADJ_MAX,
        "; Count >= ",
        MIN_TERM_COUNT
      ),
      x = NULL,
      y = "Terms",
      fill = "Database / ontology"
    ) +
    ggplot2::theme_classic(
      base_size = 10
    ) +
    ggplot2::theme(
      axis.text.x =
        ggplot2::element_text(
          angle = 30,
          hjust = 1
        )
    )

} else {

  p3b <- null_plot(
    "No significant enrichment terms were returned."
  )
}

grDevices::pdf(
  OUT_FIG3,
  width = 13,
  height = 6.5,
  useDingbats = FALSE
)

print(p3a | p3b)
grDevices::dev.off()

# ==============================================================================
# 13. FIGURE 4 | PROGRAM-LEVEL FUNCTIONAL SUMMARY
# ==============================================================================

top_program_go <- top_terms(
  program_go,
  TOP_PROGRAM_GO
)

top_program_kegg <- top_terms(
  program_kegg,
  TOP_PROGRAM_KEGG
)

if (nrow(top_program_go) > 0) {

  top_program_go <- top_program_go |>
    dplyr::mutate(
      Description_wrapped =
        wrap_term(
          Description,
          52
        )
    )

  p4a <- ggplot2::ggplot(
    top_program_go,
    ggplot2::aes(
      x = GeneRatio_numeric,
      y = stats::reorder(
        Description_wrapped,
        GeneRatio_numeric
      ),
      size = Count,
      color = neglog10_FDR
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.85
    ) +
    ggplot2::facet_wrap(
      ~group_name,
      scales = "free_y",
      ncol = 1
    ) +
    ggplot2::labs(
      title = "Program-level GO Biological Process",
      x = "Gene ratio",
      y = NULL,
      size = "Genes",
      color = "-log10(FDR)"
    ) +
    ggplot2::theme_classic(
      base_size = 9
    )

} else {

  p4a <- null_plot(
    "No program-level GO BP terms passed the predefined threshold."
  )
}

if (nrow(top_program_kegg) > 0) {

  top_program_kegg <- top_program_kegg |>
    dplyr::mutate(
      Description_wrapped =
        wrap_term(
          Description,
          52
        )
    )

  p4b <- ggplot2::ggplot(
    top_program_kegg,
    ggplot2::aes(
      x = GeneRatio_numeric,
      y = stats::reorder(
        Description_wrapped,
        GeneRatio_numeric
      ),
      size = Count,
      color = neglog10_FDR
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.85
    ) +
    ggplot2::facet_wrap(
      ~group_name,
      scales = "free_y",
      ncol = 1
    ) +
    ggplot2::labs(
      title = "Program-level KEGG",
      x = "Gene ratio",
      y = NULL,
      size = "Genes",
      color = "-log10(FDR)"
    ) +
    ggplot2::theme_classic(
      base_size = 9
    )

} else {

  p4b <- null_plot(
    "No program-level KEGG terms passed the predefined threshold."
  )
}

grDevices::pdf(
  OUT_FIG4,
  width = 13,
  height = 15,
  useDingbats = FALSE
)

print(
  p4a /
    p4b +
    patchwork::plot_annotation(
      title = "Neutral program-level functional summary"
    )
)

grDevices::dev.off()

# ==============================================================================
# 14. SAVE TABLES
# ==============================================================================

utils::write.csv(
  module_go,
  file.path(
    TABDIR,
    "01_module_GO_all_significant.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  module_kegg,
  file.path(
    TABDIR,
    "02_module_KEGG_all_significant.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  program_go,
  file.path(
    TABDIR,
    "03_program_GO_BP_significant.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  program_kegg,
  file.path(
    TABDIR,
    "04_program_KEGG_significant.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  mapping_audit,
  file.path(
    TABDIR,
    "05_symbol_entrez_mapping_audit.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 15. EXCEL
# ==============================================================================

parameter_tbl <- tibble::tibble(
  parameter = c(
    "Species",
    "OrgDb",
    "KEGG organism",
    "Enrichment background symbols",
    "Mapped background symbols",
    "Background mapping rate",
    "FDR threshold",
    "Minimum term gene count",
    "Minimum gene-set size",
    "Maximum gene-set size",
    "GO redundancy cutoff",
    "Step6B selected K",
    "Robust Arm2 modules",
    "Arm1-retained shared modules",
    "Exploratory Arm2 modules",
    "Module GO significant terms",
    "Module KEGG significant terms",
    "Program GO significant terms",
    "Program KEGG significant terms"
  ),
  value = c(
    "Mus musculus",
    "org.Mm.eg.db",
    "mmu",
    length(background_symbols),
    length(background_map$mapped_symbols),
    background_map$mapping_rate,
    PADJ_MAX,
    MIN_TERM_COUNT,
    MIN_GENESET_SIZE,
    MAX_GENESET_SIZE,
    GO_SIMPLIFY_CUTOFF,
    step6b$selected_k,
    paste(
      robust_arm2_modules,
      collapse = "; "
    ),
    paste(
      arm1_retained_modules,
      collapse = "; "
    ),
    paste(
      exploratory_arm2_modules,
      collapse = "; "
    ),
    nrow(module_go),
    nrow(module_kegg),
    nrow(program_go),
    nrow(program_kegg)
  )
)

wb <- openxlsx::createWorkbook()

safe_sheet(wb, "Parameters", parameter_tbl)
safe_sheet(wb, "Program_definition", program_definition)
safe_sheet(wb, "Mapping_audit", mapping_audit)
safe_sheet(wb, "Module_GO", module_go)
safe_sheet(wb, "Module_GO_BP_top", top_module_go_bp)
safe_sheet(wb, "Module_KEGG", module_kegg)
safe_sheet(wb, "Module_KEGG_top", top_module_kegg)
safe_sheet(wb, "Program_GO_BP", program_go)
safe_sheet(wb, "Program_KEGG", program_kegg)
safe_sheet(wb, "Term_counts", term_counts)
safe_sheet(wb, "Background_mapping", mapping_table)

openxlsx::saveWorkbook(
  wb,
  OUT_XLSX,
  overwrite = TRUE
)

# ==============================================================================
# 16. SAVE RDS
# ==============================================================================

result_object <- list(
  species = "Mus musculus",
  background_symbols = background_symbols,
  background_mapping = mapping_table,
  mapping_audit = mapping_audit,
  module_gene_lists = module_gene_lists,
  program_definition = program_definition,
  program_gene_lists = program_gene_lists,
  module_GO = module_go,
  module_KEGG = module_kegg,
  program_GO_BP = program_go,
  program_KEGG = program_kegg,
  term_counts = term_counts,
  thresholds = list(
    FDR = PADJ_MAX,
    minimum_term_count = MIN_TERM_COUNT,
    minimum_gene_set_size = MIN_GENESET_SIZE,
    maximum_gene_set_size = MAX_GENESET_SIZE,
    GO_simplify_cutoff = GO_SIMPLIFY_CUTOFF
  )
)

saveRDS(
  result_object,
  OUT_RDS,
  compress = TRUE
)

# ==============================================================================
# 17. REPORT
# ==============================================================================

top_module_go_report <- top_module_go_bp |>
  dplyr::select(
    group_name,
    ID,
    Description,
    GeneRatio,
    Count,
    p.adjust,
    geneID
  )

top_module_kegg_report <- top_module_kegg |>
  dplyr::select(
    group_name,
    ID,
    Description,
    GeneRatio,
    Count,
    p.adjust,
    geneID
  )

top_program_go_report <- top_program_go |>
  dplyr::select(
    group_name,
    ID,
    Description,
    GeneRatio,
    Count,
    p.adjust,
    geneID
  )

top_program_kegg_report <- top_program_kegg |>
  dplyr::select(
    group_name,
    ID,
    Description,
    GeneRatio,
    Count,
    p.adjust,
    geneID
  )

report_lines <- c(
  "======================================================================",
  "GSE162610 | STEP 6C REPORT",
  "Dynamic-module GO / KEGG enrichment",
  "======================================================================",
  "",
  "FORMAL POSITION:",
  "Step6C performs functional enrichment after Step6B dynamic-module discovery.",
  "Repair/scar labels were NOT used to define gene sets or enrichment input.",
  "GSEA remains Step6D.",
  "",
  "SPECIES / DATABASE:",
  "Species: Mus musculus",
  "GO OrgDb: org.Mm.eg.db",
  "KEGG organism: mmu",
  "",
  "TESTED-GENE BACKGROUND:",
  paste0(
    "Step6B genes entering clustering: ",
    length(background_symbols)
  ),
  paste0(
    "Mapped background symbols: ",
    length(background_map$mapped_symbols),
    " (",
    round(
      100 * background_map$mapping_rate,
      2
    ),
    "%)"
  ),
  "The enrichment universe was restricted to Step6B-tested genes rather than",
  "the whole genome.",
  "",
  "PROGRAM DEFINITIONS:",
  capture_table(program_definition),
  "",
  "MAPPING AUDIT:",
  capture_table(mapping_audit),
  "",
  "SIGNIFICANT TERM COUNTS:",
  capture_table(term_counts),
  "",
  "TOP MODULE GO BP TERMS:",
  capture_table(top_module_go_report),
  "",
  "TOP MODULE KEGG TERMS:",
  capture_table(top_module_kegg_report),
  "",
  "TOP PROGRAM-LEVEL GO BP TERMS:",
  capture_table(top_program_go_report),
  "",
  "TOP PROGRAM-LEVEL KEGG TERMS:",
  capture_table(top_program_kegg_report),
  "",
  "INTERPRETATION RULE:",
  "Primary interpretation should prioritize:",
  "  1. coherent GO/KEGG themes repeated across robust Arm2 modules;",
  "  2. distinct functions of the Arm1-retained shared/prebranch program;",
  "  3. terms supported by multiple genes rather than isolated single pathways;",
  "  4. biological consistency with Step5C-R and Step6A without circular",
  "     selection from repair/scar candidate lists.",
  "",
  "M4/M5 are exploratory because Step6B replicate effect size was weaker.",
  "Their pathways should not outweigh the robust M1/M2/M3 evidence.",
  "",
  "STEP6C PASS CRITERIA:",
  "  1. adequate mouse SYMBOL-to-Entrez mapping;",
  "  2. significant, interpretable terms for at least the robust M1/M2/M3 or",
  "     their corresponding neutral program-level gene sets;",
  "  3. functional distinction between Arm2 core terminal and Arm1-retained",
  "     shared programs;",
  "  4. no conclusion based only on one isolated term.",
  "",
  "NEXT:",
  "After Step6C interpretation, proceed to Step6D ranked-list GSEA.",
  "",
  "OUTPUTS:",
  OUT_FIG1,
  OUT_FIG2,
  OUT_FIG3,
  OUT_FIG4,
  OUT_XLSX,
  OUT_RDS,
  OUT_REP,
  OUT_SESS
)

writeLines(
  report_lines,
  OUT_REP,
  useBytes = TRUE
)

writeLines(
  capture.output(sessionInfo()),
  OUT_SESS,
  useBytes = TRUE
)

# ==============================================================================
# 18. FINAL CONSOLE SUMMARY
# ==============================================================================

msg("====================================================")
msg("STEP 6C COMPLETE")
msg("====================================================")
msg(
  "Background mapping rate: ",
  round(
    100 * background_map$mapping_rate,
    1
  ),
  "%"
)
msg(
  "Module GO significant terms: ",
  nrow(module_go)
)
msg(
  "Module KEGG significant terms: ",
  nrow(module_kegg)
)
msg(
  "Program GO significant terms: ",
  nrow(program_go)
)
msg(
  "Program KEGG significant terms: ",
  nrow(program_kegg)
)

msg("")
msg("请上传以下5个结果用于Step6C PASS/FAIL判断：")
msg("1) figures/Fig6C_01_module_GO_BP_dotplot.pdf")
msg("2) figures/Fig6C_02_module_KEGG_dotplot.pdf")
msg("3) figures/Fig6C_03_enrichment_quality_and_term_counts.pdf")
msg("4) figures/Fig6C_04_program_level_functional_summary.pdf")
msg("5) Step6C_report.txt")
msg("")
msg("如需精确通路基因，再补充 Table6C_dynamic_module_GO_KEGG.xlsx。")
