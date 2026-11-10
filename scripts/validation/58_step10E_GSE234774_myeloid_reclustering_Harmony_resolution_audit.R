# =============================================================================
# 58_step10E_GSE234774_myeloid_reclustering_Harmony_resolution_audit_v2.R
#
# Adult SCI activated/DAM-like microglia project
# GSE234774 Step 10E:
#   Myeloid-only reanalysis, conservative-Harmony sensitivity, and
#   higher-resolution clustering audit
#
# INPUT
#   <PROJECT_ROOT>/04_results/
#     step10D_GSE234774_marker_first_major_lineage_annotation/RDS/
#       GSE234774_step10D_major_lineage_annotation_evidence_BPCells.rds
#
# FROZEN PARENT-ATLAS SELECTION
#   Parent whole-atlas cluster 1:
#     activated microglia-dominant myeloid compartment
#   Parent whole-atlas cluster 21:
#     cycling myeloid compartment
#
# PURPOSE
#   1. Extract parent clusters 1 and 21 without adding other major lineages.
#   2. Rebuild a stable myeloid-only BPCells/Seurat object from raw counts.
#   3. Re-run LogNormalize and myeloid-only variable-feature selection.
#   4. Score cell cycle but do NOT regress it.
#   5. Construct an identity-focused PCA excluding canonical cell-cycle,
#      mitochondrial, ribosomal, and hemoglobin genes from PCA features.
#   6. Compare unintegrated myeloid PCA against conservative Harmony.
#   7. Sweep Harmony clustering resolutions 0.2, 0.4, 0.6, 0.8, and 1.0.
#   8. Export UMAP, marker-panel, mixing, cluster-size, dominance, stability,
#      and audit evidence for manual resolution selection.
#
# ANALYSIS BOUNDARIES
#   - No final myeloid subtype label is frozen.
#   - No cell is removed after extraction.
#   - No cell-cycle regression is run.
#   - Harmony changes low-dimensional embeddings only.
#   - RNA counts and LogNormalize expression remain uncorrected.
#   - No GSE162610 Arm_1/Arm_2 score is calculated or redefined.
#   - No differential expression across time, pseudobulk, trajectory,
#     regulon, CellChat, or spatial analysis is run.
#
# REQUIRED REVIEW OUTPUTS
# ROBUSTNESS FIX IN V2
#   - Uses type-safe dimension/name comparisons for the normalized BPCells layer.
#   - Automatically rebuilds a stale/incomplete/mismatched normalized directory.
#
# REQUIRED REVIEW OUTPUTS
#   Upload exactly seven files:
#     3 figures + 2 tables + 1 audit + 1 report.
# =============================================================================

options(
  stringsAsFactors = FALSE,
  warn = 1,
  Seurat.object.assay.version = "v5",
  future.globals.maxSize = 16 * 1024^3
)

# -----------------------------------------------------------------------------
# 0. User controls
# -----------------------------------------------------------------------------
SEED <- 20260803L

PARENT_CLUSTER_COLUMN <- "hcon_res_0_2"
PARENT_CLUSTERS_KEEP <- c("1", "21")
EXPECTED_MYLOID_NUCLEI <- 25205L

N_VARIABLE_FEATURES <- 3000L
N_PCS_COMPUTE <- 40L
PCS_USE <- 1:30

NORMALIZATION_SCALE_FACTOR <- 10000
NORMALIZATION_BLOCK_SIZE <- 5000L

K_GRAPH <- 30L
K_MIXING <- 30L
N_TREES <- 50L
SNN_PRUNE <- 1 / 15

RESOLUTIONS <- c(
  0.2,
  0.4,
  0.6,
  0.8,
  1.0
)

DISPLAY_RESOLUTION <- 0.6

UMAP_N_NEIGHBORS <- 30L
UMAP_MIN_DIST <- 0.25

# Conservative Harmony settings frozen from the whole-atlas sensitivity step.
HARMONY_THETA <- 1
HARMONY_LAMBDA <- 2
HARMONY_MAX_ITER <- 20L
HARMONY_NCORES <- 4L

# Rebuild controls. Leave FALSE during normal use.
FORCE_REWRITE_MYLOID_COUNTS <- FALSE
FORCE_REWRITE_MYLOID_DATA <- FALSE

# Review flags only.
LIBRARY_DOMINANCE_THRESHOLD <- 0.80
TIMEPOINT_DOMINANCE_THRESHOLD <- 0.95
PARENT_CLUSTER_DOMINANCE_THRESHOLD <- 0.90
CYCLING_PHASE_DOMINANCE_THRESHOLD <- 0.70
MIN_CLUSTER_SIZE_FOR_DOMINANCE <- 200L
MIN_ACCEPTABLE_MEAN_KNN_OVERLAP <- 0.30

expected_timepoints <- c(
  "Uninjured",
  "1 dpi",
  "4 dpi",
  "7 dpi",
  "14 dpi",
  "1 month",
  "2 months"
)

branch_specs <- data.frame(
  branch = c(
    "myeloid_pca_identity",
    "myeloid_harmony_conservative"
  ),
  reduction = c(
    "pca_myeloid_identity",
    "harmony_myeloid_conservative"
  ),
  graph_prefix = c(
    "myeloid_pca",
    "myeloid_hcon"
  ),
  cluster_prefix = c(
    "mypca",
    "myhcon"
  ),
  umap = c(
    "umap_myeloid_pca",
    "umap_myeloid_harmony"
  ),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# 1. Curated marker panels for resolution review
# -----------------------------------------------------------------------------
marker_panels_requested <- list(
  Myeloid_general = c(
    "Ptprc",
    "Aif1",
    "Tyrobp",
    "Lyz2",
    "Csf1r"
  ),
  Resident_homeostatic = c(
    "P2ry12",
    "Tmem119",
    "Sall1",
    "Fcrls",
    "Hexb",
    "Cx3cr1",
    "Gpr34",
    "Olfml3"
  ),
  Monocyte = c(
    "Ccr2",
    "Ly6c2",
    "Plac8",
    "S100a8",
    "S100a9",
    "Retnlg",
    "Fcgr3"
  ),
  Macrophage_activated = c(
    "Lgals3",
    "Ms4a7",
    "Cd68",
    "Apoe",
    "Cst7",
    "Ctsb",
    "Ctsd",
    "Trem2",
    "Itgax"
  ),
  BAM = c(
    "Mrc1",
    "Lyve1",
    "Pf4",
    "Folr2",
    "Cd163",
    "Cbr2"
  ),
  Dendritic = c(
    "Flt3",
    "Zbtb46",
    "Clec10a",
    "Cd209a",
    "Cd74",
    "H2-Ab1"
  ),
  Interferon_response = c(
    "Isg15",
    "Ifit1",
    "Ifit3",
    "Stat1",
    "Rsad2"
  ),
  Inflammatory_chemokine = c(
    "Il1b",
    "Tnf",
    "Nfkbia",
    "Ccl2",
    "Ccl3",
    "Ccl4"
  ),
  Cell_cycle = c(
    "Mki67",
    "Top2a",
    "Stmn1",
    "Tubb5",
    "Cenpf",
    "Ube2c"
  )
)

# -----------------------------------------------------------------------------
# 2. Paths
# -----------------------------------------------------------------------------
project_root <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
input_object_path <- file.path(
  project_root,
  "04_results",
  "step10D_GSE234774_marker_first_major_lineage_annotation",
  "RDS",
  "GSE234774_step10D_major_lineage_annotation_evidence_BPCells.rds"
)

output_dir <- file.path(
  project_root,
  "04_results",
  "step10E_GSE234774_myeloid_reclustering_Harmony_resolution"
)

tables_dir <- file.path(
  output_dir,
  "tables"
)

figures_dir <- file.path(
  output_dir,
  "figures"
)

rds_dir <- file.path(
  output_dir,
  "RDS"
)

bpcells_dir <- file.path(
  output_dir,
  "BPCells"
)

for (directory in c(
  output_dir,
  tables_dir,
  figures_dir,
  rds_dir,
  bpcells_dir
)) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

run_log <- file.path(
  output_dir,
  "00_run_log.txt"
)

if (file.exists(
  run_log
)) {
  file.remove(
    run_log
  )
}

myeloid_counts_dir <- file.path(
  bpcells_dir,
  "GSE234774_step10E_myeloid_counts"
)

myeloid_data_dir <- file.path(
  bpcells_dir,
  "GSE234774_step10E_myeloid_LogNormalize_data"
)

output_object_path <- file.path(
  rds_dir,
  "GSE234774_step10E_myeloid_PCA_Harmony_resolution_comparison_BPCells.rds"
)

# -----------------------------------------------------------------------------
# 3. Helpers
# -----------------------------------------------------------------------------
log_message <- function(...) {
  message_text <- paste0(...)
  timestamp <- format(
    Sys.time(),
    "%Y-%m-%d %H:%M:%S"
  )
  line <- paste0(
    "[",
    timestamp,
    "] ",
    message_text
  )

  cat(
    line,
    "\n"
  )

  cat(
    line,
    "\n",
    file = run_log,
    append = TRUE
  )
}

write_csv_safe <- function(
    x,
    filename
) {
  path <- file.path(
    tables_dir,
    filename
  )

  utils::write.csv(
    x,
    path,
    row.names = FALSE,
    na = ""
  )

  log_message(
    "Wrote: ",
    path
  )

  invisible(
    path
  )
}

stop_with_report <- function(reason) {
  report_path <- file.path(
    output_dir,
    "10_STEP10E_FINAL_REPORT.txt"
  )

  report <- c(
    "GSE234774 STEP 10E MYELOID RECLUSTERING / HARMONY / RESOLUTION",
    "============================================================",
    paste0(
      "Run time: ",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      )
    ),
    paste0(
      "Output directory: ",
      output_dir
    ),
    "",
    "FINAL DECISION: STOP",
    "",
    paste0(
      "Reason: ",
      reason
    ),
    "",
    "Do not freeze a myeloid resolution or resident-microglia object."
  )

  writeLines(
    report,
    report_path,
    useBytes = TRUE
  )

  stop(
    reason,
    call. = FALSE
  )
}

safe_package_version <- function(
    package_name
) {
  if (!requireNamespace(
    package_name,
    quietly = TRUE
  )) {
    return(
      NA_character_
    )
  }

  as.character(
    utils::packageVersion(
      package_name
    )
  )
}

open_bpcells_safe <- function(path) {
  tryCatch(
    BPCells::open_matrix_dir(
      path
    ),
    error = function(error_object) {
      error_object
    }
  )
}

resolution_tag <- function(resolution) {
  gsub(
    "\\.",
    "_",
    format(
      resolution,
      scientific = FALSE,
      trim = TRUE
    )
  )
}

cluster_column_name <- function(
    prefix,
    resolution
) {
  paste0(
    prefix,
    "_res_",
    resolution_tag(
      resolution
    )
  )
}

resolve_features <- function(
    requested_features,
    available_features
) {
  available_upper <- toupper(
    available_features
  )

  first_index <- !duplicated(
    available_upper
  )

  lookup <- setNames(
    available_features[
      first_index
    ],
    available_upper[
      first_index
    ]
  )

  resolved <- unname(
    lookup[
      toupper(
        requested_features
      )
    ]
  )

  names(
    resolved
  ) <- requested_features

  resolved
}

deduplicate_plot_panels <- function(
    panel_list
) {
  used <- character(0)
  output <- list()

  for (panel_name in names(
    panel_list
  )) {
    genes <- panel_list[[panel_name]]
    genes <- genes[
      !genes %in%
        used
    ]

    if (length(
      genes
    ) > 0L) {
      output[[panel_name]] <- genes
      used <- c(
        used,
        genes
      )
    }
  }

  output
}

adjusted_rand_index <- function(
    x,
    y
) {
  x <- as.character(
    x
  )

  y <- as.character(
    y
  )

  valid <- !is.na(
    x
  ) &
    !is.na(
      y
    )

  x <- x[
    valid
  ]

  y <- y[
    valid
  ]

  contingency <- table(
    x,
    y
  )

  n <- sum(
    contingency
  )

  if (n < 2L) {
    return(
      NA_real_
    )
  }

  choose_two <- function(z) {
    z * (
      z - 1
    ) / 2
  }

  sum_cells <- sum(
    choose_two(
      contingency
    )
  )

  sum_rows <- sum(
    choose_two(
      rowSums(
        contingency
      )
    )
  )

  sum_columns <- sum(
    choose_two(
      colSums(
        contingency
      )
    )
  )

  total_pairs <- choose_two(
    n
  )

  expected_index <- (
    sum_rows *
      sum_columns
  ) /
    total_pairs

  maximum_index <- 0.5 *
    (
      sum_rows +
        sum_columns
    )

  denominator <- maximum_index -
    expected_index

  if (abs(
    denominator
  ) < .Machine$double.eps) {
    return(
      1
    )
  }

  (
    sum_cells -
      expected_index
  ) /
    denominator
}

clean_neighbor_indices <- function(
    index_matrix,
    k
) {
  n_rows <- nrow(
    index_matrix
  )

  cleaned <- matrix(
    NA_integer_,
    nrow = n_rows,
    ncol = k
  )

  for (row_index in seq_len(
    n_rows
  )) {
    values <- index_matrix[
      row_index,
    ]

    values <- values[
      !is.na(
        values
      ) &
        values !=
          row_index
    ]

    take <- min(
      k,
      length(
        values
      )
    )

    if (take > 0L) {
      cleaned[
        row_index,
        seq_len(
          take
        )
      ] <- values[
        seq_len(
          take
        )
      ]
    }
  }

  cleaned
}

calculate_within_timepoint_mixing <- function(
    embeddings,
    metadata,
    reduction_name,
    k = 30L
) {
  summary_rows <- list()
  neighbor_cache <- list()

  for (timepoint_name in expected_timepoints) {
    cells <- rownames(
      metadata
    )[
      as.character(
        metadata$timepoint
      ) ==
        timepoint_name
    ]

    sub_metadata <- metadata[
      cells,
      ,
      drop = FALSE
    ]

    sub_embeddings <- embeddings[
      cells,
      ,
      drop = FALSE
    ]

    libraries <- as.character(
      sub_metadata$biological_library
    )

    library_levels <- unique(
      libraries
    )

    if (
      length(
        library_levels
      ) < 2L ||
        nrow(
          sub_embeddings
        ) <=
        k + 1L
    ) {
      summary_rows[[timepoint_name]] <- data.frame(
        reduction = reduction_name,
        timepoint = timepoint_name,
        n_cells = nrow(
          sub_embeddings
        ),
        n_libraries = length(
          library_levels
        ),
        mean_other_library_fraction =
          NA_real_,
        median_other_library_fraction =
          NA_real_,
        mean_normalized_library_entropy =
          NA_real_,
        median_normalized_library_entropy =
          NA_real_,
        fraction_with_other_library_neighbor =
          NA_real_,
        stringsAsFactors = FALSE
      )

      neighbor_cache[[timepoint_name]] <- NULL
      next
    }

    log_message(
      "Finding within-timepoint neighbors: ",
      reduction_name,
      " / ",
      timepoint_name,
      "."
    )

    neighbor_object <- Seurat::FindNeighbors(
      object = sub_embeddings,
      k.param = min(
        k + 1L,
        nrow(
          sub_embeddings
        ) - 1L
      ),
      return.neighbor = TRUE,
      compute.SNN = FALSE,
      nn.method = "annoy",
      n.trees = N_TREES,
      annoy.metric = "euclidean",
      verbose = FALSE
    )

    raw_indices <- SeuratObject::Indices(
      neighbor_object
    )

    clean_indices <- clean_neighbor_indices(
      raw_indices,
      k = k
    )

    neighbor_libraries <- matrix(
      libraries[
        clean_indices
      ],
      nrow = nrow(
        clean_indices
      ),
      ncol = ncol(
        clean_indices
      )
    )

    valid_neighbor_counts <- rowSums(
      !is.na(
        neighbor_libraries
      )
    )

    same_library_matrix <- sweep(
      neighbor_libraries,
      MARGIN = 1L,
      STATS = libraries,
      FUN = "=="
    )

    same_library_fraction <- rowSums(
      same_library_matrix,
      na.rm = TRUE
    ) /
      valid_neighbor_counts

    other_library_fraction <-
      1 -
      same_library_fraction

    has_other_library_neighbor <- rowSums(
      !same_library_matrix &
        !is.na(
          same_library_matrix
        ),
      na.rm = TRUE
    ) > 0L

    probability_matrix <- vapply(
      library_levels,
      function(library_name) {
        rowSums(
          neighbor_libraries ==
            library_name,
          na.rm = TRUE
        ) /
          valid_neighbor_counts
      },
      numeric(
        nrow(
          neighbor_libraries
        )
      )
    )

    if (is.null(
      dim(
        probability_matrix
      )
    )) {
      probability_matrix <- matrix(
        probability_matrix,
        ncol = 1L
      )
    }

    entropy_terms <- ifelse(
      probability_matrix > 0,
      probability_matrix *
        log(
          probability_matrix
        ),
      0
    )

    normalized_entropy <- -rowSums(
      entropy_terms
    ) /
      log(
        length(
          library_levels
        )
      )

    summary_rows[[timepoint_name]] <- data.frame(
      reduction = reduction_name,
      timepoint = timepoint_name,
      n_cells = length(
        cells
      ),
      n_libraries = length(
        library_levels
      ),
      mean_other_library_fraction =
        mean(
          other_library_fraction,
          na.rm = TRUE
        ),
      median_other_library_fraction =
        median(
          other_library_fraction,
          na.rm = TRUE
        ),
      mean_normalized_library_entropy =
        mean(
          normalized_entropy,
          na.rm = TRUE
        ),
      median_normalized_library_entropy =
        median(
          normalized_entropy,
          na.rm = TRUE
        ),
      fraction_with_other_library_neighbor =
        mean(
          has_other_library_neighbor,
          na.rm = TRUE
        ),
      stringsAsFactors = FALSE
    )

    neighbor_cache[[timepoint_name]] <- list(
      cells = cells,
      indices = clean_indices
    )

    rm(
      neighbor_object,
      raw_indices,
      clean_indices,
      neighbor_libraries,
      probability_matrix
    )

    invisible(
      gc()
    )
  }

  list(
    summaries = do.call(
      rbind,
      summary_rows
    ),
    neighbor_cache = neighbor_cache
  )
}

calculate_neighbor_overlap <- function(
    reference_cache,
    comparison_cache
) {
  rows <- list()

  for (timepoint_name in expected_timepoints) {
    reference_entry <-
      reference_cache[[timepoint_name]]

    comparison_entry <-
      comparison_cache[[timepoint_name]]

    if (
      is.null(
        reference_entry
      ) ||
        is.null(
          comparison_entry
        )
    ) {
      rows[[timepoint_name]] <- data.frame(
        timepoint = timepoint_name,
        n_cells = NA_integer_,
        mean_KNN_overlap = NA_real_,
        median_KNN_overlap = NA_real_,
        q10_KNN_overlap = NA_real_,
        q90_KNN_overlap = NA_real_,
        stringsAsFactors = FALSE
      )

      next
    }

    if (!identical(
      reference_entry$cells,
      comparison_entry$cells
    )) {
      stop(
        paste0(
          "Cell order differs in neighbor caches at ",
          timepoint_name,
          "."
        )
      )
    }

    reference_indices <-
      reference_entry$indices

    comparison_indices <-
      comparison_entry$indices

    overlap_fraction <- vapply(
      seq_len(
        nrow(
          reference_indices
        )
      ),
      function(row_index) {
        reference_values <-
          reference_indices[
            row_index,
          ]

        comparison_values <-
          comparison_indices[
            row_index,
          ]

        reference_values <-
          reference_values[
            !is.na(
              reference_values
            )
          ]

        comparison_values <-
          comparison_values[
            !is.na(
              comparison_values
            )
          ]

        denominator <- min(
          length(
            reference_values
          ),
          length(
            comparison_values
          )
        )

        if (denominator == 0L) {
          return(
            NA_real_
          )
        }

        length(
          intersect(
            reference_values,
            comparison_values
          )
        ) /
          denominator
      },
      numeric(1)
    )

    rows[[timepoint_name]] <- data.frame(
      timepoint = timepoint_name,
      n_cells = length(
        overlap_fraction
      ),
      mean_KNN_overlap = mean(
        overlap_fraction,
        na.rm = TRUE
      ),
      median_KNN_overlap = median(
        overlap_fraction,
        na.rm = TRUE
      ),
      q10_KNN_overlap = as.numeric(
        stats::quantile(
          overlap_fraction,
          0.10,
          na.rm = TRUE
        )
      ),
      q90_KNN_overlap = as.numeric(
        stats::quantile(
          overlap_fraction,
          0.90,
          na.rm = TRUE
        )
      ),
      stringsAsFactors = FALSE
    )
  }

  do.call(
    rbind,
    rows
  )
}

run_harmony_matrix_compat <- function(
    pca_embeddings,
    metadata
) {
  harmony_version <- utils::packageVersion(
    "harmony"
  )

  log_message(
    "Running conservative myeloid Harmony with theta=",
    HARMONY_THETA,
    ", lambda=",
    HARMONY_LAMBDA,
    ", harmony version ",
    harmony_version,
    "."
  )

  set.seed(
    SEED
  )

  if (
    harmony_version >=
      "2.0.0"
  ) {
    output <- harmony::RunHarmony(
      data_mat = pca_embeddings,
      meta_data = metadata,
      vars_use =
        "biological_library",
      theta = HARMONY_THETA,
      lambda = HARMONY_LAMBDA,
      sigma = 0.1,
      nclust = NULL,
      max_iter =
        HARMONY_MAX_ITER,
      early_stop = TRUE,
      ncores =
        HARMONY_NCORES,
      plot_convergence = FALSE,
      return_object = FALSE,
      verbose = TRUE
    )
  } else if (
    exists(
      "HarmonyMatrix",
      where = asNamespace(
        "harmony"
      ),
      inherits = FALSE
    )
  ) {
    output <- harmony::HarmonyMatrix(
      data_mat = pca_embeddings,
      meta_data = metadata,
      vars_use =
        "biological_library",
      do_pca = FALSE,
      theta = HARMONY_THETA,
      lambda = HARMONY_LAMBDA,
      sigma = 0.1,
      nclust = NULL,
      tau = 0,
      block.size = 0.05,
      max.iter.harmony =
        HARMONY_MAX_ITER,
      plot_convergence = FALSE,
      return_object = FALSE,
      verbose = TRUE
    )
  } else {
    stop(
      paste0(
        "Unsupported harmony version: ",
        harmony_version,
        "."
      )
    )
  }

  output <- as.matrix(
    output
  )

  if (
    nrow(
      output
    ) !=
      nrow(
        pca_embeddings
      )
  ) {
    if (
      ncol(
        output
      ) ==
        nrow(
          pca_embeddings
        )
    ) {
      output <- t(
        output
      )
    } else {
      stop(
        "Harmony returned unexpected dimensions."
      )
    }
  }

  rownames(
    output
  ) <- rownames(
    pca_embeddings
  )

  output
}

summarize_clustering <- function(
    metadata,
    cluster_column,
    branch,
    resolution
) {
  cluster_values <- as.character(
    metadata[
      ,
      cluster_column
    ]
  )

  cluster_data <- data.table::data.table(
    cluster = cluster_values,
    biological_library =
      as.character(
        metadata$biological_library
      ),
    timepoint =
      as.character(
        metadata$timepoint
      ),
    parent_atlas_cluster =
      as.character(
        metadata$parent_atlas_cluster
      ),
    Phase =
      as.character(
        metadata$Phase
      )
  )

  sizes <- cluster_data[
    ,
    .(
      n_cells = .N,
      n_libraries = uniqueN(
        biological_library
      ),
      n_timepoints = uniqueN(
        timepoint
      ),
      n_parent_clusters = uniqueN(
        parent_atlas_cluster
      )
    ),
    by = cluster
  ]

  dominant_fraction <- function(
      value_column,
      output_value_name,
      output_fraction_name
  ) {
    counts <- cluster_data[
      ,
      .(
        n_value_cells = .N
      ),
      by = c(
        "cluster",
        value_column
      )
    ]

    counts <- merge(
      counts,
      sizes[
        ,
        c(
          "cluster",
          "n_cells"
        ),
        with = FALSE
      ],
      by = "cluster",
      all.x = TRUE,
      sort = FALSE
    )

    counts$fraction <-
      counts$n_value_cells /
      counts$n_cells

    top <- counts[
      order(
        cluster,
        -fraction
      ),
      .SD[
        1L
      ],
      by = cluster
    ]

    output <- data.frame(
      cluster = top$cluster,
      dominant_value = top[[value_column]],
      dominant_fraction =
        top$fraction,
      stringsAsFactors = FALSE
    )

    names(
      output
    )[
      names(
        output
      ) ==
        "dominant_value"
    ] <- output_value_name

    names(
      output
    )[
      names(
        output
      ) ==
        "dominant_fraction"
    ] <- output_fraction_name

    output
  }

  library_top <- dominant_fraction(
    "biological_library",
    "dominant_library",
    "maximum_library_fraction"
  )

  timepoint_top <- dominant_fraction(
    "timepoint",
    "dominant_timepoint",
    "maximum_timepoint_fraction"
  )

  parent_top <- dominant_fraction(
    "parent_atlas_cluster",
    "dominant_parent_cluster",
    "maximum_parent_cluster_fraction"
  )

  phase_top <- dominant_fraction(
    "Phase",
    "dominant_phase",
    "maximum_phase_fraction"
  )

  detail <- Reduce(
    function(x, y) {
      merge(
        x,
        y,
        by = "cluster",
        all.x = TRUE,
        sort = FALSE
      )
    },
    list(
      as.data.frame(
        sizes
      ),
      library_top,
      timepoint_top,
      parent_top,
      phase_top
    )
  )

  detail$branch <- branch
  detail$resolution <- resolution
  detail$cluster_column <-
    cluster_column

  detail$library_dominance_review <-
    detail$n_cells >=
      MIN_CLUSTER_SIZE_FOR_DOMINANCE &
    detail$maximum_library_fraction >=
      LIBRARY_DOMINANCE_THRESHOLD

  detail$timepoint_dominance_review <-
    detail$n_cells >=
      MIN_CLUSTER_SIZE_FOR_DOMINANCE &
    detail$maximum_timepoint_fraction >=
      TIMEPOINT_DOMINANCE_THRESHOLD

  detail$parent_cluster_dominance_review <-
    detail$n_cells >=
      MIN_CLUSTER_SIZE_FOR_DOMINANCE &
    detail$maximum_parent_cluster_fraction >=
      PARENT_CLUSTER_DOMINANCE_THRESHOLD

  detail$cycling_phase_review <-
    detail$n_cells >=
      MIN_CLUSTER_SIZE_FOR_DOMINANCE &
    detail$dominant_phase %in%
      c(
        "S",
        "G2M"
      ) &
    detail$maximum_phase_fraction >=
      CYCLING_PHASE_DOMINANCE_THRESHOLD

  summary <- data.frame(
    record_type = "resolution_summary",
    branch = branch,
    resolution = resolution,
    cluster_column = cluster_column,
    cluster = NA_character_,
    n_cells = nrow(
      metadata
    ),
    n_clusters = nrow(
      detail
    ),
    minimum_cluster_size = min(
      detail$n_cells
    ),
    median_cluster_size = median(
      detail$n_cells
    ),
    maximum_cluster_size = max(
      detail$n_cells
    ),
    n_clusters_lt_100 = sum(
      detail$n_cells < 100L
    ),
    n_clusters_lt_500 = sum(
      detail$n_cells < 500L
    ),
    n_library_dominant_clusters = sum(
      detail$library_dominance_review
    ),
    n_timepoint_dominant_clusters = sum(
      detail$timepoint_dominance_review
    ),
    n_parent_dominant_clusters = sum(
      detail$parent_cluster_dominance_review
    ),
    n_cycling_phase_dominant_clusters = sum(
      detail$cycling_phase_review
    ),
    adjacent_resolution_ARI =
      NA_real_,
    dominant_library =
      NA_character_,
    maximum_library_fraction =
      NA_real_,
    dominant_timepoint =
      NA_character_,
    maximum_timepoint_fraction =
      NA_real_,
    dominant_parent_cluster =
      NA_character_,
    maximum_parent_cluster_fraction =
      NA_real_,
    dominant_phase =
      NA_character_,
    maximum_phase_fraction =
      NA_real_,
    stringsAsFactors = FALSE
  )

  detail_output <- data.frame(
    record_type = "cluster_detail",
    branch = detail$branch,
    resolution =
      detail$resolution,
    cluster_column =
      detail$cluster_column,
    cluster =
      detail$cluster,
    n_cells =
      detail$n_cells,
    n_clusters =
      NA_integer_,
    minimum_cluster_size =
      NA_integer_,
    median_cluster_size =
      NA_real_,
    maximum_cluster_size =
      NA_integer_,
    n_clusters_lt_100 =
      NA_integer_,
    n_clusters_lt_500 =
      NA_integer_,
    n_library_dominant_clusters =
      NA_integer_,
    n_timepoint_dominant_clusters =
      NA_integer_,
    n_parent_dominant_clusters =
      NA_integer_,
    n_cycling_phase_dominant_clusters =
      NA_integer_,
    adjacent_resolution_ARI =
      NA_real_,
    dominant_library =
      detail$dominant_library,
    maximum_library_fraction =
      detail$maximum_library_fraction,
    dominant_timepoint =
      detail$dominant_timepoint,
    maximum_timepoint_fraction =
      detail$maximum_timepoint_fraction,
    dominant_parent_cluster =
      detail$dominant_parent_cluster,
    maximum_parent_cluster_fraction =
      detail$maximum_parent_cluster_fraction,
    dominant_phase =
      detail$dominant_phase,
    maximum_phase_fraction =
      detail$maximum_phase_fraction,
    stringsAsFactors = FALSE
  )

  list(
    summary = summary,
    detail = detail_output
  )
}

# -----------------------------------------------------------------------------
# 4. Package preflight
# -----------------------------------------------------------------------------
log_message(
  "Starting package preflight."
)

required_packages <- c(
  "Matrix",
  "data.table",
  "ggplot2",
  "patchwork",
  "scales",
  "Seurat",
  "SeuratObject",
  "BPCells",
  "harmony"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (
  "harmony" %in%
    missing_packages
) {
  log_message(
    "Installing harmony from CRAN."
  )

  try(
    utils::install.packages(
      "harmony",
      repos =
        "https://cloud.r-project.org",
      dependencies = TRUE
    ),
    silent = TRUE
  )
}

missing_after_install <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (
  length(
    missing_after_install
  ) > 0L
) {
  stop_with_report(
    paste0(
      "Missing required package(s): ",
      paste(
        missing_after_install,
        collapse = ", "
      ),
      "."
    )
  )
}

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(Seurat)
  library(SeuratObject)
  library(BPCells)
  library(harmony)
})

package_versions <- data.frame(
  package = c(
    "R",
    "Matrix",
    "data.table",
    "ggplot2",
    "patchwork",
    "scales",
    "Seurat",
    "SeuratObject",
    "BPCells",
    "harmony"
  ),
  version = c(
    paste0(
      R.version$major,
      ".",
      R.version$minor
    ),
    safe_package_version(
      "Matrix"
    ),
    safe_package_version(
      "data.table"
    ),
    safe_package_version(
      "ggplot2"
    ),
    safe_package_version(
      "patchwork"
    ),
    safe_package_version(
      "scales"
    ),
    safe_package_version(
      "Seurat"
    ),
    safe_package_version(
      "SeuratObject"
    ),
    safe_package_version(
      "BPCells"
    ),
    safe_package_version(
      "harmony"
    )
  ),
  stringsAsFactors = FALSE
)

write_csv_safe(
  package_versions,
  "00_package_versions.csv"
)

# -----------------------------------------------------------------------------
# 5. Load Step 10D evidence object and extract the frozen myeloid candidates
# -----------------------------------------------------------------------------
if (!file.exists(
  input_object_path
)) {
  stop_with_report(
    paste0(
      "Input Step 10D object not found: ",
      input_object_path
    )
  )
}

log_message(
  "Loading Step 10D evidence object."
)

atlas_obj <- readRDS(
  input_object_path
)

if (!inherits(
  atlas_obj,
  "Seurat"
)) {
  stop_with_report(
    "Input RDS is not a Seurat object."
  )
}

if (
  nrow(
    atlas_obj
  ) != 26918L ||
    ncol(
      atlas_obj
    ) != 144359L
) {
  stop_with_report(
    paste0(
      "Unexpected Step 10D dimensions: ",
      nrow(
        atlas_obj
      ),
      " x ",
      ncol(
        atlas_obj
      ),
      "."
    )
  )
}

required_metadata <- c(
  "biological_library",
  "timepoint",
  "replicate_id",
  "nCount_RNA",
  "nFeature_RNA",
  "percent.mt",
  "percent.HB",
  "scDblFinder.class",
  PARENT_CLUSTER_COLUMN
)

missing_metadata <- setdiff(
  required_metadata,
  names(
    atlas_obj@meta.data
  )
)

if (
  length(
    missing_metadata
  ) > 0L
) {
  stop_with_report(
    paste0(
      "Missing required metadata column(s): ",
      paste(
        missing_metadata,
        collapse = ", "
      ),
      "."
    )
  )
}

if (!all(
  as.character(
    atlas_obj$scDblFinder.class
  ) ==
    "singlet"
)) {
  stop_with_report(
    "Step 10D object contains non-singlet nuclei."
  )
}

DefaultAssay(
  atlas_obj
) <- "RNA"

parent_cluster_values <- as.character(
  atlas_obj@meta.data[
    ,
    PARENT_CLUSTER_COLUMN
  ]
)

myeloid_cells <- rownames(
  atlas_obj@meta.data
)[
  parent_cluster_values %in%
    PARENT_CLUSTERS_KEEP
]

if (
  length(
    myeloid_cells
  ) !=
    EXPECTED_MYLOID_NUCLEI
) {
  stop_with_report(
    paste0(
      "Expected ",
      EXPECTED_MYLOID_NUCLEI,
      " nuclei from parent clusters 1 and 21; detected ",
      length(
        myeloid_cells
      ),
      "."
    )
  )
}

parent_count_table <- as.data.frame(
  table(
    parent_atlas_cluster =
      parent_cluster_values[
        parent_cluster_values %in%
          PARENT_CLUSTERS_KEEP
      ]
  ),
  stringsAsFactors = FALSE
)

names(
  parent_count_table
)[
  names(
    parent_count_table
  ) ==
    "Freq"
] <- "n_nuclei"

write_csv_safe(
  parent_count_table,
  "01_parent_cluster_extraction_check.csv"
)

# -----------------------------------------------------------------------------
# 6. Write a stable myeloid-only BPCells counts matrix
# -----------------------------------------------------------------------------
rna_layers <- SeuratObject::Layers(
  object = atlas_obj,
  assay = "RNA",
  search = NA
)

counts_layers <- rna_layers[
  grepl(
    "^counts($|\\.)",
    rna_layers
  )
]

if (
  length(
    counts_layers
  ) != 1L
) {
  stop_with_report(
    paste0(
      "Expected one RNA counts layer; detected ",
      length(
        counts_layers
      ),
      "."
    )
  )
}

counts_layer_name <- counts_layers[[1L]]

atlas_counts <- SeuratObject::LayerData(
  object = atlas_obj,
  assay = "RNA",
  layer = counts_layer_name,
  fast = FALSE
)

myeloid_cell_indices <- match(
  myeloid_cells,
  colnames(
    atlas_counts
  )
)

if (anyNA(
  myeloid_cell_indices
)) {
  stop_with_report(
    "At least one extracted myeloid barcode is absent from the RNA counts layer."
  )
}

if (
  isTRUE(
    FORCE_REWRITE_MYLOID_COUNTS
  ) &&
    dir.exists(
      myeloid_counts_dir
    )
) {
  log_message(
    "Deleting existing myeloid counts BPCells directory."
  )

  unlink(
    myeloid_counts_dir,
    recursive = TRUE,
    force = TRUE
  )
}

if (!dir.exists(
  myeloid_counts_dir
)) {
  log_message(
    "Writing the stable myeloid-only BPCells counts matrix."
  )

  myeloid_counts_view <- atlas_counts[
    ,
    myeloid_cell_indices,
    drop = FALSE
  ]

  BPCells::write_matrix_dir(
    mat = myeloid_counts_view,
    dir = myeloid_counts_dir,
    compress = TRUE,
    overwrite = FALSE
  )

  rm(
    myeloid_counts_view
  )

  invisible(
    gc()
  )
} else {
  log_message(
    "Reusing existing myeloid counts BPCells directory: ",
    myeloid_counts_dir
  )
}

myeloid_counts <- open_bpcells_safe(
  myeloid_counts_dir
)

if (inherits(
  myeloid_counts,
  "error"
)) {
  stop_with_report(
    paste0(
      "Myeloid counts BPCells directory cannot be opened: ",
      conditionMessage(
        myeloid_counts
      ),
      ". Set FORCE_REWRITE_MYLOID_COUNTS <- TRUE and rerun."
    )
  )
}

if (
  nrow(
    myeloid_counts
  ) !=
    nrow(
      atlas_obj
    ) ||
    ncol(
      myeloid_counts
    ) !=
    EXPECTED_MYLOID_NUCLEI ||
    !identical(
      colnames(
        myeloid_counts
      ),
      myeloid_cells
    )
) {
  stop_with_report(
    "Myeloid counts BPCells matrix dimensions or barcode order are incorrect."
  )
}

# -----------------------------------------------------------------------------
# 7. Create a fresh myeloid Seurat object
# -----------------------------------------------------------------------------
myeloid_metadata <- atlas_obj@meta.data[
  myeloid_cells,
  ,
  drop = FALSE
]

myeloid_metadata$parent_atlas_cluster <-
  as.character(
    myeloid_metadata[
      ,
      PARENT_CLUSTER_COLUMN
    ]
  )

myeloid_metadata$myeloid_subset_source <-
  "whole_atlas_clusters_1_and_21"

myeloid_metadata$myeloid_annotation_status <-
  "RECLUSTERING_REVIEW_NOT_FROZEN"

# Remove stale whole-atlas nCount/nFeature values so CreateSeuratObject creates
# values from the stable myeloid counts matrix.
if (
  "nCount_RNA" %in%
    names(
      myeloid_metadata
    )
) {
  myeloid_metadata$atlas_nCount_RNA <-
    myeloid_metadata$nCount_RNA

  myeloid_metadata$nCount_RNA <- NULL
}

if (
  "nFeature_RNA" %in%
    names(
      myeloid_metadata
    )
) {
  myeloid_metadata$atlas_nFeature_RNA <-
    myeloid_metadata$nFeature_RNA

  myeloid_metadata$nFeature_RNA <- NULL
}

log_message(
  "Creating the fresh myeloid-only Seurat object."
)

myeloid_obj <- Seurat::CreateSeuratObject(
  counts = myeloid_counts,
  assay = "RNA",
  project = "GSE234774_myeloid_reclustering",
  meta.data = myeloid_metadata,
  min.cells = 0,
  min.features = 0
)

DefaultAssay(
  myeloid_obj
) <- "RNA"

if (
  nrow(
    myeloid_obj
  ) != 26918L ||
    ncol(
      myeloid_obj
    ) !=
    EXPECTED_MYLOID_NUCLEI
) {
  stop_with_report(
    "Fresh myeloid Seurat object dimensions are incorrect."
  )
}

rm(
  atlas_obj,
  atlas_counts
)

invisible(
  gc()
)

# -----------------------------------------------------------------------------
# 8. Re-run LogNormalize on the myeloid subset
# -----------------------------------------------------------------------------
log_message(
  "Running myeloid-only LogNormalize."
)

myeloid_obj <- Seurat::NormalizeData(
  object = myeloid_obj,
  assay = "RNA",
  normalization.method =
    "LogNormalize",
  scale.factor =
    NORMALIZATION_SCALE_FACTOR,
  block.size =
    NORMALIZATION_BLOCK_SIZE,
  verbose = TRUE
)

myeloid_data_layers <- SeuratObject::Layers(
  object = myeloid_obj,
  assay = "RNA",
  search = "^data"
)

if (
  length(
    myeloid_data_layers
  ) != 1L
) {
  stop_with_report(
    paste0(
      "Expected one normalized data layer; detected ",
      length(
        myeloid_data_layers
      ),
      "."
    )
  )
}

myeloid_data_layer_name <-
  myeloid_data_layers[[1L]]

myeloid_data_view <- SeuratObject::LayerData(
  object = myeloid_obj,
  assay = "RNA",
  layer =
    myeloid_data_layer_name,
  fast = FALSE
)

if (
  isTRUE(
    FORCE_REWRITE_MYLOID_DATA
  ) &&
    dir.exists(
      myeloid_data_dir
    )
) {
  log_message(
    "Deleting existing myeloid normalized-data BPCells directory."
  )

  unlink(
    myeloid_data_dir,
    recursive = TRUE,
    force = TRUE
  )
}

# Explicitly preserve feature/cell names before writing to BPCells.
rownames(myeloid_data_view) <- rownames(
  myeloid_obj
)

colnames(myeloid_data_view) <- colnames(
  myeloid_obj
)

validate_myeloid_data_matrix <- function(
    matrix_object,
    seurat_object
) {
  observed_dim <- as.integer(
    dim(
      matrix_object
    )
  )

  expected_dim <- as.integer(
    dim(
      seurat_object
    )
  )

  observed_features <- as.character(
    rownames(
      matrix_object
    )
  )

  expected_features <- as.character(
    rownames(
      seurat_object
    )
  )

  observed_cells <- as.character(
    colnames(
      matrix_object
    )
  )

  expected_cells <- as.character(
    colnames(
      seurat_object
    )
  )

  same_dimensions <- identical(
    observed_dim,
    expected_dim
  )

  same_features <- identical(
    observed_features,
    expected_features
  )

  same_cells <- identical(
    observed_cells,
    expected_cells
  )

  list(
    pass =
      same_dimensions &&
      same_features &&
      same_cells,
    same_dimensions =
      same_dimensions,
    same_features =
      same_features,
    same_cells =
      same_cells,
    observed_dim =
      observed_dim,
    expected_dim =
      expected_dim
  )
}

write_myeloid_data_matrix <- function() {
  if (dir.exists(
    myeloid_data_dir
  )) {
    unlink(
      myeloid_data_dir,
      recursive = TRUE,
      force = TRUE
    )
  }

  log_message(
    "Writing the stable myeloid LogNormalize data BPCells matrix."
  )

  BPCells::write_matrix_dir(
    mat = myeloid_data_view,
    dir = myeloid_data_dir,
    compress = TRUE,
    overwrite = FALSE
  )
}

# Write when absent or explicitly requested.
if (
  !dir.exists(
    myeloid_data_dir
  ) ||
    isTRUE(
      FORCE_REWRITE_MYLOID_DATA
    )
) {
  write_myeloid_data_matrix()
} else {
  log_message(
    "Checking the existing myeloid normalized-data BPCells directory before reuse: ",
    myeloid_data_dir
  )
}

myeloid_data_disk <- open_bpcells_safe(
  myeloid_data_dir
)

# Automatically rebuild a missing, incomplete, stale, or mismatched directory.
needs_rebuild <- inherits(
  myeloid_data_disk,
  "error"
)

if (!needs_rebuild) {
  initial_validation <- validate_myeloid_data_matrix(
    matrix_object =
      myeloid_data_disk,
    seurat_object =
      myeloid_obj
  )

  needs_rebuild <- !isTRUE(
    initial_validation$pass
  )

  if (needs_rebuild) {
    log_message(
      "Existing normalized BPCells directory failed validation and will be rebuilt. ",
      "same_dimensions=",
      initial_validation$same_dimensions,
      "; same_features=",
      initial_validation$same_features,
      "; same_cells=",
      initial_validation$same_cells,
      "; observed_dim=",
      paste(
        initial_validation$observed_dim,
        collapse = "x"
      ),
      "; expected_dim=",
      paste(
        initial_validation$expected_dim,
        collapse = "x"
      ),
      "."
    )
  }
} else {
  log_message(
    "Existing normalized BPCells directory could not be opened and will be rebuilt: ",
    conditionMessage(
      myeloid_data_disk
    )
  )
}

if (needs_rebuild) {
  rm(
    myeloid_data_disk
  )

  invisible(
    gc()
  )

  write_myeloid_data_matrix()

  myeloid_data_disk <- open_bpcells_safe(
    myeloid_data_dir
  )
}

if (inherits(
  myeloid_data_disk,
  "error"
)) {
  stop_with_report(
    paste0(
      "Myeloid normalized BPCells directory still cannot be opened after automatic rebuild: ",
      conditionMessage(
        myeloid_data_disk
      ),
      "."
    )
  )
}

final_data_validation <- validate_myeloid_data_matrix(
  matrix_object =
    myeloid_data_disk,
  seurat_object =
    myeloid_obj
)

if (!isTRUE(
  final_data_validation$pass
)) {
  stop_with_report(
    paste0(
      "Myeloid normalized BPCells matrix still fails validation after automatic rebuild. ",
      "same_dimensions=",
      final_data_validation$same_dimensions,
      "; same_features=",
      final_data_validation$same_features,
      "; same_cells=",
      final_data_validation$same_cells,
      "; observed_dim=",
      paste(
        final_data_validation$observed_dim,
        collapse = "x"
      ),
      "; expected_dim=",
      paste(
        final_data_validation$expected_dim,
        collapse = "x"
      ),
      "."
    )
  )
}

log_message(
  "Myeloid normalized BPCells matrix validation PASS: ",
  paste(
    final_data_validation$observed_dim,
    collapse = " x "
  ),
  "."
)

myeloid_obj[["RNA"]]$data <-
  myeloid_data_disk

rm(
  myeloid_data_view
)

invisible(
  gc()
)

# -----------------------------------------------------------------------------
# 9. Cell-cycle scoring without regression
# -----------------------------------------------------------------------------
log_message(
  "Scoring cell cycle without regression."
)

utils::data(
  "cc.genes.updated.2019",
  package = "Seurat"
)

s_features <- Seurat::CaseMatch(
  search =
    cc.genes.updated.2019$s.genes,
  match = rownames(
    myeloid_obj
  )
)

g2m_features <- Seurat::CaseMatch(
  search =
    cc.genes.updated.2019$g2m.genes,
  match = rownames(
    myeloid_obj
  )
)

s_features <- unique(
  unname(
    s_features[
      !is.na(
        s_features
      )
    ]
  )
)

g2m_features <- unique(
  unname(
    g2m_features[
      !is.na(
        g2m_features
      )
    ]
  )
)

if (
  length(
    s_features
  ) < 20L ||
    length(
      g2m_features
    ) < 20L
) {
  stop_with_report(
    paste0(
      "Insufficient cell-cycle genes were resolved: S=",
      length(
        s_features
      ),
      "; G2M=",
      length(
        g2m_features
      ),
      "."
    )
  )
}

myeloid_obj <- Seurat::CellCycleScoring(
  object = myeloid_obj,
  s.features = s_features,
  g2m.features = g2m_features,
  set.ident = FALSE
)

phase_summary <- data.table::as.data.table(
  myeloid_obj@meta.data
)[
  ,
  .(
    n_cells = .N
  ),
  by = .(
    parent_atlas_cluster,
    timepoint,
    Phase
  )
]

write_csv_safe(
  as.data.frame(
    phase_summary
  ),
  "02_parent_cluster_timepoint_phase_summary.csv"
)

# -----------------------------------------------------------------------------
# 10. Myeloid-only variable features and identity-focused PCA
# -----------------------------------------------------------------------------
log_message(
  "Selecting myeloid-only variable features."
)

myeloid_obj <- Seurat::FindVariableFeatures(
  object = myeloid_obj,
  assay = "RNA",
  selection.method = "vst",
  nfeatures =
    N_VARIABLE_FEATURES,
  verbose = TRUE
)

variable_features <- Seurat::VariableFeatures(
  object = myeloid_obj,
  assay = "RNA"
)

if (
  length(
    variable_features
  ) !=
    N_VARIABLE_FEATURES
) {
  stop_with_report(
    paste0(
      "Expected ",
      N_VARIABLE_FEATURES,
      " variable features; detected ",
      length(
        variable_features
      ),
      "."
    )
  )
}

cell_cycle_features <- unique(
  c(
    s_features,
    g2m_features
  )
)

technical_identity_exclusions <-
  variable_features[
    grepl(
      "^(mt-|Rpl|Rps|Hb[ab])",
      variable_features,
      ignore.case = TRUE
    ) |
      variable_features %in%
        cell_cycle_features
  ]

identity_pca_features <- setdiff(
  variable_features,
  technical_identity_exclusions
)

if (
  length(
    identity_pca_features
  ) < 1800L
) {
  stop_with_report(
    paste0(
      "Only ",
      length(
        identity_pca_features
      ),
      " identity-focused PCA features remain after exclusions."
    )
  )
}

variable_feature_audit <- data.frame(
  feature = variable_features,
  excluded_from_identity_PCA =
    variable_features %in%
      technical_identity_exclusions,
  exclusion_reason = ifelse(
    grepl(
      "^mt-",
      variable_features,
      ignore.case = TRUE
    ),
    "mitochondrial",
    ifelse(
      grepl(
        "^(Rpl|Rps)",
        variable_features,
        ignore.case = TRUE
      ),
      "ribosomal",
      ifelse(
        grepl(
          "^Hb[ab]",
          variable_features,
          ignore.case = TRUE
        ),
        "hemoglobin",
        ifelse(
          variable_features %in%
            cell_cycle_features,
          "cell_cycle",
          ""
        )
      )
    )
  ),
  stringsAsFactors = FALSE
)

write_csv_safe(
  variable_feature_audit,
  "03_myeloid_variable_feature_audit.csv"
)

log_message(
  "Scaling identity-focused PCA features without regression."
)

myeloid_obj <- Seurat::ScaleData(
  object = myeloid_obj,
  assay = "RNA",
  features =
    identity_pca_features,
  vars.to.regress = NULL,
  do.center = TRUE,
  do.scale = TRUE,
  block.size = 1000L,
  verbose = TRUE
)

log_message(
  "Running identity-focused myeloid PCA."
)

set.seed(
  SEED
)

myeloid_obj <- Seurat::RunPCA(
  object = myeloid_obj,
  assay = "RNA",
  features =
    identity_pca_features,
  npcs =
    N_PCS_COMPUTE,
  reduction.name =
    "pca_myeloid_identity",
  reduction.key =
    "MYPC_",
  approx = TRUE,
  seed.use =
    SEED,
  verbose = TRUE
)

pca_embeddings <- Seurat::Embeddings(
  object = myeloid_obj,
  reduction =
    "pca_myeloid_identity"
)

if (
  nrow(
    pca_embeddings
  ) !=
    ncol(
      myeloid_obj
    ) ||
    ncol(
      pca_embeddings
    ) <
    max(
      PCS_USE
    )
) {
  stop_with_report(
    "Myeloid PCA embedding dimensions are insufficient."
  )
}

pca_embeddings_use <- pca_embeddings[
  ,
  PCS_USE,
  drop = FALSE
]

# -----------------------------------------------------------------------------
# 11. Conservative myeloid Harmony
# -----------------------------------------------------------------------------
metadata_for_harmony <-
  myeloid_obj@meta.data

harmony_embeddings <-
  run_harmony_matrix_compat(
    pca_embeddings =
      pca_embeddings_use,
    metadata =
      metadata_for_harmony
  )

if (
  nrow(
    harmony_embeddings
  ) !=
    ncol(
      myeloid_obj
    ) ||
    ncol(
      harmony_embeddings
    ) !=
    length(
      PCS_USE
    )
) {
  stop_with_report(
    "Myeloid Harmony embedding dimensions are incorrect."
  )
}

colnames(
  harmony_embeddings
) <- paste0(
  "MYHCON_",
  seq_len(
    ncol(
      harmony_embeddings
    )
  )
)

myeloid_obj[["harmony_myeloid_conservative"]] <- SeuratObject::CreateDimReducObject(
  embeddings =
    harmony_embeddings,
  key = "MYHCON_",
  assay = "RNA"
)

# -----------------------------------------------------------------------------
# 12. PCA and Harmony graphs, cluster sweep, and UMAPs
# -----------------------------------------------------------------------------
summary_rows <- list()
detail_rows <- list()

for (
  branch_index in seq_len(
    nrow(
      branch_specs
    )
  )
) {
  branch <- branch_specs$branch[
    branch_index
  ]

  reduction_name <- branch_specs$reduction[
    branch_index
  ]

  graph_prefix <- branch_specs$graph_prefix[
    branch_index
  ]

  cluster_prefix <- branch_specs$cluster_prefix[
    branch_index
  ]

  umap_name <- branch_specs$umap[
    branch_index
  ]

  graph_names <- c(
    paste0(
      graph_prefix,
      "_nn"
    ),
    paste0(
      graph_prefix,
      "_snn"
    )
  )

  log_message(
    "Building myeloid graph for ",
    branch,
    "."
  )

  myeloid_obj <- Seurat::FindNeighbors(
    object = myeloid_obj,
    reduction =
      reduction_name,
    dims =
      PCS_USE,
    k.param =
      K_GRAPH,
    compute.SNN = TRUE,
    prune.SNN =
      SNN_PRUNE,
    nn.method = "annoy",
    n.trees =
      N_TREES,
    annoy.metric =
      "euclidean",
    graph.name =
      graph_names,
    verbose = TRUE
  )

  for (
    resolution_index in seq_along(
      RESOLUTIONS
    )
  ) {
    resolution <-
      RESOLUTIONS[
        resolution_index
      ]

    cluster_column <- cluster_column_name(
      cluster_prefix,
      resolution
    )

    log_message(
      "Clustering ",
      branch,
      " at resolution ",
      resolution,
      "."
    )

    myeloid_obj <- Seurat::FindClusters(
      object = myeloid_obj,
      graph.name =
        graph_names[[2L]],
      cluster.name =
        cluster_column,
      resolution =
        resolution,
      algorithm = 1,
      n.start = 10L,
      n.iter = 10L,
      random.seed =
        SEED +
        branch_index * 100L +
        resolution_index,
      group.singletons = TRUE,
      verbose = TRUE
    )

    clustering_result <-
      summarize_clustering(
        metadata =
          myeloid_obj@meta.data,
        cluster_column =
          cluster_column,
        branch =
          branch,
        resolution =
          resolution
      )

    summary_rows[[length(summary_rows) + 1L]] <- clustering_result$summary

    detail_rows[[length(detail_rows) + 1L]] <- clustering_result$detail
  }

  log_message(
    "Running UMAP for ",
    branch,
    "."
  )

  myeloid_obj <- Seurat::RunUMAP(
    object = myeloid_obj,
    reduction =
      reduction_name,
    dims =
      PCS_USE,
    reduction.name =
      umap_name,
    reduction.key = if (
      branch ==
        "myeloid_pca_identity"
    ) {
      "MYPCUMAP_"
    } else {
      "MYHUMAP_"
    },
    umap.method = "uwot",
    n.neighbors =
      UMAP_N_NEIGHBORS,
    n.components = 2L,
    metric = "cosine",
    min.dist =
      UMAP_MIN_DIST,
    uwot.sgd = TRUE,
    seed.use =
      SEED +
      branch_index * 1000L,
    verbose = TRUE
  )
}

resolution_summary <- do.call(
  rbind,
  summary_rows
)

resolution_detail <- do.call(
  rbind,
  detail_rows
)

# Adjacent-resolution ARI within each branch.
for (
  branch_index in seq_len(
    nrow(
      branch_specs
    )
  )
) {
  branch <- branch_specs$branch[
    branch_index
  ]

  prefix <- branch_specs$cluster_prefix[
    branch_index
  ]

  for (
    resolution_index in seq_len(
      length(
        RESOLUTIONS
      ) - 1L
    )
  ) {
    resolution_a <-
      RESOLUTIONS[
        resolution_index
      ]

    resolution_b <-
      RESOLUTIONS[
        resolution_index +
          1L
      ]

    column_a <- cluster_column_name(
      prefix,
      resolution_a
    )

    column_b <- cluster_column_name(
      prefix,
      resolution_b
    )

    ari_value <- adjusted_rand_index(
      myeloid_obj@meta.data[
        ,
        column_a
      ],
      myeloid_obj@meta.data[
        ,
        column_b
      ]
    )

    target_index <- which(
      resolution_summary$branch ==
        branch &
        resolution_summary$resolution ==
        resolution_a
    )

    resolution_summary$adjacent_resolution_ARI[
      target_index
    ] <- ari_value
  }
}

resolution_output <- rbind(
  resolution_summary,
  resolution_detail
)

resolution_output <- resolution_output[
  order(
    match(
      resolution_output$branch,
      branch_specs$branch
    ),
    resolution_output$resolution,
    resolution_output$record_type,
    suppressWarnings(
      as.numeric(
        resolution_output$cluster
      )
    )
  ),
  ,
  drop = FALSE
]

write_csv_safe(
  resolution_output,
  "04_myeloid_resolution_cluster_summary_and_detail.csv"
)

# -----------------------------------------------------------------------------
# 13. Within-timepoint library mixing and local-structure preservation
# -----------------------------------------------------------------------------
log_message(
  "Calculating within-timepoint library mixing for myeloid PCA."
)

mixing_pca <- calculate_within_timepoint_mixing(
  embeddings =
    pca_embeddings_use,
  metadata =
    myeloid_obj@meta.data,
  reduction_name =
    "myeloid_pca_identity",
  k =
    K_MIXING
)

log_message(
  "Calculating within-timepoint library mixing for myeloid Harmony."
)

mixing_harmony <-
  calculate_within_timepoint_mixing(
    embeddings =
      harmony_embeddings,
    metadata =
      myeloid_obj@meta.data,
    reduction_name =
      "myeloid_harmony_conservative",
    k =
      K_MIXING
  )

mixing_summary <- rbind(
  mixing_pca$summaries,
  mixing_harmony$summaries
)

mixing_summary$timepoint <- factor(
  mixing_summary$timepoint,
  levels =
    expected_timepoints,
  ordered = TRUE
)

pca_reference <- mixing_summary[
  mixing_summary$reduction ==
    "myeloid_pca_identity",
  c(
    "timepoint",
    "mean_other_library_fraction",
    "mean_normalized_library_entropy"
  ),
  drop = FALSE
]

names(
  pca_reference
)[
  names(
    pca_reference
  ) ==
    "mean_other_library_fraction"
] <- "PCA_mean_other_library_fraction"

names(
  pca_reference
)[
  names(
    pca_reference
  ) ==
    "mean_normalized_library_entropy"
] <- "PCA_mean_normalized_library_entropy"

mixing_summary <- merge(
  mixing_summary,
  pca_reference,
  by = "timepoint",
  all.x = TRUE,
  sort = FALSE
)

mixing_summary$gain_vs_PCA_other_library_fraction <-
  mixing_summary$mean_other_library_fraction -
  mixing_summary$PCA_mean_other_library_fraction

mixing_summary$gain_vs_PCA_entropy <-
  mixing_summary$mean_normalized_library_entropy -
  mixing_summary$PCA_mean_normalized_library_entropy

neighbor_overlap <-
  calculate_neighbor_overlap(
    reference_cache =
      mixing_pca$neighbor_cache,
    comparison_cache =
      mixing_harmony$neighbor_cache
  )

neighbor_overlap$reduction <-
  "myeloid_harmony_conservative"

mixing_summary <- merge(
  mixing_summary,
  neighbor_overlap[
    ,
    c(
      "reduction",
      "timepoint",
      "mean_KNN_overlap",
      "median_KNN_overlap",
      "q10_KNN_overlap",
      "q90_KNN_overlap"
    ),
    drop = FALSE
  ],
  by = c(
    "reduction",
    "timepoint"
  ),
  all.x = TRUE,
  sort = FALSE
)

mixing_summary$mean_KNN_overlap[
  mixing_summary$reduction ==
    "myeloid_pca_identity"
] <- 1

mixing_summary$median_KNN_overlap[
  mixing_summary$reduction ==
    "myeloid_pca_identity"
] <- 1

mixing_summary$q10_KNN_overlap[
  mixing_summary$reduction ==
    "myeloid_pca_identity"
] <- 1

mixing_summary$q90_KNN_overlap[
  mixing_summary$reduction ==
    "myeloid_pca_identity"
] <- 1

mixing_summary <- mixing_summary[
  order(
    mixing_summary$timepoint,
    match(
      mixing_summary$reduction,
      branch_specs$branch
    )
  ),
  ,
  drop = FALSE
]

write_csv_safe(
  mixing_summary,
  "05_myeloid_within_timepoint_library_mixing.csv"
)

# -----------------------------------------------------------------------------
# 14. Resolve curated marker panels
# -----------------------------------------------------------------------------
all_requested_features <- unique(
  unlist(
    marker_panels_requested,
    use.names = FALSE
  )
)

resolved_features <- resolve_features(
  requested_features =
    all_requested_features,
  available_features =
    rownames(
      myeloid_obj
    )
)

marker_panels_present <- lapply(
  marker_panels_requested,
  function(panel_genes) {
    resolved <- resolved_features[
      panel_genes
    ]

    unique(
      unname(
        resolved[
          !is.na(
            resolved
          )
        ]
      )
    )
  }
)

marker_panel_availability <- data.frame(
  panel = names(
    marker_panels_requested
  ),
  requested_n = vapply(
    marker_panels_requested,
    length,
    integer(1)
  ),
  present_n = vapply(
    marker_panels_present,
    length,
    integer(1)
  ),
  present_features = vapply(
    marker_panels_present,
    function(features) {
      paste(
        features,
        collapse = ";"
      )
    },
    character(1)
  ),
  stringsAsFactors = FALSE
)

marker_panel_availability$present_fraction <-
  marker_panel_availability$present_n /
  marker_panel_availability$requested_n

write_csv_safe(
  marker_panel_availability,
  "06_myeloid_marker_panel_availability.csv"
)

if (
  any(
    marker_panel_availability$present_n <
      2L
  )
) {
  failed_panels <-
    marker_panel_availability$panel[
      marker_panel_availability$present_n <
        2L
    ]

  stop_with_report(
    paste0(
      "Fewer than two genes were resolved for marker panel(s): ",
      paste(
        failed_panels,
        collapse = ", "
      ),
      "."
    )
  )
}

marker_plot_panels <- deduplicate_plot_panels(
  marker_panels_present
)

# -----------------------------------------------------------------------------
# 15. Figures
# -----------------------------------------------------------------------------
log_message(
  "Creating myeloid PCA/Harmony/resolution review figures."
)

# Figure 1: PCA versus Harmony by parent cluster, timepoint, library, and phase.
comparison_plots <- list()

for (
  branch_index in seq_len(
    nrow(
      branch_specs
    )
  )
) {
  branch <- branch_specs$branch[
    branch_index
  ]

  umap_name <- branch_specs$umap[
    branch_index
  ]

  p_parent <- Seurat::DimPlot(
    object = myeloid_obj,
    reduction =
      umap_name,
    group.by =
      "parent_atlas_cluster",
    label = TRUE,
    repel = TRUE,
    raster = TRUE
  ) +
    ggplot2::labs(
      title = paste0(
        branch,
        ": parent atlas cluster"
      )
    )

  p_timepoint <- Seurat::DimPlot(
    object = myeloid_obj,
    reduction =
      umap_name,
    group.by =
      "timepoint",
    raster = TRUE
  ) +
    ggplot2::labs(
      title = paste0(
        branch,
        ": timepoint"
      )
    )

  p_library <- Seurat::DimPlot(
    object = myeloid_obj,
    reduction =
      umap_name,
    group.by =
      "biological_library",
    raster = TRUE
  ) +
    ggplot2::labs(
      title = paste0(
        branch,
        ": biological library"
      )
    )

  p_phase <- Seurat::DimPlot(
    object = myeloid_obj,
    reduction =
      umap_name,
    group.by =
      "Phase",
    raster = TRUE
  ) +
    ggplot2::labs(
      title = paste0(
        branch,
        ": cell-cycle phase"
      )
    )

  comparison_plots[[branch]] <- (
    p_parent |
      p_timepoint
  ) / (
    p_library |
      p_phase
  )
}

ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "01_myeloid_PCA_Harmony_UMAP_comparison.pdf"
  ),
  plot = comparison_plots[["myeloid_pca_identity"]] |
    comparison_plots[["myeloid_harmony_conservative"]],
  width = 22,
  height = 14,
  limitsize = FALSE
)

# Figure 2: Harmony UMAP resolution sweep.
resolution_plots <- lapply(
  RESOLUTIONS,
  function(resolution) {
    cluster_column <- cluster_column_name(
      "myhcon",
      resolution
    )

    Seurat::DimPlot(
      object = myeloid_obj,
      reduction =
        "umap_myeloid_harmony",
      group.by =
        cluster_column,
      label = TRUE,
      repel = TRUE,
      raster = TRUE
    ) +
      ggplot2::labs(
        title = paste0(
          "Myeloid conservative Harmony | resolution ",
          resolution
        )
      ) +
      Seurat::NoLegend()
  }
)

ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "02_myeloid_Harmony_resolution_0.2_to_1.0_UMAP.pdf"
  ),
  plot = patchwork::wrap_plots(
    resolution_plots,
    ncol = 2
  ),
  width = 16,
  height = 18,
  limitsize = FALSE
)

# Figure 3: one marker DotPlot page per Harmony resolution.
marker_dotplot_path <- file.path(
  figures_dir,
  "03_myeloid_marker_DotPlots_by_resolution.pdf"
)

grDevices::pdf(
  marker_dotplot_path,
  width = 20,
  height = 10,
  onefile = TRUE
)

for (resolution in RESOLUTIONS) {
  cluster_column <- cluster_column_name(
    "myhcon",
    resolution
  )

  print(
    Seurat::DotPlot(
      object = myeloid_obj,
      features =
        marker_plot_panels,
      assay = "RNA",
      group.by =
        cluster_column,
      cluster.idents = FALSE,
      dot.scale = 6
    ) +
      Seurat::RotatedAxis() +
      ggplot2::labs(
        title = paste0(
          "Myeloid identity/state marker evidence | Harmony resolution ",
          resolution
        ),
        x = NULL,
        y = "Myeloid Harmony cluster"
      ) +
      ggplot2::theme(
        axis.text.x =
          ggplot2::element_text(
            size = 7
          )
      )
  )
}

grDevices::dev.off()

log_message(
  "Wrote: ",
  marker_dotplot_path
)

# -----------------------------------------------------------------------------
# 16. Save the comparison object
# -----------------------------------------------------------------------------
# PCA/graphs/UMAPs are already stored. Remove the dense scale.data layer before
# saving to reduce object size. RNA counts and LogNormalize data remain.
myeloid_obj[["RNA"]]$scale.data <-
  NULL

myeloid_obj@misc$GSE234774_step10E <- list(
  input_object = normalizePath(
    input_object_path,
    winslash = "/",
    mustWork = TRUE
  ),
  parent_cluster_column =
    PARENT_CLUSTER_COLUMN,
  parent_clusters_kept =
    PARENT_CLUSTERS_KEEP,
  expected_myeloid_nuclei =
    EXPECTED_MYLOID_NUCLEI,
  normalization =
    "LogNormalize",
  variable_features =
    N_VARIABLE_FEATURES,
  identity_PCA_features =
    length(
      identity_pca_features
    ),
  excluded_cell_cycle_features =
    length(
      intersect(
        variable_features,
        cell_cycle_features
      )
    ),
  cell_cycle_scored = TRUE,
  cell_cycle_regressed = FALSE,
  PCA_reduction =
    "pca_myeloid_identity",
  Harmony_reduction =
    "harmony_myeloid_conservative",
  Harmony_group =
    "biological_library",
  Harmony_theta =
    HARMONY_THETA,
  Harmony_lambda =
    HARMONY_LAMBDA,
  UMAP_PCA =
    "umap_myeloid_pca",
  UMAP_Harmony =
    "umap_myeloid_harmony",
  resolutions =
    RESOLUTIONS,
  final_resolution_frozen = FALSE,
  final_myeloid_labels_frozen = FALSE,
  resident_microglia_object_frozen = FALSE,
  Arm_1_Arm_2_scoring_run = FALSE,
  completed_at =
    as.character(
      Sys.time()
    )
)

log_message(
  "Saving the Step 10E myeloid comparison object."
)

saveRDS(
  myeloid_obj,
  output_object_path,
  compress = FALSE
)

reload_obj <- tryCatch(
  readRDS(
    output_object_path
  ),
  error = function(error_object) {
    error_object
  }
)

required_reductions <- c(
  "pca_myeloid_identity",
  "harmony_myeloid_conservative",
  "umap_myeloid_pca",
  "umap_myeloid_harmony"
)

required_cluster_columns <- unlist(
  lapply(
    branch_specs$cluster_prefix,
    function(prefix) {
      vapply(
        RESOLUTIONS,
        function(resolution) {
          cluster_column_name(
            prefix,
            resolution
          )
        },
        character(1)
      )
    }
  ),
  use.names = FALSE
)

reload_pass <-
  !inherits(
    reload_obj,
    "error"
  ) &&
  inherits(
    reload_obj,
    "Seurat"
  ) &&
  nrow(
    reload_obj
  ) == 26918L &&
  ncol(
    reload_obj
  ) ==
    EXPECTED_MYLOID_NUCLEI &&
  all(
    required_reductions %in%
      names(
        reload_obj@reductions
      )
  ) &&
  all(
    required_cluster_columns %in%
      names(
        reload_obj@meta.data
      )
  )

rm(
  reload_obj
)

invisible(
  gc()
)

# -----------------------------------------------------------------------------
# 17. Audit and final decision
# -----------------------------------------------------------------------------
valid_mixing_rows <- mixing_summary[
  mixing_summary$reduction ==
    "myeloid_harmony_conservative" &
    is.finite(
      mixing_summary$mean_KNN_overlap
    ),
  ,
  drop = FALSE
]

mean_harmony_KNN_overlap <- if (
  nrow(
    valid_mixing_rows
  ) > 0L
) {
  mean(
    valid_mixing_rows$mean_KNN_overlap,
    na.rm = TRUE
  )
} else {
  NA_real_
}

display_cluster_column <- cluster_column_name(
  "myhcon",
  DISPLAY_RESOLUTION
)

display_cluster_count <- length(
  unique(
    as.character(
      myeloid_obj@meta.data[
        ,
        display_cluster_column
      ]
    )
  )
)

required_figure_paths <- file.path(
  figures_dir,
  c(
    "01_myeloid_PCA_Harmony_UMAP_comparison.pdf",
    "02_myeloid_Harmony_resolution_0.2_to_1.0_UMAP.pdf",
    "03_myeloid_marker_DotPlots_by_resolution.pdf"
  )
)

audit_checks <- data.frame(
  check = c(
    "Step 10D input contains 26,918 features and 144,359 singlets",
    "Exactly 25,205 nuclei were extracted from parent clusters 1 and 21",
    "Stable myeloid counts and LogNormalize BPCells matrices are readable",
    "Cell-cycle scores were calculated without regression",
    "Exactly 3,000 myeloid variable features were selected",
    "At least 1,800 identity-focused PCA features remain after exclusions",
    "Identity-focused PCA and conservative Harmony contain all myeloid nuclei",
    "PCA and Harmony graphs, UMAPs, and all requested resolutions were created",
    "Within-timepoint mixing and KNN-overlap metrics were calculated",
    "All three required review PDFs were created",
    "Saved myeloid comparison object reloads successfully",
    "A final myeloid resolution has been frozen",
    "Final resident/infiltrating subtype labels have been frozen",
    "A clean resident-microglia object has been frozen",
    "GSE162610 Arm_1/Arm_2 scoring was run"
  ),
  required_for_pass = c(
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    FALSE,
    FALSE,
    FALSE,
    FALSE
  ),
  status = c(
    "PASS",
    ifelse(
      ncol(
        myeloid_obj
      ) ==
        EXPECTED_MYLOID_NUCLEI,
      "PASS",
      "STOP"
    ),
    ifelse(
      ncol(
        myeloid_counts
      ) ==
        EXPECTED_MYLOID_NUCLEI &&
        ncol(
          myeloid_data_disk
        ) ==
        EXPECTED_MYLOID_NUCLEI,
      "PASS",
      "STOP"
    ),
    ifelse(
      all(
        c(
          "S.Score",
          "G2M.Score",
          "Phase"
        ) %in%
          names(
            myeloid_obj@meta.data
          )
      ),
      "PASS",
      "STOP"
    ),
    ifelse(
      length(
        variable_features
      ) ==
        N_VARIABLE_FEATURES,
      "PASS",
      "STOP"
    ),
    ifelse(
      length(
        identity_pca_features
      ) >=
        1800L,
      "PASS",
      "STOP"
    ),
    ifelse(
      nrow(
        pca_embeddings_use
      ) ==
        EXPECTED_MYLOID_NUCLEI &&
        nrow(
          harmony_embeddings
        ) ==
        EXPECTED_MYLOID_NUCLEI,
      "PASS",
      "STOP"
    ),
    ifelse(
      all(
        required_reductions %in%
          names(
            myeloid_obj@reductions
          )
      ) &&
        all(
          required_cluster_columns %in%
            names(
              myeloid_obj@meta.data
            )
        ),
      "PASS",
      "STOP"
    ),
    ifelse(
      nrow(
        mixing_summary
      ) ==
        2L *
        length(
          expected_timepoints
        ),
      "PASS",
      "STOP"
    ),
    ifelse(
      all(
        file.exists(
          required_figure_paths
        )
      ),
      "PASS",
      "STOP"
    ),
    ifelse(
      reload_pass,
      "PASS",
      "STOP"
    ),
    "NOT_FROZEN",
    "NOT_FROZEN",
    "NOT_FROZEN",
    "CONFIRMED_NOT_RUN"
  ),
  detail = c(
    "Input object dimensions verified.",
    paste0(
      "myeloid nuclei=",
      ncol(
        myeloid_obj
      ),
      "; parent clusters=",
      paste(
        PARENT_CLUSTERS_KEEP,
        collapse = ","
      )
    ),
    paste(
      normalizePath(
        myeloid_counts_dir,
        winslash = "/",
        mustWork = TRUE
      ),
      normalizePath(
        myeloid_data_dir,
        winslash = "/",
        mustWork = TRUE
      ),
      sep = " ; "
    ),
    paste0(
      "S genes=",
      length(
        s_features
      ),
      "; G2M genes=",
      length(
        g2m_features
      ),
      "; vars.to.regress=NULL"
    ),
    paste0(
      "selected=",
      length(
        variable_features
      )
    ),
    paste0(
      "identity PCA features=",
      length(
        identity_pca_features
      ),
      "; technical/cycle exclusions=",
      length(
        technical_identity_exclusions
      )
    ),
    paste0(
      "PCA=",
      nrow(
        pca_embeddings_use
      ),
      " x ",
      ncol(
        pca_embeddings_use
      ),
      "; Harmony=",
      nrow(
        harmony_embeddings
      ),
      " x ",
      ncol(
        harmony_embeddings
      )
    ),
    paste0(
      "cluster columns=",
      length(
        required_cluster_columns
      ),
      "; display resolution ",
      DISPLAY_RESOLUTION,
      " has ",
      display_cluster_count,
      " clusters"
    ),
    paste0(
      "mixing rows=",
      nrow(
        mixing_summary
      ),
      "; mean Harmony-vs-PCA KNN overlap=",
      round(
        mean_harmony_KNN_overlap,
        3
      )
    ),
    paste(
      basename(
        required_figure_paths
      ),
      collapse = "; "
    ),
    ifelse(
      reload_pass,
      "RDS reload PASS",
      "RDS reload failed"
    ),
    "Pending review of figures, marker patterns, resolution stability, and cluster composition.",
    "Pending marker-first myeloid subtype annotation.",
    "Pending explicit resident-versus-infiltrating decision.",
    "Arm_1/Arm_2 remain frozen from GSE162610."
  ),
  stringsAsFactors = FALSE
)

write_csv_safe(
  audit_checks,
  "09_STEP10E_audit_checks.csv"
)

required_failures <-
  audit_checks$required_for_pass &
  audit_checks$status !=
    "PASS"

final_decision <- if (
  any(
    required_failures
  )
) {
  "STOP"
} else {
  "PASS_WITH_MYLOID_RESOLUTION_REVIEW"
}

# -----------------------------------------------------------------------------
# 18. Final report
# -----------------------------------------------------------------------------
summary_only <- resolution_summary[
  order(
    match(
      resolution_summary$branch,
      branch_specs$branch
    ),
    resolution_summary$resolution
  ),
  ,
  drop = FALSE
]

resolution_lines <- paste0(
  summary_only$branch,
  " | resolution ",
  summary_only$resolution,
  ": ",
  summary_only$n_clusters,
  " clusters; min/median/max=",
  summary_only$minimum_cluster_size,
  "/",
  round(
    summary_only$median_cluster_size,
    1
  ),
  "/",
  summary_only$maximum_cluster_size,
  "; <500=",
  summary_only$n_clusters_lt_500,
  "; library-dominant=",
  summary_only$n_library_dominant_clusters,
  "; timepoint-dominant=",
  summary_only$n_timepoint_dominant_clusters,
  "; parent-dominant=",
  summary_only$n_parent_dominant_clusters,
  "; cycling-phase-dominant=",
  summary_only$n_cycling_phase_dominant_clusters,
  "; adjacent ARI=",
  ifelse(
    is.na(
      summary_only$adjacent_resolution_ARI
    ),
    "NA",
    round(
      summary_only$adjacent_resolution_ARI,
      3
    )
  )
)

mixing_lines <- apply(
  mixing_summary,
  1L,
  function(row) {
    paste0(
      row[["timepoint"]],
      " | ",
      row[["reduction"]],
      ": other-library fraction=",
      ifelse(
        is.na(
          as.numeric(
            row[["mean_other_library_fraction"]]
          )
        ),
        "NA",
        round(
          as.numeric(
            row[["mean_other_library_fraction"]]
          ),
          3
        )
      ),
      "; entropy=",
      ifelse(
        is.na(
          as.numeric(
            row[["mean_normalized_library_entropy"]]
          )
        ),
        "NA",
        round(
          as.numeric(
            row[["mean_normalized_library_entropy"]]
          ),
          3
        )
      ),
      "; KNN overlap with PCA=",
      ifelse(
        is.na(
          as.numeric(
            row[["mean_KNN_overlap"]]
          )
        ),
        "NA",
        round(
          as.numeric(
            row[["mean_KNN_overlap"]]
          ),
          3
        )
      )
    )
  }
)

report_lines <- c(
  "GSE234774 STEP 10E MYELOID RECLUSTERING / HARMONY / RESOLUTION",
  "============================================================",
  paste0(
    "Run time: ",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  paste0(
    "Input object: ",
    input_object_path
  ),
  paste0(
    "Output directory: ",
    output_dir
  ),
  "",
  paste0(
    "FINAL DECISION: ",
    final_decision
  ),
  "",
  "FROZEN EXTRACTION",
  paste0(
    "Parent cluster column: ",
    PARENT_CLUSTER_COLUMN
  ),
  paste0(
    "Parent clusters retained: ",
    paste(
      PARENT_CLUSTERS_KEEP,
      collapse = ", "
    )
  ),
  paste0(
    "Myeloid nuclei: ",
    ncol(
      myeloid_obj
    )
  ),
  "",
  "PREPROCESSING",
  "Normalization: myeloid-only LogNormalize",
  paste0(
    "Variable features: ",
    length(
      variable_features
    )
  ),
  paste0(
    "Identity-focused PCA features: ",
    length(
      identity_pca_features
    )
  ),
  "Canonical S/G2M genes were scored but excluded from identity-focused PCA features.",
  "Cell-cycle scores were not regressed.",
  "Mitochondrial, ribosomal, and hemoglobin genes were excluded from identity-focused PCA features.",
  "",
  "DIMENSIONAL REDUCTION",
  paste0(
    "PCA dimensions computed: ",
    N_PCS_COMPUTE
  ),
  paste0(
    "Dimensions used: ",
    min(
      PCS_USE
    ),
    "-",
    max(
      PCS_USE
    )
  ),
  paste0(
    "Harmony: theta=",
    HARMONY_THETA,
    ", lambda=",
    HARMONY_LAMBDA,
    ", group=biological_library"
  ),
  paste0(
    "Mean Harmony-versus-PCA within-timepoint KNN overlap: ",
    round(
      mean_harmony_KNN_overlap,
      3
    )
  ),
  paste0(
    "Review threshold: mean overlap below ",
    MIN_ACCEPTABLE_MEAN_KNN_OVERLAP,
    " indicates strong structural alteration."
  ),
  "",
  "RESOLUTION SUMMARY",
  resolution_lines,
  "",
  "WITHIN-TIMEPOINT MIXING",
  mixing_lines,
  "",
  "INTERPRETATION BOUNDARIES",
  "1. Parent clusters 1 and 21 are a myeloid candidate compartment, not a clean resident-microglia object.",
  "2. Excluding canonical cell-cycle genes from PCA prevents proliferation from defining the primary identity axis; Phase remains available for interpretation.",
  "3. Harmony is used only for embeddings, neighbors, clustering, and visualization.",
  "4. Better library mixing alone is not sufficient to select the final myeloid reduction or resolution.",
  "5. Resolution selection must consider marker separation, cluster size, adjacent-resolution stability, library/timepoint dominance, and cycling composition.",
  "6. Activated resident microglia can lose P2ry12/Tmem119 and gain Apoe/Cst7/Lgals3; such cells must not be automatically called infiltrating macrophages.",
  "7. No final myeloid subtype, clean resident-microglia set, or Arm_1/Arm_2 result was frozen.",
  "",
  "PRIMARY COMPARISON OBJECT",
  output_object_path,
  "",
  "NEXT DECISION AFTER SEVEN-FILE REVIEW",
  "Select the myeloid reduction and working resolution.",
  "Then run marker-first myeloid subtype annotation to separate resident, activated/DAM-like, cycling, monocyte/macrophage-like, BAM-like, and ambiguous myeloid states."
)

if (
  final_decision ==
    "STOP"
) {
  failed <- audit_checks[
    audit_checks$required_for_pass &
      audit_checks$status !=
        "PASS",
    ,
    drop = FALSE
  ]

  report_lines <- c(
    report_lines,
    "",
    "STOP REASONS"
  )

  for (
    index in seq_len(
      nrow(
        failed
      )
    )
  ) {
    report_lines <- c(
      report_lines,
      paste0(
        index,
        ". ",
        failed$check[
          index
        ],
        " -- ",
        failed$detail[
          index
        ]
      )
    )
  }
}

report_path <- file.path(
  output_dir,
  "10_STEP10E_FINAL_REPORT.txt"
)

writeLines(
  report_lines,
  report_path,
  useBytes = TRUE
)

log_message(
  "Final decision: ",
  final_decision
)

log_message(
  "Final report: ",
  report_path
)

cat("\n")
cat("============================================================\n")
cat("GSE234774 STEP 10E COMPLETED\n")
cat(
  "FINAL DECISION:",
  final_decision,
  "\n"
)
cat(
  "Output:",
  output_dir,
  "\n"
)
cat("Primary comparison object:\n")
cat(
  "  ",
  output_object_path,
  "\n",
  sep = ""
)
cat("Upload exactly these seven files for review:\n")
cat("  1) figures/01_myeloid_PCA_Harmony_UMAP_comparison.pdf\n")
cat("  2) figures/02_myeloid_Harmony_resolution_0.2_to_1.0_UMAP.pdf\n")
cat("  3) figures/03_myeloid_marker_DotPlots_by_resolution.pdf\n")
cat("  4) tables/04_myeloid_resolution_cluster_summary_and_detail.csv\n")
cat("  5) tables/05_myeloid_within_timepoint_library_mixing.csv\n")
cat("  6) tables/09_STEP10E_audit_checks.csv\n")
cat("  7) 10_STEP10E_FINAL_REPORT.txt\n")
cat("============================================================\n")

if (
  final_decision ==
    "STOP"
) {
  stop(
    "GSE234774 Step 10E decision is STOP. Review the final report.",
    call. = FALSE
  )
}
