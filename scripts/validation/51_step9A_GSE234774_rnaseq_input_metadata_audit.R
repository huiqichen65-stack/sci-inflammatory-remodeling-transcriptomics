# =============================================================================
# 51_step9A_GSE234774_rnaseq_input_metadata_audit.R
#
# Adult SCI activated/DAM-like microglia project
# Step 9A (redefined): GSE234774 long-term snRNA-seq input and metadata audit
#
# This script ONLY audits the local processed inputs and metadata.
# It does NOT create a Seurat object, perform QC, integration, annotation,
# signature scoring, label transfer, pseudotime, or spatial analysis.
#
# Main goals
#   1. Confirm the processed snRNA-seq matrix, barcode, feature, and metadata
#      files are present and structurally consistent.
#   2. Recover the metadata columns needed to define the adult standard-SCI
#      validation cohort: biological sample, timepoint, age, injury/model,
#      severity, treatment, and cell annotation.
#   3. Confirm that metadata rows correspond to expression-matrix barcodes.
#   4. Audit the availability and dimensions of the 2D/3D spatial inputs,
#      without reconstructing spatial objects.
#   5. Return a strict PASS or STOP decision for entry into Step 9B.
#
# Project root
#   <PROJECT_ROOT>
#
# Input directory
#   01_raw_data/GSE234774
#
# Output directory
#   04_results/step9A_GSE234774_rnaseq_input_metadata_audit
#
# No internet download is performed.
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

# -----------------------------------------------------------------------------
# 0. Fixed paths
# -----------------------------------------------------------------------------
project_root <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
input_dir <- file.path(project_root, "01_raw_data", "GSE234774")
output_dir <- file.path(
  project_root,
  "04_results",
  "step9A_GSE234774_rnaseq_input_metadata_audit"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

run_log <- file.path(output_dir, "00_run_log.txt")
if (file.exists(run_log)) {
  file.remove(run_log)
}

log_message <- function(...) {
  msg <- paste0(...)
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", stamp, "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = run_log, append = TRUE)
}

write_csv_safe <- function(x, filename) {
  out <- file.path(output_dir, filename)
  utils::write.csv(x, out, row.names = FALSE, na = "")
  log_message("Wrote: ", out)
  invisible(out)
}

trim_na <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

first_nonempty <- function(x) {
  x <- trim_na(x)
  x <- x[!is.na(x)]
  if (length(x) == 0L) {
    return(NA_character_)
  }
  x[[1L]]
}

collapse_unique <- function(x, sep = " | ", max_values = 30L) {
  x <- trim_na(x)
  x <- unique(x[!is.na(x)])
  if (length(x) == 0L) {
    return(NA_character_)
  }
  if (length(x) > max_values) {
    x <- c(x[seq_len(max_values)], paste0("... +", length(x) - max_values, " more"))
  }
  paste(x, collapse = sep)
}

normalise_name <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

safe_fraction <- function(x) {
  if (length(x) == 0L) {
    return(0)
  }
  mean(!is.na(trim_na(x)))
}

open_text_connection <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
}

read_text_lines <- function(path, n = -1L) {
  con <- open_text_connection(path)
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  readLines(con, n = n, warn = FALSE, encoding = "UTF-8")
}

count_text_lines <- function(path, chunk_size = 100000L) {
  con <- open_text_connection(path)
  on.exit(try(close(con), silent = TRUE), add = TRUE)

  total <- 0L
  repeat {
    block <- readLines(con, n = chunk_size, warn = FALSE)
    n_block <- length(block)
    total <- total + n_block
    if (n_block < chunk_size) {
      break
    }
  }
  total
}

probe_file <- function(path) {
  if (!file.exists(path)) {
    return(list(ok = FALSE, detail = "file does not exist"))
  }

  info <- file.info(path)
  if (is.na(info$size) || info$size <= 0) {
    return(list(ok = FALSE, detail = "file size is zero or unavailable"))
  }

  if (file.access(path, mode = 4L) != 0L) {
    return(list(ok = FALSE, detail = "file is not readable"))
  }

  con <- NULL
  result <- tryCatch(
    {
      if (grepl("\\.gz$", path, ignore.case = TRUE)) {
        con <- gzfile(path, open = "rb")
      } else {
        con <- file(path, open = "rb")
      }
      x <- readBin(con, what = "raw", n = 32L)
      list(
        ok = length(x) > 0L,
        detail = paste0("readable; first_bytes=", length(x))
      )
    },
    error = function(e) {
      list(ok = FALSE, detail = conditionMessage(e))
    },
    finally = {
      if (!is.null(con)) {
        try(close(con), silent = TRUE)
      }
    }
  )
  result
}

read_mtx_header <- function(path) {
  result <- list(
    valid = FALSE,
    matrix_market_header = NA_character_,
    n_features = NA_real_,
    n_barcodes = NA_real_,
    n_nonzero = NA_real_,
    detail = NA_character_
  )

  if (!file.exists(path)) {
    result$detail <- "matrix file does not exist"
    return(result)
  }

  con <- open_text_connection(path)
  on.exit(try(close(con), silent = TRUE), add = TRUE)

  first_line <- tryCatch(
    readLines(con, n = 1L, warn = FALSE),
    error = function(e) character(0)
  )

  if (length(first_line) != 1L) {
    result$detail <- "could not read Matrix Market header"
    return(result)
  }

  result$matrix_market_header <- first_line[[1L]]

  if (!grepl("^%%MatrixMarket\\s+matrix", first_line[[1L]], ignore.case = TRUE)) {
    result$detail <- paste0(
      "first line is not a Matrix Market header: ",
      first_line[[1L]]
    )
    return(result)
  }

  dimension_line <- NA_character_

  repeat {
    x <- readLines(con, n = 1L, warn = FALSE)
    if (length(x) == 0L) {
      break
    }
    x <- trimws(x[[1L]])
    if (x == "" || startsWith(x, "%")) {
      next
    }
    dimension_line <- x
    break
  }

  if (is.na(dimension_line)) {
    result$detail <- "Matrix Market dimension line was not found"
    return(result)
  }

  parts <- strsplit(dimension_line, "\\s+", perl = TRUE)[[1L]]
  numbers <- suppressWarnings(as.numeric(parts))

  if (length(numbers) < 3L || any(is.na(numbers[1:3]))) {
    result$detail <- paste0("invalid dimension line: ", dimension_line)
    return(result)
  }

  result$valid <- TRUE
  result$n_features <- numbers[[1L]]
  result$n_barcodes <- numbers[[2L]]
  result$n_nonzero <- numbers[[3L]]
  result$detail <- "Matrix Market header and dimensions parsed successfully"
  result
}

read_delimited_gz <- function(path, nrows = -1L) {
  con <- open_text_connection(path)
  on.exit(try(close(con), silent = TRUE), add = TRUE)

  x <- utils::read.delim(
    con,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    nrows = nrows
  )

  names(x) <- make.unique(names(x), sep = "_dup")
  x
}

read_feature_table <- function(path) {
  con <- open_text_connection(path)
  on.exit(try(close(con), silent = TRUE), add = TRUE)

  x <- utils::read.delim(
    con,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (ncol(x) == 1L) {
    names(x) <- "feature_id"
  } else if (ncol(x) == 2L) {
    names(x) <- c("feature_id", "feature_name")
  } else {
    names(x)[1:3] <- c("feature_id", "feature_name", "feature_type")
  }

  x
}

format_examples <- function(x, n = 8L) {
  x <- unique(trim_na(x))
  x <- x[!is.na(x)]
  if (length(x) == 0L) {
    return(NA_character_)
  }
  paste(head(x, n), collapse = " | ")
}

summarise_metadata_columns <- function(meta) {
  if (ncol(meta) == 0L) {
    return(data.frame())
  }

  out <- lapply(seq_along(meta), function(i) {
    x <- meta[[i]]
    x_chr <- trim_na(x)
    nonmissing <- x_chr[!is.na(x_chr)]
    n_unique <- length(unique(nonmissing))

    data.frame(
      column_index = i,
      column_name = names(meta)[[i]],
      normalised_name = normalise_name(names(meta)[[i]]),
      class = paste(class(x), collapse = "/"),
      n_rows = length(x),
      n_nonmissing = length(nonmissing),
      missing_fraction = if (length(x) == 0L) NA_real_ else mean(is.na(x_chr)),
      n_unique_nonmissing = n_unique,
      unique_fraction_nonmissing = if (length(nonmissing) == 0L) {
        NA_real_
      } else {
        n_unique / length(nonmissing)
      },
      example_values = format_examples(nonmissing),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

field_patterns <- list(
  barcode = c(
    "^barcode$", "cell_barcode", "cell_id", "^cell$", "nucleus_id",
    "barcode_id", "cb"
  ),
  sample = c(
    "^sample$", "^sample_id$", "sample_name", "orig_ident", "orig.ident",
    "library", "library_id", "replicate", "mouse", "animal", "subject",
    "donor", "biological_sample", "specimen"
  ),
  timepoint = c(
    "^time$", "timepoint", "time_point", "dpi", "days_post", "day_post",
    "post_injury", "postinjury", "week", "month", "stage"
  ),
  age = c(
    "^age$", "age_group", "animal_age", "mouse_age", "young", "adult",
    "aged", "old"
  ),
  injury = c(
    "^injury$", "injury_model", "injury_type", "condition", "disease",
    "sci", "contusion", "compression", "lesion", "sham"
  ),
  severity = c(
    "severity", "injury_severity", "mild", "moderate", "severe",
    "impact_force", "force"
  ),
  treatment = c(
    "^treatment$", "treatment_group", "drug", "therapy", "intervention",
    "vehicle", "compound", "dose"
  ),
  sex = c(
    "^sex$", "^gender$"
  ),
  celltype = c(
    "cell_type", "celltype", "cell_annotation", "annotation",
    "broad_type", "major_type", "cluster", "class", "subtype"
  ),
  region = c(
    "region", "location", "anatomical", "lesion_site", "spinal_level"
  )
)

value_regex <- list(
  timepoint = paste0(
    "(^|[^a-z0-9])(",
    "uninjured|naive|sham|intact|",
    "[0-9]+\\s*(h|hr|hour|dpi|d|day|w|wk|week|month|mo)",
    ")([^a-z0-9]|$)"
  ),
  age = "(young|adult|aged|old|juvenile|[0-9]+\\s*(week|wk|month|mo))",
  injury = "(spinal cord|sci|contusion|compression|injur|lesion|sham|control)",
  severity = "(mild|moderate|severe|severity)",
  treatment = "(vehicle|saline|drug|treat|therapy|dose|compound|control)",
  sex = "(^|[^a-z])(male|female|m|f)([^a-z]|$)",
  celltype = paste0(
    "(microglia|macrophage|monocyte|astrocyte|oligodendro|",
    "endothelial|fibroblast|pericyte|neuron|t cell|b cell|myeloid)"
  )
)

score_column_for_field <- function(meta, column_name, field) {
  x <- trim_na(meta[[column_name]])
  x_nonmissing <- x[!is.na(x)]
  n <- length(x)
  n_nonmissing <- length(x_nonmissing)
  n_unique <- length(unique(x_nonmissing))
  name_norm <- normalise_name(column_name)

  patterns <- field_patterns[[field]]
  name_score <- 0

  if (!is.null(patterns)) {
    exact_tokens <- gsub("^\\^|\\$$", "", patterns)
    exact_hits <- name_norm %in% exact_tokens
    contain_hits <- vapply(
      patterns,
      function(p) grepl(p, name_norm, perl = TRUE, ignore.case = TRUE),
      logical(1)
    )

    if (any(exact_hits)) {
      name_score <- 100
    } else if (any(contain_hits)) {
      name_score <- 60
    }
  }

  value_score <- 0
  regex <- value_regex[[field]]

  if (!is.null(regex) && n_nonmissing > 0L) {
    sample_values <- unique(x_nonmissing)
    if (length(sample_values) > 5000L) {
      sample_values <- sample_values[seq_len(5000L)]
    }

    hit_fraction <- mean(
      grepl(regex, tolower(sample_values), perl = TRUE, ignore.case = TRUE)
    )
    value_score <- 40 * hit_fraction
  }

  structure_score <- 0

  if (field == "sample" && n_nonmissing > 0L) {
    repeated_fraction <- 1 - (n_unique / n_nonmissing)
    reasonable_cardinality <- n_unique >= 2L &&
      n_unique <= max(500L, floor(n_nonmissing / 5L))

    if (reasonable_cardinality) {
      structure_score <- 25 * max(0, repeated_fraction)
    } else if (n_unique == n_nonmissing) {
      structure_score <- -30
    }
  }

  if (field == "barcode" && n_nonmissing > 0L) {
    if (n_unique == n_nonmissing) {
      structure_score <- 20
    }
  }

  data.frame(
    field = field,
    column_name = column_name,
    name_score = name_score,
    value_score = value_score,
    structure_score = structure_score,
    total_score = name_score + value_score + structure_score,
    n_nonmissing = n_nonmissing,
    n_unique = n_unique,
    example_values = format_examples(x_nonmissing),
    stringsAsFactors = FALSE
  )
}

detect_metadata_fields <- function(meta) {
  fields <- names(field_patterns)
  all_scores <- list()
  counter <- 0L

  for (field in fields) {
    for (column_name in names(meta)) {
      counter <- counter + 1L
      all_scores[[counter]] <- score_column_for_field(
        meta = meta,
        column_name = column_name,
        field = field
      )
    }
  }

  score_table <- do.call(rbind, all_scores)

  selected <- lapply(fields, function(field) {
    z <- score_table[score_table$field == field, , drop = FALSE]
    z <- z[order(-z$total_score, -z$name_score, z$column_name), , drop = FALSE]
    best <- z[1L, , drop = FALSE]

    minimum_score <- if (field %in% c("sample", "barcode")) 50 else 35
    selected_column <- if (
      nrow(best) == 1L &&
      is.finite(best$total_score[[1L]]) &&
      best$total_score[[1L]] >= minimum_score
    ) {
      best$column_name[[1L]]
    } else {
      NA_character_
    }

    data.frame(
      field = field,
      selected_column = selected_column,
      selected_score = best$total_score[[1L]],
      second_best_column = if (nrow(z) >= 2L) z$column_name[[2L]] else NA_character_,
      second_best_score = if (nrow(z) >= 2L) z$total_score[[2L]] else NA_real_,
      confidence = if (is.na(selected_column)) {
        "UNRESOLVED"
      } else if (
        nrow(z) == 1L ||
        (best$total_score[[1L]] - z$total_score[[2L]]) >= 20
      ) {
        "HIGH"
      } else {
        "REVIEW"
      },
      stringsAsFactors = FALSE
    )
  })

  list(
    selected = do.call(rbind, selected),
    scores = score_table
  )
}

get_selected_column <- function(field_detection, field) {
  x <- field_detection$selected
  z <- x$selected_column[x$field == field]
  first_nonempty(z)
}

summarise_selected_values <- function(meta, field_detection) {
  rows <- list()
  counter <- 0L

  for (field in field_detection$selected$field) {
    col <- get_selected_column(field_detection, field)

    if (is.na(col)) {
      counter <- counter + 1L
      rows[[counter]] <- data.frame(
        field = field,
        selected_column = NA_character_,
        value = NA_character_,
        n_rows = NA_integer_,
        proportion = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    x <- trim_na(meta[[col]])
    x[is.na(x)] <- "NA"
    tab <- sort(table(x), decreasing = TRUE)

    for (i in seq_along(tab)) {
      counter <- counter + 1L
      rows[[counter]] <- data.frame(
        field = field,
        selected_column = col,
        value = names(tab)[[i]],
        n_rows = as.integer(tab[[i]]),
        proportion = as.integer(tab[[i]]) / length(x),
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, rows)
}

make_sample_manifest <- function(meta, field_detection) {
  sample_col <- get_selected_column(field_detection, "sample")

  if (is.na(sample_col)) {
    return(data.frame(
      sample_id = character(0),
      n_cells_or_nuclei = integer(0),
      stringsAsFactors = FALSE
    ))
  }

  sample_id <- trim_na(meta[[sample_col]])
  sample_id[is.na(sample_id)] <- "UNRESOLVED_SAMPLE"

  fields_to_include <- c(
    "timepoint", "age", "injury", "severity",
    "treatment", "sex", "region", "celltype"
  )

  split_idx <- split(seq_len(nrow(meta)), sample_id)

  out <- lapply(names(split_idx), function(s) {
    idx <- split_idx[[s]]
    row <- list(
      sample_id = s,
      n_cells_or_nuclei = length(idx)
    )

    for (field in fields_to_include) {
      col <- get_selected_column(field_detection, field)
      value <- if (is.na(col)) {
        NA_character_
      } else {
        collapse_unique(meta[[col]][idx], max_values = 20L)
      }
      row[[paste0(field, "_values")]] <- value
    }

    as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
  })

  result <- do.call(rbind, out)
  rownames(result) <- NULL
  result
}

make_sample_timepoint_counts <- function(meta, field_detection) {
  sample_col <- get_selected_column(field_detection, "sample")
  time_col <- get_selected_column(field_detection, "timepoint")

  if (is.na(sample_col) || is.na(time_col)) {
    return(data.frame(
      sample_id = character(0),
      timepoint = character(0),
      n_cells_or_nuclei = integer(0),
      stringsAsFactors = FALSE
    ))
  }

  sample_id <- trim_na(meta[[sample_col]])
  timepoint <- trim_na(meta[[time_col]])
  sample_id[is.na(sample_id)] <- "UNRESOLVED_SAMPLE"
  timepoint[is.na(timepoint)] <- "UNRESOLVED_TIMEPOINT"

  tab <- as.data.frame(
    table(sample_id = sample_id, timepoint = timepoint),
    stringsAsFactors = FALSE
  )
  tab <- tab[tab$Freq > 0L, , drop = FALSE]
  names(tab)[names(tab) == "Freq"] <- "n_cells_or_nuclei"
  tab[order(tab$timepoint, tab$sample_id), , drop = FALSE]
}

audit_barcode_mapping <- function(meta, barcodes) {
  if (nrow(meta) == 0L || length(barcodes) == 0L) {
    return(list(
      summary = data.frame(),
      best_column = NA_character_,
      best_overlap_fraction = 0
    ))
  }

  rows <- lapply(names(meta), function(col) {
    x <- trim_na(meta[[col]])
    x_nonmissing <- x[!is.na(x)]

    exact_order_match <- length(x) == length(barcodes) &&
      all(!is.na(x)) &&
      identical(as.character(x), as.character(barcodes))

    metadata_in_barcode_fraction <- if (length(x_nonmissing) == 0L) {
      0
    } else {
      mean(x_nonmissing %in% barcodes)
    }

    barcode_in_metadata_fraction <- if (length(x_nonmissing) == 0L) {
      0
    } else {
      mean(barcodes %in% x_nonmissing)
    }

    data.frame(
      metadata_column = col,
      n_nonmissing = length(x_nonmissing),
      n_unique_nonmissing = length(unique(x_nonmissing)),
      exact_order_match = exact_order_match,
      metadata_values_in_barcode_fraction = metadata_in_barcode_fraction,
      barcodes_in_metadata_values_fraction = barcode_in_metadata_fraction,
      mean_bidirectional_overlap = mean(c(
        metadata_in_barcode_fraction,
        barcode_in_metadata_fraction
      )),
      stringsAsFactors = FALSE
    )
  })

  summary <- do.call(rbind, rows)
  summary <- summary[
    order(
      -summary$exact_order_match,
      -summary$mean_bidirectional_overlap,
      summary$metadata_column
    ),
    ,
    drop = FALSE
  ]

  best_column <- summary$metadata_column[[1L]]
  best_overlap <- summary$mean_bidirectional_overlap[[1L]]

  list(
    summary = summary,
    best_column = best_column,
    best_overlap_fraction = best_overlap
  )
}

audit_processed_modality <- function(
  modality,
  matrix_path,
  barcode_path,
  feature_path,
  metadata_path,
  read_full_metadata = FALSE
) {
  mtx <- read_mtx_header(matrix_path)

  barcode_count <- if (file.exists(barcode_path)) {
    count_text_lines(barcode_path)
  } else {
    NA_integer_
  }

  feature_count <- if (file.exists(feature_path)) {
    count_text_lines(feature_path)
  } else {
    NA_integer_
  }

  metadata_rows <- NA_integer_
  metadata_columns <- NA_integer_
  metadata_readable <- FALSE
  metadata_error <- NA_character_

  if (file.exists(metadata_path)) {
    meta_probe <- tryCatch(
      {
        if (isTRUE(read_full_metadata)) {
          z <- read_delimited_gz(metadata_path)
        } else {
          z <- read_delimited_gz(metadata_path, nrows = 50L)
        }
        metadata_rows <- if (isTRUE(read_full_metadata)) {
          nrow(z)
        } else {
          count_text_lines(metadata_path) - 1L
        }
        metadata_columns <- ncol(z)
        metadata_readable <- TRUE
        TRUE
      },
      error = function(e) {
        metadata_error <<- conditionMessage(e)
        FALSE
      }
    )
    invisible(meta_probe)
  }

  data.frame(
    modality = modality,
    matrix_file = basename(matrix_path),
    barcode_file = basename(barcode_path),
    feature_file = basename(feature_path),
    metadata_file = basename(metadata_path),
    matrix_header_valid = isTRUE(mtx$valid),
    matrix_n_features = mtx$n_features,
    matrix_n_barcodes = mtx$n_barcodes,
    matrix_n_nonzero = mtx$n_nonzero,
    barcode_line_count = barcode_count,
    feature_line_count = feature_count,
    metadata_row_count = metadata_rows,
    metadata_column_count = metadata_columns,
    metadata_readable = metadata_readable,
    feature_dimension_match = isTRUE(mtx$valid) &&
      !is.na(feature_count) &&
      identical(as.numeric(feature_count), as.numeric(mtx$n_features)),
    barcode_dimension_match = isTRUE(mtx$valid) &&
      !is.na(barcode_count) &&
      identical(as.numeric(barcode_count), as.numeric(mtx$n_barcodes)),
    metadata_row_match = isTRUE(mtx$valid) &&
      !is.na(metadata_rows) &&
      identical(as.numeric(metadata_rows), as.numeric(mtx$n_barcodes)),
    detail = first_nonempty(c(metadata_error, mtx$detail)),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# 1. Input directory and expected files
# -----------------------------------------------------------------------------
log_message("Project root: ", project_root)
log_message("Input directory: ", input_dir)
log_message("Output directory: ", output_dir)

if (!dir.exists(input_dir)) {
  report <- c(
    "GSE234774 STEP 9A AUDIT",
    "FINAL DECISION: STOP",
    "",
    paste0("Input directory does not exist: ", input_dir)
  )
  writeLines(report, file.path(output_dir, "11_STEP9A_FINAL_REPORT.txt"))
  stop(paste(report, collapse = "\n"), call. = FALSE)
}

expected_files <- c(
  "GSE234774_RAW.tar",
  "GSE234774_rnaseq_barcodes.txt.gz",
  "GSE234774_rnaseq_features.txt.gz",
  "GSE234774_rnaseq_filtered_scRNA.mtx.gz",
  "GSE234774_rnaseq_meta.txt.gz",
  "GSE234774_spatial_2d_barcodes.txt.gz",
  "GSE234774_spatial_2d_features.txt.gz",
  "GSE234774_spatial_2d_filtered_spatial_2d.mtx.gz",
  "GSE234774_spatial_2d_meta.txt.gz",
  "GSE234774_spatial_3d_barcodes.txt.gz",
  "GSE234774_spatial_3d_features.txt.gz",
  "GSE234774_spatial_3d_filtered_spatial_3d.mtx.gz",
  "GSE234774_spatial_3d_meta.txt.gz"
)

expected_paths <- file.path(input_dir, expected_files)
probe_results <- lapply(expected_paths, probe_file)
file_info <- file.info(expected_paths)

inventory <- data.frame(
  file_name = expected_files,
  full_path = normalizePath(
    expected_paths,
    winslash = "/",
    mustWork = FALSE
  ),
  exists = file.exists(expected_paths),
  size_bytes = as.numeric(file_info$size),
  size_MB = round(as.numeric(file_info$size) / 1024^2, 3),
  readable = vapply(probe_results, function(z) isTRUE(z$ok), logical(1)),
  readability_detail = vapply(
    probe_results,
    function(z) z$detail,
    character(1)
  ),
  stringsAsFactors = FALSE
)

write_csv_safe(inventory, "01_input_file_inventory.csv")

# -----------------------------------------------------------------------------
# 2. Main snRNA-seq matrix / barcode / feature dimension audit
# -----------------------------------------------------------------------------
rnaseq_matrix <- file.path(
  input_dir,
  "GSE234774_rnaseq_filtered_scRNA.mtx.gz"
)
rnaseq_barcodes_path <- file.path(
  input_dir,
  "GSE234774_rnaseq_barcodes.txt.gz"
)
rnaseq_features_path <- file.path(
  input_dir,
  "GSE234774_rnaseq_features.txt.gz"
)
rnaseq_meta_path <- file.path(
  input_dir,
  "GSE234774_rnaseq_meta.txt.gz"
)

log_message("Reading snRNA-seq Matrix Market header only; full matrix is not loaded.")
rnaseq_dimension_audit <- audit_processed_modality(
  modality = "rnaseq",
  matrix_path = rnaseq_matrix,
  barcode_path = rnaseq_barcodes_path,
  feature_path = rnaseq_features_path,
  metadata_path = rnaseq_meta_path,
  read_full_metadata = FALSE
)

write_csv_safe(
  rnaseq_dimension_audit,
  "02_RNAseq_dimension_consistency.csv"
)

# -----------------------------------------------------------------------------
# 3. Read barcodes, features, and the full snRNA-seq metadata
# -----------------------------------------------------------------------------
log_message("Reading snRNA-seq barcodes.")
rnaseq_barcodes <- tryCatch(
  read_text_lines(rnaseq_barcodes_path),
  error = function(e) {
    log_message("Barcode read error: ", conditionMessage(e))
    character(0)
  }
)
rnaseq_barcodes <- trimws(rnaseq_barcodes)
rnaseq_barcodes <- rnaseq_barcodes[rnaseq_barcodes != ""]

log_message("Reading snRNA-seq feature table.")
rnaseq_features <- tryCatch(
  read_feature_table(rnaseq_features_path),
  error = function(e) {
    log_message("Feature read error: ", conditionMessage(e))
    data.frame()
  }
)

log_message("Reading complete snRNA-seq metadata.")
rnaseq_meta <- tryCatch(
  read_delimited_gz(rnaseq_meta_path),
  error = function(e) {
    log_message("Metadata read error: ", conditionMessage(e))
    data.frame()
  }
)

if (nrow(rnaseq_meta) > 0L) {
  metadata_preview <- head(rnaseq_meta, 50L)
} else {
  metadata_preview <- data.frame()
}
write_csv_safe(metadata_preview, "03_RNAseq_metadata_preview50.csv")

# -----------------------------------------------------------------------------
# 4. Basic barcode and feature integrity
# -----------------------------------------------------------------------------
feature_id <- if (
  nrow(rnaseq_features) > 0L &&
  "feature_id" %in% names(rnaseq_features)
) {
  trim_na(rnaseq_features$feature_id)
} else {
  character(0)
}

feature_name <- if (
  nrow(rnaseq_features) > 0L &&
  "feature_name" %in% names(rnaseq_features)
) {
  trim_na(rnaseq_features$feature_name)
} else {
  feature_id
}

basic_integrity <- data.frame(
  item = c(
    "RNAseq barcode count",
    "RNAseq duplicated barcodes",
    "RNAseq feature count",
    "RNAseq duplicated feature IDs",
    "RNAseq duplicated feature names",
    "RNAseq metadata row count",
    "RNAseq metadata column count"
  ),
  value = c(
    length(rnaseq_barcodes),
    sum(duplicated(rnaseq_barcodes)),
    nrow(rnaseq_features),
    sum(duplicated(feature_id[!is.na(feature_id)])),
    sum(duplicated(feature_name[!is.na(feature_name)])),
    nrow(rnaseq_meta),
    ncol(rnaseq_meta)
  ),
  stringsAsFactors = FALSE
)

write_csv_safe(basic_integrity, "04_RNAseq_basic_integrity.csv")

# -----------------------------------------------------------------------------
# 5. Metadata column audit and automatic field detection
# -----------------------------------------------------------------------------
metadata_column_audit <- summarise_metadata_columns(rnaseq_meta)
write_csv_safe(
  metadata_column_audit,
  "05_RNAseq_metadata_column_audit.csv"
)

if (ncol(rnaseq_meta) > 0L) {
  field_detection <- detect_metadata_fields(rnaseq_meta)
} else {
  field_detection <- list(
    selected = data.frame(
      field = names(field_patterns),
      selected_column = NA_character_,
      selected_score = NA_real_,
      second_best_column = NA_character_,
      second_best_score = NA_real_,
      confidence = "UNRESOLVED",
      stringsAsFactors = FALSE
    ),
    scores = data.frame()
  )
}

write_csv_safe(
  field_detection$selected,
  "06_RNAseq_selected_metadata_fields.csv"
)

write_csv_safe(
  field_detection$scores,
  "06B_RNAseq_metadata_field_candidate_scores.csv"
)

metadata_value_summary <- if (nrow(rnaseq_meta) > 0L) {
  summarise_selected_values(rnaseq_meta, field_detection)
} else {
  data.frame()
}

write_csv_safe(
  metadata_value_summary,
  "07_RNAseq_selected_field_value_counts.csv"
)

# -----------------------------------------------------------------------------
# 6. Barcode-to-metadata correspondence
# -----------------------------------------------------------------------------
barcode_mapping <- audit_barcode_mapping(
  meta = rnaseq_meta,
  barcodes = rnaseq_barcodes
)

write_csv_safe(
  barcode_mapping$summary,
  "08_RNAseq_barcode_metadata_mapping.csv"
)

# -----------------------------------------------------------------------------
# 7. Automatic sample manifest and timepoint structure
# -----------------------------------------------------------------------------
sample_manifest <- if (nrow(rnaseq_meta) > 0L) {
  make_sample_manifest(rnaseq_meta, field_detection)
} else {
  data.frame()
}

write_csv_safe(
  sample_manifest,
  "09_RNAseq_sample_manifest_auto.csv"
)

sample_timepoint_counts <- if (nrow(rnaseq_meta) > 0L) {
  make_sample_timepoint_counts(rnaseq_meta, field_detection)
} else {
  data.frame()
}

write_csv_safe(
  sample_timepoint_counts,
  "09B_RNAseq_sample_timepoint_counts.csv"
)

# -----------------------------------------------------------------------------
# 8. Spatial 2D and 3D processed input audit
# -----------------------------------------------------------------------------
log_message("Auditing spatial 2D and 3D processed inputs without object reconstruction.")

spatial_2d_audit <- audit_processed_modality(
  modality = "spatial_2d",
  matrix_path = file.path(
    input_dir,
    "GSE234774_spatial_2d_filtered_spatial_2d.mtx.gz"
  ),
  barcode_path = file.path(
    input_dir,
    "GSE234774_spatial_2d_barcodes.txt.gz"
  ),
  feature_path = file.path(
    input_dir,
    "GSE234774_spatial_2d_features.txt.gz"
  ),
  metadata_path = file.path(
    input_dir,
    "GSE234774_spatial_2d_meta.txt.gz"
  ),
  read_full_metadata = FALSE
)

spatial_3d_audit <- audit_processed_modality(
  modality = "spatial_3d",
  matrix_path = file.path(
    input_dir,
    "GSE234774_spatial_3d_filtered_spatial_3d.mtx.gz"
  ),
  barcode_path = file.path(
    input_dir,
    "GSE234774_spatial_3d_barcodes.txt.gz"
  ),
  feature_path = file.path(
    input_dir,
    "GSE234774_spatial_3d_features.txt.gz"
  ),
  metadata_path = file.path(
    input_dir,
    "GSE234774_spatial_3d_meta.txt.gz"
  ),
  read_full_metadata = FALSE
)

spatial_audit <- rbind(spatial_2d_audit, spatial_3d_audit)
write_csv_safe(
  spatial_audit,
  "10_spatial_processed_input_audit.csv"
)

# -----------------------------------------------------------------------------
# 9. RAW.tar member audit
# -----------------------------------------------------------------------------
raw_tar <- file.path(input_dir, "GSE234774_RAW.tar")

raw_archive_summary <- data.frame(
  archive_file = basename(raw_tar),
  exists = file.exists(raw_tar),
  listable = FALSE,
  n_members = 0L,
  n_spatial_members = 0L,
  n_image_like_members = 0L,
  n_position_like_members = 0L,
  member_examples = NA_character_,
  detail = NA_character_,
  stringsAsFactors = FALSE
)

if (file.exists(raw_tar)) {
  raw_members <- tryCatch(
    utils::untar(raw_tar, list = TRUE),
    error = function(e) e
  )

  if (inherits(raw_members, "error")) {
    raw_archive_summary$detail <- conditionMessage(raw_members)
  } else {
    raw_archive_summary$listable <- length(raw_members) > 0L
    raw_archive_summary$n_members <- length(raw_members)
    raw_archive_summary$n_spatial_members <- sum(
      grepl("spatial|visium|section|slice", raw_members, ignore.case = TRUE)
    )
    raw_archive_summary$n_image_like_members <- sum(
      grepl("\\.(png|jpg|jpeg|tif|tiff)$", raw_members, ignore.case = TRUE)
    )
    raw_archive_summary$n_position_like_members <- sum(
      grepl(
        "position|tissue_positions|scalefactor|coordinate",
        raw_members,
        ignore.case = TRUE
      )
    )
    raw_archive_summary$member_examples <- collapse_unique(
      head(raw_members, 30L),
      sep = " | ",
      max_values = 30L
    )
    raw_archive_summary$detail <- "archive listed successfully"
  }
}

write_csv_safe(
  raw_archive_summary,
  "10B_RAW_archive_summary.csv"
)

# -----------------------------------------------------------------------------
# 10. Strict decision logic
# -----------------------------------------------------------------------------
required_rnaseq_files <- c(
  "GSE234774_rnaseq_barcodes.txt.gz",
  "GSE234774_rnaseq_features.txt.gz",
  "GSE234774_rnaseq_filtered_scRNA.mtx.gz",
  "GSE234774_rnaseq_meta.txt.gz"
)

required_file_rows <- inventory$file_name %in% required_rnaseq_files
required_files_pass <- all(
  inventory$exists[required_file_rows] &
    inventory$readable[required_file_rows]
)

dimension_pass <- nrow(rnaseq_dimension_audit) == 1L &&
  isTRUE(rnaseq_dimension_audit$matrix_header_valid[[1L]]) &&
  isTRUE(rnaseq_dimension_audit$feature_dimension_match[[1L]]) &&
  isTRUE(rnaseq_dimension_audit$barcode_dimension_match[[1L]])

metadata_read_pass <- nrow(rnaseq_meta) > 0L && ncol(rnaseq_meta) > 0L

metadata_row_dimension_pass <- isTRUE(
  nrow(rnaseq_meta) == rnaseq_dimension_audit$matrix_n_barcodes[[1L]]
)

barcode_overlap_pass <- is.finite(
  barcode_mapping$best_overlap_fraction
) &&
  barcode_mapping$best_overlap_fraction >= 0.99

sample_col <- get_selected_column(field_detection, "sample")
timepoint_col <- get_selected_column(field_detection, "timepoint")
age_col <- get_selected_column(field_detection, "age")
injury_col <- get_selected_column(field_detection, "injury")
severity_col <- get_selected_column(field_detection, "severity")
treatment_col <- get_selected_column(field_detection, "treatment")
celltype_col <- get_selected_column(field_detection, "celltype")

sample_pass <- !is.na(sample_col) &&
  length(unique(trim_na(rnaseq_meta[[sample_col]]))) >= 2L

timepoint_pass <- !is.na(timepoint_col) &&
  length(unique(trim_na(rnaseq_meta[[timepoint_col]]))) >= 2L

adult_cohort_field_pass <- !is.na(age_col)
injury_field_pass <- !is.na(injury_col)
condition_stratification_pass <- !is.na(severity_col) || !is.na(treatment_col)

celltype_field_available <- !is.na(celltype_col)

spatial_2d_rows <- spatial_audit$modality == "spatial_2d"
spatial_2d_technical_pass <- isTRUE(
  nrow(spatial_audit) >= 1L &&
    any(spatial_2d_rows) &&
    all(c(
      spatial_audit$matrix_header_valid[spatial_2d_rows],
      spatial_audit$feature_dimension_match[spatial_2d_rows],
      spatial_audit$barcode_dimension_match[spatial_2d_rows]
    ), na.rm = FALSE)
)

raw_archive_pass <- isTRUE(raw_archive_summary$exists[[1L]]) &&
  isTRUE(raw_archive_summary$listable[[1L]])

audit_checks <- data.frame(
  check = c(
    "Required snRNA-seq processed files exist and are readable",
    "snRNA-seq Matrix Market dimensions match feature and barcode files",
    "snRNA-seq metadata is readable",
    "snRNA-seq metadata row count matches matrix barcode dimension",
    "A metadata column maps to >=99% of expression barcodes",
    "Biological sample column is recoverable with >=2 sample levels",
    "Timepoint column is recoverable with >=2 time levels",
    "Age field is recoverable for adult-cohort selection",
    "Injury/model field is recoverable for SCI-cohort selection",
    "Severity or treatment field is recoverable for standard-cohort selection",
    "Cell-type/annotation field is available for later microglia review",
    "Spatial 2D processed matrix/barcode/feature structure is technically consistent",
    "RAW.tar is present and listable"
  ),
  required_for_step9B = c(
    TRUE, TRUE, TRUE, TRUE, TRUE,
    TRUE, TRUE, TRUE, TRUE, TRUE,
    FALSE, FALSE, FALSE
  ),
  status = c(
    ifelse(required_files_pass, "PASS", "STOP"),
    ifelse(dimension_pass, "PASS", "STOP"),
    ifelse(metadata_read_pass, "PASS", "STOP"),
    ifelse(metadata_row_dimension_pass, "PASS", "STOP"),
    ifelse(barcode_overlap_pass, "PASS", "STOP"),
    ifelse(sample_pass, "PASS", "STOP"),
    ifelse(timepoint_pass, "PASS", "STOP"),
    ifelse(adult_cohort_field_pass, "PASS", "STOP"),
    ifelse(injury_field_pass, "PASS", "STOP"),
    ifelse(condition_stratification_pass, "PASS", "STOP"),
    ifelse(celltype_field_available, "PASS", "REVIEW"),
    ifelse(spatial_2d_technical_pass, "PASS", "REVIEW"),
    ifelse(raw_archive_pass, "PASS", "REVIEW")
  ),
  detail = c(
    paste0(
      sum(inventory$exists[required_file_rows] & inventory$readable[required_file_rows]),
      "/",
      sum(required_file_rows),
      " required files pass"
    ),
    paste0(
      "matrix=",
      rnaseq_dimension_audit$matrix_n_features[[1L]],
      " x ",
      rnaseq_dimension_audit$matrix_n_barcodes[[1L]],
      "; feature lines=",
      rnaseq_dimension_audit$feature_line_count[[1L]],
      "; barcode lines=",
      rnaseq_dimension_audit$barcode_line_count[[1L]]
    ),
    paste0(
      "metadata rows=",
      nrow(rnaseq_meta),
      "; columns=",
      ncol(rnaseq_meta)
    ),
    paste0(
      "metadata rows=",
      nrow(rnaseq_meta),
      "; matrix barcodes=",
      rnaseq_dimension_audit$matrix_n_barcodes[[1L]]
    ),
    paste0(
      "best barcode column=",
      barcode_mapping$best_column,
      "; mean bidirectional overlap=",
      round(100 * barcode_mapping$best_overlap_fraction, 3),
      "%"
    ),
    paste0(
      "selected sample column=",
      first_nonempty(sample_col)
    ),
    paste0(
      "selected timepoint column=",
      first_nonempty(timepoint_col)
    ),
    paste0(
      "selected age column=",
      first_nonempty(age_col)
    ),
    paste0(
      "selected injury/model column=",
      first_nonempty(injury_col)
    ),
    paste0(
      "selected severity column=",
      first_nonempty(severity_col),
      "; selected treatment column=",
      first_nonempty(treatment_col)
    ),
    paste0(
      "selected cell-type column=",
      first_nonempty(celltype_col)
    ),
    paste0(
      "2D matrix=",
      spatial_2d_audit$matrix_n_features[[1L]],
      " x ",
      spatial_2d_audit$matrix_n_barcodes[[1L]]
    ),
    paste0(
      "archive members=",
      raw_archive_summary$n_members[[1L]],
      "; image-like members=",
      raw_archive_summary$n_image_like_members[[1L]],
      "; position-like members=",
      raw_archive_summary$n_position_like_members[[1L]]
    )
  ),
  stringsAsFactors = FALSE
)

write_csv_safe(
  audit_checks,
  "11_STEP9A_audit_checks.csv"
)

required_fail <- audit_checks$required_for_step9B &
  audit_checks$status != "PASS"
required_fail[is.na(required_fail)] <- TRUE

final_decision <- if (any(required_fail)) "STOP" else "PASS"

# -----------------------------------------------------------------------------
# 11. Final report
# -----------------------------------------------------------------------------
selected_field_lines <- vapply(
  seq_len(nrow(field_detection$selected)),
  function(i) {
    paste0(
      field_detection$selected$field[[i]],
      ": ",
      first_nonempty(field_detection$selected$selected_column[[i]]),
      " [",
      field_detection$selected$confidence[[i]],
      "]"
    )
  },
  character(1)
)

failed_rows <- audit_checks[
  audit_checks$required_for_step9B &
    audit_checks$status != "PASS",
  ,
  drop = FALSE
]

report_lines <- c(
  "GSE234774 STEP 9A snRNA-seq INPUT AND METADATA AUDIT",
  "============================================================",
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Project root: ", project_root),
  paste0("Input directory: ", input_dir),
  paste0("Output directory: ", output_dir),
  "",
  paste0("FINAL DECISION: ", final_decision),
  "",
  "SCOPE",
  "This step audited files, dimensions, barcode correspondence, and metadata only.",
  "The 3.1-GB expression matrix was not loaded into memory; only its Matrix Market header was read.",
  "No Seurat object, QC filtering, integration, annotation, signature transfer, pseudotime, or spatial analysis was run.",
  "No GSE162610 Step 1-Step 8D result was recomputed.",
  "",
  "MAIN snRNA-seq INPUT",
  paste0(
    "Matrix dimensions: ",
    rnaseq_dimension_audit$matrix_n_features[[1L]],
    " features x ",
    rnaseq_dimension_audit$matrix_n_barcodes[[1L]],
    " barcodes"
  ),
  paste0(
    "Nonzero entries reported in Matrix Market header: ",
    rnaseq_dimension_audit$matrix_n_nonzero[[1L]]
  ),
  paste0("Barcode lines: ", length(rnaseq_barcodes)),
  paste0("Feature rows: ", nrow(rnaseq_features)),
  paste0("Metadata: ", nrow(rnaseq_meta), " rows x ", ncol(rnaseq_meta), " columns"),
  paste0("Duplicated barcodes: ", sum(duplicated(rnaseq_barcodes))),
  paste0(
    "Best metadata-barcode column: ",
    barcode_mapping$best_column,
    "; overlap=",
    round(100 * barcode_mapping$best_overlap_fraction, 3),
    "%"
  ),
  "",
  "AUTOMATIC METADATA FIELD SELECTION",
  selected_field_lines,
  "",
  "INTERPRETATION RULES",
  "1. The adult standard-SCI cohort must be selected only from actual metadata values.",
  "2. Age, injury/model, severity, treatment, and timepoint must not be guessed from the paper or filenames.",
  "3. Biological sample, not individual nuclei/cells, is the statistical replicate.",
  "4. Spatial 2D/3D inputs are only audited here; spatial reconstruction is deferred.",
  "5. GSE234774 will validate frozen GSE162610 signatures and will not redefine Arm_1/Arm_2."
)

if (final_decision == "STOP") {
  report_lines <- c(
    report_lines,
    "",
    "STOP REASONS"
  )

  if (nrow(failed_rows) > 0L) {
    for (i in seq_len(nrow(failed_rows))) {
      report_lines <- c(
        report_lines,
        paste0(
          i,
          ". ",
          failed_rows$check[[i]],
          " -- ",
          failed_rows$detail[[i]]
        )
      )
    }
  }

  report_lines <- c(
    report_lines,
    "",
    "Do not start Step 9B yet.",
    "Upload the five key output files listed at the end of the R console.",
    "The detected metadata columns and their values must be reviewed before cohort filtering."
  )
} else {
  report_lines <- c(
    report_lines,
    "",
    "STEP 9A PASS",
    "The snRNA-seq processed inputs and required metadata structure are sufficient to design Step 9B.",
    "Before object reconstruction, manually confirm the exact adult standard-SCI values in",
    "07_RNAseq_selected_field_value_counts.csv and 09_RNAseq_sample_manifest_auto.csv."
  )
}

writeLines(
  report_lines,
  con = file.path(output_dir, "12_STEP9A_FINAL_REPORT.txt"),
  useBytes = TRUE
)

log_message("Final decision: ", final_decision)
log_message(
  "Report: ",
  file.path(output_dir, "12_STEP9A_FINAL_REPORT.txt")
)

cat("\n")
cat("============================================================\n")
cat("GSE234774 STEP 9A COMPLETED\n")
cat("FINAL DECISION:", final_decision, "\n")
cat("Output:", output_dir, "\n")
cat("Please upload these five key files for review:\n")
cat("  1) 02_RNAseq_dimension_consistency.csv\n")
cat("  2) 06_RNAseq_selected_metadata_fields.csv\n")
cat("  3) 07_RNAseq_selected_field_value_counts.csv\n")
cat("  4) 09_RNAseq_sample_manifest_auto.csv\n")
cat("  5) 11_STEP9A_audit_checks.csv\n")
cat("Also keep 12_STEP9A_FINAL_REPORT.txt locally; upload it if the decision is STOP.\n")
cat("============================================================\n")

if (final_decision == "STOP") {
  stop(
    paste(
      "GSE234774 Step 9A decision is STOP.",
      "Review 11_STEP9A_audit_checks.csv and 12_STEP9A_FINAL_REPORT.txt."
    ),
    call. = FALSE
  )
}
