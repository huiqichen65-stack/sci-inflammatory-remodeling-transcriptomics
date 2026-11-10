# ==============================================================================
# 63_step12A_GSE234774_spatial_input_section_metadata_audit.R
#
# Project: Adult SCI activated/DAM-like microglia
# Step 12A: GSE234774 spatial 2D input and section metadata audit
#
# Purpose:
#   1) Audit local spatial 2D matrix / barcode / feature / metadata files.
#   2) Inspect GSE234774_RAW.tar member names without manually extracting it.
#   3) Recover evidence-supported section / animal / timepoint / condition fields.
#   4) Test matrix-barcode-feature-metadata consistency.
#   5) Test whether expression / barcode / position / image / scalefactor assets
#      can be linked to every section.
#
# IMPORTANT:
#   - This script DOES NOT create a Seurat object.
#   - This script DOES NOT run RCTD, signature scoring, spatial statistics,
#     normalization, clustering, or plotting.
#   - This script DOES NOT extract GSE234774_RAW.tar.
#   - Missing or ambiguous critical metadata causes STOP, not silent guessing.
#
# Run in RStudio:
#   source("<PROJECT_ROOT>/63_step12A_GSE234774_spatial_input_section_metadata_audit.R")
# ==============================================================================

options(stringsAsFactors = FALSE)
options(warn = 1)

# ------------------------------- 0. Paths --------------------------------------
PROJECT_ROOT <- Sys.getenv("GLIA_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
RAW_DIR <- file.path(PROJECT_ROOT, "01_raw_data", "GSE234774")
OUT_DIR <- file.path(
  PROJECT_ROOT,
  "04_results",
  "step12A_GSE234774_spatial_input_section_metadata_audit"
)

EXPECTED_FILES <- c(
  matrix_2d   = "GSE234774_spatial_2d_filtered_spatial_2d.mtx.gz",
  barcodes_2d = "GSE234774_spatial_2d_barcodes.txt.gz",
  features_2d = "GSE234774_spatial_2d_features.txt.gz",
  metadata_2d = "GSE234774_spatial_2d_meta.txt.gz",
  raw_tar     = "GSE234774_RAW.tar"
)

EXPECTED_PATHS <- file.path(RAW_DIR, unname(EXPECTED_FILES))
names(EXPECTED_PATHS) <- names(EXPECTED_FILES)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------- 1. Package check ---------------------------------
required_packages <- c("data.table")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    paste0(
      "Missing required R package(s): ",
      paste(missing_packages, collapse = ", "),
      "\nInstall once in RStudio with:\ninstall.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "),
      "))"
    ),
    call. = FALSE
  )
}

# ----------------------------- 2. Utilities ------------------------------------
trim_na <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL", "null", ".")] <- NA_character_
  x
}

collapse_unique <- function(x, max_values = 20L) {
  x <- trim_na(x)
  x <- unique(x[!is.na(x)])
  if (length(x) == 0L) {
    return(NA_character_)
  }
  shown <- head(x, max_values)
  out <- paste(shown, collapse = " | ")
  if (length(x) > max_values) {
    out <- paste0(out, " | ... [", length(x), " unique values]")
  }
  out
}

normalize_field_name <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

normalize_token <- function(x) {
  tolower(gsub("[^a-z0-9]", "", as.character(x)))
}

relative_to_raw <- function(path) {
  path_norm <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root_norm <- normalizePath(RAW_DIR, winslash = "/", mustWork = FALSE)
  prefix <- paste0(root_norm, "/")
  out <- ifelse(
    startsWith(path_norm, prefix),
    substring(path_norm, nchar(prefix) + 1L),
    path_norm
  )
  out
}

file_extension_full <- function(path) {
  b <- tolower(basename(path))
  known <- c(
    ".tar.gz", ".nii.gz", ".mtx.gz", ".tsv.gz", ".txt.gz", ".csv.gz",
    ".json.gz", ".h5ad.gz", ".tar", ".tgz", ".zip", ".gz"
  )
  hit <- known[endsWith(b, known)]
  if (length(hit) > 0L) {
    return(sub("^\\.", "", hit[[1L]]))
  }
  ext <- tools::file_ext(b)
  ifelse(nzchar(ext), ext, NA_character_)
}

classify_asset <- function(path) {
  p <- tolower(gsub("\\\\", "/", path))
  b <- basename(p)

  if (grepl("\\.(tar|tar\\.gz|tgz|zip)$", p)) {
    return("archive")
  }
  if (grepl("scalefactors.*\\.json(\\.gz)?$", b)) {
    return("scalefactor")
  }
  if (grepl(
    "(tissue_positions|positions_list|position|coordinates|coords).*(\\.csv|\\.tsv|\\.txt)(\\.gz)?$",
    b
  )) {
    return("position")
  }
  if (grepl(
    "(tissue|histology|hires|lowres).*(\\.png|\\.jpg|\\.jpeg|\\.tif|\\.tiff)$",
    b
  )) {
    return("tissue_image")
  }
  if (grepl("\\.(png|jpg|jpeg|tif|tiff)$", b)) {
    return("image_other")
  }
  if (grepl(
    "(matrix\\.mtx|filtered_feature_bc_matrix.*\\.h5|filtered_spatial.*\\.mtx)(\\.gz)?$",
    b
  )) {
    return("expression_matrix")
  }
  if (grepl("barcode.*(\\.tsv|\\.txt|\\.csv)(\\.gz)?$", b)) {
    return("barcode")
  }
  if (grepl("(feature|gene).*(\\.tsv|\\.txt|\\.csv)(\\.gz)?$", b)) {
    return("feature")
  }
  if (grepl("(meta|metadata).*(\\.tsv|\\.txt|\\.csv)(\\.gz)?$", b)) {
    return("metadata")
  }
  "other"
}

classify_scope <- function(path) {
  p <- tolower(path)
  if (grepl("spatial[_-]?2d", p)) {
    return("spatial_2d")
  }
  if (grepl("spatial[_-]?3d", p)) {
    return("spatial_3d")
  }
  if (grepl("rnaseq|scrna|snrna", p)) {
    return("rnaseq")
  }
  "unspecified"
}

extract_section_token <- function(path) {
  p <- gsub("\\\\", "/", as.character(path))

  gsm_hit <- regmatches(
    p,
    regexpr("GSM[0-9]+", p, ignore.case = TRUE, perl = TRUE)
  )
  if (length(gsm_hit) == 1L && nzchar(gsm_hit)) {
    return(toupper(gsm_hit))
  }

  parts <- unlist(strsplit(p, "/", fixed = TRUE), use.names = FALSE)
  parts <- parts[nzchar(parts)]

  section_like <- parts[
    grepl(
      "(^|[_-])(section|slice|sample|library|lib)[_-]?[a-z0-9.-]+",
      parts,
      ignore.case = TRUE,
      perl = TRUE
    )
  ]
  if (length(section_like) > 0L) {
    return(section_like[[length(section_like)]])
  }

  NA_character_
}

safe_write_csv <- function(x, filename) {
  if (!is.data.frame(x)) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  }
  data.table::fwrite(
    x,
    file = file.path(OUT_DIR, filename),
    sep = ",",
    na = "",
    quote = TRUE,
    row.names = FALSE
  )
}

safe_read_single_column <- function(path) {
  if (!file.exists(path)) {
    return(character())
  }
  out <- tryCatch(
    data.table::fread(
      path,
      header = FALSE,
      select = 1L,
      data.table = FALSE,
      showProgress = FALSE
    )[[1L]],
    error = function(e) {
      warning("Failed to read: ", path, "\n", conditionMessage(e))
      character()
    }
  )
  trim_na(out)
}

read_mtx_header <- function(path) {
  result <- list(
    success = FALSE,
    n_features = NA_real_,
    n_barcodes = NA_real_,
    n_nonzero = NA_real_,
    message = "File missing or unreadable"
  )

  if (!file.exists(path)) {
    return(result)
  }

  con <- NULL
  tryCatch(
    {
      con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
        gzfile(path, open = "rt")
      } else {
        file(path, open = "rt")
      }

      first_line <- readLines(con, n = 1L, warn = FALSE)
      if (length(first_line) != 1L ||
          !startsWith(trimws(first_line), "%%MatrixMarket")) {
        stop("First line is not a MatrixMarket header.")
      }

      dimension_line <- character()
      repeat {
        one_line <- readLines(con, n = 1L, warn = FALSE)
        if (length(one_line) == 0L) {
          break
        }
        one_line <- trimws(one_line)
        if (!nzchar(one_line) || startsWith(one_line, "%")) {
          next
        }
        dimension_line <- one_line
        break
      }

      if (length(dimension_line) != 1L) {
        stop("Matrix dimension line was not found.")
      }

      fields <- strsplit(dimension_line, "\\s+", perl = TRUE)[[1L]]
      if (length(fields) < 3L) {
        stop("Matrix dimension line has fewer than three values.")
      }

      dims <- suppressWarnings(as.numeric(fields[1:3]))
      if (anyNA(dims)) {
        stop("Matrix dimensions are not numeric.")
      }

      result$success <- TRUE
      result$n_features <- dims[[1L]]
      result$n_barcodes <- dims[[2L]]
      result$n_nonzero <- dims[[3L]]
      result$message <- "MatrixMarket header parsed successfully."
    },
    error = function(e) {
      result$message <- conditionMessage(e)
    },
    finally = {
      if (!is.null(con)) {
        try(close(con), silent = TRUE)
      }
    }
  )

  result
}

barcode_overlap_summary <- function(dt, barcodes) {
  if (is.null(dt) || ncol(dt) == 0L || length(barcodes) == 0L) {
    return(data.frame(
      column = character(),
      overlap_n = integer(),
      overlap_fraction = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  barcodes <- unique(trim_na(barcodes))
  barcodes <- barcodes[!is.na(barcodes)]

  rows <- lapply(names(dt), function(nm) {
    values <- unique(trim_na(dt[[nm]]))
    values <- values[!is.na(values)]
    overlap_n <- sum(barcodes %in% values)
    data.frame(
      column = nm,
      overlap_n = overlap_n,
      overlap_fraction = if (length(barcodes) > 0L) {
        overlap_n / length(barcodes)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

read_metadata_robust <- function(path, expected_barcodes) {
  empty <- list(
    data = data.frame(),
    mode = NA_character_,
    message = "Metadata file missing or unreadable",
    overlap = data.frame()
  )

  if (!file.exists(path)) {
    return(empty)
  }

  read_one <- function(header_value) {
    tryCatch(
      data.table::fread(
        path,
        header = header_value,
        data.table = FALSE,
        showProgress = FALSE,
        check.names = FALSE
      ),
      error = function(e) NULL
    )
  }

  dt_auto <- read_one("auto")
  dt_noheader <- read_one(FALSE)

  choices <- list()
  if (!is.null(dt_auto)) {
    choices[["auto_header"]] <- dt_auto
  }
  if (!is.null(dt_noheader)) {
    choices[["no_header"]] <- dt_noheader
  }

  if (length(choices) == 0L) {
    return(empty)
  }

  n_expected <- length(expected_barcodes)

  choice_metrics <- lapply(names(choices), function(mode_name) {
    one <- choices[[mode_name]]
    overlap <- barcode_overlap_summary(one, expected_barcodes)
    best_overlap <- if (nrow(overlap) > 0L) {
      max(overlap$overlap_fraction, na.rm = TRUE)
    } else {
      0
    }
    if (!is.finite(best_overlap)) {
      best_overlap <- 0
    }

    row_distance <- if (n_expected > 0L) {
      abs(nrow(one) - n_expected)
    } else {
      Inf
    }

    exact_row_bonus <- if (n_expected > 0L && nrow(one) == n_expected) 10 else 0
    auto_tie_bonus <- if (mode_name == "auto_header") 0.001 else 0

    data.frame(
      mode = mode_name,
      n_rows = nrow(one),
      n_columns = ncol(one),
      row_distance = row_distance,
      best_barcode_overlap = best_overlap,
      selection_score = exact_row_bonus + best_overlap + auto_tie_bonus,
      stringsAsFactors = FALSE
    )
  })

  metrics <- do.call(rbind, choice_metrics)
  best_idx <- order(
    -metrics$selection_score,
    metrics$row_distance,
    -metrics$best_barcode_overlap
  )[[1L]]
  selected_mode <- metrics$mode[[best_idx]]
  selected <- choices[[selected_mode]]
  selected_overlap <- barcode_overlap_summary(selected, expected_barcodes)

  list(
    data = selected,
    mode = selected_mode,
    message = paste0(
      "Selected ", selected_mode,
      "; rows=", nrow(selected),
      "; columns=", ncol(selected),
      "; expected barcode rows=", n_expected, "."
    ),
    overlap = selected_overlap,
    comparison = metrics
  )
}

role_priority <- list(
  section = c(
    "section_id", "section", "slice_id", "slice",
    "library_id", "library", "sample_id", "sample",
    "orig_ident", "origident"
  ),
  animal = c(
    "animal_id", "animal", "mouse_id", "mouse",
    "subject_id", "subject", "donor_id", "donor",
    "biological_replicate", "replicate"
  ),
  timepoint = c(
    "timepoint", "time_point", "dpi", "days_post_injury",
    "day_post_injury", "post_injury_time", "time"
  ),
  condition = c(
    "injury_condition", "condition", "injury", "group", "status"
  ),
  age = c("age", "age_group", "agegroup"),
  severity = c(
    "severity", "injury_severity", "contusion_severity",
    "lesion_severity", "force"
  ),
  treatment = c(
    "treatment", "drug", "therapy", "intervention", "compound"
  ),
  orientation = c(
    "orientation", "section_orientation", "slice_orientation",
    "cord_level", "spinal_level", "level"
  )
)

role_regex <- list(
  section = "(section|slice|library|sample|orig[._]?ident)",
  animal = "(animal|mouse|subject|donor|biological.*replicate)",
  timepoint = "(timepoint|time_point|dpi|post.*injury|day|week|month)",
  condition = "(condition|injury|group|status|sham|control)",
  age = "(^|_)age($|_)",
  severity = "(severity|contusion|lesion|force)",
  treatment = "(treatment|drug|therapy|intervention|compound)",
  orientation = "(orientation|cord.*level|spinal.*level|section.*level)"
)

role_name_score <- function(column_name, role) {
  nm <- normalize_field_name(column_name)
  priorities <- role_priority[[role]]

  if (nm %in% priorities) {
    position <- match(nm, priorities)
    return(200 - position)
  }

  if (grepl(role_regex[[role]], nm, ignore.case = TRUE, perl = TRUE)) {
    return(50)
  }

  0
}

role_plausibility <- function(n_unique, n_nonmissing, n_rows, role) {
  if (is.na(n_unique) || is.na(n_nonmissing) || n_nonmissing == 0L) {
    return(FALSE)
  }

  if (role %in% c("section", "animal")) {
    return(
      n_unique >= 1L &&
      n_unique <= 500L &&
      n_unique <= max(1L, floor(n_rows / 5L))
    )
  }

  if (role %in% c("timepoint", "condition", "age", "severity", "treatment")) {
    return(n_unique >= 1L && n_unique <= 100L)
  }

  if (role == "orientation") {
    return(n_unique >= 1L && n_unique <= 500L)
  }

  TRUE
}

build_metadata_field_audit <- function(meta, overlap_table) {
  if (ncol(meta) == 0L) {
    return(data.frame(
      column = character(),
      normalized_column = character(),
      storage_class = character(),
      n_rows = integer(),
      n_nonmissing = integer(),
      n_unique = integer(),
      examples = character(),
      barcode_overlap_n = integer(),
      barcode_overlap_fraction = numeric(),
      candidate_roles = character(),
      stringsAsFactors = FALSE
    ))
  }

  rows <- lapply(names(meta), function(nm) {
    values <- trim_na(meta[[nm]])
    nonmissing <- values[!is.na(values)]
    role_hits <- names(role_priority)[
      vapply(
        names(role_priority),
        function(role) role_name_score(nm, role) > 0,
        logical(1)
      )
    ]

    overlap_row <- overlap_table[overlap_table$column == nm, , drop = FALSE]
    overlap_n <- if (nrow(overlap_row) == 1L) {
      overlap_row$overlap_n[[1L]]
    } else {
      0L
    }
    overlap_fraction <- if (nrow(overlap_row) == 1L) {
      overlap_row$overlap_fraction[[1L]]
    } else {
      0
    }

    data.frame(
      column = nm,
      normalized_column = normalize_field_name(nm),
      storage_class = paste(class(meta[[nm]]), collapse = "/"),
      n_rows = length(values),
      n_nonmissing = length(nonmissing),
      n_unique = length(unique(nonmissing)),
      examples = collapse_unique(nonmissing, max_values = 8L),
      barcode_overlap_n = overlap_n,
      barcode_overlap_fraction = overlap_fraction,
      candidate_roles = if (length(role_hits) > 0L) {
        paste(role_hits, collapse = " | ")
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

select_barcode_column <- function(field_audit) {
  if (nrow(field_audit) == 0L) {
    return(list(column = NA_character_, confidence = "none", reason = "No metadata columns."))
  }

  best_fraction <- max(field_audit$barcode_overlap_fraction, na.rm = TRUE)
  if (!is.finite(best_fraction) || best_fraction <= 0) {
    return(list(
      column = NA_character_,
      confidence = "none",
      reason = "No metadata column overlaps the processed barcode list."
    ))
  }

  candidates <- field_audit$column[
    field_audit$barcode_overlap_fraction == best_fraction
  ]

  if (length(candidates) != 1L) {
    return(list(
      column = NA_character_,
      confidence = "ambiguous",
      reason = paste0(
        "Multiple columns share the best barcode overlap: ",
        paste(candidates, collapse = ", ")
      )
    ))
  }

  list(
    column = candidates[[1L]],
    confidence = if (best_fraction == 1) "exact_full_overlap" else "partial_overlap",
    reason = paste0(
      "Unique best exact barcode overlap fraction = ",
      format(best_fraction, digits = 6)
    )
  )
}

select_role_column <- function(field_audit, role) {
  if (nrow(field_audit) == 0L) {
    return(list(
      column = NA_character_,
      confidence = "none",
      reason = "No metadata columns."
    ))
  }

  scores <- vapply(
    field_audit$column,
    function(nm) role_name_score(nm, role),
    numeric(1)
  )

  plausible <- mapply(
    role_plausibility,
    n_unique = field_audit$n_unique,
    n_nonmissing = field_audit$n_nonmissing,
    MoreArgs = list(
      n_rows = if (nrow(field_audit) > 0L) field_audit$n_rows[[1L]] else 0L,
      role = role
    ),
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  )

  scores[!plausible] <- 0
  best_score <- max(scores, na.rm = TRUE)

  if (!is.finite(best_score) || best_score <= 0) {
    return(list(
      column = NA_character_,
      confidence = "none",
      reason = paste0("No plausible ", role, " column identified from field names and cardinality.")
    ))
  }

  candidates <- field_audit$column[scores == best_score]
  if (length(candidates) != 1L) {
    return(list(
      column = NA_character_,
      confidence = "ambiguous",
      reason = paste0(
        "Multiple equally ranked ", role, " columns: ",
        paste(candidates, collapse = ", ")
      )
    ))
  }

  confidence <- if (best_score >= 100) {
    "exact_name_and_cardinality"
  } else {
    "regex_candidate_requires_review"
  }

  list(
    column = candidates[[1L]],
    confidence = confidence,
    reason = paste0(
      "Selected by field-name score ", best_score,
      " plus cardinality plausibility."
    )
  )
}

make_local_inventory <- function() {
  local_paths <- if (dir.exists(RAW_DIR)) {
    list.files(
      RAW_DIR,
      recursive = TRUE,
      full.names = TRUE,
      all.files = TRUE,
      include.dirs = FALSE,
      no.. = TRUE
    )
  } else {
    character()
  }

  if (length(local_paths) > 0L) {
    info <- file.info(local_paths)
    inventory <- data.frame(
      source = "local",
      relative_path = relative_to_raw(local_paths),
      basename = basename(local_paths),
      extension = vapply(local_paths, file_extension_full, character(1)),
      size_bytes = as.numeric(info$size),
      modified_time = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
      scope = vapply(local_paths, classify_scope, character(1)),
      asset_category = vapply(local_paths, classify_asset, character(1)),
      candidate_section_token = vapply(
        local_paths,
        extract_section_token,
        character(1)
      ),
      expected_role = NA_character_,
      exists = TRUE,
      stringsAsFactors = FALSE
    )
  } else {
    inventory <- data.frame(
      source = character(),
      relative_path = character(),
      basename = character(),
      extension = character(),
      size_bytes = numeric(),
      modified_time = character(),
      scope = character(),
      asset_category = character(),
      candidate_section_token = character(),
      expected_role = character(),
      exists = logical(),
      stringsAsFactors = FALSE
    )
  }

  for (role in names(EXPECTED_FILES)) {
    file_name <- EXPECTED_FILES[[role]]
    idx <- which(tolower(inventory$basename) == tolower(file_name))

    if (length(idx) > 0L) {
      inventory$expected_role[idx] <- role
    } else {
      inventory <- rbind(
        inventory,
        data.frame(
          source = "expected_missing",
          relative_path = file_name,
          basename = file_name,
          extension = file_extension_full(file_name),
          size_bytes = NA_real_,
          modified_time = NA_character_,
          scope = classify_scope(file_name),
          asset_category = classify_asset(file_name),
          candidate_section_token = NA_character_,
          expected_role = role,
          exists = FALSE,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  inventory[order(inventory$relative_path), , drop = FALSE]
}

make_tar_inventory <- function(tar_path) {
  empty <- data.frame(
    outer_archive = character(),
    member_path = character(),
    member_directory = character(),
    member_basename = character(),
    extension = character(),
    scope = character(),
    asset_category = character(),
    candidate_section_token = character(),
    nested_archive = logical(),
    unsafe_member_path = logical(),
    stringsAsFactors = FALSE
  )

  if (!file.exists(tar_path)) {
    attr(empty, "listing_success") <- FALSE
    attr(empty, "listing_message") <- "GSE234774_RAW.tar is missing."
    return(empty)
  }

  members <- tryCatch(
    suppressWarnings(utils::untar(tar_path, list = TRUE)),
    error = function(e) structure(character(), error_message = conditionMessage(e))
  )

  error_message <- attr(members, "error_message")
  if (length(members) == 0L) {
    attr(empty, "listing_success") <- FALSE
    attr(empty, "listing_message") <- if (!is.null(error_message)) {
      error_message
    } else {
      "RAW.tar listing returned zero members."
    }
    return(empty)
  }

  member_path_norm <- gsub("\\\\", "/", members)
  unsafe <- grepl("(^/|^[A-Za-z]:|(^|/)\\.\\.(/|$))", member_path_norm, perl = TRUE)

  out <- data.frame(
    outer_archive = basename(tar_path),
    member_path = member_path_norm,
    member_directory = dirname(member_path_norm),
    member_basename = basename(member_path_norm),
    extension = vapply(member_path_norm, file_extension_full, character(1)),
    scope = vapply(member_path_norm, classify_scope, character(1)),
    asset_category = vapply(member_path_norm, classify_asset, character(1)),
    candidate_section_token = vapply(
      member_path_norm,
      extract_section_token,
      character(1)
    ),
    nested_archive = grepl(
      "\\.(tar|tar\\.gz|tgz|zip)$",
      member_path_norm,
      ignore.case = TRUE
    ),
    unsafe_member_path = unsafe,
    stringsAsFactors = FALSE
  )

  attr(out, "listing_success") <- TRUE
  attr(out, "listing_message") <- paste0(
    "Listed ", nrow(out),
    " outer TAR members. No files were extracted."
  )
  out
}

build_section_manifest <- function(meta, selections) {
  empty <- data.frame(
    section_id = character(),
    n_spots = integer(),
    n_unique_barcodes = integer(),
    animal_id = character(),
    animal_n_unique = integer(),
    timepoint = character(),
    timepoint_n_unique = integer(),
    injury_condition = character(),
    condition_n_unique = integer(),
    age = character(),
    age_n_unique = integer(),
    severity = character(),
    severity_n_unique = integer(),
    treatment = character(),
    treatment_n_unique = integer(),
    orientation = character(),
    orientation_n_unique = integer(),
    section_field = character(),
    barcode_field = character(),
    animal_field = character(),
    timepoint_field = character(),
    condition_field = character(),
    stringsAsFactors = FALSE
  )

  section_col <- selections$section$column
  barcode_col <- selections$barcode$column

  if (nrow(meta) == 0L || is.na(section_col) || !section_col %in% names(meta)) {
    return(empty)
  }

  section_values <- trim_na(meta[[section_col]])
  section_values[is.na(section_values)] <- "<MISSING_SECTION>"
  split_index <- split(seq_len(nrow(meta)), section_values, drop = TRUE)

  get_column_values <- function(idx, role) {
    col <- selections[[role]]$column
    if (is.na(col) || !col %in% names(meta)) {
      return(rep(NA_character_, length(idx)))
    }
    meta[[col]][idx]
  }

  unique_count <- function(x) {
    x <- trim_na(x)
    length(unique(x[!is.na(x)]))
  }

  rows <- lapply(names(split_index), function(section_name) {
    idx <- split_index[[section_name]]

    barcode_values <- if (!is.na(barcode_col) && barcode_col %in% names(meta)) {
      trim_na(meta[[barcode_col]][idx])
    } else {
      rep(NA_character_, length(idx))
    }

    animal_values <- get_column_values(idx, "animal")
    time_values <- get_column_values(idx, "timepoint")
    condition_values <- get_column_values(idx, "condition")
    age_values <- get_column_values(idx, "age")
    severity_values <- get_column_values(idx, "severity")
    treatment_values <- get_column_values(idx, "treatment")
    orientation_values <- get_column_values(idx, "orientation")

    data.frame(
      section_id = section_name,
      n_spots = length(idx),
      n_unique_barcodes = unique_count(barcode_values),
      animal_id = collapse_unique(animal_values),
      animal_n_unique = unique_count(animal_values),
      timepoint = collapse_unique(time_values),
      timepoint_n_unique = unique_count(time_values),
      injury_condition = collapse_unique(condition_values),
      condition_n_unique = unique_count(condition_values),
      age = collapse_unique(age_values),
      age_n_unique = unique_count(age_values),
      severity = collapse_unique(severity_values),
      severity_n_unique = unique_count(severity_values),
      treatment = collapse_unique(treatment_values),
      treatment_n_unique = unique_count(treatment_values),
      orientation = collapse_unique(orientation_values),
      orientation_n_unique = unique_count(orientation_values),
      section_field = section_col,
      barcode_field = barcode_col,
      animal_field = selections$animal$column,
      timepoint_field = selections$timepoint$column,
      condition_field = selections$condition$column,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$section_id), , drop = FALSE]
}

build_asset_table <- function(local_inventory, tar_inventory) {
  local_assets <- data.frame(
    source = "local",
    asset_path = local_inventory$relative_path,
    asset_category = local_inventory$asset_category,
    candidate_section_token = local_inventory$candidate_section_token,
    stringsAsFactors = FALSE
  )

  tar_assets <- data.frame(
    source = "RAW.tar",
    asset_path = tar_inventory$member_path,
    asset_category = tar_inventory$asset_category,
    candidate_section_token = tar_inventory$candidate_section_token,
    stringsAsFactors = FALSE
  )

  out <- rbind(local_assets, tar_assets)
  out <- out[
    out$asset_category %in% c(
      "expression_matrix", "barcode", "feature", "metadata",
      "position", "scalefactor", "tissue_image", "image_other", "archive"
    ),
    ,
    drop = FALSE
  ]
  rownames(out) <- NULL
  out
}

match_assets_to_section <- function(asset_table, section_id) {
  if (nrow(asset_table) == 0L || is.na(section_id) || !nzchar(section_id)) {
    return(rep(FALSE, nrow(asset_table)))
  }

  exact_path <- grepl(
    section_id,
    asset_table$asset_path,
    fixed = TRUE,
    ignore.case = TRUE
  )

  normalized_section <- normalize_token(section_id)
  normalized_candidates <- normalize_token(asset_table$candidate_section_token)
  token_match <- !is.na(asset_table$candidate_section_token) &
    nzchar(normalized_section) &
    normalized_candidates == normalized_section

  exact_path | token_match
}

build_completeness <- function(
  section_manifest,
  asset_table,
  matrix_header,
  n_barcodes,
  n_features
) {
  empty <- data.frame(
    section_id = character(),
    animal_id = character(),
    timepoint = character(),
    injury_condition = character(),
    n_spots_in_metadata = integer(),
    processed_expression_available = logical(),
    processed_barcode_available = logical(),
    processed_feature_available = logical(),
    position_file_count = integer(),
    tissue_image_count = integer(),
    other_image_count = integer(),
    scalefactor_file_count = integer(),
    section_specific_asset_match = logical(),
    expression_barcode_position_image_complete = logical(),
    missing_assets = character(),
    matched_asset_examples = character(),
    stringsAsFactors = FALSE
  )

  if (nrow(section_manifest) == 0L) {
    return(empty)
  }

  processed_expression_ok <- isTRUE(matrix_header$success)
  processed_barcode_ok <- n_barcodes > 0L
  processed_feature_ok <- n_features > 0L

  rows <- lapply(seq_len(nrow(section_manifest)), function(i) {
    one <- section_manifest[i, , drop = FALSE]
    section_id <- one$section_id[[1L]]
    matched <- match_assets_to_section(asset_table, section_id)
    matched_assets <- asset_table[matched, , drop = FALSE]

    count_category <- function(category) {
      sum(matched_assets$asset_category == category, na.rm = TRUE)
    }

    position_n <- count_category("position")
    tissue_image_n <- count_category("tissue_image")
    other_image_n <- count_category("image_other")
    scale_n <- count_category("scalefactor")

    section_match <- nrow(matched_assets) > 0L
    complete <- processed_expression_ok &&
      processed_barcode_ok &&
      processed_feature_ok &&
      one$n_spots[[1L]] > 0L &&
      position_n > 0L &&
      (tissue_image_n + other_image_n) > 0L &&
      scale_n > 0L

    missing <- character()
    if (!processed_expression_ok) missing <- c(missing, "processed expression")
    if (!processed_barcode_ok) missing <- c(missing, "processed barcode")
    if (!processed_feature_ok) missing <- c(missing, "processed feature")
    if (one$n_spots[[1L]] <= 0L) missing <- c(missing, "section spots")
    if (position_n == 0L) missing <- c(missing, "position")
    if ((tissue_image_n + other_image_n) == 0L) missing <- c(missing, "image")
    if (scale_n == 0L) missing <- c(missing, "scalefactor")

    data.frame(
      section_id = section_id,
      animal_id = one$animal_id[[1L]],
      timepoint = one$timepoint[[1L]],
      injury_condition = one$injury_condition[[1L]],
      n_spots_in_metadata = one$n_spots[[1L]],
      processed_expression_available = processed_expression_ok,
      processed_barcode_available = processed_barcode_ok,
      processed_feature_available = processed_feature_ok,
      position_file_count = position_n,
      tissue_image_count = tissue_image_n,
      other_image_count = other_image_n,
      scalefactor_file_count = scale_n,
      section_specific_asset_match = section_match,
      expression_barcode_position_image_complete = complete,
      missing_assets = if (length(missing) == 0L) {
        NA_character_
      } else {
        paste(missing, collapse = " | ")
      },
      matched_asset_examples = collapse_unique(
        matched_assets$asset_path,
        max_values = 8L
      ),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# ---------------------------- 3. Input inventory -------------------------------
local_inventory <- make_local_inventory()
safe_write_csv(local_inventory, "01_spatial_file_inventory.csv")

# ----------------------------- 4. RAW.tar audit --------------------------------
tar_inventory <- make_tar_inventory(EXPECTED_PATHS[["raw_tar"]])
tar_listing_success <- isTRUE(attr(tar_inventory, "listing_success"))
tar_listing_message <- attr(tar_inventory, "listing_message")
safe_write_csv(tar_inventory, "02_RAW_tar_section_member_inventory.csv")

# ------------------------- 5. Processed file audit -----------------------------
matrix_header <- read_mtx_header(EXPECTED_PATHS[["matrix_2d"]])
barcodes <- safe_read_single_column(EXPECTED_PATHS[["barcodes_2d"]])
features <- safe_read_single_column(EXPECTED_PATHS[["features_2d"]])

metadata_result <- read_metadata_robust(
  EXPECTED_PATHS[["metadata_2d"]],
  barcodes
)
meta <- metadata_result$data

field_audit <- build_metadata_field_audit(meta, metadata_result$overlap)

barcode_selection <- select_barcode_column(field_audit)
role_selections <- lapply(
  names(role_priority),
  function(role) select_role_column(field_audit, role)
)
names(role_selections) <- names(role_priority)
selections <- c(list(barcode = barcode_selection), role_selections)

selected_role_by_column <- rep(NA_character_, nrow(field_audit))
if (nrow(field_audit) > 0L) {
  for (role in names(selections)) {
    selected_column <- selections[[role]]$column
    if (!is.na(selected_column)) {
      idx <- which(field_audit$column == selected_column)
      selected_role_by_column[idx] <- ifelse(
        is.na(selected_role_by_column[idx]),
        role,
        paste(selected_role_by_column[idx], role, sep = " | ")
      )
    }
  }
}
field_audit$selected_for_role <- selected_role_by_column
field_audit$metadata_read_mode <- metadata_result$mode
safe_write_csv(field_audit, "03_spatial_metadata_field_audit.csv")

# ---------------------------- 6. Section manifest ------------------------------
section_manifest <- build_section_manifest(meta, selections)
safe_write_csv(
  section_manifest,
  "04_section_animal_timepoint_condition_manifest.csv"
)

# --------------------------- 7. Asset completeness -----------------------------
asset_table <- build_asset_table(local_inventory, tar_inventory)
completeness <- build_completeness(
  section_manifest = section_manifest,
  asset_table = asset_table,
  matrix_header = matrix_header,
  n_barcodes = length(barcodes),
  n_features = length(features)
)
safe_write_csv(
  completeness,
  "05_section_expression_position_image_completeness.csv"
)

# ------------------------------ 8. Audit checks --------------------------------
checks <- data.frame(
  check_id = character(),
  category = character(),
  description = character(),
  status = character(),
  observed = character(),
  required = character(),
  critical = logical(),
  action_if_not_pass = character(),
  stringsAsFactors = FALSE
)

add_check <- function(
  check_id,
  category,
  description,
  status,
  observed,
  required,
  critical,
  action_if_not_pass
) {
  checks <<- rbind(
    checks,
    data.frame(
      check_id = check_id,
      category = category,
      description = description,
      status = status,
      observed = as.character(observed),
      required = as.character(required),
      critical = isTRUE(critical),
      action_if_not_pass = action_if_not_pass,
      stringsAsFactors = FALSE
    )
  )
}

expected_exists <- file.exists(EXPECTED_PATHS)
add_check(
  "C01", "input",
  "All five required Step 12A input files exist.",
  if (all(expected_exists)) "PASS" else "FAIL",
  paste(
    paste(names(expected_exists), expected_exists, sep = "="),
    collapse = "; "
  ),
  "matrix/barcodes/features/meta/RAW.tar all TRUE",
  TRUE,
  "Restore the missing local file(s); do not download silently."
)

add_check(
  "C02", "matrix",
  "MatrixMarket header is readable.",
  if (isTRUE(matrix_header$success)) "PASS" else "FAIL",
  matrix_header$message,
  "Readable MatrixMarket dimensions",
  TRUE,
  "Check whether the spatial 2D MTX file is complete and not corrupted."
)

add_check(
  "C03", "matrix",
  "Matrix feature dimension equals feature-file row count.",
  if (
    isTRUE(matrix_header$success) &&
    matrix_header$n_features == length(features)
  ) "PASS" else "FAIL",
  paste0(
    "matrix rows=", matrix_header$n_features,
    "; feature rows=", length(features)
  ),
  "Equal",
  TRUE,
  "Do not create a spatial object until feature dimensions match."
)

add_check(
  "C04", "matrix",
  "Matrix barcode dimension equals barcode-file row count.",
  if (
    isTRUE(matrix_header$success) &&
    matrix_header$n_barcodes == length(barcodes)
  ) "PASS" else "FAIL",
  paste0(
    "matrix columns=", matrix_header$n_barcodes,
    "; barcode rows=", length(barcodes)
  ),
  "Equal",
  TRUE,
  "Do not create a spatial object until barcode dimensions match."
)

add_check(
  "C05", "metadata",
  "Spatial metadata is readable.",
  if (nrow(meta) > 0L && ncol(meta) > 0L) "PASS" else "FAIL",
  paste0(
    metadata_result$message,
    " Metadata dimensions=", nrow(meta), " x ", ncol(meta), "."
  ),
  "Non-empty metadata table",
  TRUE,
  "Inspect the metadata file delimiter/header and file integrity."
)

best_barcode_overlap <- if (
  !is.na(barcode_selection$column) &&
  barcode_selection$column %in% field_audit$column
) {
  field_audit$barcode_overlap_fraction[
    match(barcode_selection$column, field_audit$column)
  ]
} else {
  0
}

add_check(
  "C06", "metadata",
  "A unique barcode column is identified by exact overlap.",
  if (!is.na(barcode_selection$column)) "PASS" else "FAIL",
  paste0(
    "column=", barcode_selection$column,
    "; confidence=", barcode_selection$confidence,
    "; ", barcode_selection$reason
  ),
  "One unique metadata barcode column",
  TRUE,
  "Do not guess the merge key; inspect metadata fields."
)

add_check(
  "C07", "metadata",
  "Every processed barcode is represented exactly in metadata.",
  if (
    !is.na(barcode_selection$column) &&
    isTRUE(all.equal(best_barcode_overlap, 1))
  ) "PASS" else "FAIL",
  paste0(
    "exact barcode overlap fraction=",
    format(best_barcode_overlap, digits = 8)
  ),
  "1.0",
  TRUE,
  "Resolve barcode formatting or missing metadata rows before merging."
)

metadata_barcode_duplicates <- if (
  !is.na(barcode_selection$column) &&
  barcode_selection$column %in% names(meta)
) {
  sum(duplicated(trim_na(meta[[barcode_selection$column]])), na.rm = TRUE)
} else {
  NA_integer_
}

add_check(
  "C08", "metadata",
  "Metadata barcode key is unique.",
  if (!is.na(metadata_barcode_duplicates) &&
      metadata_barcode_duplicates == 0L) "PASS" else "FAIL",
  paste0("duplicated metadata barcodes=", metadata_barcode_duplicates),
  "0",
  TRUE,
  "Determine why one barcode maps to multiple metadata rows."
)

for (role in c("section", "animal", "timepoint", "condition")) {
  selection <- selections[[role]]
  pass <- !is.na(selection$column) &&
    selection$confidence == "exact_name_and_cardinality"

  add_check(
    paste0(
      "C",
      sprintf("%02d", match(role, c("section", "animal", "timepoint", "condition")) + 8L)
    ),
    "metadata",
    paste0("A reliable ", role, " field is recovered."),
    if (pass) "PASS" else "FAIL",
    paste0(
      "column=", selection$column,
      "; confidence=", selection$confidence,
      "; ", selection$reason
    ),
    "Exact-name field with plausible cardinality",
    TRUE,
    paste0(
      "Do not infer ", role,
      " from memory; review metadata and source naming evidence."
    )
  )
}

add_check(
  "C13", "archive",
  "RAW.tar member listing succeeds without extraction.",
  if (tar_listing_success) "PASS" else "FAIL",
  tar_listing_message,
  "Non-empty TAR member inventory",
  TRUE,
  "Verify RAW.tar integrity and local path."
)

global_position_n <- sum(asset_table$asset_category == "position", na.rm = TRUE)
global_image_n <- sum(
  asset_table$asset_category %in% c("tissue_image", "image_other"),
  na.rm = TRUE
)
global_scalefactor_n <- sum(
  asset_table$asset_category == "scalefactor",
  na.rm = TRUE
)

add_check(
  "C14", "spatial_asset",
  "At least one position/coordinate file is visible locally or in RAW.tar.",
  if (global_position_n > 0L) "PASS" else "FAIL",
  paste0("position files=", global_position_n),
  ">0",
  TRUE,
  "If RAW.tar contains nested archives, a later controlled nested-archive audit may be required."
)

add_check(
  "C15", "spatial_asset",
  "At least one tissue/image file is visible locally or in RAW.tar.",
  if (global_image_n > 0L) "PASS" else "FAIL",
  paste0("image files=", global_image_n),
  ">0",
  TRUE,
  "Do not proceed to Seurat spatial loading without section images."
)

add_check(
  "C16", "spatial_asset",
  "At least one scalefactor JSON is visible locally or in RAW.tar.",
  if (global_scalefactor_n > 0L) "PASS" else "FAIL",
  paste0("scalefactor files=", global_scalefactor_n),
  ">0",
  TRUE,
  "Locate the correct scalefactor file for every section."
)

all_sections_complete <- nrow(completeness) > 0L &&
  all(completeness$expression_barcode_position_image_complete)

add_check(
  "C17", "section_linkage",
  "Every recovered section has expression/barcode/position/image/scalefactor support.",
  if (all_sections_complete) "PASS" else "FAIL",
  if (nrow(completeness) == 0L) {
    "No section-level completeness rows were generated."
  } else {
    paste0(
      sum(completeness$expression_barcode_position_image_complete),
      "/", nrow(completeness), " sections complete."
    )
  },
  "All sections complete",
  TRUE,
  "Review section identifiers and section-specific RAW.tar member matching."
)

section_consistent <- nrow(section_manifest) > 0L &&
  all(
    section_manifest$animal_n_unique <= 1L &
    section_manifest$timepoint_n_unique <= 1L &
    section_manifest$condition_n_unique <= 1L
  )

add_check(
  "C18", "section_linkage",
  "Each section maps to no more than one animal/timepoint/condition.",
  if (section_consistent) "PASS" else "FAIL",
  if (nrow(section_manifest) == 0L) {
    "No section manifest generated."
  } else {
    paste0(
      "sections with conflicts=",
      sum(
        section_manifest$animal_n_unique > 1L |
        section_manifest$timepoint_n_unique > 1L |
        section_manifest$condition_n_unique > 1L
      )
    )
  },
  "Zero conflicting sections",
  TRUE,
  "Resolve metadata conflicts before any section-level analysis."
)

nested_archive_n <- sum(tar_inventory$nested_archive, na.rm = TRUE)
add_check(
  "C19", "archive",
  "Nested archives inside RAW.tar are documented.",
  if (nested_archive_n == 0L) "PASS" else "REVIEW",
  paste0("nested archive members=", nested_archive_n),
  "0, or separately audited later",
  FALSE,
  "Do not manually unpack them; use a controlled script only if direct assets are not visible."
)

unsafe_member_n <- sum(tar_inventory$unsafe_member_path, na.rm = TRUE)
add_check(
  "C20", "archive",
  "RAW.tar member paths contain no absolute or parent-traversal paths.",
  if (tar_listing_success && unsafe_member_n == 0L) "PASS" else "FAIL",
  paste0("unsafe member paths=", unsafe_member_n),
  "0",
  TRUE,
  "Do not extract an archive containing unsafe member paths."
)

for (role in c("age", "severity", "treatment", "orientation")) {
  selection <- selections[[role]]
  status <- if (
    !is.na(selection$column) &&
    selection$confidence == "exact_name_and_cardinality"
  ) {
    "PASS"
  } else {
    "REVIEW"
  }

  add_check(
    paste0("R_", toupper(role)),
    "optional_metadata_review",
    paste0("Audit ", role, " field for adult standard SCI filtering."),
    status,
    paste0(
      "column=", selection$column,
      "; confidence=", selection$confidence,
      "; ", selection$reason
    ),
    "Recover when present; otherwise document the limitation",
    FALSE,
    "Review source metadata before defining the adult standard SCI cohort."
  )
}

safe_write_csv(checks, "09_STEP12A_audit_checks.csv")

# ----------------------------- 9. Final report ---------------------------------
critical_failures <- checks[
  checks$critical & checks$status == "FAIL",
  ,
  drop = FALSE
]
review_items <- checks[checks$status == "REVIEW", , drop = FALSE]
FINAL_STATUS <- if (nrow(critical_failures) == 0L) "PASS" else "STOP"

selected_field_line <- function(role) {
  one <- selections[[role]]
  paste0(
    role, "=",
    ifelse(is.na(one$column), "<UNRESOLVED>", one$column),
    " [", one$confidence, "]"
  )
}

category_counts <- if (nrow(tar_inventory) > 0L) {
  sort(table(tar_inventory$asset_category), decreasing = TRUE)
} else {
  integer()
}

category_count_text <- if (length(category_counts) > 0L) {
  paste(
    paste(names(category_counts), as.integer(category_counts), sep = "="),
    collapse = "; "
  )
} else {
  "<none>"
}

report_lines <- c(
  "============================================================",
  "GSE234774 Step 12A spatial input / section metadata audit",
  "============================================================",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Project root: ", PROJECT_ROOT),
  paste0("Input directory: ", RAW_DIR),
  paste0("Output directory: ", OUT_DIR),
  "",
  paste0("FINAL_STATUS: ", FINAL_STATUS),
  "",
  "SCOPE",
  "- Audit only.",
  "- No Seurat object was created.",
  "- No RCTD, signature scoring, normalization, clustering, or spatial statistics were run.",
  "- GSE234774_RAW.tar was listed but not extracted.",
  "",
  "PROCESSED SPATIAL 2D INPUT",
  paste0(
    "- Matrix dimensions: features=", matrix_header$n_features,
    "; barcodes=", matrix_header$n_barcodes,
    "; nonzero=", matrix_header$n_nonzero
  ),
  paste0("- Feature-file rows: ", length(features)),
  paste0("- Barcode-file rows: ", length(barcodes)),
  paste0("- Metadata dimensions: ", nrow(meta), " x ", ncol(meta)),
  paste0("- Metadata read mode: ", metadata_result$mode),
  "",
  "SELECTED METADATA FIELDS",
  paste0("- ", selected_field_line("barcode")),
  paste0("- ", selected_field_line("section")),
  paste0("- ", selected_field_line("animal")),
  paste0("- ", selected_field_line("timepoint")),
  paste0("- ", selected_field_line("condition")),
  paste0("- ", selected_field_line("age")),
  paste0("- ", selected_field_line("severity")),
  paste0("- ", selected_field_line("treatment")),
  paste0("- ", selected_field_line("orientation")),
  "",
  "RAW.TAR",
  paste0("- Listing: ", tar_listing_message),
  paste0("- Category counts: ", category_count_text),
  paste0("- Nested archive members: ", nested_archive_n),
  paste0("- Unsafe member paths: ", unsafe_member_n),
  "",
  "SECTION SUMMARY",
  paste0("- Recovered sections: ", nrow(section_manifest)),
  paste0(
    "- Fully linked sections: ",
    if (nrow(completeness) > 0L) {
      sum(completeness$expression_barcode_position_image_complete)
    } else {
      0L
    },
    "/", nrow(completeness)
  ),
  "",
  "CRITICAL CHECKS",
  paste0(
    "- PASS: ",
    sum(checks$critical & checks$status == "PASS"),
    "; FAIL: ",
    sum(checks$critical & checks$status == "FAIL")
  )
)

if (nrow(critical_failures) > 0L) {
  report_lines <- c(
    report_lines,
    "",
    "STOP REASONS"
  )
  for (i in seq_len(nrow(critical_failures))) {
    report_lines <- c(
      report_lines,
      paste0(
        "- ", critical_failures$check_id[[i]], ": ",
        critical_failures$description[[i]],
        " Observed: ", critical_failures$observed[[i]],
        " Action: ", critical_failures$action_if_not_pass[[i]]
      )
    )
  }
}

if (nrow(review_items) > 0L) {
  report_lines <- c(
    report_lines,
    "",
    "REVIEW ITEMS"
  )
  for (i in seq_len(nrow(review_items))) {
    report_lines <- c(
      report_lines,
      paste0(
        "- ", review_items$check_id[[i]], ": ",
        review_items$description[[i]],
        " Observed: ", review_items$observed[[i]]
      )
    )
  }
}

report_lines <- c(
  report_lines,
  "",
  "NEXT ACTION",
  if (FINAL_STATUS == "PASS") {
    paste(
      "Step 12A PASS. Upload the seven key outputs for review.",
      "Only after review should Step 12B read and QC each section separately."
    )
  } else {
    paste(
      "Step 12A STOP. Upload the seven key outputs for diagnosis.",
      "Do not create a Seurat spatial object and do not run RCTD or target scoring."
    )
  },
  "",
  "KEY OUTPUTS",
  "- 01_spatial_file_inventory.csv",
  "- 02_RAW_tar_section_member_inventory.csv",
  "- 03_spatial_metadata_field_audit.csv",
  "- 04_section_animal_timepoint_condition_manifest.csv",
  "- 05_section_expression_position_image_completeness.csv",
  "- 09_STEP12A_audit_checks.csv",
  "- 10_STEP12A_FINAL_REPORT.txt",
  "============================================================"
)

writeLines(
  report_lines,
  con = file.path(OUT_DIR, "10_STEP12A_FINAL_REPORT.txt"),
  useBytes = TRUE
)

message("============================================================")
message("Step 12A completed.")
message("FINAL_STATUS: ", FINAL_STATUS)
message("Output directory: ", OUT_DIR)
message("Please upload the seven listed key outputs for review.")
message("============================================================")
