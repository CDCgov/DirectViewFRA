#!/usr/bin/env Rscript
# Run from the repository root: Rscript tests/check_serum_aliases.R
# Checks actual readers; no model fitting or output-folder changes.

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
stopifnot(length(script_arg) == 1L)
test_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
repo_dir <- dirname(dirname(test_path))
source(file.path(repo_dir, "directview_fra_analysis.R"))
stopifnot(requireNamespace("readr", quietly = TRUE))

local({
  work_dir <- tempfile("serum_alias_check_")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  # Reader metadata retains the original CSV headers. Compare canonical column
  # names, row counts, types and values without comparing that source metadata.
  same_columns <- function(a, b) {
    identical(names(a), names(b)) && nrow(a) == nrow(b) &&
      all(vapply(seq_along(a), function(i) identical(a[[i]], b[[i]]), logical(1)))
  }

  check_reader <- function(path, reader, label) {
    original <- read_csv_minimal(path, label)
    stopifnot("Serum" %in% names(original))
    expected <- reader(path)
    stopifnot("sample_id" %in% names(expected), !"serum" %in% names(expected))

    # Every existing alias remains equivalent, including legacy Global.Name.
    for (header in c("Global.Name", "sample_id", "sample", "sera", "antibody_id", " SERUM ")) {
      variant <- original
      names(variant)[names(variant) == "Serum"] <- header
      candidate <- file.path(work_dir, paste0(label, "_variant.csv"))
      readr::write_csv(variant, candidate, na = "")
      stopifnot(same_columns(expected, reader(candidate)))
    }

    # Equal contents still constitute ambiguous metadata: never pick silently.
    ambiguous <- original
    ambiguous[["Global.Name"]] <- ambiguous[["Serum"]]
    candidate <- file.path(work_dir, paste0(label, "_ambiguous.csv"))
    readr::write_csv(ambiguous, candidate, na = "")
    failure <- tryCatch({ reader(candidate); NULL }, error = identity)
    stopifnot(inherits(failure, "error"))
    stopifnot(grepl("multiple columns that could map to 'sample_id'", conditionMessage(failure), fixed = TRUE))
    message(label, ": alias equivalence and ambiguity checks passed.")
    expected
  }

  raw <- check_reader(
    file.path(repo_dir, "example-data", "mixed_quality_example_raw_counts.csv"),
    standardize_raw_counts, "raw"
  )
  fits <- check_reader(
    file.path(repo_dir, "example-data", "mixed_quality_example_ll4_fit_curve.csv"),
    standardize_gen5_fits, "fits"
  )

  # A minimal synthetic lookup exercises the third reader independently.
  lookup_path <- file.path(work_dir, "lookup.csv")
  readr::write_csv(data.frame(
    well_id = "S1", plate_id = "P1", Serum = "example", virus = "example-virus"
  ), lookup_path)
  check_reader(lookup_path, standardize_lookup, "lookup")

  keys <- c("plate_number", "plate_id", "virus", "sample_id")
  samples <- classify_controls(raw)
  raw_keys <- unique(as.data.frame(samples[samples$control_type == "sample", keys]))
  fit_keys <- unique(as.data.frame(fits[, keys]))
  stopifnot(nrow(raw) == 1152L, nrow(fits) == 24L, nrow(raw_keys) == 24L)
  sort_keys <- function(x) {
    x <- x[do.call(order, x), , drop = FALSE]
    rownames(x) <- NULL
    x
  }
  stopifnot(same_columns(sort_keys(raw_keys), sort_keys(fit_keys)))
  message("All Serum input checks passed; all 24 raw curves match an external-fit key.")
})
