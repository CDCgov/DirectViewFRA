#!/usr/bin/env Rscript

# DirectView-FRA neutralization analysis
# --------------------------------------
#
# This script analyzes per-well influenza focus counts, normalizes each sample
# to a plate- and virus-matched control mean by default, fits a four-parameter
# log-logistic (LL4) model, estimates NT50 by both LL4 inverse prediction and
# observed-value interpolation with uniroot(), applies quality-control checks,
# and writes each enabled plot branch to its own output folder.
#
# Design principles:
#   * Do not report extrapolated NT50 values. The observed dilution means must
#     bracket the requested threshold, and each reported estimate must be inside
#     the tested dilution range.
#   * Fit an unconstrained LL4 first. If its Hill estimate is outside the
#     prespecified range, refit only that curve with Hill fixed at 1.0.
#   * Support more than one virus on a physical plate. Virus controls are matched
#     by explicit metadata by default. A row/metadata conflict is not silently
#     reassigned; optional row inference can be enabled for validated layouts.
#   * Calculate LL4 and uniroot estimates independently when their respective
#     inputs allow it. The uniroot estimate uses observed dilution-level mean
#     percent neutralization values.
#   * Treat confidence intervals as method-specific. The LL4 ED estimate can
#     receive a model-based interval. The observed-value uniroot estimate is
#     exported without borrowing the LL4 interval.
#   * Run only explicitly enabled plot branches. Each branch has its own folder,
#     and any plot fallback is local to the failed plot.
#   * Build paths with base R path functions so the same script works on Windows
#     and Linux. Forward slashes are accepted on Windows, and relative paths are
#     resolved against the configured input folder.

SCRIPT_VERSION <- "1.0.0"
MIN_R_VERSION <- "4.1.0"

# The following libraries are required: drc, dplyr, ggplot2, readr, 
# scales, and tibble. The lack of aatached packages is an intentional design 
# choice, and instead fully qualified package calls are used throughout.  This
# was to make it simplier for handling certain errors for novice users.  Look
# to the function `check_required_packages` if you would like to edit this behavior
# in any local deployment. 

# =============================================================================
# CONFIGURATION SOURCE
# =============================================================================
# TRUE (default): use the complete HARDCODED_SETTINGS block below. This is the
# easiest route in Jupyter, RStudio, and other interactive environments.
# FALSE: ignore HARDCODED_SETTINGS and begin from the built-in defaults.
# Command-line arguments can override either route and always win.
USE_HARDCODED_SETTINGS <- TRUE

# =============================================================================
# HARDCODED SETTINGS
# Edit values inside this one bracket when USE_HARDCODED_SETTINGS is TRUE.
# All paths and optional analysis, QC, and plotting choices are exposed here.
# =============================================================================
HARDCODED_SETTINGS <- list(
  paths = list(
    # Reset getwd() to the folder containing this script when that folder can be
    # detected. In Jupyter, a fully specified input_folder remains authoritative.
    # This can be useful for more novice users, and helpped in beta testing.
    set_working_directory_to_script_folder = TRUE,

    # Folder containing the input files. Use forward slashes on Windows to avoid
    # escaping backslashes, for example "C:/DirectViewFRA/Run_01".
    input_folder = "example-data"
  ),
  # Relative paths are created inside input_folder.
  output_dir = "directview_fra_results",
  

  input = list(
    # Required input. May be a filename inside input_folder or an absolute path.
    raw_counts_file = "mixed_quality_example_raw_counts.csv",

    # Optional inputs. Leave as NULL when unused.
     
    # Can be used to help with mapping, but likely will need edits to the code
    # based on the data layout used by each lab.
    lookup_file = NULL,
    
    # Gen5 as a placeholder for the fits that may or may not be provided from 
    # any commercial software. This is currently setup to take files in the 
    # format of the provided example CSV file for comparison.
    gen5_fits_file = NULL
    #gen5_fits_file = "mixed_quality_example_ll4_fit_curve.csv"
  ),

  normalization = list(
    # Multiple viruses may occur on one physical plate. Keep FALSE unless the
    # row-block layout has been validated and row identity should override
    # conflicting or missing control-virus metadata.
    infer_control_virus_from_well_row = FALSE,

    # Cross-plate same-virus VCs are an explicit opt-in. Matches strings exactly.
    allow_cross_plate_virus_controls = FALSE,

    # Optional restriction on contributing source plates. NULL permits all
    # otherwise eligible same-virus VC plates in the input dataset.
    # Here you could limit the specific plates that are allowed to use a pooled
    # virus control (for and grouping preference such as temperature testing if
    # not captured elsewhere).
    
    cross_plate_vc_plate_numbers = NULL,

    # Separate permission for a cross-plate denominator to pass the reporting
    # gate. Keep FALSE unless this was prespecified and validated study-wide.
    
    allow_cross_plate_vc_for_reporting = FALSE,

    # TRUE gives each contributing plate mean equal weight when VCs are pooled.
    pool_equal_weight_by_plate = TRUE
  ),

  analysis = list(
    threshold_pct = 50,
    reporting_method = "ll4",              # "ll4" or "uniroot"
    uniroot_interpolation_scale = "linear", # "linear" or "log2"
    ll4_ci_method = "delta",                # "none", "delta", or "inv"
    confidence_level = 0.95,
    min_unique_dilutions = 4L,
    min_total_observations = 8L,
    prediction_grid_points = 400L
  ),

  qc = list(
    vc_count_range = c(100, 1600),
    hill_range = c(0.5, 2.0),
    fixed_hill = 1.0,

    # Optional review thresholds. NULL calculates and exports the diagnostic but
    # does not use it to gate reportability.
    vc_cv_max_pct = NULL,
    max_replicate_sd_pp = NULL,
    cell_control_count_max = NULL,
    lack_of_fit_alpha = NULL,
    min_ll4_pseudo_r2 = NULL,
    max_ll4_rmse_pp = NULL,
    max_upward_mean_steps = NULL,
    max_ll4_uniroot_disagreement_fold = NULL,
    max_ll4_gen5_disagreement_fold = NULL,

    min_vc_replicates = 2L,
    require_min_vc_replicates_for_reporting = FALSE
  ),

  plot = list(
    clear_enabled_branch_directories = TRUE,

    # Readability controls for primary QC flags below individual plots.
    flag_wrap_width = 72L,
    blank_lines_between_flags = 1L,
    individual_height_per_flag_line = 0.30,
    max_individual_plot_height = 48,

    # Multiple branches may be enabled together. Every enabled branch writes to
    # its own plots/<branch_name>/ folder.
    branches = list(
      individual_review_png = list(
        enabled = TRUE,
        layout = "individual",
        format = "png",
        point_display = "both",            # "raw", "summary", or "both"
        nt50_display = "both",              # "primary", "ll4", "uniroot", "both", "none"
        show_comparison_details = TRUE,     # top-right LL4/uniroot detail box
        show_ll4_curve = TRUE,
        show_gen5_overlay = FALSE,
        show_model_confidence_band = FALSE,
        show_nt50_ci = TRUE,
        y_axis = "data",                    # "data" shows the full observed range
        y_limits = c(-10, 110),
        width = 8,
        height = 8,
        dpi = 300L,
        facet_columns = 2L,
        max_x_breaks = 8L,
        fallback_on_plot_error = TRUE
      ),
      faceted_summary_pdf = list(
        enabled = FALSE,
        layout = "faceted",
        format = "pdf",
        point_display = "summary",
        nt50_display = "primary",
        show_comparison_details = TRUE,
        show_ll4_curve = TRUE,
        show_gen5_overlay = FALSE,
        show_model_confidence_band = FALSE,
        show_nt50_ci = FALSE,
        y_axis = "data",
        y_limits = c(-10, 110),
        width = 8,
        height = 8,
        dpi = 300L,
        facet_columns = 2L,
        max_x_breaks = 8L,
        fallback_on_plot_error = TRUE
      )
    )
  )
)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    candidate <- sub("^--file=", "", file_arg[[1]])
    return(normalizePath(candidate, winslash = "/", mustWork = FALSE))
  }

  # When the file is run with source(), R commonly records the source path in an
  # `ofile` entry on one of the active frames.
  frame_paths <- vapply(
    sys.frames(),
    function(frame) {
      candidate <- frame$ofile %||% ""
      if (length(candidate) != 1L || is.na(candidate)) "" else as.character(candidate)
    },
    character(1)
  )
  frame_paths <- frame_paths[nzchar(frame_paths)]
  if (length(frame_paths) > 0L) {
    return(normalizePath(tail(frame_paths, 1L), winslash = "/", mustWork = FALSE))
  }

  # RStudio support is optional and does not add a required package dependency.
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      isTRUE(tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE))) {
    context <- tryCatch(
      rstudioapi::getActiveDocumentContext(),
      error = function(e) NULL
    )
    if (!is.null(context) && length(context$path) == 1L && nzchar(context$path)) {
      return(normalizePath(context$path, winslash = "/", mustWork = FALSE))
    }
  }

  NA_character_
}

SCRIPT_PATH <- get_script_path()
SCRIPT_DIRECTORY <- if (length(SCRIPT_PATH) == 1L && !is.na(SCRIPT_PATH) &&
                            nzchar(SCRIPT_PATH)) {
  normalizePath(dirname(SCRIPT_PATH), winslash = "/", mustWork = FALSE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

is_absolute_path <- function(path) {
  path <- trimws(as.character(path))
  grepl("^/", path) ||
    grepl("^[A-Za-z]:[/\\\\]", path) ||
    grepl("^[/\\\\]{2}", path)
}

resolve_portable_path <- function(path, base_dir = getwd(), must_work = FALSE) {
  if (is.null(path)) {
    return(NULL)
  }

  path <- path.expand(trimws(as.character(path)))
  if (!nzchar(path)) {
    return(NULL)
  }

  # A relative path copied from Windows may contain backslashes. On Linux or
  # macOS, treat those as directory separators rather than literal characters.
  if (.Platform$OS.type != "windows" && !is_absolute_path(path)) {
    path <- gsub("\\\\", "/", path)
  }

  if (!is_absolute_path(path)) {
    path <- file.path(base_dir, path)
  }

  # Do not reinterpret a Windows drive path as a relative Linux path. It will
  # subsequently fail with a clear file-not-found message if used on Linux.
  if (.Platform$OS.type != "windows" && grepl("^[A-Za-z]:[/\\\\]", path)) {
    return(path)
  }

  normalizePath(path, winslash = "/", mustWork = must_work)
}

ANALYSIS_RESULT_REQUIRED_FIELDS <- c(
  "plate_number",
  "plate_id",
  "sample_id",
  "virus",
  "primary_reporting_method",
  "reported_nt50",
  "ll4_calculated_nt50",
  "uniroot_calculated_nt50",
  "hard_qc_flags",
  "review_qc_flags",
  "analysis_script_version"
)

make_plot_branch <- function(
  enabled = FALSE,
  layout = "individual",
  format = "png",
  point_display = "both",
  nt50_display = "both",
  show_comparison_details = TRUE,
  show_ll4_curve = TRUE,
  show_gen5_overlay = FALSE,
  show_model_confidence_band = FALSE,
  show_nt50_ci = TRUE,
  y_axis = "data",
  y_limits = c(-10, 110),
  width = 8,
  height = 8,
  dpi = 300L,
  facet_columns = 2L,
  max_x_breaks = 8L,
  fallback_on_plot_error = TRUE
) {
  list(
    enabled = enabled,
    layout = layout,
    format = format,
    point_display = point_display,
    nt50_display = nt50_display,
    show_comparison_details = show_comparison_details,
    show_ll4_curve = show_ll4_curve,
    show_gen5_overlay = show_gen5_overlay,
    show_model_confidence_band = show_model_confidence_band,
    show_nt50_ci = show_nt50_ci,
    y_axis = y_axis,
    y_limits = y_limits,
    width = width,
    height = height,
    dpi = dpi,
    facet_columns = facet_columns,
    max_x_breaks = max_x_breaks,
    fallback_on_plot_error = fallback_on_plot_error
  )
}

BASE_CONFIG <- list(
  paths = list(
    set_working_directory_to_script_folder = FALSE,
    input_folder = "."
  ),

  input = list(
    raw_counts_file = NULL,
    lookup_file = NULL,
    gen5_fits_file = NULL
  ),

  normalization = list(
    infer_control_virus_from_well_row = FALSE,

    # Strict built-in defaults for command-line mode. Cross-plate controls are
    # used only when explicitly requested with --allow-cross-plate-vc true.
    allow_cross_plate_virus_controls = FALSE,
    cross_plate_vc_plate_numbers = NULL,
    allow_cross_plate_vc_for_reporting = FALSE,
    pool_equal_weight_by_plate = TRUE
  ),

  output_dir = "directview_fra_results",

  analysis = list(
    threshold_pct = 50,
    reporting_method = "ll4",
    uniroot_interpolation_scale = "linear",
    ll4_ci_method = "delta",
    confidence_level = 0.95,
    min_unique_dilutions = 4L,
    min_total_observations = 8L,
    prediction_grid_points = 400L
  ),

  qc = list(
    vc_count_range = c(100, 1600),
    hill_range = c(0.5, 2.0),
    fixed_hill = 1.0,
    vc_cv_max_pct = NULL,
    max_replicate_sd_pp = NULL,
    cell_control_count_max = NULL,
    lack_of_fit_alpha = NULL,
    min_ll4_pseudo_r2 = NULL,
    max_ll4_rmse_pp = NULL,
    max_upward_mean_steps = NULL,
    max_ll4_uniroot_disagreement_fold = NULL,
    max_ll4_gen5_disagreement_fold = NULL,
    min_vc_replicates = 2L,
    require_min_vc_replicates_for_reporting = FALSE
  ),

  plot = list(
    clear_enabled_branch_directories = TRUE,
    flag_wrap_width = 72L,
    blank_lines_between_flags = 1L,
    individual_height_per_flag_line = 0.30,
    max_individual_plot_height = 48,
    branches = list(
      individual_review_png = make_plot_branch(
        enabled = TRUE,
        layout = "individual",
        format = "png",
        point_display = "both",
        nt50_display = "both",
        show_comparison_details = TRUE,
        show_ll4_curve = TRUE,
        show_gen5_overlay = FALSE,
        show_model_confidence_band = FALSE,
        show_nt50_ci = TRUE
      ),
      faceted_summary_pdf = make_plot_branch(
        enabled = FALSE,
        layout = "faceted",
        format = "pdf",
        point_display = "summary",
        nt50_display = "primary",
        show_comparison_details = TRUE,
        show_ll4_curve = TRUE,
        show_gen5_overlay = FALSE,
        show_model_confidence_band = FALSE,
        show_nt50_ci = FALSE
      )
    )
  )
)

build_initial_config <- function(use_hardcoded_settings = USE_HARDCODED_SETTINGS) {
  if (!is.logical(use_hardcoded_settings) || length(use_hardcoded_settings) != 1L ||
      is.na(use_hardcoded_settings)) {
    stop("use_hardcoded_settings must be TRUE or FALSE.", call. = FALSE)
  }

  config <- if (isTRUE(use_hardcoded_settings)) {
    utils::modifyList(BASE_CONFIG, HARDCODED_SETTINGS, keep.null = TRUE)
  } else {
    BASE_CONFIG
  }
  config
}


print_usage <- function() {
  cat(
    paste0(
      "DirectView-FRA analysis script v", SCRIPT_VERSION, "\n\n",
      "Usage:\n",
      "  Rscript directview_fra_analysis.R [options]\n\n",
      "Run mode:\n",
      "  --use-hardcoded-settings BOOL  TRUE uses HARDCODED_SETTINGS; FALSE ignores\n",
      "                                 that block and starts from built-in defaults.\n\n",
      "Input and output:\n",
      "  --input-folder PATH            Folder containing input CSV files.\n",
      "  --set-working-directory BOOL   Reset to the script folder.\n",
      "  --raw-counts FILE_OR_PATH      Required per-well raw-count CSV.\n",
      "  --lookup FILE_OR_PATH|none     Optional plate-layout lookup CSV.\n",
      "  --gen5-fits FILE_OR_PATH|none  Optional Gen5/imaging-system LL4 CSV.\n",
      "  --output-dir PATH              Output folder; relative paths use input-folder.\n\n",
      "Analysis and normalization:\n",
      "  --reporting-method VALUE       ll4 or uniroot. Both estimates remain in CSV.\n",
      "  --uniroot-scale VALUE          linear or log2 interpolation of observed means.\n",
      "  --infer-control-virus-from-row BOOL  Match VCs to the sample virus in the same row.\n",
      "  --allow-cross-plate-vc BOOL    Allow same-virus VCs from other plates.\n",
      "  --cross-plate-vc-plates LIST|none  Optional comma-separated source plates.\n",
      "  --allow-cross-plate-vc-reporting BOOL  Permit validated cross-plate normalization for reporting.\n",
      "  --ll4-ci VALUE                 none, delta, or inv.\n",
      "  --plot-branches VALUE          Comma-separated configured branch names,\n",
      "                                 or none to disable all plots.\n",
      "  --show-comparison-details BOOL Show or hide the top-right method detail box\n",
      "                                 for all configured plot branches.\n",
      "  --help                         Show this help text.\n\n",
      "The script defaults to USE_HARDCODED_SETTINGS <- TRUE for Jupyter/RStudio.\n",
      "Explicit command-line values override the selected starting configuration.\n",
      "Set the top switch to FALSE, or pass --use-hardcoded-settings false, to\n",
      "ignore the hardcoded bracket. In that mode, pass --raw-counts explicitly.\n",
      "Relative inputs are resolved inside --input-folder. Base R path functions\n",
      "provide Windows and Linux compatibility.\n"
    )
  )
}

parse_cli_args <- function(args) {
  parsed <- list()
  i <- 1L

  while (i <= length(args)) {
    arg <- args[[i]]

    if (identical(arg, "--help")) {
      parsed$help <- TRUE
      i <- i + 1L
      next
    }

    if (!startsWith(arg, "--")) {
      stop("Unexpected positional argument: ", arg, call. = FALSE)
    }

    arg_body <- substring(arg, 3L)
    if (grepl("=", arg_body, fixed = TRUE)) {
      pieces <- strsplit(arg_body, "=", fixed = TRUE)[[1]]
      key <- pieces[[1]]
      value <- paste(pieces[-1], collapse = "=")
      parsed[[gsub("-", "_", key, fixed = TRUE)]] <- value
      i <- i + 1L
      next
    }

    key <- gsub("-", "_", arg_body, fixed = TRUE)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      parsed[[key]] <- TRUE
      i <- i + 1L
    } else {
      parsed[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }

  parsed
}

parse_boolean <- function(x, option_name) {
  if (is.logical(x) && length(x) == 1L) {
    return(x)
  }

  value <- tolower(trimws(as.character(x)))
  if (value %in% c("true", "t", "yes", "y", "1")) {
    return(TRUE)
  }
  if (value %in% c("false", "f", "no", "n", "0")) {
    return(FALSE)
  }

  stop("Invalid boolean for --", option_name, ": ", x, call. = FALSE)
}

null_if_none <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (tolower(trimws(as.character(x))) %in% c("none", "null", "na", "")) {
    return(NULL)
  }
  as.character(x)
}

apply_cli_overrides <- function(config, cli) {
  if (isTRUE(cli$help %||% FALSE)) {
    print_usage()
    quit(save = "no", status = 0L)
  }

  if (!is.null(cli$input_folder)) {
    config$paths$input_folder <- as.character(cli$input_folder)
  }
  if (!is.null(cli$set_working_directory)) {
    config$paths$set_working_directory_to_script_folder <- parse_boolean(
      cli$set_working_directory,
      "set-working-directory"
    )
  }
  if (!is.null(cli$raw_counts)) {
    config$input$raw_counts_file <- as.character(cli$raw_counts)
  }
  if (!is.null(cli$lookup)) {
    config$input$lookup_file <- null_if_none(cli$lookup)
  }
  if (!is.null(cli$gen5_fits)) {
    config$input$gen5_fits_file <- null_if_none(cli$gen5_fits)
  }
  if (!is.null(cli$output_dir)) {
    config$output_dir <- as.character(cli$output_dir)
  }
  if (!is.null(cli$reporting_method)) {
    config$analysis$reporting_method <- as.character(cli$reporting_method)
  }
  if (!is.null(cli$uniroot_scale)) {
    config$analysis$uniroot_interpolation_scale <- as.character(cli$uniroot_scale)
  }
  if (!is.null(cli$infer_control_virus_from_row)) {
    config$normalization$infer_control_virus_from_well_row <- parse_boolean(
      cli$infer_control_virus_from_row,
      "infer-control-virus-from-row"
    )
  }
  if (!is.null(cli$allow_cross_plate_vc)) {
    config$normalization$allow_cross_plate_virus_controls <- parse_boolean(
      cli$allow_cross_plate_vc,
      "allow-cross-plate-vc"
    )
  }
  if (!is.null(cli$cross_plate_vc_plates)) {
    cross_plate_text <- trimws(as.character(cli$cross_plate_vc_plates))
    config$normalization$cross_plate_vc_plate_numbers <-
      if (tolower(cross_plate_text) %in% c("none", "null", "")) {
        NULL
      } else {
        unique(trimws(strsplit(cross_plate_text, ",", fixed = TRUE)[[1]]))
      }
  }
  if (!is.null(cli$allow_cross_plate_vc_reporting)) {
    config$normalization$allow_cross_plate_vc_for_reporting <- parse_boolean(
      cli$allow_cross_plate_vc_reporting,
      "allow-cross-plate-vc-reporting"
    )
  }
  if (!is.null(cli$ll4_ci)) {
    config$analysis$ll4_ci_method <- as.character(cli$ll4_ci)
  }

  if (!is.null(cli$show_comparison_details)) {
    comparison_value <- parse_boolean(
      cli$show_comparison_details,
      "show-comparison-details"
    )
    for (branch_name in names(config$plot$branches)) {
      config$plot$branches[[branch_name]]$show_comparison_details <- comparison_value
    }
  }

  if (!is.null(cli$plot_branches)) {
    requested_text <- trimws(as.character(cli$plot_branches))
    requested <- if (tolower(requested_text) %in% c("none", "null", "")) {
      character()
    } else {
      unique(trimws(strsplit(requested_text, ",", fixed = TRUE)[[1]]))
    }

    available <- names(config$plot$branches)
    unknown_branches <- setdiff(requested, available)
    if (length(unknown_branches) > 0L) {
      stop(
        "Unknown plot branch(es): ",
        paste(unknown_branches, collapse = ", "),
        ". Configured branches are: ",
        paste(available, collapse = ", "),
        call. = FALSE
      )
    }

    for (branch_name in available) {
      config$plot$branches[[branch_name]]$enabled <- branch_name %in% requested
    }
  }

  known <- c(
    "help", "use_hardcoded_settings", "input_folder",
    "set_working_directory", "raw_counts", "lookup", "gen5_fits",
    "output_dir", "reporting_method", "uniroot_scale",
    "infer_control_virus_from_row", "allow_cross_plate_vc",
    "cross_plate_vc_plates", "allow_cross_plate_vc_reporting",
    "ll4_ci", "plot_branches", "show_comparison_details"
  )
  unknown <- setdiff(names(cli), known)
  if (length(unknown) > 0L) {
    stop(
      "Unknown option(s): ",
      paste0("--", gsub("_", "-", unknown, fixed = TRUE), collapse = ", "),
      call. = FALSE
    )
  }

  config
}


build_run_config <- function(cli = list()) {
  use_hardcoded <- USE_HARDCODED_SETTINGS
  if (!is.null(cli$use_hardcoded_settings)) {
    use_hardcoded <- parse_boolean(
      cli$use_hardcoded_settings,
      "use-hardcoded-settings"
    )
  }

  config <- build_initial_config(use_hardcoded)

  # Explicit command-line values override the selected starting configuration.
  # This permits a CLI user to pass only the values that differ from the top
  # block, or to ignore that block entirely with --use-hardcoded-settings false.
  config <- apply_cli_overrides(config, cli)
  override_names <- setdiff(names(cli), c("help", "use_hardcoded_settings"))
  config$run_mode <- list(
    use_hardcoded_settings = isTRUE(use_hardcoded),
    configuration_source = if (isTRUE(use_hardcoded)) {
      "hardcoded_settings"
    } else {
      "command_line_defaults"
    },
    cli_overrides = paste(override_names, collapse = ", ")
  )
  config
}

validate_choice <- function(value, choices, name) {
  value <- tolower(trimws(as.character(value)))
  if (!value %in% choices) {
    stop(
      name,
      " must be one of: ",
      paste(choices, collapse = ", "),
      ". Received: ",
      value,
      call. = FALSE
    )
  }
  value
}

validate_config <- function(config) {
  config$analysis$reporting_method <- validate_choice(
    config$analysis$reporting_method,
    c("ll4", "uniroot"),
    "analysis$reporting_method"
  )
  config$analysis$uniroot_interpolation_scale <- validate_choice(
    config$analysis$uniroot_interpolation_scale,
    c("linear", "log2"),
    "analysis$uniroot_interpolation_scale"
  )
  config$analysis$ll4_ci_method <- validate_choice(
    config$analysis$ll4_ci_method,
    c("none", "delta", "inv"),
    "analysis$ll4_ci_method"
  )

  validate_finite_scalar <- function(value, name) {
    value <- suppressWarnings(as.numeric(value))
    if (length(value) != 1L || !is.finite(value)) {
      stop(name, " must be one finite numeric value.", call. = FALSE)
    }
    value
  }

  validate_positive_integer <- function(value, name) {
    value <- validate_finite_scalar(value, name)
    if (value < 1 || value != floor(value)) {
      stop(name, " must be a positive integer.", call. = FALSE)
    }
    as.integer(value)
  }

  validate_nonnegative_integer <- function(value, name) {
    value <- validate_finite_scalar(value, name)
    if (value < 0 || value != floor(value)) {
      stop(name, " must be a nonnegative integer.", call. = FALSE)
    }
    as.integer(value)
  }

  validate_optional_positive <- function(value, name, allow_zero = FALSE) {
    if (is.null(value)) {
      return(NULL)
    }
    value <- validate_finite_scalar(value, name)
    invalid <- if (allow_zero) value < 0 else value <= 0
    if (invalid) {
      stop(
        name,
        " must be NULL or a ",
        if (allow_zero) "nonnegative" else "positive",
        " number.",
        call. = FALSE
      )
    }
    value
  }

  validate_boolean_scalar <- function(value, name) {
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(name, " must be TRUE or FALSE.", call. = FALSE)
    }
    value
  }

  validate_nonblank_path_setting <- function(value, name, optional = FALSE) {
    if (is.null(value)) {
      if (optional) return(NULL)
      stop(name, " may not be NULL.", call. = FALSE)
    }
    value <- trimws(as.character(value))
    if (length(value) != 1L || is.na(value) || !nzchar(value)) {
      stop(name, " must be one nonblank path or filename.", call. = FALSE)
    }
    value
  }

  config$paths$set_working_directory_to_script_folder <- validate_boolean_scalar(
    config$paths$set_working_directory_to_script_folder,
    "paths$set_working_directory_to_script_folder"
  )
  config$paths$input_folder <- validate_nonblank_path_setting(
    config$paths$input_folder,
    "paths$input_folder"
  )
  config$input$raw_counts_file <- validate_nonblank_path_setting(
    config$input$raw_counts_file,
    "input$raw_counts_file"
  )
  config$input$lookup_file <- validate_nonblank_path_setting(
    config$input$lookup_file,
    "input$lookup_file",
    optional = TRUE
  )
  config$input$gen5_fits_file <- validate_nonblank_path_setting(
    config$input$gen5_fits_file,
    "input$gen5_fits_file",
    optional = TRUE
  )
  config$output_dir <- validate_nonblank_path_setting(
    config$output_dir,
    "output_dir"
  )

  threshold <- validate_finite_scalar(
    config$analysis$threshold_pct,
    "analysis$threshold_pct"
  )
  if (threshold <= 0 || threshold >= 100) {
    stop("analysis$threshold_pct must be between 0 and 100.", call. = FALSE)
  }
  config$analysis$threshold_pct <- threshold

  confidence_level <- validate_finite_scalar(
    config$analysis$confidence_level,
    "analysis$confidence_level"
  )
  if (confidence_level <= 0 || confidence_level >= 1) {
    stop("analysis$confidence_level must be between 0 and 1.", call. = FALSE)
  }
  config$analysis$confidence_level <- confidence_level

  config$normalization$infer_control_virus_from_well_row <-
    validate_boolean_scalar(
      config$normalization$infer_control_virus_from_well_row,
      "normalization$infer_control_virus_from_well_row"
    )
  config$normalization$allow_cross_plate_virus_controls <-
    validate_boolean_scalar(
      config$normalization$allow_cross_plate_virus_controls,
      "normalization$allow_cross_plate_virus_controls"
    )
  config$normalization$allow_cross_plate_vc_for_reporting <-
    validate_boolean_scalar(
      config$normalization$allow_cross_plate_vc_for_reporting,
      "normalization$allow_cross_plate_vc_for_reporting"
    )
  config$normalization$pool_equal_weight_by_plate <- validate_boolean_scalar(
    config$normalization$pool_equal_weight_by_plate,
    "normalization$pool_equal_weight_by_plate"
  )
  cross_plate_plates <- config$normalization$cross_plate_vc_plate_numbers
  if (!is.null(cross_plate_plates)) {
    cross_plate_plates <- unique(trimws(as.character(cross_plate_plates)))
    cross_plate_plates <- cross_plate_plates[
      nzchar(cross_plate_plates) & !is.na(cross_plate_plates)
    ]
    if (length(cross_plate_plates) == 0L) cross_plate_plates <- NULL
  }
  config$normalization$cross_plate_vc_plate_numbers <- cross_plate_plates

  config$analysis$min_unique_dilutions <- validate_positive_integer(
    config$analysis$min_unique_dilutions,
    "analysis$min_unique_dilutions"
  )
  config$analysis$min_total_observations <- validate_positive_integer(
    config$analysis$min_total_observations,
    "analysis$min_total_observations"
  )
  config$analysis$prediction_grid_points <- validate_positive_integer(
    config$analysis$prediction_grid_points,
    "analysis$prediction_grid_points"
  )
  if (config$analysis$prediction_grid_points < 25L) {
    stop("analysis$prediction_grid_points must be at least 25.", call. = FALSE)
  }

  config$qc$vc_count_range <- suppressWarnings(as.numeric(config$qc$vc_count_range))
  if (length(config$qc$vc_count_range) != 2L ||
      any(!is.finite(config$qc$vc_count_range)) ||
      config$qc$vc_count_range[[1]] < 0 ||
      config$qc$vc_count_range[[1]] >= config$qc$vc_count_range[[2]]) {
    stop(
      "qc$vc_count_range must contain two increasing nonnegative finite values.",
      call. = FALSE
    )
  }

  config$qc$hill_range <- suppressWarnings(as.numeric(config$qc$hill_range))
  if (length(config$qc$hill_range) != 2L ||
      any(!is.finite(config$qc$hill_range)) ||
      config$qc$hill_range[[1]] <= 0 ||
      config$qc$hill_range[[1]] >= config$qc$hill_range[[2]]) {
    stop(
      "qc$hill_range must contain two increasing positive finite values.",
      call. = FALSE
    )
  }

  config$qc$fixed_hill <- validate_finite_scalar(config$qc$fixed_hill, "qc$fixed_hill")
  if (config$qc$fixed_hill <= 0) {
    stop("qc$fixed_hill must be positive.", call. = FALSE)
  }

  config$qc$vc_cv_max_pct <- validate_optional_positive(
    config$qc$vc_cv_max_pct,
    "qc$vc_cv_max_pct"
  )
  config$qc$max_replicate_sd_pp <- validate_optional_positive(
    config$qc$max_replicate_sd_pp,
    "qc$max_replicate_sd_pp"
  )
  config$qc$cell_control_count_max <- validate_optional_positive(
    config$qc$cell_control_count_max,
    "qc$cell_control_count_max",
    allow_zero = TRUE
  )
  config$qc$max_ll4_rmse_pp <- validate_optional_positive(
    config$qc$max_ll4_rmse_pp,
    "qc$max_ll4_rmse_pp"
  )
  config$qc$max_ll4_uniroot_disagreement_fold <- validate_optional_positive(
    config$qc$max_ll4_uniroot_disagreement_fold,
    "qc$max_ll4_uniroot_disagreement_fold"
  )
  config$qc$max_ll4_gen5_disagreement_fold <- validate_optional_positive(
    config$qc$max_ll4_gen5_disagreement_fold,
    "qc$max_ll4_gen5_disagreement_fold"
  )
  for (field in c(
    "max_ll4_uniroot_disagreement_fold",
    "max_ll4_gen5_disagreement_fold"
  )) {
    value <- config$qc[[field]]
    if (!is.null(value) && value < 1) {
      stop("qc$", field, " must be NULL or at least 1.", call. = FALSE)
    }
  }

  if (!is.null(config$qc$min_ll4_pseudo_r2)) {
    config$qc$min_ll4_pseudo_r2 <- validate_finite_scalar(
      config$qc$min_ll4_pseudo_r2,
      "qc$min_ll4_pseudo_r2"
    )
    if (config$qc$min_ll4_pseudo_r2 < -1 || config$qc$min_ll4_pseudo_r2 > 1) {
      stop("qc$min_ll4_pseudo_r2 must be NULL or between -1 and 1.", call. = FALSE)
    }
  }

  if (!is.null(config$qc$max_upward_mean_steps)) {
    config$qc$max_upward_mean_steps <- validate_nonnegative_integer(
      config$qc$max_upward_mean_steps,
      "qc$max_upward_mean_steps"
    )
  }

  if (!is.null(config$qc$lack_of_fit_alpha)) {
    config$qc$lack_of_fit_alpha <- validate_finite_scalar(
      config$qc$lack_of_fit_alpha,
      "qc$lack_of_fit_alpha"
    )
    if (config$qc$lack_of_fit_alpha <= 0 ||
        config$qc$lack_of_fit_alpha >= 1) {
      stop(
        "qc$lack_of_fit_alpha must be NULL or between 0 and 1.",
        call. = FALSE
      )
    }
  }
  config$qc$min_vc_replicates <- validate_positive_integer(
    config$qc$min_vc_replicates,
    "qc$min_vc_replicates"
  )
  config$qc$require_min_vc_replicates_for_reporting <- validate_boolean_scalar(
    config$qc$require_min_vc_replicates_for_reporting,
    "qc$require_min_vc_replicates_for_reporting"
  )

  config$plot$clear_enabled_branch_directories <- validate_boolean_scalar(
    config$plot$clear_enabled_branch_directories,
    "plot$clear_enabled_branch_directories"
  )
  config$plot$flag_wrap_width <- validate_positive_integer(
    config$plot$flag_wrap_width,
    "plot$flag_wrap_width"
  )
  config$plot$blank_lines_between_flags <- validate_nonnegative_integer(
    config$plot$blank_lines_between_flags,
    "plot$blank_lines_between_flags"
  )
  config$plot$individual_height_per_flag_line <- validate_optional_positive(
    config$plot$individual_height_per_flag_line,
    "plot$individual_height_per_flag_line"
  )
  config$plot$max_individual_plot_height <- validate_optional_positive(
    config$plot$max_individual_plot_height,
    "plot$max_individual_plot_height"
  )

  branches <- config$plot$branches
  if (!is.list(branches) || length(branches) == 0L || is.null(names(branches)) ||
      any(names(branches) == "") || anyDuplicated(names(branches))) {
    stop("plot$branches must be a nonempty, uniquely named list.", call. = FALSE)
  }
  invalid_branch_names <- !grepl("^[A-Za-z0-9._-]+$", names(branches))
  if (any(invalid_branch_names)) {
    stop(
      "Plot branch names may contain only letters, numbers, dot, underscore, and hyphen: ",
      paste(names(branches)[invalid_branch_names], collapse = ", "),
      call. = FALSE
    )
  }

  for (branch_name in names(branches)) {
    branch <- branches[[branch_name]]
    prefix <- paste0("plot$branches$", branch_name, "$")

    branch$enabled <- validate_boolean_scalar(branch$enabled, paste0(prefix, "enabled"))
    branch$layout <- validate_choice(
      branch$layout,
      c("individual", "faceted"),
      paste0(prefix, "layout")
    )
    branch$format <- validate_choice(
      branch$format,
      c("png", "pdf"),
      paste0(prefix, "format")
    )
    branch$point_display <- validate_choice(
      branch$point_display,
      c("raw", "summary", "both"),
      paste0(prefix, "point_display")
    )
    branch$nt50_display <- validate_choice(
      branch$nt50_display,
      c("primary", "ll4", "uniroot", "both", "none"),
      paste0(prefix, "nt50_display")
    )
    branch$y_axis <- validate_choice(
      branch$y_axis,
      c("fixed", "data"),
      paste0(prefix, "y_axis")
    )

    for (field in c(
      "show_comparison_details", "show_ll4_curve", "show_gen5_overlay",
      "show_model_confidence_band", "show_nt50_ci", "fallback_on_plot_error"
    )) {
      branch[[field]] <- validate_boolean_scalar(
        branch[[field]],
        paste0(prefix, field)
      )
    }

    branch$y_limits <- suppressWarnings(as.numeric(branch$y_limits))
    if (length(branch$y_limits) != 2L ||
        any(!is.finite(branch$y_limits)) ||
        branch$y_limits[[1]] >= branch$y_limits[[2]]) {
      stop(paste0(prefix, "y_limits must contain two increasing finite values."), call. = FALSE)
    }

    branch$width <- validate_optional_positive(branch$width, paste0(prefix, "width"))
    branch$height <- validate_optional_positive(branch$height, paste0(prefix, "height"))
    branch$dpi <- validate_positive_integer(branch$dpi, paste0(prefix, "dpi"))
    branch$facet_columns <- validate_positive_integer(
      branch$facet_columns,
      paste0(prefix, "facet_columns")
    )
    branch$max_x_breaks <- validate_positive_integer(
      branch$max_x_breaks,
      paste0(prefix, "max_x_breaks")
    )

    branches[[branch_name]] <- branch
  }
  config$plot$branches <- branches

  individual_heights <- vapply(
    branches[vapply(branches, function(branch) branch$layout == "individual", logical(1))],
    function(branch) branch$height,
    numeric(1)
  )
  if (length(individual_heights) > 0L &&
      config$plot$max_individual_plot_height < max(individual_heights)) {
    stop(
      "plot$max_individual_plot_height must be at least as large as every individual branch height.",
      call. = FALSE
    )
  }

  enabled_names <- names(branches)[vapply(
    branches,
    function(branch) isTRUE(branch$enabled),
    logical(1)
  )]
  if (length(enabled_names) > 1L) {
    signatures <- vapply(
      branches[enabled_names],
      function(branch) {
        branch$enabled <- NULL
        branch <- branch[sort(names(branch))]
        paste(capture.output(dput(branch)), collapse = "")
      },
      character(1)
    )
    duplicated_signatures <- unique(signatures[duplicated(signatures)])
    if (length(duplicated_signatures) > 0L) {
      duplicate_groups <- vapply(
        duplicated_signatures,
        function(signature) paste(enabled_names[signatures == signature], collapse = ", "),
        character(1)
      )
      stop(
        "Enabled plot branches may not request exact duplicate outputs. Duplicate configuration group(s): ",
        paste(duplicate_groups, collapse = " | "),
        call. = FALSE
      )
    }
  }

  config
}

prepare_runtime_paths <- function(config) {
  original_working_directory <- normalizePath(
    getwd(),
    winslash = "/",
    mustWork = FALSE
  )

  if (isTRUE(config$paths$set_working_directory_to_script_folder)) {
    if (!dir.exists(SCRIPT_DIRECTORY)) {
      stop("Script directory does not exist: ", SCRIPT_DIRECTORY, call. = FALSE)
    }
    setwd(SCRIPT_DIRECTORY)
  }

  active_working_directory <- normalizePath(
    getwd(),
    winslash = "/",
    mustWork = FALSE
  )
  input_folder <- resolve_portable_path(
    config$paths$input_folder,
    base_dir = active_working_directory,
    must_work = FALSE
  )
  if (is.null(input_folder) || !dir.exists(input_folder)) {
    stop(
      "Input folder does not exist: ",
      input_folder %||% "<blank>",
      call. = FALSE
    )
  }
  input_folder <- normalizePath(input_folder, winslash = "/", mustWork = TRUE)

  raw_counts_value <- config$input$raw_counts_file
  if (is.null(raw_counts_value) ||
      !nzchar(trimws(as.character(raw_counts_value)))) {
    stop(
      "A raw-count CSV is required. Edit HARDCODED_SETTINGS$input$raw_counts_file ",
      "or, in command-line mode, pass --raw-counts FILE_OR_PATH.",
      call. = FALSE
    )
  }

  resolve_input_file <- function(value) {
    resolve_portable_path(value, base_dir = input_folder, must_work = FALSE)
  }

  config$paths$input_folder <- input_folder
  config$input$raw_counts_file <- resolve_input_file(config$input$raw_counts_file)
  config$input$lookup_file <- resolve_input_file(config$input$lookup_file)
  config$input$gen5_fits_file <- resolve_input_file(config$input$gen5_fits_file)
  config$output_dir <- resolve_portable_path(
    config$output_dir,
    base_dir = input_folder,
    must_work = FALSE
  )

  config$runtime <- list(
    script_path = SCRIPT_PATH,
    script_directory = SCRIPT_DIRECTORY,
    original_working_directory = original_working_directory,
    working_directory = active_working_directory,
    working_directory_reset_to_script_folder =
      config$paths$set_working_directory_to_script_folder,
    use_hardcoded_settings = isTRUE(
      config$run_mode$use_hardcoded_settings %||%
        identical(config$settings_source %||% "", "hardcoded")
    ),
    configuration_source = config$run_mode$configuration_source %||%
      config$settings_source %||% "custom_config",
    cli_overrides = config$run_mode$cli_overrides %||% ""
  )

  config
}

check_r_version <- function() {
  if (utils::compareVersion(as.character(getRversion()), MIN_R_VERSION) < 0L) {
    stop(
      "This script requires R ", MIN_R_VERSION, " or newer because it uses the base R pipe (|>). ",
      "Detected R ", as.character(getRversion()), ".",
      call. = FALSE
    )
  }
}

check_required_packages <- function() {
  required <- c("drc", "dplyr", "ggplot2", "readr", "scales", "tibble")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing) > 0L) {
    stop(
      "Missing required R package(s): ",
      paste(missing, collapse = ", "),
      ". Install them before running this script.",
      call. = FALSE
    )
  }
}

normalize_header <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

rename_by_alias <- function(data, alias_map, required, context) {
  original_names <- names(data)
  normalized_names <- normalize_header(original_names)

  if (anyDuplicated(normalized_names)) {
    duplicates <- unique(normalized_names[duplicated(normalized_names)])
    stop(
      context,
      " contains column names that become duplicates after normalization: ",
      paste(duplicates, collapse = ", "),
      call. = FALSE
    )
  }

  names(data) <- normalized_names

  for (canonical in names(alias_map)) {
    candidates <- unique(normalize_header(alias_map[[canonical]]))
    hits <- intersect(candidates, names(data))

    if (length(hits) > 1L) {
      stop(
        context,
        " has multiple columns that could map to '",
        canonical,
        "': ",
        paste(hits, collapse = ", "),
        call. = FALSE
      )
    }

    if (length(hits) == 1L && !identical(hits, canonical)) {
      names(data)[names(data) == hits] <- canonical
    }
  }

  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(
      context,
      " is missing required column(s): ",
      paste(missing, collapse = ", "),
      ". Original columns were: ",
      paste(original_names, collapse = ", "),
      call. = FALSE
    )
  }

  data
}

parse_numeric_strict <- function(x, field_name, row_ids = seq_along(x)) {
  text <- trimws(as.character(x))
  missing_text <- is.na(x) | text == "" | toupper(text) %in% c("NA", "N/A", "NULL")
  parsed <- suppressWarnings(as.numeric(text))
  invalid <- !missing_text & !is.finite(parsed)

  if (any(invalid)) {
    bad_rows <- head(row_ids[invalid], 10L)
    stop(
      "Column '",
      field_name,
      "' contains nonnumeric value(s) at source row(s): ",
      paste(bad_rows, collapse = ", "),
      if (sum(invalid) > 10L) " ..." else "",
      call. = FALSE
    )
  }

  parsed[missing_text] <- NA_real_
  parsed
}

trim_character_columns <- function(data) {
  is_character <- vapply(data, is.character, logical(1))
  data[is_character] <- lapply(data[is_character], trimws)
  data
}

read_csv_minimal <- function(path, context) {
  if (is.null(path)) {
    return(NULL)
  }
  if (!file.exists(path)) {
    stop(context, " file does not exist: ", path, call. = FALSE)
  }

  readr::read_csv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    na = c("", "NA", "N/A", "NULL"),
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal"
  )
}

raw_aliases <- list(
  well_id = c("well_id", "wellid", "local_name", "localname"),
  well = c("well", "well_position", "well_name"),
  dilution = c(
    "conc_dil", "concentration_dilution", "dilution", "dilution_factor",
    "reciprocal_dilution"
  ),
  focus_count = c("cell_count", "focus_count", "foci_count", "spot_count"),
  plate_id = c("plate_id", "plateid"),
  plate_number = c("plate_number", "plate_no", "platenumber"),
  virus = c("virus", "virus_name"),
  sample_id = c("global_name", "sample_id", "sample", "sera", "antibody_id"),
  control_type = c(
    "control_type", "control_designation", "well_type", "sample_type"
  ),
  reading_datetime = c("reading_date_time", "reading_datetime", "read_datetime")
)

lookup_aliases <- list(
  well_id = c("well_id", "local_name", "localname"),
  plate_id = c("plate_id", "plateid"),
  plate_number = c("plate_number", "plate_no", "platenumber"),
  sample_id = c("sample_id", "global_name", "sera", "sample", "antibody_id"),
  virus = c("virus", "virus_name"),
  control_type = c(
    "control_type", "control_designation", "well_type", "sample_type"
  )
)

gen5_aliases <- list(
  curve_name = c("curve_name"),
  local_name = c("name", "local_name"),
  curve_formula = c("curve_formula", "y_formula"),
  a = c("a"),
  b = c("b"),
  c = c("c"),
  d = c("d"),
  r2 = c("r2", "r_squared"),
  fit_f_prob = c("fit_f_prob", "fit_f_probability"),
  plate_number = c("plate_number", "plate_no", "platenumber"),
  plate_id = c("plate_id", "plateid"),
  virus = c("virus", "virus_name"),
  sample_id = c("global_name", "sample_id", "sample", "sera", "antibody_id")
)

standardize_raw_counts <- function(path) {
  data <- read_csv_minimal(path, "Raw-count")
  if (nrow(data) == 0L) {
    stop("Raw-count CSV has no data rows.", call. = FALSE)
  }
  data$source_row <- seq_len(nrow(data)) + 1L

  # Virus and sample metadata can be supplied later by the lookup table, so they
  # are not required at this stage.
  data <- rename_by_alias(
    data,
    raw_aliases,
    required = c("well_id", "well", "dilution", "focus_count", "plate_id", "plate_number"),
    context = "Raw-count CSV"
  )
  data <- trim_character_columns(data)

  data$dilution <- parse_numeric_strict(data$dilution, "dilution", data$source_row)
  data$focus_count <- parse_numeric_strict(data$focus_count, "focus_count", data$source_row)

  if (!"virus" %in% names(data)) {
    data$virus <- NA_character_
  }
  if (!"sample_id" %in% names(data)) {
    data$sample_id <- NA_character_
  }
  if (!"control_type" %in% names(data)) {
    data$control_type <- NA_character_
  }

  data$virus <- as.character(data$virus)
  data$sample_id <- as.character(data$sample_id)
  data$control_type <- as.character(data$control_type)
  data
}

standardize_lookup <- function(path) {
  if (is.null(path)) {
    return(NULL)
  }

  data <- read_csv_minimal(path, "Lookup")
  if (nrow(data) == 0L) {
    stop("Lookup CSV has no data rows.", call. = FALSE)
  }
  data <- rename_by_alias(
    data,
    lookup_aliases,
    required = c("well_id", "plate_id", "sample_id", "virus"),
    context = "Lookup CSV"
  )
  data <- trim_character_columns(data)

  if (!"plate_number" %in% names(data)) {
    data$plate_number <- NA_character_
  }
  if (!"control_type" %in% names(data)) {
    data$control_type <- NA_character_
  }

  required_text <- c("well_id", "plate_id", "sample_id", "virus")
  for (field in required_text) {
    value <- trimws(as.character(data[[field]]))
    missing <- is.na(data[[field]]) | value == ""
    if (any(missing)) {
      stop(
        "Lookup CSV has missing '", field, "' value(s) at data row(s): ",
        paste(head(which(missing) + 1L, 10L), collapse = ", "),
        call. = FALSE
      )
    }
  }

  key <- paste(data$plate_id, data$well_id, sep = "\r")
  if (anyDuplicated(key)) {
    duplicate_keys <- unique(key[duplicated(key)])
    stop(
      "Lookup CSV has duplicate plate_id + well_id mapping(s): ",
      paste(head(gsub("\r", " / ", duplicate_keys, fixed = TRUE), 10L), collapse = ", "),
      call. = FALSE
    )
  }

  data
}

resolve_metadata_field <- function(raw_value, lookup_value, field_name, row_ids) {
  raw_text <- trimws(as.character(raw_value))
  lookup_text <- trimws(as.character(lookup_value))

  raw_missing <- is.na(raw_value) | raw_text == ""
  lookup_missing <- is.na(lookup_value) | lookup_text == ""
  conflict <- !raw_missing & !lookup_missing & raw_text != lookup_text

  if (any(conflict)) {
    rows <- head(row_ids[conflict], 10L)
    stop(
      "Raw-count and lookup metadata disagree for '",
      field_name,
      "' at source row(s): ",
      paste(rows, collapse = ", "),
      call. = FALSE
    )
  }

  result <- raw_text
  result[raw_missing] <- lookup_text[raw_missing]
  result[result == ""] <- NA_character_
  result
}

merge_lookup_metadata <- function(raw_data, lookup_data) {
  if (is.null(lookup_data)) {
    # Sample-level metadata requirements are checked after control classification,
    # so explicitly designated controls may omit sample or virus text.
    return(raw_data)
  }

  lookup_for_join <- lookup_data[, c(
    "plate_id", "well_id", "plate_number", "sample_id", "virus", "control_type"
  )]
  names(lookup_for_join) <- c(
    "plate_id", "well_id", "lookup_plate_number", "lookup_sample_id", "lookup_virus",
    "lookup_control_type"
  )

  merged <- dplyr::left_join(
    raw_data,
    lookup_for_join,
    by = c("plate_id", "well_id")
  )

  missing_mapping <- is.na(merged$lookup_sample_id) | is.na(merged$lookup_virus)
  if (any(missing_mapping)) {
    examples <- unique(paste(
      merged$plate_id[missing_mapping],
      merged$well_id[missing_mapping],
      sep = " / "
    ))
    stop(
      "Lookup CSV does not map every raw-count plate_id + well_id. Example(s): ",
      paste(head(examples, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  merged$sample_id <- resolve_metadata_field(
    merged$sample_id,
    merged$lookup_sample_id,
    "sample_id",
    merged$source_row
  )
  merged$virus <- resolve_metadata_field(
    merged$virus,
    merged$lookup_virus,
    "virus",
    merged$source_row
  )
  merged$control_type <- resolve_metadata_field(
    canonicalize_control_type_values(
      merged$control_type,
      merged$source_row,
      "Raw-count control designation"
    ),
    canonicalize_control_type_values(
      merged$lookup_control_type,
      merged$source_row,
      "Lookup control designation"
    ),
    "control_type",
    merged$source_row
  )

  plate_lookup_available <- !is.na(merged$lookup_plate_number) &
    trimws(merged$lookup_plate_number) != ""
  plate_conflict <- plate_lookup_available &
    trimws(merged$plate_number) != trimws(merged$lookup_plate_number)
  if (any(plate_conflict)) {
    rows <- head(merged$source_row[plate_conflict], 10L)
    stop(
      "Raw-count and lookup plate_number values disagree at source row(s): ",
      paste(rows, collapse = ", "),
      call. = FALSE
    )
  }

  merged$lookup_plate_number <- NULL
  merged$lookup_sample_id <- NULL
  merged$lookup_virus <- NULL
  merged$lookup_control_type <- NULL
  merged
}

canonicalize_control_type_values <- function(values, row_ids, context) {
  text <- tolower(trimws(as.character(values)))
  missing <- is.na(values) | text == ""
  key <- gsub("[^a-z0-9]+", "_", text)
  key <- gsub("^_+|_+$", "", key)

  result <- rep(NA_character_, length(key))
  result[key %in% c(
    "sample", "test", "test_sample", "serum", "sera", "antibody", "specimen"
  )] <- "sample"
  result[key %in% c("virus_control", "viruscontrol", "vc")] <- "virus_control"
  result[key %in% c("cell_control", "cellcontrol", "cc")] <- "cell_control"
  result[missing] <- NA_character_

  invalid <- !missing & is.na(result)
  if (any(invalid)) {
    stop(
      context,
      " contains unsupported value(s) at source row(s): ",
      paste(head(row_ids[invalid], 10L), collapse = ", "),
      ". Accepted categories are sample, virus_control/VC, and cell_control/CC.",
      call. = FALSE
    )
  }

  result
}

classify_controls <- function(data) {
  explicit <- canonicalize_control_type_values(
    data$control_type,
    data$source_row,
    "Control designation"
  )
  well_upper <- toupper(trimws(data$well_id))
  sample_lower <- tolower(trimws(data$sample_id))

  inferred <- rep("sample", nrow(data))
  inferred[
    well_upper %in% c("VC", "VIRUS CONTROL", "VIRUS_CONTROL") |
      grepl("^virus[ _-]*control$", sample_lower)
  ] <- "virus_control"
  inferred[
    well_upper %in% c("CC", "CELL CONTROL", "CELL_CONTROL") |
      grepl("^cell[ _-]*control$", sample_lower)
  ] <- "cell_control"

  explicit_available <- !is.na(explicit)
  inferred_control <- inferred %in% c("virus_control", "cell_control")
  conflict <- explicit_available & inferred_control & explicit != inferred
  if (any(conflict)) {
    stop(
      "Explicit control designation conflicts with the Well ID or sample name at source row(s): ",
      paste(head(data$source_row[conflict], 10L), collapse = ", "),
      call. = FALSE
    )
  }

  data$control_type <- ifelse(explicit_available, explicit, inferred)
  data$control_type_source <- ifelse(explicit_available, "explicit", "inferred")
  data
}

well_row_from_position <- function(well) {
  text <- toupper(trimws(as.character(well)))
  row <- sub("[0-9].*$", "", text)
  valid <- nzchar(row) & grepl("^[A-Z]+$", row)
  row[!valid] <- NA_character_
  row
}

assign_control_reference_viruses <- function(data, config) {
  data$well_row <- well_row_from_position(data$well)
  data$normalization_virus <- NA_character_
  data$control_virus_source <- "not_applicable"
  data$control_virus_mapping_flag <- ""

  sample_rows <- data$control_type == "sample"
  data$normalization_virus[sample_rows] <- trimws(as.character(data$virus[sample_rows]))
  data$control_virus_source[sample_rows] <- "sample_metadata"

  plate_ids <- unique(data$plate_id)
  for (plate_id in plate_ids) {
    on_plate <- data$plate_id == plate_id
    sample_on_plate <- on_plate & sample_rows
    sample_viruses <- sort(unique(trimws(as.character(data$virus[sample_on_plate]))))
    sample_viruses <- sample_viruses[!is.na(sample_viruses) & nzchar(sample_viruses)]

    vc_indices <- which(on_plate & data$control_type == "virus_control")
    for (index in vc_indices) {
      explicit_virus <- trimws(as.character(data$virus[[index]] %||% ""))
      if (is.na(explicit_virus)) explicit_virus <- ""

      row_candidates <- sort(unique(trimws(as.character(data$virus[
        sample_on_plate & data$well_row == data$well_row[[index]]
      ]))))
      row_candidates <- row_candidates[
        !is.na(row_candidates) & nzchar(row_candidates)
      ]
      unique_row_virus <- if (length(row_candidates) == 1L) {
        row_candidates[[1]]
      } else {
        ""
      }

      # When explicit metadata and the unique sample virus in that row agree,
      # the control can be matched without inference.
      if (nzchar(explicit_virus) && nzchar(unique_row_virus) &&
          identical(explicit_virus, unique_row_virus)) {
        data$normalization_virus[[index]] <- explicit_virus
        data$control_virus_source[[index]] <- "explicit_matches_row"
        next
      }

      # A row/metadata conflict is deliberately conservative. By default it is
      # left unresolved; an explicitly enabled row-inference mode may override
      # the metadata for a validated row-block plate layout.
      if (nzchar(explicit_virus) && nzchar(unique_row_virus) &&
          !identical(explicit_virus, unique_row_virus)) {
        if (isTRUE(config$normalization$infer_control_virus_from_well_row)) {
          data$normalization_virus[[index]] <- unique_row_virus
          data$control_virus_source[[index]] <-
            "row_inferred_overrode_control_metadata"
          data$control_virus_mapping_flag[[index]] <-
            "virus_control_metadata_overridden_by_row_inference"
        } else {
          data$control_virus_source[[index]] <-
            "control_metadata_conflicts_with_row_sample_virus"
          data$control_virus_mapping_flag[[index]] <-
            "virus_control_metadata_conflicts_with_row_sample_virus"
        }
        next
      }

      if (nzchar(explicit_virus) && explicit_virus %in% sample_viruses) {
        data$normalization_virus[[index]] <- explicit_virus
        data$control_virus_source[[index]] <- "explicit_control_metadata"
      } else if (isTRUE(config$normalization$infer_control_virus_from_well_row) &&
                 nzchar(unique_row_virus)) {
        data$normalization_virus[[index]] <- unique_row_virus
        data$control_virus_source[[index]] <- "row_inferred"
        data$control_virus_mapping_flag[[index]] <-
          "virus_control_virus_inferred_from_row"
      } else if (length(sample_viruses) == 1L) {
        data$normalization_virus[[index]] <- sample_viruses[[1]]
        data$control_virus_source[[index]] <- "single_plate_virus_inferred"
        data$control_virus_mapping_flag[[index]] <-
          "virus_control_virus_inferred_from_single_plate_virus"
      } else {
        data$control_virus_source[[index]] <- "unresolved"
        data$control_virus_mapping_flag[[index]] <-
          "virus_control_virus_unresolved"
      }
    }
  }

  data
}

validate_raw_counts <- function(data) {
  required_all_rows <- c("well_id", "well", "plate_id", "plate_number")
  for (field in required_all_rows) {
    value <- trimws(as.character(data[[field]]))
    missing <- is.na(data[[field]]) | value == ""
    if (any(missing)) {
      rows <- head(data$source_row[missing], 10L)
      stop(
        "Required field '", field, "' is missing at source row(s): ",
        paste(rows, collapse = ", "),
        call. = FALSE
      )
    }
  }

  sample_rows <- data$control_type == "sample"
  for (field in c("virus", "sample_id")) {
    value <- trimws(as.character(data[[field]]))
    missing <- sample_rows & (is.na(data[[field]]) | value == "")
    if (any(missing)) {
      rows <- head(data$source_row[missing], 10L)
      stop(
        "Sample field '", field, "' is missing at source row(s): ",
        paste(rows, collapse = ", "),
        call. = FALSE
      )
    }
  }

  invalid_count <- !is.finite(data$focus_count) | data$focus_count < 0
  if (any(invalid_count)) {
    rows <- head(data$source_row[invalid_count], 10L)
    stop(
      "focus_count must be finite and nonnegative at source row(s): ",
      paste(rows, collapse = ", "),
      call. = FALSE
    )
  }

  invalid_dilution <- sample_rows & (!is.finite(data$dilution) | data$dilution <= 0)
  if (any(invalid_dilution)) {
    rows <- head(data$source_row[invalid_dilution], 10L)
    stop(
      "Sample dilution must be finite and greater than zero at source row(s): ",
      paste(rows, collapse = ", "),
      call. = FALSE
    )
  }

  well_key <- paste(data$plate_id, data$well, sep = "\r")
  if (anyDuplicated(well_key)) {
    duplicate_keys <- unique(well_key[duplicated(well_key)])
    stop(
      "Duplicate physical well(s) detected within a plate: ",
      paste(head(gsub("\r", " / ", duplicate_keys, fixed = TRUE), 10L), collapse = ", "),
      call. = FALSE
    )
  }

  plate_map <- unique(data[, c("plate_id", "plate_number")])
  if (anyDuplicated(plate_map$plate_id)) {
    stop("A plate_id maps to more than one plate_number.", call. = FALSE)
  }

  data
}

# Natural ordering helper for labels such as "Plate 2" and "Plate 10" without
# adding another package dependency.
first_numeric_token <- function(x) {
  text <- as.character(x)
  match <- regexpr("[0-9]+(?:\\.[0-9]+)?", text, perl = TRUE)
  value <- rep(NA_real_, length(text))
  found <- match > 0L
  value[found] <- suppressWarnings(as.numeric(regmatches(text, match)[found]))
  value
}

natural_order <- function(x) {
  text <- trimws(as.character(x))
  number <- first_numeric_token(text)
  order(is.na(number), number, tolower(text), na.last = TRUE)
}

# Deterministic within-group replicate label. The physical well remains the
# authoritative trace-back field; replicate_index is only a convenience label
# within each plate/virus/sample/dilution group.
derive_replicate_index <- function(data) {
  well_text <- toupper(trimws(as.character(data$well)))
  row_label <- sub("[0-9]+$", "", well_text)
  column_text <- sub("^[A-Z]+", "", well_text)
  row_order <- match(row_label, LETTERS)
  column_order <- suppressWarnings(as.integer(column_text))
  row_order[is.na(row_order)] <- .Machine$integer.max
  column_order[is.na(column_order)] <- .Machine$integer.max

  data$.well_row_order <- row_order
  data$.well_column_order <- column_order

  data <- data |>
    dplyr::arrange(
      .data$plate_id,
      .data$normalization_virus,
      .data$virus,
      .data$sample_id,
      .data$control_type,
      .data$dilution,
      .data$.well_row_order,
      .data$.well_column_order,
      .data$well
    ) |>
    dplyr::group_by(
      .data$plate_number,
      .data$plate_id,
      .data$normalization_virus,
      .data$virus,
      .data$sample_id,
      .data$control_type,
      .data$dilution
    ) |>
    dplyr::mutate(replicate_index = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::arrange(.data$source_row)

  data$.well_row_order <- NULL
  data$.well_column_order <- NULL
  data
}

standardize_gen5_fits <- function(path) {
  if (is.null(path)) {
    return(NULL)
  }

  data <- read_csv_minimal(path, "Gen5 fit")
  if (nrow(data) == 0L) {
    warning("Gen5 fit CSV has no data rows; continuing without Gen5 overlays.")
    return(NULL)
  }
  data <- rename_by_alias(
    data,
    gen5_aliases,
    required = c("a", "b", "c", "d", "plate_number", "plate_id", "virus", "sample_id"),
    context = "Gen5 fit CSV"
  )
  data <- trim_character_columns(data)

  for (field in c("a", "b", "c", "d")) {
    data[[field]] <- parse_numeric_strict(data[[field]], field, seq_len(nrow(data)) + 1L)
  }
  if ("r2" %in% names(data)) {
    data$r2 <- parse_numeric_strict(data$r2, "r2", seq_len(nrow(data)) + 1L)
  } else {
    data$r2 <- NA_real_
  }
  if ("fit_f_prob" %in% names(data)) {
    data$fit_f_prob <- parse_numeric_strict(
      data$fit_f_prob,
      "fit_f_prob",
      seq_len(nrow(data)) + 1L
    )
  } else {
    data$fit_f_prob <- NA_real_
  }

  key_fields <- c("plate_number", "plate_id", "virus", "sample_id")
  for (field in key_fields) {
    value <- trimws(as.character(data[[field]]))
    missing <- is.na(data[[field]]) | value == ""
    if (any(missing)) {
      stop(
        "Gen5 fit CSV has missing '", field, "' value(s) at data row(s): ",
        paste(head(which(missing) + 1L, 10L), collapse = ", "),
        call. = FALSE
      )
    }
  }

  key <- do.call(paste, c(data[key_fields], sep = "\r"))
  if (anyDuplicated(key)) {
    stop(
      "Gen5 fit CSV has duplicate plate/sample/virus fit rows.",
      call. = FALSE
    )
  }

  data
}

# ---- Small NA-tolerant summary helpers -------------------------------------
# Defined together so every helper is available before any analysis function
# calls it.

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else stats::sd(x)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else mean(x)
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else min(x)
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else max(x)
}

safe_cv_pct <- function(x) {
  mean_value <- safe_mean(x)
  sd_value <- safe_sd(x)
  if (!is.finite(mean_value) || mean_value == 0 || !is.finite(sd_value)) {
    return(NA_real_)
  }
  100 * sd_value / mean_value
}

collapse_flags <- function(flags) {
  flags <- trimws(as.character(flags))
  flags <- unique(flags[!is.na(flags) & nzchar(flags)])
  if (length(flags) == 0L) "" else paste(flags, collapse = "; ")
}

split_flags <- function(value) {
  value <- as.character(value %||% "")
  if (length(value) == 0L || is.na(value[[1]]) || !nzchar(value[[1]])) {
    return(character())
  }
  unique(trimws(strsplit(value[[1]], ";", fixed = TRUE)[[1]]))
}

fold_difference <- function(x, y) {
  if (!is.finite(x) || !is.finite(y) || x <= 0 || y <= 0) {
    return(NA_real_)
  }
  max(x, y) / min(x, y)
}

log2_ratio <- function(x, y) {
  if (!is.finite(x) || !is.finite(y) || x <= 0 || y <= 0) {
    return(NA_real_)
  }
  log2(x / y)
}

model_warnings_requiring_review <- function(warnings) {
  warnings <- unique(as.character(warnings))
  warnings <- warnings[!is.na(warnings) & nzchar(warnings)]
  if (length(warnings) == 0L) return(character())

  # drc may emit transient warnings while searching for starting values even
  # when the final model is finite and usable. Preserve every warning in the
  # CSV, but elevate only messages that indicate convergence or structural
  # numerical problems. Final parameters and predictions are checked separately.
  review_pattern <- paste(
    c(
      "converg", "singular", "rank[- ]?deficient",
      "not positive definite", "iteration limit", "maximum.*iteration",
      "failed", "failure", "non[- ]?finite", "computationally singular",
      "ill[- ]?conditioned", "gradient", "step factor"
    ),
    collapse = "|"
  )
  warnings[grepl(review_pattern, tolower(warnings), perl = TRUE)]
}

validate_output_table <- function(
  data,
  output_name,
  required_fields = character()
) {
  duplicate_names <- unique(names(data)[duplicated(names(data))])
  if (length(duplicate_names) > 0L) {
    stop(
      output_name,
      " contains duplicate column name(s): ",
      paste(duplicate_names, collapse = ", "),
      call. = FALSE
    )
  }

  missing_required <- setdiff(required_fields, names(data))
  if (length(missing_required) > 0L) {
    stop(
      output_name,
      " is missing required canonical field(s): ",
      paste(missing_required, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

compute_normalization_qc <- function(data, config) {
  sample_data <- data[data$control_type == "sample", , drop = FALSE]
  curve_keys <- unique(sample_data[, c("plate_number", "plate_id", "virus")])
  plate_key <- first_numeric_token(curve_keys$plate_number)
  curve_keys <- curve_keys[order(
    is.na(plate_key),
    plate_key,
    tolower(curve_keys$plate_number),
    curve_keys$plate_id,
    curve_keys$virus
  ), , drop = FALSE]

  all_vc <- data[data$control_type == "virus_control", , drop = FALSE]
  rows <- vector("list", nrow(curve_keys))
  vc_lower <- config$qc$vc_count_range[[1]]
  vc_upper <- config$qc$vc_count_range[[2]]

  for (i in seq_len(nrow(curve_keys))) {
    plate_number <- curve_keys$plate_number[[i]]
    plate_id <- curve_keys$plate_id[[i]]
    virus <- curve_keys$virus[[i]]

    on_plate <- data$plate_number == plate_number & data$plate_id == plate_id
    plate_vc <- all_vc[
      all_vc$plate_number == plate_number & all_vc$plate_id == plate_id,
      ,
      drop = FALSE
    ]
    matched_vc <- plate_vc[
      !is.na(plate_vc$normalization_virus) &
        plate_vc$normalization_virus == virus,
      ,
      drop = FALSE
    ]

    plate_matched_vc_available <- nrow(matched_vc) > 0L
    reference_controls <- matched_vc
    reference_source <- if (plate_matched_vc_available) {
      if (any(grepl("row_inferred", matched_vc$control_virus_source, fixed = TRUE))) {
        "plate_matched_row_inferred"
      } else {
        "plate_matched_explicit"
      }
    } else {
      "unavailable"
    }
    reference_is_plate_matched <- plate_matched_vc_available
    denominator <- if (plate_matched_vc_available) {
      safe_mean(matched_vc$focus_count)
    } else {
      NA_real_
    }

    hard_flags <- character()
    review_flags <- character()
    if (!plate_matched_vc_available) {
      review_flags <- c(review_flags, "no_plate_matched_virus_control")

      if (isTRUE(config$normalization$allow_cross_plate_virus_controls)) {
        candidates <- all_vc[
          !is.na(all_vc$normalization_virus) &
            all_vc$normalization_virus == virus &
            all_vc$plate_id != plate_id,
          ,
          drop = FALSE
        ]
        selected_plates <- config$normalization$cross_plate_vc_plate_numbers
        if (!is.null(selected_plates)) {
          candidates <- candidates[
            candidates$plate_number %in% selected_plates,
            ,
            drop = FALSE
          ]
        }

        if (nrow(candidates) > 0L) {
          reference_controls <- candidates
          if (isTRUE(config$normalization$pool_equal_weight_by_plate)) {
            plate_means <- tapply(
              candidates$focus_count,
              candidates$plate_id,
              safe_mean
            )
            denominator <- safe_mean(as.numeric(plate_means))
          } else {
            denominator <- safe_mean(candidates$focus_count)
          }
          reference_source <- "cross_plate_same_virus_fallback"
          reference_is_plate_matched <- FALSE
          review_flags <- c(review_flags, "cross_plate_virus_control_used")
        } else {
          review_flags <- c(
            review_flags,
            "cross_plate_virus_control_enabled_but_no_eligible_controls_found"
          )
        }
      }
    }

    curve_sample_rows <- unique(sample_data$well_row[
      sample_data$plate_number == plate_number &
        sample_data$plate_id == plate_id &
        sample_data$virus == virus
    ])
    relevant_current_plate_vc <- plate_vc[
      !is.na(plate_vc$well_row) & plate_vc$well_row %in% curve_sample_rows,
      ,
      drop = FALSE
    ]
    mapping_flags <- unique(c(
      reference_controls$control_virus_mapping_flag,
      relevant_current_plate_vc$control_virus_mapping_flag
    ))
    mapping_flags <- mapping_flags[!is.na(mapping_flags) & nzchar(mapping_flags)]
    review_flags <- c(review_flags, mapping_flags)
    plate_matched_vc_required_metadata_inference <-
      plate_matched_vc_available && any(grepl(
        "inferred",
        reference_controls$control_virus_source,
        fixed = TRUE
      ))

    vc <- reference_controls$focus_count
    vc_mean <- denominator
    vc_sd <- safe_sd(vc)
    vc_cv <- if (is.finite(vc_mean) && vc_mean != 0 && is.finite(vc_sd)) {
      100 * vc_sd / vc_mean
    } else {
      NA_real_
    }
    vc_replicate_requirement_met <- length(vc) >= config$qc$min_vc_replicates
    vc_mean_in_range <- is.finite(vc_mean) && vc_mean >= vc_lower && vc_mean <= vc_upper
    vc_well_in_range <- is.finite(vc) & vc >= vc_lower & vc <= vc_upper
    all_vc_wells_in_range <- length(vc) > 0L && all(vc_well_in_range)
    n_vc_outside_range <- if (length(vc) == 0L) 0L else sum(!vc_well_in_range)
    vc_overall_in_range <- vc_mean_in_range && all_vc_wells_in_range

    current_cc <- data$focus_count[on_plate & data$control_type == "cell_control"]
    cell_control_max <- safe_max(current_cc)
    cell_control_background_high <- if (is.null(config$qc$cell_control_count_max)) {
      NA
    } else {
      is.finite(cell_control_max) &&
        cell_control_max > config$qc$cell_control_count_max
    }
    vc_cv_high <- if (is.null(config$qc$vc_cv_max_pct)) {
      NA
    } else {
      is.finite(vc_cv) && vc_cv > config$qc$vc_cv_max_pct
    }

    normalization_available <- is.finite(vc_mean) && vc_mean > 0
    if (!normalization_available) {
      hard_flags <- c(hard_flags, "virus_control_reference_unavailable")
    }
    if (normalization_available && !vc_overall_in_range) {
      hard_flags <- c(hard_flags, "virus_control_reference_outside_working_range")
    }
    if (!reference_is_plate_matched && normalization_available &&
        !isTRUE(config$normalization$allow_cross_plate_vc_for_reporting)) {
      hard_flags <- c(
        hard_flags,
        "cross_plate_virus_control_not_allowed_for_reporting"
      )
    }
    if (!vc_replicate_requirement_met) {
      if (isTRUE(config$qc$require_min_vc_replicates_for_reporting)) {
        hard_flags <- c(hard_flags, "too_few_virus_control_replicates")
      } else {
        review_flags <- c(review_flags, "too_few_virus_control_replicates")
      }
    }
    if (isTRUE(vc_cv_high)) {
      review_flags <- c(review_flags, "high_virus_control_cv")
    }
    if (isTRUE(cell_control_background_high)) {
      review_flags <- c(review_flags, "high_cell_control_background")
    }

    hard_flags <- unique(hard_flags)
    review_flags <- unique(review_flags)
    normalization_reportable <- length(hard_flags) == 0L
    normalization_status <- if (!normalization_reportable) {
      "FAIL"
    } else if (length(review_flags) > 0L) {
      "REVIEW"
    } else {
      "PASS"
    }

    reference_plate_numbers <- unique(as.character(reference_controls$plate_number))
    if (length(reference_plate_numbers) > 0L) {
      reference_plate_numbers <- reference_plate_numbers[
        natural_order(reference_plate_numbers)
      ]
    }
    reference_plate_ids <- unique(as.character(reference_controls$plate_id))

    rows[[i]] <- data.frame(
      plate_number = plate_number,
      plate_id = plate_id,
      virus = virus,
      normalization_source = reference_source,
      normalization_available = normalization_available,
      normalization_reportable = normalization_reportable,
      normalization_status = normalization_status,
      normalization_reference_is_plate_matched = reference_is_plate_matched,
      plate_matched_virus_control_available = plate_matched_vc_available,
      plate_matched_vc_required_metadata_inference =
        plate_matched_vc_required_metadata_inference,
      cross_plate_virus_controls_allowed =
        config$normalization$allow_cross_plate_virus_controls,
      cross_plate_vc_reporting_allowed =
        config$normalization$allow_cross_plate_vc_for_reporting,
      cross_plate_vc_plate_numbers_requested = paste(
        config$normalization$cross_plate_vc_plate_numbers %||% character(),
        collapse = " | "
      ),
      normalization_reference_plate_numbers = paste(
        reference_plate_numbers,
        collapse = " | "
      ),
      normalization_reference_plate_ids = paste(
        reference_plate_ids,
        collapse = " | "
      ),
      n_normalization_reference_plates = length(reference_plate_ids),
      n_virus_control = length(vc),
      minimum_required_virus_controls = config$qc$min_vc_replicates,
      virus_control_replicate_requirement_met = vc_replicate_requirement_met,
      virus_control_mean = vc_mean,
      virus_control_sd = vc_sd,
      virus_control_cv_pct = vc_cv,
      virus_control_min = safe_min(vc),
      virus_control_max = safe_max(vc),
      virus_control_mean_in_working_range = vc_mean_in_range,
      all_virus_control_wells_in_working_range = all_vc_wells_in_range,
      n_virus_control_outside_working_range = n_vc_outside_range,
      virus_control_in_working_range = vc_overall_in_range,
      virus_control_cv_threshold_pct = config$qc$vc_cv_max_pct %||% NA_real_,
      virus_control_cv_high = vc_cv_high,
      n_cell_control = length(current_cc),
      cell_control_mean = safe_mean(current_cc),
      cell_control_sd = safe_sd(current_cc),
      cell_control_max = cell_control_max,
      cell_control_count_threshold = config$qc$cell_control_count_max %||% NA_real_,
      cell_control_background_high = cell_control_background_high,
      control_virus_sources = paste(
        unique(reference_controls$control_virus_source),
        collapse = " | "
      ),
      normalization_hard_qc_flags = collapse_flags(hard_flags),
      normalization_review_qc_flags = collapse_flags(review_flags),
      stringsAsFactors = FALSE
    )
  }

  dplyr::bind_rows(rows)
}

normalize_to_controls <- function(data, normalization_qc) {
  qc_for_join <- normalization_qc
  names(qc_for_join)[names(qc_for_join) == "virus"] <- "normalization_virus"

  normalized <- dplyr::left_join(
    data,
    qc_for_join,
    by = c("plate_number", "plate_id", "normalization_virus")
  )

  valid_denominator <- is.finite(normalized$virus_control_mean) &
    normalized$virus_control_mean > 0
  eligible_rows <- normalized$control_type != "cell_control"
  normalized$percent_neutralization <- NA_real_
  normalized$percent_neutralization[valid_denominator & eligible_rows] <-
    (1 - normalized$focus_count[valid_denominator & eligible_rows] /
      normalized$virus_control_mean[valid_denominator & eligible_rows]) * 100

  normalized
}

summarize_dilutions <- function(normalized_data) {
  sample_data <- normalized_data[normalized_data$control_type == "sample", , drop = FALSE]
  if (nrow(sample_data) == 0L) {
    stop("No sample wells were found after control classification.", call. = FALSE)
  }

  summary_data <- sample_data |>
    dplyr::group_by(
      .data$plate_number,
      .data$plate_id,
      .data$virus,
      .data$sample_id,
      .data$dilution
    ) |>
    dplyr::summarise(
      n_replicates = dplyr::n(),
      mean_focus_count = safe_mean(.data$focus_count),
      sd_focus_count = safe_sd(.data$focus_count),
      mean_percent_neutralization = safe_mean(.data$percent_neutralization),
      sd_percent_neutralization = safe_sd(.data$percent_neutralization),
      se_percent_neutralization = if (dplyr::n() > 1L) {
        safe_sd(.data$percent_neutralization) / sqrt(dplyr::n())
      } else {
        NA_real_
      },
      min_percent_neutralization = safe_min(.data$percent_neutralization),
      max_percent_neutralization = safe_max(.data$percent_neutralization),
      .groups = "drop"
    )

  plate_key <- first_numeric_token(summary_data$plate_number)
  summary_data[order(
    is.na(plate_key),
    plate_key,
    tolower(summary_data$plate_number),
    summary_data$plate_id,
    summary_data$virus,
    summary_data$sample_id,
    summary_data$dilution
  ), , drop = FALSE]
}

count_threshold_crossings <- function(x, y, threshold, tolerance = 1e-8) {
  if (length(x) != length(y) || length(x) < 2L) {
    return(list(
      exact_indices = integer(),
      sign_crossing_indices = integer(),
      total = 0L
    ))
  }

  ordering <- order(x)
  x <- x[ordering]
  y <- y[ordering]
  z <- y - threshold

  exact <- which(abs(z) <= tolerance)
  sign_crossings <- which(z[-length(z)] * z[-1L] < 0)

  list(
    exact_indices = exact,
    sign_crossing_indices = sign_crossings,
    total = length(exact) + length(sign_crossings)
  )
}

# Secondary NT50 estimate from observed values.
#
# This method first evaluates the COMPLETE observed dilution series to determine
# whether there is exactly one threshold crossing. If and only if the full series
# contains one unique crossing, the numerical root is obtained by linear
# interpolation inside the two adjacent observed dilution means that contain that
# crossing. Using a piecewise-linear interpolant over the full series would give
# the same numerical root for a unique crossing; restricting the final root call
# to the local bracket avoids an arbitrary root when irregular data cross the
# threshold more than once. It does not use LL4 or Gen5 predictions and does not
# extrapolate beyond the tested dilution range.
estimate_uniroot_nt50 <- function(
  summary_data,
  threshold,
  interpolation_scale = "linear"
) {
  method_name <- if (interpolation_scale == "log2") {
    "observed_mean_linear_interpolation_on_log2_dilution_scale"
  } else {
    "observed_mean_linear_interpolation_on_dilution_scale"
  }

  empty_result <- function(
    status,
    warnings = character(),
    error = "",
    n_dilutions = 0L,
    n_crossings = 0L
  ) {
    list(
      estimate = NA_real_,
      estimate_in_range = FALSE,
      status = status,
      method = method_name,
      interpolation_scale = interpolation_scale,
      n_dilutions = as.integer(n_dilutions),
      n_crossings = as.integer(n_crossings),
      warnings = unique(warnings),
      error = error %||% ""
    )
  }

  finite <- is.finite(summary_data$dilution) & summary_data$dilution > 0 &
    is.finite(summary_data$mean_percent_neutralization)
  x <- as.numeric(summary_data$dilution[finite])
  y <- as.numeric(summary_data$mean_percent_neutralization[finite])
  if (length(x) < 2L) {
    return(empty_result(
      "insufficient_observed_dilution_means",
      n_dilutions = length(x)
    ))
  }

  # Defensive aggregation in case a nonstandard input reaches this function with
  # more than one summary row for the same dilution.
  unique_x <- sort(unique(x))
  y <- vapply(
    unique_x,
    function(value) safe_mean(y[x == value]),
    numeric(1)
  )
  x <- unique_x
  finite <- is.finite(x) & is.finite(y)
  x <- x[finite]
  y <- y[finite]
  if (length(x) < 2L) {
    return(empty_result(
      "insufficient_finite_observed_dilution_means",
      n_dilutions = length(x)
    ))
  }

  crossings <- count_threshold_crossings(x, y, threshold)
  if (length(crossings$exact_indices) == 1L &&
      length(crossings$sign_crossing_indices) == 0L) {
    estimate <- x[crossings$exact_indices[[1]]]
    return(list(
      estimate = estimate,
      estimate_in_range = TRUE,
      status = "exact_observed_mean_at_tested_dilution",
      method = method_name,
      interpolation_scale = interpolation_scale,
      n_dilutions = as.integer(length(x)),
      n_crossings = 1L,
      warnings = character(),
      error = ""
    ))
  }

  if (crossings$total == 0L) {
    return(empty_result(
      "observed_means_do_not_bracket_threshold",
      n_dilutions = length(x),
      n_crossings = 0L
    ))
  }
  if (crossings$total > 1L) {
    return(empty_result(
      "multiple_observed_mean_threshold_crossings",
      n_dilutions = length(x),
      n_crossings = crossings$total
    ))
  }

  index <- crossings$sign_crossing_indices[[1]]
  x_pair <- x[c(index, index + 1L)]
  y_pair <- y[c(index, index + 1L)]
  working_x <- if (interpolation_scale == "log2") log2(x_pair) else x_pair
  interpolator <- stats::approxfun(
    x = working_x,
    y = y_pair,
    method = "linear",
    rule = 1L,
    ties = "ordered"
  )

  root_attempt <- capture_model_fit(
    stats::uniroot(
      function(value) interpolator(value) - threshold,
      interval = range(working_x),
      extendInt = "no"
    )$root
  )
  if (is.null(root_attempt$value) || !is.finite(root_attempt$value)) {
    return(empty_result(
      "uniroot_failed_inside_observed_bracket",
      warnings = root_attempt$warnings,
      error = root_attempt$error,
      n_dilutions = length(x),
      n_crossings = 1L
    ))
  }

  estimate <- if (interpolation_scale == "log2") {
    2^as.numeric(root_attempt$value)
  } else {
    as.numeric(root_attempt$value)
  }
  in_range <- estimate >= min(x) && estimate <= max(x)
  list(
    estimate = estimate,
    estimate_in_range = in_range,
    status = if (in_range) "in_range" else "estimate_outside_tested_range",
    method = method_name,
    interpolation_scale = interpolation_scale,
    n_dilutions = as.integer(length(x)),
    n_crossings = 1L,
    warnings = unique(root_attempt$warnings),
    error = root_attempt$error %||% ""
  )
}

capture_model_fit <- function(expression) {
  warnings <- character()
  error_message <- NULL

  value <- withCallingHandlers(
    tryCatch(
      expression,
      error = function(e) {
        error_message <<- conditionMessage(e)
        NULL
      }
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  list(
    value = value,
    warnings = unique(warnings),
    error = error_message
  )
}

extract_named_parameter <- function(coefficients, prefix) {
  coefficient_names <- names(coefficients)
  if (is.null(coefficient_names)) {
    return(NA_real_)
  }

  hit <- grep(paste0("^", prefix, "($|:)"), coefficient_names)
  if (length(hit) == 0L) {
    return(NA_real_)
  }
  as.numeric(coefficients[[hit[[1]]]])
}

extract_ll4_parameters <- function(model, forced_hill_value = NA_real_) {
  coefficients <- stats::coef(model)

  hill <- extract_named_parameter(coefficients, "hill")
  lower <- extract_named_parameter(coefficients, "lower")
  upper <- extract_named_parameter(coefficients, "upper")
  midpoint <- extract_named_parameter(coefficients, "midpoint")

  # Defensive fallback for older drc naming conventions.
  if (!is.finite(lower)) lower <- extract_named_parameter(coefficients, "c")
  if (!is.finite(upper)) upper <- extract_named_parameter(coefficients, "d")
  if (!is.finite(midpoint)) midpoint <- extract_named_parameter(coefficients, "e")
  if (!is.finite(hill)) hill <- extract_named_parameter(coefficients, "b")

  if (!is.finite(hill) && is.finite(forced_hill_value)) {
    hill <- forced_hill_value
  }

  # Final positional fallback. For a fixed-Hill LL4, the estimated coefficients
  # are lower, upper, midpoint; otherwise they are hill, lower, upper, midpoint.
  if (any(!is.finite(c(lower, upper, midpoint)))) {
    values <- as.numeric(coefficients)
    if (length(values) == 4L) {
      hill <- values[[1]]
      lower <- values[[2]]
      upper <- values[[3]]
      midpoint <- values[[4]]
    } else if (length(values) == 3L && is.finite(forced_hill_value)) {
      hill <- forced_hill_value
      lower <- values[[1]]
      upper <- values[[2]]
      midpoint <- values[[3]]
    }
  }

  list(
    hill = hill,
    lower = lower,
    upper = upper,
    midpoint = midpoint
  )
}

fit_ll4_curve <- function(curve_data, config) {
  n_unique <- length(unique(curve_data$dilution[is.finite(curve_data$dilution)]))
  n_observations <- sum(is.finite(curve_data$dilution) &
    is.finite(curve_data$percent_neutralization))

  if (n_unique < config$analysis$min_unique_dilutions ||
      n_observations < config$analysis$min_total_observations) {
    return(list(
      model = NULL,
      forced_hill = FALSE,
      initial_hill = NA_real_,
      parameters = list(
        hill = NA_real_, lower = NA_real_, upper = NA_real_, midpoint = NA_real_
      ),
      warnings = character(),
      error = "insufficient_data_for_ll4"
    ))
  }

  unconstrained <- capture_model_fit(
    drc::drm(
      percent_neutralization ~ dilution,
      data = curve_data,
      fct = drc::LL.4(names = c("hill", "lower", "upper", "midpoint"))
    )
  )

  if (is.null(unconstrained$value)) {
    return(list(
      model = NULL,
      forced_hill = FALSE,
      initial_hill = NA_real_,
      parameters = list(
        hill = NA_real_, lower = NA_real_, upper = NA_real_, midpoint = NA_real_
      ),
      warnings = unconstrained$warnings,
      error = unconstrained$error %||% "unconstrained_ll4_failed"
    ))
  }

  unconstrained_parameters <- extract_ll4_parameters(unconstrained$value)
  initial_hill <- unconstrained_parameters$hill
  hill_outside_range <- !is.finite(initial_hill) ||
    initial_hill < config$qc$hill_range[[1]] ||
    initial_hill > config$qc$hill_range[[2]]

  if (!hill_outside_range) {
    return(list(
      model = unconstrained$value,
      forced_hill = FALSE,
      initial_hill = initial_hill,
      parameters = unconstrained_parameters,
      warnings = unconstrained$warnings,
      error = NULL
    ))
  }

  fixed_hill_value <- config$qc$fixed_hill
  fixed_fit <- capture_model_fit(
    drc::drm(
      percent_neutralization ~ dilution,
      data = curve_data,
      fct = drc::LL.4(
        fixed = c(fixed_hill_value, NA, NA, NA),
        names = c("hill", "lower", "upper", "midpoint")
      )
    )
  )

  if (is.null(fixed_fit$value)) {
    return(list(
      model = NULL,
      forced_hill = TRUE,
      initial_hill = initial_hill,
      parameters = list(
        hill = NA_real_, lower = NA_real_, upper = NA_real_, midpoint = NA_real_
      ),
      warnings = unique(c(unconstrained$warnings, fixed_fit$warnings)),
      error = paste0(
        "hill_outside_range_and_fixed_hill_refit_failed: ",
        fixed_fit$error %||% "unknown error"
      )
    ))
  }

  list(
    model = fixed_fit$value,
    forced_hill = TRUE,
    initial_hill = initial_hill,
    parameters = extract_ll4_parameters(fixed_fit$value, fixed_hill_value),
    warnings = unique(c(unconstrained$warnings, fixed_fit$warnings)),
    error = NULL
  )
}

extract_ed_output <- function(ed_object) {
  matrix_value <- as.matrix(ed_object)
  if (nrow(matrix_value) < 1L || ncol(matrix_value) < 1L) {
    return(list(
      estimate = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_
    ))
  }

  columns <- tolower(colnames(matrix_value) %||% rep("", ncol(matrix_value)))
  named_value <- function(pattern) {
    hit <- grep(pattern, columns)
    if (length(hit) == 0L) NA_real_ else as.numeric(matrix_value[1L, hit[[1]]])
  }

  estimate <- named_value("estimate")
  se <- named_value("std|se")
  lower <- named_value("lower")
  upper <- named_value("upper")

  # Positional fallbacks cover older drc outputs. A three-column result is
  # interpreted as Estimate/Lower/Upper; a four-column result is interpreted as
  # Estimate/Std. Error/Lower/Upper. This avoids mislabeling a lower limit as SE.
  if (!is.finite(estimate)) estimate <- as.numeric(matrix_value[1L, 1L])
  if (ncol(matrix_value) >= 4L) {
    if (!is.finite(se)) se <- as.numeric(matrix_value[1L, 2L])
    if (!is.finite(lower)) lower <- as.numeric(matrix_value[1L, 3L])
    if (!is.finite(upper)) upper <- as.numeric(matrix_value[1L, 4L])
  } else if (ncol(matrix_value) == 3L) {
    if (!is.finite(lower)) lower <- as.numeric(matrix_value[1L, 2L])
    if (!is.finite(upper)) upper <- as.numeric(matrix_value[1L, 3L])
  }

  list(estimate = estimate, se = se, lower = lower, upper = upper)
}

drc_ed_supported_intervals <- function() {
  fallback <- c("none", "delta", "fls", "tfls", "inv")
  method <- tryCatch(
    getFromNamespace("ED.drc", "drc"),
    error = function(e) NULL
  )
  if (is.null(method)) {
    return(fallback)
  }

  supported <- tryCatch(
    eval(formals(method)$interval, envir = environment(method)),
    error = function(e) NULL
  )
  if (is.null(supported)) fallback else unique(as.character(supported))
}

estimate_ll4_nt50 <- function(model, parameters, summary_data, config) {
  threshold <- config$analysis$threshold_pct
  min_dilution <- safe_min(summary_data$dilution)
  max_dilution <- safe_max(summary_data$dilution)

  empty_result <- function(status, ci_status = "not_attempted", warnings = character(), error = "") {
    list(
      estimate = NA_real_,
      se = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      estimate_in_range = FALSE,
      ci_available = FALSE,
      ci_in_range = FALSE,
      status = status,
      ci_status = ci_status,
      warnings = unique(warnings),
      error = error %||% ""
    )
  }

  if (is.null(model)) {
    return(empty_result("ll4_fit_unavailable"))
  }
  if (!is.finite(min_dilution) || !is.finite(max_dilution) ||
      min_dilution <= 0 || max_dilution <= min_dilution) {
    return(empty_result("invalid_tested_dilution_range"))
  }

  fitted_bounds <- sort(c(parameters$lower, parameters$upper))
  if (any(!is.finite(fitted_bounds)) ||
      threshold <= fitted_bounds[[1]] || threshold >= fitted_bounds[[2]]) {
    return(empty_result("threshold_outside_fitted_asymptotes"))
  }

  requested_interval <- config$analysis$ll4_ci_method
  supported_intervals <- drc_ed_supported_intervals()
  interval_for_call <- requested_interval
  ci_status <- if (requested_interval == "none") "not_requested" else requested_interval
  if (requested_interval == "inv" && !"inv" %in% supported_intervals) {
    interval_for_call <- "none"
    ci_status <- "inv_not_supported_by_installed_drc_estimate_only"
  }

  ed_attempt <- capture_model_fit(
    drc::ED(
      model,
      respLev = threshold,
      type = "absolute",
      interval = interval_for_call,
      level = config$analysis$confidence_level,
      display = FALSE
    )
  )

  all_warnings <- ed_attempt$warnings
  all_errors <- ed_attempt$error %||% ""

  if (is.null(ed_attempt$value) && interval_for_call != "none") {
    estimate_only_attempt <- capture_model_fit(
      drc::ED(
        model,
        respLev = threshold,
        type = "absolute",
        interval = "none",
        display = FALSE
      )
    )
    all_warnings <- unique(c(all_warnings, estimate_only_attempt$warnings))
    all_errors <- paste(
      c(all_errors, estimate_only_attempt$error %||% "")[nzchar(c(all_errors, estimate_only_attempt$error %||% ""))],
      collapse = " | "
    )
    ed_attempt <- estimate_only_attempt
    ci_status <- "requested_ci_failed_estimate_only"
  }

  if (is.null(ed_attempt$value)) {
    return(empty_result(
      "ll4_inverse_prediction_failed",
      ci_status,
      all_warnings,
      all_errors
    ))
  }

  extracted <- extract_ed_output(ed_attempt$value)
  ci_limits <- sort(c(extracted$lower, extracted$upper))
  lower <- if (length(ci_limits) == 2L) ci_limits[[1]] else NA_real_
  upper <- if (length(ci_limits) == 2L) ci_limits[[2]] else NA_real_

  estimate_in_range <- is.finite(extracted$estimate) &&
    extracted$estimate >= min_dilution && extracted$estimate <= max_dilution
  ci_available <- is.finite(lower) && is.finite(upper) && lower > 0 && upper > 0
  ci_in_range <- ci_available && lower >= min_dilution && upper <= max_dilution

  list(
    estimate = extracted$estimate,
    se = extracted$se,
    lower = lower,
    upper = upper,
    estimate_in_range = estimate_in_range,
    ci_available = ci_available,
    ci_in_range = ci_in_range,
    status = if (estimate_in_range) "in_range" else "estimate_outside_tested_range",
    ci_status = if (!ci_available) {
      ci_status
    } else if (ci_in_range) {
      paste0(ci_status, "_in_range")
    } else {
      paste0(ci_status, "_extends_outside_tested_range")
    },
    warnings = unique(all_warnings),
    error = all_errors
  )
}

enabled_plot_branch_names <- function(config) {
  branches <- config$plot$branches %||% list()
  names(branches)[vapply(
    branches,
    function(branch) isTRUE(branch$enabled),
    logical(1)
  )]
}

any_enabled_plot_branch_requests <- function(config, field) {
  enabled_names <- enabled_plot_branch_names(config)
  if (length(enabled_names) == 0L) {
    return(FALSE)
  }
  any(vapply(
    config$plot$branches[enabled_names],
    function(branch) isTRUE(branch[[field]]),
    logical(1)
  ))
}

predict_ll4_curve <- function(model, min_dilution, max_dilution, config) {
  if (is.null(model) || !is.finite(min_dilution) || !is.finite(max_dilution) ||
      min_dilution <= 0 || max_dilution <= min_dilution) {
    return(tibble::tibble())
  }

  grid <- 2^seq(
    log2(min_dilution),
    log2(max_dilution),
    length.out = config$analysis$prediction_grid_points
  )
  new_data <- data.frame(dilution = grid)
  request_band <- any_enabled_plot_branch_requests(
    config,
    "show_model_confidence_band"
  )

  if (request_band) {
    prediction_attempt <- capture_model_fit(
      stats::predict(
        model,
        newdata = new_data,
        interval = "confidence",
        level = config$analysis$confidence_level
      )
    )

    if (!is.null(prediction_attempt$value)) {
      prediction_matrix <- as.matrix(prediction_attempt$value)
      columns <- tolower(colnames(prediction_matrix) %||% rep("", ncol(prediction_matrix)))
      fitted_index <- grep("prediction|fit", columns)
      lower_index <- grep("lower", columns)
      upper_index <- grep("upper", columns)

      fitted <- if (length(fitted_index) > 0L) {
        prediction_matrix[, fitted_index[[1]]]
      } else {
        prediction_matrix[, 1L]
      }
      lower <- if (length(lower_index) > 0L) {
        prediction_matrix[, lower_index[[1]]]
      } else if (ncol(prediction_matrix) >= 2L) {
        prediction_matrix[, 2L]
      } else {
        rep(NA_real_, length(grid))
      }
      upper <- if (length(upper_index) > 0L) {
        prediction_matrix[, upper_index[[1]]]
      } else if (ncol(prediction_matrix) >= 3L) {
        prediction_matrix[, 3L]
      } else {
        rep(NA_real_, length(grid))
      }

      return(tibble::tibble(
        dilution = grid,
        fitted = as.numeric(fitted),
        lower = as.numeric(lower),
        upper = as.numeric(upper),
        model_band_available = TRUE,
        prediction_warnings = paste(prediction_attempt$warnings, collapse = " | ")
      ))
    }
  }

  prediction_attempt <- capture_model_fit(
    stats::predict(model, newdata = new_data)
  )
  if (is.null(prediction_attempt$value)) {
    return(tibble::tibble())
  }

  tibble::tibble(
    dilution = grid,
    fitted = as.numeric(prediction_attempt$value),
    lower = NA_real_,
    upper = NA_real_,
    model_band_available = FALSE,
    prediction_warnings = paste(prediction_attempt$warnings, collapse = " | ")
  )
}

extract_lack_of_fit_p <- function(model) {
  if (is.null(model)) {
    return(NA_real_)
  }

  fit_test <- tryCatch(
    drc::modelFit(model),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  if (is.null(fit_test)) {
    return(NA_real_)
  }

  fit_matrix <- as.matrix(fit_test)
  p_columns <- grep("pr\\(|p.value|p-value", tolower(colnames(fit_matrix) %||% ""))
  if (length(p_columns) == 0L) {
    return(NA_real_)
  }

  candidates <- suppressWarnings(as.numeric(fit_matrix[, p_columns[[1]]]))
  candidates <- candidates[is.finite(candidates)]
  if (length(candidates) == 0L) NA_real_ else tail(candidates, 1L)
}

compute_fit_metrics <- function(model, curve_data) {
  if (is.null(model)) {
    return(list(rmse = NA_real_, pseudo_r2 = NA_real_, residual_df = NA_real_))
  }

  fitted_values <- tryCatch(
    as.numeric(stats::predict(model, newdata = curve_data)),
    error = function(e) rep(NA_real_, nrow(curve_data))
  )
  observed <- curve_data$percent_neutralization
  finite <- is.finite(fitted_values) & is.finite(observed)

  if (sum(finite) < 2L) {
    return(list(rmse = NA_real_, pseudo_r2 = NA_real_, residual_df = NA_real_))
  }

  residuals <- observed[finite] - fitted_values[finite]
  sse <- sum(residuals^2)
  sst <- sum((observed[finite] - mean(observed[finite]))^2)

  list(
    rmse = sqrt(mean(residuals^2)),
    pseudo_r2 = if (sst > 0) 1 - sse / sst else NA_real_,
    residual_df = tryCatch(stats::df.residual(model), error = function(e) NA_real_)
  )
}

gen5_count_equation <- function(x, a, b, c, d) {
  (a - d) / (1 + (x / c)^b) + d
}

extract_gen5_coefficients <- function(gen5_row) {
  if (is.null(gen5_row) || nrow(gen5_row) != 1L) {
    return(list(valid = FALSE, a = NA_real_, b = NA_real_, c = NA_real_, d = NA_real_))
  }

  values <- vapply(
    c("a", "b", "c", "d"),
    function(field) suppressWarnings(as.numeric(gen5_row[[field]][[1]])),
    numeric(1)
  )

  list(
    valid = all(is.finite(values)) && values[["c"]] > 0,
    a = values[["a"]],
    b = values[["b"]],
    c = values[["c"]],
    d = values[["d"]]
  )
}

predict_gen5_curve <- function(gen5_row, vc_mean, min_dilution, max_dilution, config) {
  coefficients <- extract_gen5_coefficients(gen5_row)
  if (!coefficients$valid || !is.finite(vc_mean) || vc_mean <= 0 ||
      !is.finite(min_dilution) || !is.finite(max_dilution) ||
      min_dilution <= 0 || max_dilution <= min_dilution) {
    return(tibble::tibble())
  }

  grid <- 2^seq(
    log2(min_dilution),
    log2(max_dilution),
    length.out = config$analysis$prediction_grid_points
  )
  fitted_counts <- gen5_count_equation(
    grid,
    coefficients$a,
    coefficients$b,
    coefficients$c,
    coefficients$d
  )

  tibble::tibble(
    dilution = grid,
    fitted = (1 - fitted_counts / vc_mean) * 100
  )
}

# Reference continuous inversion (Gen5/imaging-system fit). Unlike the LL4
# uniroot cross-check above, this inverts the continuous fitted equation on the
# log2 dilution scale. It exists to compare against the imaging-system fit, and
# also serves as the template if a true continuous LL4 inversion is ever desired.
estimate_gen5_nt50 <- function(gen5_row, vc_mean, min_dilution, max_dilution, threshold) {
  coefficients <- extract_gen5_coefficients(gen5_row)
  if (!coefficients$valid || !is.finite(vc_mean) || vc_mean <= 0) {
    return(list(estimate = NA_real_, status = "gen5_fit_unavailable_or_invalid"))
  }
  if (!is.finite(min_dilution) || !is.finite(max_dilution) ||
      min_dilution <= 0 || max_dilution <= min_dilution) {
    return(list(estimate = NA_real_, status = "invalid_tested_dilution_range"))
  }

  function_on_log_scale <- function(log_dilution) {
    dilution <- 2^log_dilution
    fitted_count <- gen5_count_equation(
      dilution,
      coefficients$a,
      coefficients$b,
      coefficients$c,
      coefficients$d
    )
    (1 - fitted_count / vc_mean) * 100 - threshold
  }

  lower_value <- function_on_log_scale(log2(min_dilution))
  upper_value <- function_on_log_scale(log2(max_dilution))
  tolerance <- 1e-8
  if (!is.finite(lower_value) || !is.finite(upper_value)) {
    return(list(estimate = NA_real_, status = "gen5_prediction_nonfinite"))
  }
  if (abs(lower_value) <= tolerance && abs(upper_value) <= tolerance) {
    return(list(estimate = NA_real_, status = "gen5_threshold_crossing_not_unique"))
  }
  if (abs(lower_value) <= tolerance) {
    return(list(estimate = min_dilution, status = "in_range_at_lower_endpoint"))
  }
  if (abs(upper_value) <= tolerance) {
    return(list(estimate = max_dilution, status = "in_range_at_upper_endpoint"))
  }
  if (lower_value * upper_value > 0) {
    return(list(estimate = NA_real_, status = "gen5_threshold_not_bracketed"))
  }

  root <- tryCatch(
    stats::uniroot(
      function_on_log_scale,
      interval = c(log2(min_dilution), log2(max_dilution)),
      extendInt = "no"
    )$root,
    error = function(e) NA_real_
  )

  list(
    estimate = if (is.finite(root)) 2^root else NA_real_,
    status = if (is.finite(root)) "in_range" else "gen5_root_failure"
  )
}

curve_descriptors <- function(summary_data, threshold) {
  finite <- is.finite(summary_data$dilution) & summary_data$dilution > 0 &
    is.finite(summary_data$mean_percent_neutralization)
  summary_data <- summary_data[finite, , drop = FALSE]
  summary_data <- summary_data[order(summary_data$dilution), , drop = FALSE]

  if (nrow(summary_data) == 0L) {
    return(list(
      threshold_bracketed = FALSE,
      min_observed_mean = NA_real_,
      max_observed_mean = NA_real_,
      spearman_rho = NA_real_,
      n_upward_mean_steps = NA_integer_,
      max_replicate_sd_pp = NA_real_,
      min_dilution = NA_real_,
      max_dilution = NA_real_,
      n_unique_dilutions = 0L,
      min_replicates_per_dilution = NA_integer_,
      max_replicates_per_dilution = NA_integer_
    ))
  }

  x <- summary_data$dilution
  y <- summary_data$mean_percent_neutralization
  min_y <- safe_min(y)
  max_y <- safe_max(y)
  bracketed <- is.finite(min_y) && is.finite(max_y) &&
    min_y <= threshold && max_y >= threshold

  spearman <- if (length(unique(x)) >= 3L && length(unique(y)) >= 2L) {
    suppressWarnings(stats::cor(log2(x), y, method = "spearman"))
  } else {
    NA_real_
  }

  list(
    threshold_bracketed = bracketed,
    min_observed_mean = min_y,
    max_observed_mean = max_y,
    spearman_rho = spearman,
    n_upward_mean_steps = if (length(y) > 1L) sum(diff(y) > 0, na.rm = TRUE) else 0L,
    max_replicate_sd_pp = safe_max(summary_data$sd_percent_neutralization),
    min_dilution = safe_min(x),
    max_dilution = safe_max(x),
    n_unique_dilutions = length(unique(x)),
    min_replicates_per_dilution = as.integer(safe_min(summary_data$n_replicates)),
    max_replicates_per_dilution = as.integer(safe_max(summary_data$n_replicates))
  )
}

format_number <- function(x, digits = 3L) {
  if (!is.finite(x)) {
    return("NA")
  }
  formatC(x, digits = digits, format = "fg", flag = "#", big.mark = ",")
}

pretty_qc_flag <- function(flag) {
  text <- gsub("_", " ", as.character(flag), fixed = TRUE)
  paste0(toupper(substr(text, 1L, 1L)), substring(text, 2L))
}

format_qc_flags_for_plot <- function(
  flags,
  width = 66L,
  max_flags = Inf,
  bullet = TRUE,
  blank_lines_between_flags = 1L
) {
  values <- split_flags(flags)
  if (length(values) == 0L) return("")

  omitted <- 0L
  if (is.finite(max_flags) && length(values) > max_flags) {
    omitted <- length(values) - as.integer(max_flags)
    values <- head(values, as.integer(max_flags))
  }

  blocks <- lapply(values, function(value) {
    wrapped <- strwrap(pretty_qc_flag(value), width = width)
    if (length(wrapped) == 0L) return("")
    prefix <- if (bullet) "- " else ""
    continuation <- if (bullet) "  " else ""
    paste(
      c(paste0(prefix, wrapped[[1]]), paste0(continuation, wrapped[-1L])),
      collapse = "\n"
    )
  })
  blocks <- unlist(blocks, use.names = FALSE)
  blocks <- blocks[nzchar(blocks)]

  if (omitted > 0L) {
    blocks <- c(
      blocks,
      paste0(if (bullet) "- " else "", omitted, " additional flag(s) in CSV")
    )
  }

  # One configured blank line means two newline characters between flag blocks.
  separator <- strrep("\n", as.integer(blank_lines_between_flags) + 1L)
  paste(blocks, collapse = separator)
}

safe_filename <- function(...) {
  value <- paste(..., sep = "__")
  converted <- suppressWarnings(iconv(value, from = "", to = "ASCII//TRANSLIT", sub = "_"))
  if (length(converted) == 0L || is.na(converted)) {
    converted <- value
  }
  value <- converted
  value <- gsub("[^A-Za-z0-9._-]+", "_", value)
  value <- gsub("_+", "_", value)
  value <- gsub("^_+|_+$", "", value)
  if (nchar(value) == 0L) "plot" else substr(value, 1L, 180L)
}

ensure_plot_extension <- function(path, format) {
  path <- as.character(path)
  format <- tolower(trimws(as.character(format)))
  if (!format %in% c("png", "pdf")) {
    stop("Unsupported plot format: ", format, call. = FALSE)
  }

  expected_suffix <- paste0(".", format)
  if (!endsWith(tolower(path), expected_suffix)) {
    path <- paste0(path, expected_suffix)
  }
  path
}

build_plot_output_path <- function(directory, index, stem, format) {
  filename <- paste0(sprintf("%03d__", as.integer(index)), stem)
  ensure_plot_extension(file.path(directory, filename), format)
}

extensionless_plot_path <- function(path, format) {
  path <- ensure_plot_extension(path, format)
  suffix <- paste0(".", tolower(trimws(as.character(format))))
  substr(path, 1L, nchar(path) - nchar(suffix))
}

recover_extensionless_plot <- function(expected_path, format) {
  expected_path <- ensure_plot_extension(expected_path, format)
  if (file.exists(expected_path)) {
    return(TRUE)
  }

  extensionless_path <- extensionless_plot_path(expected_path, format)
  if (file.exists(extensionless_path)) {
    renamed <- file.rename(extensionless_path, expected_path)
    if (!isTRUE(renamed)) {
      copied <- file.copy(extensionless_path, expected_path, overwrite = TRUE)
      if (isTRUE(copied)) {
        unlink(extensionless_path, force = TRUE)
      } else {
        stop(
          "Plot was created without its expected extension and could not be renamed: ",
          extensionless_path,
          call. = FALSE
        )
      }
    }
  }

  file.exists(expected_path)
}

write_plot_file <- function(path, plot, branch, width, height, limitsize) {
  path <- ensure_plot_extension(path, branch$format)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  # Remove an extensionless file left by an older run before saving. The current
  # save always requests both an explicit extension and an explicit device.
  extensionless_path <- extensionless_plot_path(path, branch$format)
  if (file.exists(extensionless_path)) {
    unlink(extensionless_path, force = TRUE)
  }

  # Supplying both an explicit filename extension and an explicit graphics device
  # prevents platform-dependent device guessing and guarantees Windows-recognized
  # .png or .pdf filenames.
  
  device_fun <- switch(
    branch$format,
    png = grDevices::png,
    pdf = grDevices::pdf,
    stop("Unsupported plot format: ", branch$format, call. = FALSE)
  )
  
  ggplot2::ggsave(
    filename = basename(path),
    path = dirname(path),
    plot = plot,
    device = branch$format,
    width = width,
    height = height,
    dpi = branch$dpi,
    bg = "white",
    limitsize = limitsize
  )

  if (!recover_extensionless_plot(path, branch$format)) {
    stop("Plot device completed but the expected file was not created: ", path, call. = FALSE)
  }

  invisible(path)
}

make_curve_id <- function(plate_number, plate_id, virus, sample_id) {
  paste(plate_number, plate_id, virus, sample_id, sep = " | ")
}

match_gen5_row <- function(gen5_data, curve_data) {
  if (is.null(gen5_data)) {
    return(NULL)
  }

  matches <- gen5_data[
    gen5_data$plate_number == curve_data$plate_number[[1]] &
      gen5_data$plate_id == curve_data$plate_id[[1]] &
      gen5_data$virus == curve_data$virus[[1]] &
      gen5_data$sample_id == curve_data$sample_id[[1]],
    ,
    drop = FALSE
  ]

  if (nrow(matches) == 0L) NULL else matches[1L, , drop = FALSE]
}

method_report_status <- function(reportable, review_flags) {
  if (!isTRUE(reportable)) {
    "NOT_REPORTABLE"
  } else if (length(unique(review_flags[nzchar(review_flags)])) > 0L) {
    "REPORTABLE_REVIEW_REQUIRED"
  } else {
    "REPORTABLE"
  }
}

make_primary_annotation <- function(result, confidence_level) {
  if (!isTRUE(result$reportable) || !is.finite(result$reported_nt50)) {
    reasons <- format_qc_flags_for_plot(
      result$hard_qc_flags %||% "",
      width = 36L,
      max_flags = 3L,
      bullet = FALSE
    )
    if (!nzchar(reasons)) reasons <- "Primary method unavailable"
    return(paste0("NT50 not reported\n", reasons))
  }

  method_label <- toupper(result$primary_reporting_method)
  if (is.finite(result$reported_ci_lower) && is.finite(result$reported_ci_upper)) {
    paste0(
      method_label, " NT50 = ", format_number(result$reported_nt50),
      "\n", round(100 * confidence_level), "% CI: ",
      format_number(result$reported_ci_lower), "-",
      format_number(result$reported_ci_upper)
    )
  } else {
    paste0(method_label, " NT50 = ", format_number(result$reported_nt50))
  }
}

analyze_one_curve <- function(
  curve_data,
  summary_data,
  normalization_qc_row,
  gen5_row,
  config
) {
  if (nrow(normalization_qc_row) != 1L) {
    stop("Expected exactly one normalization-QC row for the curve.", call. = FALSE)
  }

  threshold <- config$analysis$threshold_pct
  descriptor <- curve_descriptors(summary_data, threshold)
  ll4_fit <- fit_ll4_curve(curve_data, config)
  ll4_metrics <- compute_fit_metrics(ll4_fit$model, curve_data)
  lack_of_fit_p <- extract_lack_of_fit_p(ll4_fit$model)

  ll4_nt50 <- estimate_ll4_nt50(
    ll4_fit$model,
    ll4_fit$parameters,
    summary_data,
    config
  )
  uniroot_nt50 <- estimate_uniroot_nt50(
    summary_data,
    threshold,
    config$analysis$uniroot_interpolation_scale
  )

  ll4_prediction <- predict_ll4_curve(
    ll4_fit$model,
    descriptor$min_dilution,
    descriptor$max_dilution,
    config
  )

  vc_mean <- normalization_qc_row$virus_control_mean[[1]]
  gen5_coefficients <- extract_gen5_coefficients(gen5_row)
  gen5_nt50 <- estimate_gen5_nt50(
    gen5_row,
    vc_mean,
    descriptor$min_dilution,
    descriptor$max_dilution,
    threshold
  )
  gen5_prediction <- predict_gen5_curve(
    gen5_row,
    vc_mean,
    descriptor$min_dilution,
    descriptor$max_dilution,
    config
  )

  ll4_vs_uniroot_fold <- fold_difference(
    ll4_nt50$estimate,
    uniroot_nt50$estimate
  )
  ll4_vs_uniroot_log2_ratio <- log2_ratio(
    ll4_nt50$estimate,
    uniroot_nt50$estimate
  )
  ll4_vs_gen5_fold <- fold_difference(ll4_nt50$estimate, gen5_nt50$estimate)

  shared_hard_flags <- split_flags(
    normalization_qc_row$normalization_hard_qc_flags[[1]]
  )
  shared_review_flags <- split_flags(
    normalization_qc_row$normalization_review_qc_flags[[1]]
  )
  ll4_hard_flags <- character()
  ll4_review_flags <- character()
  uniroot_hard_flags <- character()
  uniroot_review_flags <- character()
  cross_method_review_flags <- character()
  gen5_review_flags <- character()

  if (isTRUE(normalization_qc_row$normalization_available[[1]]) &&
      !descriptor$threshold_bracketed) {
    shared_hard_flags <- c(
      shared_hard_flags,
      "threshold_not_bracketed_by_observed_means"
    )
  }
  if (!is.null(config$qc$max_replicate_sd_pp) &&
      is.finite(descriptor$max_replicate_sd_pp) &&
      descriptor$max_replicate_sd_pp > config$qc$max_replicate_sd_pp) {
    shared_review_flags <- c(shared_review_flags, "high_replicate_sd")
  }
  if (!is.null(config$qc$max_upward_mean_steps) &&
      is.finite(descriptor$n_upward_mean_steps) &&
      descriptor$n_upward_mean_steps > config$qc$max_upward_mean_steps) {
    shared_review_flags <- c(
      shared_review_flags,
      "irregular_nonmonotonic_mean_curve"
    )
  }

  if (is.null(ll4_fit$model)) {
    ll4_hard_flags <- c(ll4_hard_flags, "ll4_fit_unavailable")
  }
  if (!is.finite(ll4_nt50$estimate)) {
    ll4_hard_flags <- c(ll4_hard_flags, "ll4_nt50_unavailable")
  } else if (!isTRUE(ll4_nt50$estimate_in_range)) {
    ll4_hard_flags <- c(ll4_hard_flags, "ll4_nt50_outside_tested_range")
  }
  if (isTRUE(ll4_fit$forced_hill)) {
    ll4_review_flags <- c(
      ll4_review_flags,
      "hill_fixed_to_prespecified_value"
    )
  }
  if (!is.null(config$qc$lack_of_fit_alpha) &&
      is.finite(lack_of_fit_p) &&
      lack_of_fit_p < config$qc$lack_of_fit_alpha) {
    ll4_review_flags <- c(ll4_review_flags, "significant_ll4_lack_of_fit")
  }
  if (!is.null(config$qc$min_ll4_pseudo_r2) &&
      is.finite(ll4_metrics$pseudo_r2) &&
      ll4_metrics$pseudo_r2 < config$qc$min_ll4_pseudo_r2) {
    ll4_review_flags <- c(ll4_review_flags, "low_ll4_pseudo_r2")
  }
  if (!is.null(config$qc$max_ll4_rmse_pp) &&
      is.finite(ll4_metrics$rmse) &&
      ll4_metrics$rmse > config$qc$max_ll4_rmse_pp) {
    ll4_review_flags <- c(ll4_review_flags, "high_ll4_rmse")
  }

  ll4_fit_review_warnings <- model_warnings_requiring_review(ll4_fit$warnings)
  ll4_nt50_review_warnings <- model_warnings_requiring_review(ll4_nt50$warnings)
  # Warning text is exported for manual review but is not promoted to a primary
  # QC flag by itself. Some drc warnings are emitted during starting-value
  # searches even when the final model is finite and visually appropriate.
  if (grepl("requested_ci_failed", ll4_nt50$ci_status, fixed = TRUE)) {
    ll4_review_flags <- c(ll4_review_flags, "ll4_nt50_ci_failed")
  }
  if (isTRUE(ll4_nt50$ci_available) && !isTRUE(ll4_nt50$ci_in_range)) {
    ll4_review_flags <- c(
      ll4_review_flags,
      "ll4_nt50_ci_extends_outside_tested_range"
    )
  }

  if (!is.finite(uniroot_nt50$estimate)) {
    uniroot_hard_flags <- c(uniroot_hard_flags, "uniroot_nt50_unavailable")
  } else if (!isTRUE(uniroot_nt50$estimate_in_range)) {
    uniroot_hard_flags <- c(
      uniroot_hard_flags,
      "uniroot_nt50_outside_tested_range"
    )
  }
  if (identical(
    uniroot_nt50$status,
    "multiple_observed_mean_threshold_crossings"
  )) {
    uniroot_hard_flags <- c(
      uniroot_hard_flags,
      "multiple_observed_mean_threshold_crossings"
    )
  }
  if (length(uniroot_nt50$warnings) > 0L) {
    uniroot_review_flags <- c(uniroot_review_flags, "uniroot_numerical_warning")
  }

  if (!is.null(config$qc$max_ll4_uniroot_disagreement_fold) &&
      is.finite(ll4_vs_uniroot_fold) &&
      ll4_vs_uniroot_fold > config$qc$max_ll4_uniroot_disagreement_fold) {
    cross_method_review_flags <- c(
      cross_method_review_flags,
      "ll4_uniroot_disagreement"
    )
  }
  if (!is.null(config$qc$max_ll4_gen5_disagreement_fold) &&
      is.finite(ll4_vs_gen5_fold) &&
      ll4_vs_gen5_fold > config$qc$max_ll4_gen5_disagreement_fold) {
    gen5_review_flags <- c(gen5_review_flags, "ll4_gen5_disagreement")
  }
  if (!is.null(gen5_row) && !gen5_coefficients$valid) {
    gen5_review_flags <- c(gen5_review_flags, "gen5_fit_coefficients_invalid")
  }

  shared_hard_flags <- unique(shared_hard_flags)
  shared_review_flags <- unique(shared_review_flags)
  ll4_hard_flags <- unique(ll4_hard_flags)
  ll4_review_flags <- unique(ll4_review_flags)
  uniroot_hard_flags <- unique(uniroot_hard_flags)
  uniroot_review_flags <- unique(uniroot_review_flags)
  cross_method_review_flags <- unique(cross_method_review_flags)
  gen5_review_flags <- unique(gen5_review_flags)

  ll4_reportable <- length(c(shared_hard_flags, ll4_hard_flags)) == 0L
  uniroot_reportable <- length(c(shared_hard_flags, uniroot_hard_flags)) == 0L

  ll4_status_review <- c(
    shared_review_flags,
    ll4_review_flags,
    cross_method_review_flags,
    gen5_review_flags
  )
  uniroot_status_review <- c(
    shared_review_flags,
    uniroot_review_flags,
    cross_method_review_flags,
    gen5_review_flags
  )
  ll4_report_status <- method_report_status(ll4_reportable, ll4_status_review)
  uniroot_report_status <- method_report_status(
    uniroot_reportable,
    uniroot_status_review
  )

  ll4_reported_nt50 <- if (ll4_reportable) ll4_nt50$estimate else NA_real_
  uniroot_reported_nt50 <- if (uniroot_reportable) {
    uniroot_nt50$estimate
  } else {
    NA_real_
  }

  primary_reporting_method <- config$analysis$reporting_method
  primary_reportable <- if (primary_reporting_method == "ll4") {
    ll4_reportable
  } else {
    uniroot_reportable
  }
  primary_method_hard_flags <- if (primary_reporting_method == "ll4") {
    ll4_hard_flags
  } else {
    uniroot_hard_flags
  }
  primary_method_review_flags <- if (primary_reporting_method == "ll4") {
    ll4_review_flags
  } else {
    uniroot_review_flags
  }

  primary_hard_flags <- unique(c(shared_hard_flags, primary_method_hard_flags))
  primary_review_flags <- unique(c(
    shared_review_flags,
    primary_method_review_flags,
    cross_method_review_flags,
    gen5_review_flags
  ))
  report_status <- method_report_status(primary_reportable, primary_review_flags)

  reported_nt50 <- if (!primary_reportable) {
    NA_real_
  } else if (primary_reporting_method == "ll4") {
    ll4_nt50$estimate
  } else {
    uniroot_nt50$estimate
  }

  reported_ci_lower <- NA_real_
  reported_ci_upper <- NA_real_
  reported_ci_method <- if (primary_reporting_method == "uniroot") {
    "not_available_for_observed_value_uniroot"
  } else if (config$analysis$ll4_ci_method == "none") {
    "not_requested"
  } else {
    "unavailable_or_outside_tested_range"
  }
  if (primary_reportable && primary_reporting_method == "ll4" &&
      isTRUE(ll4_nt50$ci_available) && isTRUE(ll4_nt50$ci_in_range)) {
    reported_ci_lower <- ll4_nt50$lower
    reported_ci_upper <- ll4_nt50$upper
    reported_ci_method <- config$analysis$ll4_ci_method
  }

  ll4_available <- is.finite(ll4_nt50$estimate)
  uniroot_available <- is.finite(uniroot_nt50$estimate)
  method_availability_status <- if (ll4_available && uniroot_available) {
    "both_estimates_available"
  } else if (ll4_available) {
    "ll4_estimate_only"
  } else if (uniroot_available) {
    "uniroot_estimate_only"
  } else {
    "neither_estimate_available"
  }

  curve_id <- make_curve_id(
    curve_data$plate_number[[1]],
    curve_data$plate_id[[1]],
    curve_data$virus[[1]],
    curve_data$sample_id[[1]]
  )
  ll4_fit_strategy <- if (is.null(ll4_fit$model)) {
    "failed"
  } else if (isTRUE(ll4_fit$forced_hill)) {
    "fixed_hill_refit"
  } else {
    "unconstrained"
  }

  result <- data.frame(
    analysis_script_version = SCRIPT_VERSION,
    curve_id = curve_id,
    plate_number = curve_data$plate_number[[1]],
    plate_id = curve_data$plate_id[[1]],
    virus = curve_data$virus[[1]],
    sample_id = curve_data$sample_id[[1]],
    primary_reporting_method = primary_reporting_method,
    threshold_pct = threshold,
    report_status = report_status,
    reportable = primary_reportable,
    reported_nt50 = reported_nt50,
    reported_ci_lower = reported_ci_lower,
    reported_ci_upper = reported_ci_upper,
    reported_ci_method = reported_ci_method,
    method_availability_status = method_availability_status,
    both_methods_reportable = ll4_reportable && uniroot_reportable,
    any_method_reportable = ll4_reportable || uniroot_reportable,
    alternative_method_reportable = if (primary_reporting_method == "ll4") {
      uniroot_reportable
    } else {
      ll4_reportable
    },
    tested_min_dilution = descriptor$min_dilution,
    tested_max_dilution = descriptor$max_dilution,
    n_unique_dilutions = descriptor$n_unique_dilutions,
    min_replicates_per_dilution = descriptor$min_replicates_per_dilution,
    max_replicates_per_dilution = descriptor$max_replicates_per_dilution,
    threshold_bracketed_by_observed_means = descriptor$threshold_bracketed,
    min_observed_mean_percent = descriptor$min_observed_mean,
    max_observed_mean_percent = descriptor$max_observed_mean,
    mean_curve_spearman_rho = descriptor$spearman_rho,
    n_upward_mean_steps = descriptor$n_upward_mean_steps,
    max_replicate_sd_pp = descriptor$max_replicate_sd_pp,
    normalization_source = normalization_qc_row$normalization_source[[1]],
    normalization_available = normalization_qc_row$normalization_available[[1]],
    normalization_reportable = normalization_qc_row$normalization_reportable[[1]],
    normalization_reference_is_plate_matched =
      normalization_qc_row$normalization_reference_is_plate_matched[[1]],
    plate_matched_virus_control_available =
      normalization_qc_row$plate_matched_virus_control_available[[1]],
    plate_matched_vc_required_metadata_inference =
      normalization_qc_row$plate_matched_vc_required_metadata_inference[[1]],
    normalization_reference_plate_numbers =
      normalization_qc_row$normalization_reference_plate_numbers[[1]],
    normalization_reference_plate_ids =
      normalization_qc_row$normalization_reference_plate_ids[[1]],
    n_normalization_reference_plates =
      normalization_qc_row$n_normalization_reference_plates[[1]],
    n_virus_control = normalization_qc_row$n_virus_control[[1]],
    minimum_required_virus_controls =
      normalization_qc_row$minimum_required_virus_controls[[1]],
    virus_control_replicate_requirement_met =
      normalization_qc_row$virus_control_replicate_requirement_met[[1]],
    virus_control_mean = vc_mean,
    virus_control_sd = normalization_qc_row$virus_control_sd[[1]],
    virus_control_cv_pct = normalization_qc_row$virus_control_cv_pct[[1]],
    virus_control_mean_in_working_range =
      normalization_qc_row$virus_control_mean_in_working_range[[1]],
    all_virus_control_wells_in_working_range =
      normalization_qc_row$all_virus_control_wells_in_working_range[[1]],
    n_virus_control_outside_working_range =
      normalization_qc_row$n_virus_control_outside_working_range[[1]],
    virus_control_in_working_range =
      normalization_qc_row$virus_control_in_working_range[[1]],
    cell_control_mean = normalization_qc_row$cell_control_mean[[1]],
    cell_control_max = normalization_qc_row$cell_control_max[[1]],
    cell_control_background_high =
      normalization_qc_row$cell_control_background_high[[1]],
    control_virus_sources = normalization_qc_row$control_virus_sources[[1]],
    normalization_status = normalization_qc_row$normalization_status[[1]],
    normalization_hard_qc_flags =
      normalization_qc_row$normalization_hard_qc_flags[[1]],
    normalization_review_qc_flags =
      normalization_qc_row$normalization_review_qc_flags[[1]],
    ll4_report_status = ll4_report_status,
    ll4_reportable = ll4_reportable,
    ll4_reported_nt50 = ll4_reported_nt50,
    ll4_fit_available = !is.null(ll4_fit$model),
    ll4_estimate_available = ll4_available,
    ll4_calculated_nt50 = ll4_nt50$estimate,
    ll4_fit_strategy = ll4_fit_strategy,
    ll4_forced_hill = ll4_fit$forced_hill,
    ll4_initial_hill = ll4_fit$initial_hill,
    ll4_hill = ll4_fit$parameters$hill,
    ll4_lower_asymptote = ll4_fit$parameters$lower,
    ll4_upper_asymptote = ll4_fit$parameters$upper,
    ll4_midpoint_parameter = ll4_fit$parameters$midpoint,
    ll4_rmse = ll4_metrics$rmse,
    ll4_pseudo_r2 = ll4_metrics$pseudo_r2,
    ll4_residual_df = ll4_metrics$residual_df,
    ll4_lack_of_fit_p = lack_of_fit_p,
    ll4_nt50_se = ll4_nt50$se,
    ll4_nt50_ci_lower = ll4_nt50$lower,
    ll4_nt50_ci_upper = ll4_nt50$upper,
    ll4_nt50_estimate_in_range = ll4_nt50$estimate_in_range,
    ll4_nt50_ci_available = ll4_nt50$ci_available,
    ll4_nt50_ci_in_range = ll4_nt50$ci_in_range,
    ll4_nt50_status = ll4_nt50$status,
    ll4_nt50_ci_status = ll4_nt50$ci_status,
    ll4_fit_warning_count = length(ll4_fit$warnings),
    ll4_fit_review_warning_count = length(ll4_fit_review_warnings),
    ll4_fit_warnings = paste(ll4_fit$warnings, collapse = " | "),
    ll4_fit_review_warnings = paste(ll4_fit_review_warnings, collapse = " | "),
    ll4_fit_error = ll4_fit$error %||% "",
    ll4_nt50_warning_count = length(ll4_nt50$warnings),
    ll4_nt50_review_warning_count = length(ll4_nt50_review_warnings),
    ll4_nt50_warnings = paste(ll4_nt50$warnings, collapse = " | "),
    ll4_nt50_review_warnings = paste(
      ll4_nt50_review_warnings,
      collapse = " | "
    ),
    ll4_nt50_error = ll4_nt50$error %||% "",
    uniroot_report_status = uniroot_report_status,
    uniroot_reportable = uniroot_reportable,
    uniroot_reported_nt50 = uniroot_reported_nt50,
    uniroot_estimate_available = uniroot_available,
    uniroot_calculated_nt50 = uniroot_nt50$estimate,
    uniroot_estimate_in_range = uniroot_nt50$estimate_in_range,
    uniroot_status = uniroot_nt50$status,
    uniroot_method = uniroot_nt50$method,
    uniroot_interpolation_scale = uniroot_nt50$interpolation_scale,
    uniroot_n_tested_dilutions = uniroot_nt50$n_dilutions,
    uniroot_n_crossings = uniroot_nt50$n_crossings,
    uniroot_warning_count = length(uniroot_nt50$warnings),
    uniroot_warnings = paste(uniroot_nt50$warnings, collapse = " | "),
    uniroot_error = uniroot_nt50$error %||% "",
    ll4_vs_uniroot_fold_difference = ll4_vs_uniroot_fold,
    ll4_vs_uniroot_log2_ratio = ll4_vs_uniroot_log2_ratio,
    ll4_vs_uniroot_abs_log2_difference = abs(ll4_vs_uniroot_log2_ratio),
    gen5_fit_available = !is.null(gen5_row),
    gen5_fit_valid = gen5_coefficients$valid,
    gen5_a = gen5_coefficients$a,
    gen5_b = gen5_coefficients$b,
    gen5_c = gen5_coefficients$c,
    gen5_d = gen5_coefficients$d,
    gen5_r2 = if (is.null(gen5_row)) NA_real_ else gen5_row$r2[[1]],
    gen5_fit_f_prob = if (is.null(gen5_row)) NA_real_ else gen5_row$fit_f_prob[[1]],
    gen5_calculated_nt50 = gen5_nt50$estimate,
    gen5_nt50_status = gen5_nt50$status,
    ll4_vs_gen5_fold_difference = ll4_vs_gen5_fold,
    qc_vc_count_lower = config$qc$vc_count_range[[1]],
    qc_vc_count_upper = config$qc$vc_count_range[[2]],
    qc_hill_lower = config$qc$hill_range[[1]],
    qc_hill_upper = config$qc$hill_range[[2]],
    qc_fixed_hill = config$qc$fixed_hill,
    shared_hard_qc_flags = collapse_flags(shared_hard_flags),
    shared_review_qc_flags = collapse_flags(shared_review_flags),
    ll4_hard_qc_flags = collapse_flags(ll4_hard_flags),
    ll4_review_qc_flags = collapse_flags(ll4_review_flags),
    ll4_all_qc_flags = collapse_flags(c(
      shared_hard_flags,
      ll4_hard_flags,
      shared_review_flags,
      ll4_review_flags,
      cross_method_review_flags,
      gen5_review_flags
    )),
    uniroot_hard_qc_flags = collapse_flags(uniroot_hard_flags),
    uniroot_review_qc_flags = collapse_flags(uniroot_review_flags),
    uniroot_all_qc_flags = collapse_flags(c(
      shared_hard_flags,
      uniroot_hard_flags,
      shared_review_flags,
      uniroot_review_flags,
      cross_method_review_flags,
      gen5_review_flags
    )),
    cross_method_review_qc_flags = collapse_flags(cross_method_review_flags),
    gen5_review_qc_flags = collapse_flags(gen5_review_flags),
    hard_qc_flags = collapse_flags(primary_hard_flags),
    review_qc_flags = collapse_flags(primary_review_flags),
    all_qc_flags = collapse_flags(c(primary_hard_flags, primary_review_flags)),
    stringsAsFactors = FALSE
  )
  result$plot_annotation <- make_primary_annotation(
    result,
    config$analysis$confidence_level
  )

  list(
    result = result,
    raw = curve_data,
    summary = summary_data,
    ll4_prediction = ll4_prediction,
    gen5_prediction = gen5_prediction
  )
}

analyze_all_curves <- function(
  normalized_data,
  dilution_summary,
  normalization_qc,
  gen5_data,
  config
) {
  sample_data <- normalized_data[
    normalized_data$control_type == "sample",
    ,
    drop = FALSE
  ]
  curve_keys <- unique(
    sample_data[, c("plate_number", "plate_id", "virus", "sample_id")]
  )
  plate_key <- first_numeric_token(curve_keys$plate_number)
  curve_keys <- curve_keys[order(
    curve_keys$virus,
    curve_keys$sample_id,
    is.na(plate_key),
    plate_key,
    tolower(curve_keys$plate_number),
    curve_keys$plate_id
  ), , drop = FALSE]

  if (nrow(curve_keys) == 0L) {
    stop("No sample curves were available for analysis.", call. = FALSE)
  }

  objects <- vector("list", nrow(curve_keys))
  for (i in seq_len(nrow(curve_keys))) {
    key <- curve_keys[i, , drop = FALSE]
    message(
      "Analyzing ", key$plate_number, " | ", key$virus, " | ", key$sample_id
    )

    curve_data <- sample_data[
      sample_data$plate_number == key$plate_number &
        sample_data$plate_id == key$plate_id &
        sample_data$virus == key$virus &
        sample_data$sample_id == key$sample_id,
      ,
      drop = FALSE
    ]
    summary_data <- dilution_summary[
      dilution_summary$plate_number == key$plate_number &
        dilution_summary$plate_id == key$plate_id &
        dilution_summary$virus == key$virus &
        dilution_summary$sample_id == key$sample_id,
      ,
      drop = FALSE
    ]
    normalization_qc_row <- normalization_qc[
      normalization_qc$plate_number == key$plate_number &
        normalization_qc$plate_id == key$plate_id &
        normalization_qc$virus == key$virus,
      ,
      drop = FALSE
    ]
    gen5_row <- match_gen5_row(gen5_data, curve_data)

    # Expected fit failures are returned as method-specific NA/status values by
    # the fitting functions. An unexpected programming or schema error should
    # stop the run rather than silently writing a reduced result schema.
    objects[[i]] <- analyze_one_curve(
      curve_data,
      summary_data,
      normalization_qc_row,
      gen5_row,
      config
    )
  }

  list(
    objects = objects,
    results = dplyr::bind_rows(lapply(objects, `[[`, "result"))
  )
}

select_dilution_breaks <- function(values, max_breaks) {
  values <- sort(unique(values[is.finite(values) & values > 0]))
  if (length(values) <= max_breaks) {
    return(values)
  }

  indices <- unique(as.integer(round(seq(1, length(values), length.out = max_breaks))))
  values[indices]
}

result_value <- function(result, field, default = NA) {
  if (!field %in% names(result) || length(result[[field]]) == 0L) {
    return(default)
  }
  value <- result[[field]][[1]]
  if (length(value) == 0L) default else value
}

add_branch_y_axis <- function(plot, branch) {
  if (branch$y_axis == "fixed") {
    plot + ggplot2::coord_cartesian(ylim = branch$y_limits)
  } else {
    plot
  }
}

display_requests_ll4 <- function(result, branch) {
  display <- branch$nt50_display
  if (display %in% c("ll4", "both")) {
    return(TRUE)
  }
  display == "primary" && identical(
    as.character(result_value(result, "primary_reporting_method", "")),
    "ll4"
  )
}

empty_nt50_line_data <- function() {
  data.frame(
    method = character(),
    xintercept = numeric(),
    line_color = character(),
    line_type = character(),
    stringsAsFactors = FALSE
  )
}

make_nt50_line_data <- function(result, branch) {
  display <- branch$nt50_display
  rows <- list()

  add_line <- function(method, estimate, reportable, color, linetype) {
    if (!isTRUE(reportable) || !is.finite(estimate) || estimate <= 0) {
      return(NULL)
    }
    data.frame(
      method = method,
      xintercept = estimate,
      line_color = color,
      line_type = linetype,
      stringsAsFactors = FALSE
    )
  }

  if (display == "none") {
    return(empty_nt50_line_data())
  }

  if (display == "primary") {
    rows[[1L]] <- add_line(
      paste0("primary_", result_value(result, "primary_reporting_method", "unknown")),
      suppressWarnings(as.numeric(result_value(result, "reported_nt50", NA_real_))),
      isTRUE(result_value(result, "reportable", FALSE)),
      "#333333",
      "dotted"
    )
  } else {
    if (display %in% c("ll4", "both")) {
      rows[[length(rows) + 1L]] <- add_line(
        "ll4",
        suppressWarnings(as.numeric(result_value(result, "ll4_reported_nt50", NA_real_))),
        isTRUE(result_value(result, "ll4_reportable", FALSE)),
        "#0072B2",
        "dotted"
      )
    }
    if (display %in% c("uniroot", "both")) {
      rows[[length(rows) + 1L]] <- add_line(
        "uniroot",
        suppressWarnings(as.numeric(result_value(result, "uniroot_reported_nt50", NA_real_))),
        isTRUE(result_value(result, "uniroot_reportable", FALSE)),
        "#CC79A7",
        "longdash"
      )
    }
  }

  line_data <- dplyr::bind_rows(rows)
  if (nrow(line_data) == 0L || !"method" %in% names(line_data)) {
    return(empty_nt50_line_data())
  }
  if (nrow(line_data) == 2L && display == "both") {
    ratio <- abs(log2(line_data$xintercept[[1]] / line_data$xintercept[[2]]))
    if (is.finite(ratio) && ratio <= 0.01) {
      line_data <- data.frame(
        method = "ll4_and_uniroot_coincident",
        xintercept = line_data$xintercept[[1]],
        line_color = "#333333",
        line_type = "dotdash",
        stringsAsFactors = FALSE
      )
    }
  }

  line_data
}

make_nt50_ci_data <- function(result, branch) {
  if (!isTRUE(branch$show_nt50_ci) || !display_requests_ll4(result, branch)) {
    return(data.frame())
  }

  lower <- suppressWarnings(as.numeric(result_value(
    result,
    "ll4_nt50_ci_lower",
    NA_real_
  )))
  upper <- suppressWarnings(as.numeric(result_value(
    result,
    "ll4_nt50_ci_upper",
    NA_real_
  )))
  reportable <- isTRUE(result_value(result, "ll4_reportable", FALSE))
  in_range <- isTRUE(result_value(result, "ll4_nt50_ci_in_range", FALSE))
  if (!reportable || !in_range || !is.finite(lower) || !is.finite(upper)) {
    return(data.frame())
  }

  data.frame(xmin = lower, xmax = upper, stringsAsFactors = FALSE)
}

method_annotation_line <- function(result, method, include_ci = TRUE) {
  if (method == "ll4") {
    reportable <- isTRUE(result_value(result, "ll4_reportable", FALSE))
    value <- suppressWarnings(as.numeric(result_value(
      result,
      "ll4_reported_nt50",
      NA_real_
    )))
    status <- as.character(result_value(result, "ll4_report_status", "NOT_REPORTABLE"))
    if (!reportable || !is.finite(value)) {
      return(paste0("LL4: ", status))
    }
    label <- paste0("LL4 NT50 = ", format_number(value))
    lower <- suppressWarnings(as.numeric(result_value(
      result,
      "ll4_nt50_ci_lower",
      NA_real_
    )))
    upper <- suppressWarnings(as.numeric(result_value(
      result,
      "ll4_nt50_ci_upper",
      NA_real_
    )))
    ci_in_range <- isTRUE(result_value(result, "ll4_nt50_ci_in_range", FALSE))
    if (include_ci && ci_in_range && is.finite(lower) && is.finite(upper)) {
      label <- paste0(
        label,
        " (CI ", format_number(lower), "-", format_number(upper), ")"
      )
    }
    return(label)
  }

  reportable <- isTRUE(result_value(result, "uniroot_reportable", FALSE))
  value <- suppressWarnings(as.numeric(result_value(
    result,
    "uniroot_reported_nt50",
    NA_real_
  )))
  status <- as.character(result_value(result, "uniroot_report_status", "NOT_REPORTABLE"))
  if (!reportable || !is.finite(value)) {
    return(paste0("Uniroot: ", status))
  }
  paste0("Uniroot NT50 = ", format_number(value))
}

make_plot_annotation <- function(result, branch, concise = FALSE) {
  display <- branch$nt50_display
  if (display == "none") {
    return(paste0(
      "Primary status: ",
      result_value(result, "report_status", "UNKNOWN")
    ))
  }
  if (display == "primary") {
    annotation <- as.character(result_value(result, "plot_annotation", ""))
    if (nzchar(annotation)) return(annotation)
    return(paste0(
      toupper(result_value(result, "primary_reporting_method", "primary")),
      ": ", result_value(result, "report_status", "UNKNOWN")
    ))
  }
  if (display == "ll4") {
    return(method_annotation_line(result, "ll4", include_ci = !concise))
  }
  if (display == "uniroot") {
    return(method_annotation_line(result, "uniroot", include_ci = FALSE))
  }

  lines <- c(
    method_annotation_line(result, "ll4", include_ci = !concise),
    method_annotation_line(result, "uniroot", include_ci = FALSE)
  )
  fold <- suppressWarnings(as.numeric(result_value(
    result,
    "ll4_vs_uniroot_fold_difference",
    NA_real_
  )))
  if (!concise && is.finite(fold)) {
    lines <- c(lines, paste0("Method fold difference = ", format_number(fold)))
  }
  paste(lines, collapse = "\n")
}

build_raw_focus_fallback_plot <- function(
    result,
    raw,
    config,
    branch,
    breaks
) {
  raw_finite <- raw[
    is.finite(raw$dilution) &
      raw$dilution > 0 &
      is.finite(raw$focus_count),
    ,
    drop = FALSE
  ]
  
  if (length(breaks) == 0L && nrow(raw_finite) > 0L) {
    breaks <- select_dilution_breaks(
      raw_finite$dilution,
      branch$max_x_breaks
    )
  }
  
  plot <- ggplot2::ggplot()
  
  if (nrow(raw_finite) > 0L) {
    dilution_means <- stats::aggregate(
      focus_count ~ dilution,
      data = raw_finite,
      FUN = function(values) {
        values <- values[is.finite(values)]
        if (length(values) == 0L) NA_real_ else mean(values)
      }
    )
    
    plot <- plot +
      ggplot2::geom_point(
        data = raw_finite,
        ggplot2::aes(
          x = .data$dilution,
          y = .data$focus_count
        ),
        inherit.aes = FALSE,
        color = "grey30",
        alpha = 0.70,
        size = 1.9,
        na.rm = TRUE
      ) +
      ggplot2::geom_point(
        data = dilution_means,
        ggplot2::aes(
          x = .data$dilution,
          y = .data$focus_count
        ),
        inherit.aes = FALSE,
        shape = 21,
        fill = "white",
        color = "black",
        stroke = 0.7,
        size = 2.5,
        na.rm = TRUE
      )
  }
  
  flag_caption <- format_qc_flags_for_plot(
    result_value(result, "all_qc_flags", ""),
    width = config$plot$flag_wrap_width,
    max_flags = Inf,
    bullet = TRUE,
    blank_lines_between_flags =
      config$plot$blank_lines_between_flags
  )
  
  caption_parts <- c(
    paste0(
      "Primary status: ",
      result_value(result, "report_status", "UNKNOWN")
    ),
    if (nzchar(flag_caption)) {
      paste0("Primary flags:\n\n", flag_caption)
    } else {
      "Percent neutralization was not calculated because no usable virus-control reference was selected."
    }
  )
  
  if (isTRUE(branch$show_comparison_details)) {
    plot <- plot + ggplot2::annotate(
      "label",
      x = Inf,
      y = Inf,
      label = "Raw focus counts shown\nNT50 not reported",
      hjust = 1.04,
      vjust = 1.08,
      size = 3.0,
      linewidth = 0.25,
      fill = "white"
    )
  }
  
  plot +
    ggplot2::scale_x_continuous(
      trans = "log2",
      breaks = breaks,
      labels = scales::label_comma(),
      guide = ggplot2::guide_axis(
        angle = if (length(breaks) > 8L) 45 else 0
      )
    ) +
    ggplot2::labs(
      title = result_value(result, "sample_id", "Sample"),
      subtitle = paste(
        result_value(result, "virus", ""),
        result_value(result, "plate_number", ""),
        result_value(result, "plate_id", ""),
        sep = " | "
      ),
      x = "Reciprocal dilution",
      y = "Raw focus count (normalization unavailable)",
      caption = paste(caption_parts, collapse = "\n")
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.caption = ggplot2::element_text(
        hjust = 0,
        size = 8.5,
        lineheight = 1.08
      ),
      plot.margin = ggplot2::margin(
        t = 8,
        r = 10,
        b = 18,
        l = 10
      ),
      panel.grid.minor = ggplot2::element_blank()
    )
}


build_individual_plot <- function(object, config, branch) {
  result <- object$result[1L, , drop = FALSE]
  raw <- object$raw
  summary_data <- object$summary
  ll4_prediction <- object$ll4_prediction
  gen5_prediction <- object$gen5_prediction
  threshold <- config$analysis$threshold_pct
  breaks <- select_dilution_breaks(summary_data$dilution, branch$max_x_breaks)

  plot <- ggplot2::ggplot()
  has_normalized_data <-
    any(is.finite(raw$percent_neutralization)) ||
    any(is.finite(summary_data$mean_percent_neutralization))

  if (!has_normalized_data) {
    return(build_raw_focus_fallback_plot(
      result = result,
      raw = raw,
      config = config,
      branch = branch,
      breaks = breaks
    ))
  }
  
  if (isTRUE(branch$show_model_confidence_band) &&
      nrow(ll4_prediction) > 0L &&
      any(ll4_prediction$model_band_available %in% TRUE)) {
    plot <- plot + ggplot2::geom_ribbon(
      data = ll4_prediction,
      ggplot2::aes(x = .data$dilution, ymin = .data$lower, ymax = .data$upper),
      inherit.aes = FALSE,
      fill = "#56B4E9",
      alpha = 0.18
    )
  }

  ci_data <- make_nt50_ci_data(result, branch)
  if (nrow(ci_data) > 0L) {
    plot <- plot + ggplot2::annotate(
      "rect",
      xmin = ci_data$xmin[[1]],
      xmax = ci_data$xmax[[1]],
      ymin = -Inf,
      ymax = Inf,
      fill = "#0072B2",
      alpha = 0.08
    )
  }

  if (isTRUE(branch$show_ll4_curve) && nrow(ll4_prediction) > 0L) {
    plot <- plot + ggplot2::geom_line(
      data = ll4_prediction,
      ggplot2::aes(x = .data$dilution, y = .data$fitted),
      inherit.aes = FALSE,
      color = "#0072B2",
      linewidth = 0.9
    )
  }
  if (isTRUE(branch$show_gen5_overlay) && nrow(gen5_prediction) > 0L) {
    plot <- plot + ggplot2::geom_line(
      data = gen5_prediction,
      ggplot2::aes(x = .data$dilution, y = .data$fitted),
      inherit.aes = FALSE,
      color = "#D55E00",
      linetype = "dotdash",
      linewidth = 0.8
    )
  }

  if (branch$point_display %in% c("raw", "both")) {
    plot <- plot + ggplot2::geom_point(
      data = raw,
      ggplot2::aes(x = .data$dilution, y = .data$percent_neutralization),
      inherit.aes = FALSE,
      color = "grey40",
      alpha = 0.55,
      size = 1.7,
      na.rm = TRUE
    )
  }
  if (branch$point_display %in% c("summary", "both")) {
    plot <- plot +
      ggplot2::geom_errorbar(
        data = summary_data,
        ggplot2::aes(
          x = .data$dilution,
          ymin = .data$mean_percent_neutralization - .data$sd_percent_neutralization,
          ymax = .data$mean_percent_neutralization + .data$sd_percent_neutralization
        ),
        inherit.aes = FALSE,
        width = 0,
        linewidth = 0.45,
        na.rm = TRUE
      ) +
      ggplot2::geom_point(
        data = summary_data,
        ggplot2::aes(x = .data$dilution, y = .data$mean_percent_neutralization),
        inherit.aes = FALSE,
        shape = 21,
        fill = "white",
        color = "black",
        stroke = 0.7,
        size = 2.5,
        na.rm = TRUE
      )
  }

  plot <- plot + ggplot2::geom_hline(
    yintercept = threshold,
    linetype = "dashed",
    color = "black",
    linewidth = 0.55
  )

  nt50_lines <- make_nt50_line_data(result, branch)
  if (nrow(nt50_lines) > 0L) {
    plot <- plot + ggplot2::geom_vline(
      data = nt50_lines,
      ggplot2::aes(
        xintercept = .data$xintercept,
        color = .data$line_color,
        linetype = .data$line_type
      ),
      linewidth = 0.7,
      show.legend = FALSE
    ) +
      ggplot2::scale_color_identity() +
      ggplot2::scale_linetype_identity()
  }

  flag_caption <- format_qc_flags_for_plot(
    result_value(result, "all_qc_flags", ""),
    width = config$plot$flag_wrap_width,
    max_flags = Inf,
    bullet = TRUE,
    blank_lines_between_flags = config$plot$blank_lines_between_flags
  )
  
  caption_parts <- c(
    paste0("Primary status: ", result_value(result, "report_status", "UNKNOWN")),
    if (nzchar(flag_caption)) paste0("Primary flags:\n\n", flag_caption) else NULL,
    if (isTRUE(branch$show_ll4_curve) && nrow(ll4_prediction) > 0L) {
      "Solid blue: R LL4 fit"
    } else NULL,
    if (any(nt50_lines$method == "ll4")) {
      "Dotted blue: LL4 NT50"
    } else NULL,
    if (any(nt50_lines$method == "uniroot")) {
      "Longdash pink: uniroot NT50"
    } else NULL,
    if (any(nt50_lines$method == "ll4_and_uniroot_coincident")) {
      "Dot-dash grey: LL4 and uniroot NT50 (coincident)"
    } else NULL,
    if (isTRUE(branch$show_gen5_overlay) && nrow(gen5_prediction) > 0L) {
      "Orange dot-dash: Gen5 fit"
    } else NULL
  )


  if (isTRUE(branch$show_comparison_details)) {
    plot <- plot + ggplot2::annotate(
      "label",
      x = Inf,
      y = Inf,
      label = make_plot_annotation(result, branch, concise = FALSE),
      hjust = 1.04,
      vjust = 1.08,
      size = 3.0,
      linewidth = 0.25,
      fill = "white"
    )
  }

  plot <- plot +
    ggplot2::scale_x_continuous(
      trans = "log2",
      breaks = breaks,
      labels = scales::label_comma(),
      guide = ggplot2::guide_axis(angle = if (length(breaks) > 8L) 45 else 0)
    ) +
    ggplot2::labs(
      title = result_value(result, "sample_id", "Sample"),
      subtitle = paste(
        result_value(result, "virus", ""),
        result_value(result, "plate_number", ""),
        result_value(result, "plate_id", ""),
        sep = " | "
      ),
      x = "Reciprocal dilution",
      y = "Percent neutralization",
      caption = paste(caption_parts, collapse = "\n")
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.caption = ggplot2::element_text(
        hjust = 0,
        size = 8.5,
        lineheight = 1.08
      ),
      plot.margin = ggplot2::margin(t = 8, r = 10, b = 18, l = 10),
      panel.grid.minor = ggplot2::element_blank()
    )

  add_branch_y_axis(plot, branch)
}

prepare_faceted_data <- function(objects, branch) {
  raw_rows <- list()
  summary_rows <- list()
  ll4_rows <- list()
  gen5_rows <- list()
  annotation_rows <- list()
  nt50_rows <- list()
  ci_rows <- list()

  for (i in seq_along(objects)) {
    object <- objects[[i]]
    result <- object$result[1L, , drop = FALSE]
    facet_label <- paste(
      result_value(result, "sample_id", "Sample"),
      result_value(result, "plate_number", ""),
      result_value(result, "plate_id", ""),
      sep = "\n"
    )

    raw <- object$raw
    raw$facet_label <- facet_label
    raw_rows[[i]] <- raw

    summary_data <- object$summary
    summary_data$facet_label <- facet_label
    summary_rows[[i]] <- summary_data

    ll4_prediction <- object$ll4_prediction
    if (nrow(ll4_prediction) > 0L) {
      ll4_prediction$facet_label <- facet_label
      ll4_rows[[i]] <- ll4_prediction
    }

    gen5_prediction <- object$gen5_prediction
    if (nrow(gen5_prediction) > 0L) {
      gen5_prediction$facet_label <- facet_label
      gen5_rows[[i]] <- gen5_prediction
    }

    if (isTRUE(branch$show_comparison_details)) {
      annotation_rows[[i]] <- data.frame(
        facet_label = facet_label,
        label = make_plot_annotation(result, branch, concise = TRUE),
        stringsAsFactors = FALSE
      )
    }

    lines <- make_nt50_line_data(result, branch)
    if (nrow(lines) > 0L) {
      lines$facet_label <- facet_label
      nt50_rows[[i]] <- lines
    }

    ci <- make_nt50_ci_data(result, branch)
    if (nrow(ci) > 0L) {
      ci$facet_label <- facet_label
      ci_rows[[i]] <- ci
    }
  }

  list(
    raw = dplyr::bind_rows(raw_rows),
    summary = dplyr::bind_rows(summary_rows),
    ll4 = dplyr::bind_rows(ll4_rows),
    gen5 = dplyr::bind_rows(gen5_rows),
    annotation = dplyr::bind_rows(annotation_rows),
    nt50 = dplyr::bind_rows(nt50_rows),
    ci = dplyr::bind_rows(ci_rows)
  )
}

build_faceted_plot <- function(objects, virus_name, config, branch) {
  data <- prepare_faceted_data(objects, branch)
  threshold <- config$analysis$threshold_pct
  plot <- ggplot2::ggplot()

  if (nrow(data$ci) > 0L) {
    plot <- plot + ggplot2::geom_rect(
      data = data$ci,
      ggplot2::aes(
        xmin = .data$xmin,
        xmax = .data$xmax,
        ymin = -Inf,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = "#0072B2",
      alpha = 0.08
    )
  }
  if (isTRUE(branch$show_model_confidence_band) &&
      nrow(data$ll4) > 0L &&
      any(data$ll4$model_band_available %in% TRUE)) {
    plot <- plot + ggplot2::geom_ribbon(
      data = data$ll4,
      ggplot2::aes(x = .data$dilution, ymin = .data$lower, ymax = .data$upper),
      inherit.aes = FALSE,
      fill = "#56B4E9",
      alpha = 0.18
    )
  }
  if (isTRUE(branch$show_ll4_curve) && nrow(data$ll4) > 0L) {
    plot <- plot + ggplot2::geom_line(
      data = data$ll4,
      ggplot2::aes(x = .data$dilution, y = .data$fitted),
      inherit.aes = FALSE,
      color = "#0072B2",
      linewidth = 0.8
    )
  }
  if (isTRUE(branch$show_gen5_overlay) && nrow(data$gen5) > 0L) {
    plot <- plot + ggplot2::geom_line(
      data = data$gen5,
      ggplot2::aes(x = .data$dilution, y = .data$fitted),
      inherit.aes = FALSE,
      color = "#D55E00",
      linetype = "dotdash",
      linewidth = 0.7
    )
  }
  if (branch$point_display %in% c("raw", "both") && nrow(data$raw) > 0L) {
    plot <- plot + ggplot2::geom_point(
      data = data$raw,
      ggplot2::aes(x = .data$dilution, y = .data$percent_neutralization),
      inherit.aes = FALSE,
      color = "grey45",
      alpha = 0.5,
      size = 1.1,
      na.rm = TRUE
    )
  }
  if (branch$point_display %in% c("summary", "both") && nrow(data$summary) > 0L) {
    plot <- plot +
      ggplot2::geom_errorbar(
        data = data$summary,
        ggplot2::aes(
          x = .data$dilution,
          ymin = .data$mean_percent_neutralization - .data$sd_percent_neutralization,
          ymax = .data$mean_percent_neutralization + .data$sd_percent_neutralization
        ),
        inherit.aes = FALSE,
        width = 0,
        linewidth = 0.35,
        na.rm = TRUE
      ) +
      ggplot2::geom_point(
        data = data$summary,
        ggplot2::aes(x = .data$dilution, y = .data$mean_percent_neutralization),
        inherit.aes = FALSE,
        shape = 21,
        fill = "white",
        color = "black",
        stroke = 0.55,
        size = 2.0,
        na.rm = TRUE
      )
  }

  plot <- plot + ggplot2::geom_hline(
    yintercept = threshold,
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  )
  if (nrow(data$nt50) > 0L) {
    plot <- plot + ggplot2::geom_vline(
      data = data$nt50,
      ggplot2::aes(
        xintercept = .data$xintercept,
        color = .data$line_color,
        linetype = .data$line_type
      ),
      linewidth = 0.6,
      show.legend = FALSE
    ) +
      ggplot2::scale_color_identity() +
      ggplot2::scale_linetype_identity()
  }
  if (nrow(data$annotation) > 0L) {
    plot <- plot + ggplot2::geom_label(
      data = data$annotation,
      ggplot2::aes(x = Inf, y = Inf, label = .data$label),
      inherit.aes = FALSE,
      hjust = 1.04,
      vjust = 1.08,
      size = 2.4,
      linewidth = 0.2,
      fill = "white"
    )
  }

  plot <- plot +
    ggplot2::scale_x_continuous(
      trans = "log2",
      labels = scales::label_comma()
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(facet_label),
      scales = "free_x",
      ncol = branch$facet_columns
    ) +
    ggplot2::labs(
      title = virus_name,
      subtitle = paste0(
        "Threshold ", threshold, "% | primary method: ",
        config$analysis$reporting_method,
        " | displayed NT50: ", branch$nt50_display
      ),
      x = "Reciprocal dilution",
      y = "Percent neutralization",
      caption = paste(
        c(
          if (branch$point_display %in% c("raw", "both")) {
            "Small points: individual wells"
          } else {
            NULL
          },
          if (branch$point_display %in% c("summary", "both")) {
            "Open circles/error bars: mean +/- SD"
          } else {
            NULL
          },
          if (isTRUE(branch$show_ll4_curve) && nrow(data$ll4) > 0L) {
            "Solid blue: R LL4 fit"
          } else {
            NULL
          },
          if (isTRUE(branch$show_gen5_overlay) && nrow(data$gen5) > 0L) {
            "Orange dot-dash: Gen5 fit"
          } else {
            NULL
          }
        ),
        collapse = "; "
      )
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      plot.caption = ggplot2::element_text(hjust = 0)
    )

  add_branch_y_axis(plot, branch)
}

empty_plot_manifest <- function() {
  data.frame(
    branch_name = character(),
    plot_layout = character(),
    plot_format = character(),
    plot_scope = character(),
    curve_id = character(),
    virus = character(),
    plot_file = character(),
    plot_width_in = numeric(),
    plot_height_in = numeric(),
    save_status = character(),
    fallback_used = logical(),
    fallback_reason = character(),
    exists = logical(),
    size_bytes = numeric(),
    stringsAsFactors = FALSE
  )
}

save_plot_attempt <- function(
  path,
  build_plot,
  branch,
  width,
  height,
  limitsize = TRUE
) {
  first_error <- NULL
  success <- tryCatch(
    {
      plot <- build_plot(branch)
      write_plot_file(
        path = path,
        plot = plot,
        branch = branch,
        width = width,
        height = height,
        limitsize = limitsize
      )
      TRUE
    },
    error = function(e) {
      first_error <<- conditionMessage(e)
      FALSE
    }
  )
  if (success) {
    return(list(status = "saved", fallback_used = FALSE, reason = ""))
  }

  # Remove any incomplete file left by a failed graphics device before either
  # returning or trying the local fallback. This cleanup is path-local and cannot
  # affect output from another curve or branch.
  if (file.exists(path)) {
    unlink(path, force = TRUE)
  }

  if (!isTRUE(branch$fallback_on_plot_error)) {
    return(list(status = "failed", fallback_used = FALSE, reason = first_error %||% ""))
  }

  # Local plot-only fallback. The original branch list and every subsequent plot
  # keep their configured settings.
  fallback_branch <- branch
  fallback_branch$show_model_confidence_band <- FALSE
  fallback_branch$show_gen5_overlay <- FALSE
  fallback_branch$show_nt50_ci <- FALSE
  fallback_error <- NULL
  success <- tryCatch(
    {
      if (file.exists(path)) unlink(path, force = TRUE)
      plot <- build_plot(fallback_branch)
      write_plot_file(
        path = path,
        plot = plot,
        branch = fallback_branch,
        width = width,
        height = height,
        limitsize = limitsize
      )
      TRUE
    },
    error = function(e) {
      fallback_error <<- conditionMessage(e)
      FALSE
    }
  )

  if (success) {
    list(
      status = "saved",
      fallback_used = TRUE,
      reason = first_error %||% ""
    )
  } else {
    if (file.exists(path)) {
      unlink(path, force = TRUE)
    }
    list(
      status = "failed",
      fallback_used = TRUE,
      reason = paste(
        c(first_error, fallback_error)[nzchar(c(first_error, fallback_error))],
        collapse = " | fallback also failed: "
      )
    )
  }
}

save_plots <- function(analysis_objects, output_dir, config) {
  enabled_names <- enabled_plot_branch_names(config)
  if (length(enabled_names) == 0L) {
    return(empty_plot_manifest())
  }

  plots_dir <- file.path(output_dir, "plots")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  manifest_rows <- list()
  row_index <- 0L
  generated_paths <- character()

  for (branch_name in enabled_names) {
    branch <- config$plot$branches[[branch_name]]
    branch_dir <- file.path(plots_dir, branch_name)
    if (isTRUE(config$plot$clear_enabled_branch_directories) && dir.exists(branch_dir)) {
      unlink(branch_dir, recursive = TRUE, force = TRUE)
    }
    dir.create(branch_dir, recursive = TRUE, showWarnings = FALSE)

    if (branch$layout == "individual") {
      for (i in seq_along(analysis_objects)) {
        object <- analysis_objects[[i]]
        result <- object$result[1L, , drop = FALSE]
        path <- build_plot_output_path(
          directory = branch_dir,
          index = i,
          stem = safe_filename(
            result_value(result, "plate_number", ""),
            result_value(result, "plate_id", ""),
            result_value(result, "virus", ""),
            result_value(result, "sample_id", "")
          ),
          format = branch$format
        )
        if (path %in% generated_paths) {
          stop("Duplicate plot output path generated: ", path, call. = FALSE)
        }
        generated_paths <- c(generated_paths, path)

        flag_text <- format_qc_flags_for_plot(
          result_value(result, "all_qc_flags", ""),
          width = config$plot$flag_wrap_width,
          max_flags = Inf,
          bullet = TRUE,
          blank_lines_between_flags = config$plot$blank_lines_between_flags
        )
        flag_lines <- if (nzchar(flag_text)) {
          length(strsplit(flag_text, "\n", fixed = TRUE)[[1]])
        } else {
          0L
        }
        individual_height <- min(
          config$plot$max_individual_plot_height,
          branch$height + config$plot$individual_height_per_flag_line * flag_lines
        )

        attempt <- save_plot_attempt(
          path = path,
          build_plot = function(local_branch) {
            build_individual_plot(object, config, local_branch)
          },
          branch = branch,
          width = branch$width,
          height = individual_height,
          limitsize = TRUE
        )
        if (attempt$status == "failed") {
          warning(
            "Plot failed for branch ", branch_name, " and curve ",
            result_value(result, "curve_id", "unknown"), ": ", attempt$reason
          )
        }

        row_index <- row_index + 1L
        manifest_rows[[row_index]] <- data.frame(
          branch_name = branch_name,
          plot_layout = branch$layout,
          plot_format = branch$format,
          plot_scope = "curve",
          curve_id = as.character(result_value(result, "curve_id", "")),
          virus = as.character(result_value(result, "virus", "")),
          plot_file = path,
          plot_width_in = branch$width,
          plot_height_in = individual_height,
          save_status = attempt$status,
          fallback_used = attempt$fallback_used,
          fallback_reason = attempt$reason,
          exists = file.exists(path),
          size_bytes = if (file.exists(path)) as.numeric(file.info(path)$size[[1]]) else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    } else if (branch$layout == "faceted") {
      viruses <- sort(unique(vapply(
        analysis_objects,
        function(object) as.character(result_value(object$result, "virus", "")),
        character(1)
      )))

      for (i in seq_along(viruses)) {
        virus_name <- viruses[[i]]
        virus_objects <- analysis_objects[vapply(
          analysis_objects,
          function(object) identical(
            as.character(result_value(object$result, "virus", "")),
            virus_name
          ),
          logical(1)
        )]
        n_facets <- length(virus_objects)
        n_columns <- min(branch$facet_columns, n_facets)
        n_rows <- ceiling(n_facets / n_columns)
        width <- max(branch$width, 4.5 * n_columns)
        height <- max(branch$height, 3.6 * n_rows)
        path <- build_plot_output_path(
          directory = branch_dir,
          index = i,
          stem = safe_filename(virus_name, "faceted"),
          format = branch$format
        )
        if (path %in% generated_paths) {
          stop("Duplicate plot output path generated: ", path, call. = FALSE)
        }
        generated_paths <- c(generated_paths, path)

        attempt <- save_plot_attempt(
          path = path,
          build_plot = function(local_branch) {
            build_faceted_plot(virus_objects, virus_name, config, local_branch)
          },
          branch = branch,
          width = width,
          height = height,
          limitsize = FALSE
        )
        if (attempt$status == "failed") {
          warning(
            "Faceted plot failed for branch ", branch_name, " and virus ",
            virus_name, ": ", attempt$reason
          )
        }

        row_index <- row_index + 1L
        manifest_rows[[row_index]] <- data.frame(
          branch_name = branch_name,
          plot_layout = branch$layout,
          plot_format = branch$format,
          plot_scope = "virus",
          curve_id = "",
          virus = virus_name,
          plot_file = path,
          plot_width_in = width,
          plot_height_in = height,
          save_status = attempt$status,
          fallback_used = attempt$fallback_used,
          fallback_reason = attempt$reason,
          exists = file.exists(path),
          size_bytes = if (file.exists(path)) as.numeric(file.info(path)$size[[1]]) else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(manifest_rows) == 0L) {
    empty_plot_manifest()
  } else {
    dplyr::bind_rows(manifest_rows)
  }
}

# Keep one comprehensive curve-level result table and order it for manual review.
# The first columns mirror the notebook-era workflow: pair identifiers, then all
# NT50 values, then reportability/status and flags. Normalization details, fit
# coefficients, warnings, QC settings, and version metadata follow afterward.
# Permanently duplicated aliases remain excluded.
ANALYSIS_RESULTS_OUTPUT_FIELDS <- c(
  # Pair identifiers.
  "plate_number", "sample_id", "virus", "plate_id", "curve_id",

  # NT50 values, intervals, and direct method comparisons.
  "reported_nt50", "reported_ci_lower", "reported_ci_upper",
  "reported_ci_method", "ll4_calculated_nt50", "ll4_nt50_ci_lower",
  "ll4_nt50_ci_upper", "ll4_nt50_se", "uniroot_calculated_nt50",
  "gen5_calculated_nt50", "ll4_vs_uniroot_fold_difference",
  "ll4_vs_uniroot_log2_ratio", "ll4_vs_uniroot_abs_log2_difference",
  "ll4_vs_gen5_fold_difference",

  # Actionable report status and flags immediately follow the NT50 fields.
  "report_status", "reportable", "hard_qc_flags", "review_qc_flags",
  "normalization_hard_qc_flags", "normalization_review_qc_flags",
  "shared_hard_qc_flags", "shared_review_qc_flags", "ll4_hard_qc_flags",
  "ll4_review_qc_flags", "uniroot_hard_qc_flags",
  "uniroot_review_qc_flags", "cross_method_review_qc_flags",
  "gen5_review_qc_flags",

  # Remaining reportability and method context.
  "primary_reporting_method", "threshold_pct", "method_availability_status",
  "both_methods_reportable", "any_method_reportable",
  "alternative_method_reportable", "ll4_report_status", "ll4_reportable",
  "ll4_nt50_status", "ll4_nt50_ci_status", "uniroot_report_status",
  "uniroot_reportable", "uniroot_status", "gen5_nt50_status",
  "normalization_status", "normalization_reportable",
  "threshold_bracketed_by_observed_means",

  # Tested range and observed-curve descriptors.
  "tested_min_dilution", "tested_max_dilution", "n_unique_dilutions",
  "min_replicates_per_dilution", "max_replicates_per_dilution",
  "min_observed_mean_percent", "max_observed_mean_percent",
  "mean_curve_spearman_rho", "n_upward_mean_steps",
  "max_replicate_sd_pp",

  # Normalization details.
  "normalization_source", "normalization_available",
  "normalization_reference_is_plate_matched",
  "plate_matched_virus_control_available",
  "plate_matched_vc_required_metadata_inference",
  "normalization_reference_plate_numbers",
  "normalization_reference_plate_ids", "n_normalization_reference_plates",
  "n_virus_control", "minimum_required_virus_controls",
  "virus_control_replicate_requirement_met", "virus_control_mean",
  "virus_control_sd", "virus_control_cv_pct",
  "virus_control_mean_in_working_range",
  "all_virus_control_wells_in_working_range",
  "n_virus_control_outside_working_range", "virus_control_in_working_range",
  "cell_control_mean", "cell_control_max", "cell_control_background_high",
  "control_virus_sources",

  # R LL4 fit details, coefficients, uncertainty, warnings, and errors.
  "ll4_fit_available", "ll4_estimate_available", "ll4_fit_strategy",
  "ll4_forced_hill", "ll4_initial_hill", "ll4_hill",
  "ll4_lower_asymptote", "ll4_upper_asymptote",
  "ll4_midpoint_parameter", "ll4_rmse", "ll4_pseudo_r2",
  "ll4_residual_df", "ll4_lack_of_fit_p",
  "ll4_nt50_estimate_in_range", "ll4_nt50_ci_available",
  "ll4_nt50_ci_in_range", "ll4_fit_warning_count",
  "ll4_fit_review_warning_count", "ll4_fit_warnings",
  "ll4_fit_review_warnings", "ll4_fit_error", "ll4_nt50_warning_count",
  "ll4_nt50_review_warning_count", "ll4_nt50_warnings",
  "ll4_nt50_review_warnings", "ll4_nt50_error",

  # Observed-value uniroot details.
  "uniroot_estimate_available", "uniroot_estimate_in_range",
  "uniroot_method", "uniroot_interpolation_scale",
  "uniroot_n_tested_dilutions", "uniroot_n_crossings",
  "uniroot_warning_count", "uniroot_warnings", "uniroot_error",

  # Optional Gen5 fit coefficients and diagnostics.
  "gen5_fit_available", "gen5_fit_valid", "gen5_a", "gen5_b", "gen5_c",
  "gen5_d", "gen5_r2", "gen5_fit_f_prob",

  # Prespecified QC settings and reproducibility metadata.
  "qc_vc_count_lower", "qc_vc_count_upper", "qc_hill_lower",
  "qc_hill_upper", "qc_fixed_hill", "analysis_script_version"
)

select_output_fields <- function(data, fields, output_name) {
  missing <- setdiff(fields, names(data))
  if (length(missing) > 0L) {
    stop(
      output_name,
      " could not be constructed because field(s) are missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  data[, fields, drop = FALSE]
}

file_manifest_row <- function(role, path) {
  path_text <- path %||% ""
  exists <- nzchar(path_text) && file.exists(path_text)
  info <- if (exists) file.info(path_text) else NULL

  data.frame(
    role = role,
    path = if (exists) {
      normalizePath(path_text, winslash = "/", mustWork = FALSE)
    } else {
      path_text
    },
    exists = exists,
    size_bytes = if (exists) as.numeric(info$size[[1]]) else NA_real_,
    modified_utc = if (exists) {
      format(info$mtime[[1]], tz = "UTC", usetz = TRUE)
    } else {
      ""
    },
    md5 = if (exists) unname(tools::md5sum(path_text)[[1]]) else "",
    stringsAsFactors = FALSE
  )
}

write_run_metadata <- function(config, output_dir, input_paths, generated_at_utc) {
  config_path <- file.path(output_dir, "analysis_config.txt")
  config_lines <- capture.output(dput(config))
  writeLines(config_lines, config_path)

  session_path <- file.path(output_dir, "session_info.txt")
  writeLines(capture.output(utils::sessionInfo()), session_path)

  manifest <- dplyr::bind_rows(Map(file_manifest_row, names(input_paths), input_paths))
  manifest$analysis_script_version <- SCRIPT_VERSION
  manifest$manifest_generated_utc <- generated_at_utc
  validate_output_table(manifest, "input_manifest.csv")
  readr::write_csv(manifest, file.path(output_dir, "input_manifest.csv"), na = "")
}

relative_output_path <- function(path, output_dir) {
  normalized_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  normalized_output <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  prefix <- paste0(normalized_output, "/")

  if (startsWith(normalized_path, prefix)) {
    substring(normalized_path, nchar(prefix) + 1L)
  } else {
    normalized_path
  }
}

write_plot_manifest <- function(plot_manifest, output_dir, config) {
  if (nrow(plot_manifest) == 0L) {
    manifest <- empty_plot_manifest()
  } else {
    manifest <- plot_manifest
    manifest$plot_file <- vapply(
      manifest$plot_file,
      relative_output_path,
      character(1),
      output_dir = output_dir
    )
  }

  validate_output_table(manifest, "plot_manifest.csv")
  readr::write_csv(manifest, file.path(output_dir, "plot_manifest.csv"), na = "")
}

write_run_summary <- function(results, plot_manifest, output_dir, config, generated_at_utc) {
  enabled_names <- enabled_plot_branch_names(config)
  n_saved <- if (nrow(plot_manifest) == 0L) 0L else {
    sum(plot_manifest$save_status == "saved", na.rm = TRUE)
  }
  n_failed <- if (nrow(plot_manifest) == 0L) 0L else {
    sum(plot_manifest$save_status == "failed", na.rm = TRUE)
  }
  n_plot_fallback <- if (nrow(plot_manifest) == 0L) 0L else {
    sum(plot_manifest$fallback_used %in% TRUE, na.rm = TRUE)
  }

  summary <- data.frame(
    metric = c(
      "analysis_script_version",
      "generated_utc",
      "configuration_source",
      "use_hardcoded_settings",
      "cli_overrides",
      "working_directory",
      "input_folder",
      "reporting_method_configured",
      "uniroot_interpolation_scale",
      "infer_control_virus_from_well_row",
      "allow_cross_plate_virus_controls",
      "allow_cross_plate_vc_for_reporting",
      "threshold_pct",
      "enabled_plot_branches",
      "n_curves",
      "n_primary_reportable",
      "n_primary_reportable_review_required",
      "n_primary_not_reportable",
      "n_ll4_estimates_available",
      "n_uniroot_estimates_available",
      "n_ll4_reportable",
      "n_uniroot_reportable",
      "n_both_methods_reportable",
      "n_plate_matched_normalizations",
      "n_non_plate_matched_normalizations",
      "n_row_inferred_normalizations",
      "n_plot_attempts",
      "n_plot_files_saved",
      "n_plot_failures",
      "n_plot_local_fallbacks"
    ),
    value = as.character(c(
      SCRIPT_VERSION,
      generated_at_utc,
      config$runtime$configuration_source,
      config$runtime$use_hardcoded_settings,
      config$runtime$cli_overrides,
      config$runtime$working_directory,
      config$paths$input_folder,
      config$analysis$reporting_method,
      config$analysis$uniroot_interpolation_scale,
      config$normalization$infer_control_virus_from_well_row,
      config$normalization$allow_cross_plate_virus_controls,
      config$normalization$allow_cross_plate_vc_for_reporting,
      config$analysis$threshold_pct,
      paste(enabled_names, collapse = ","),
      nrow(results),
      sum(results$reportable %in% TRUE, na.rm = TRUE),
      sum(results$report_status == "REPORTABLE_REVIEW_REQUIRED", na.rm = TRUE),
      sum(results$report_status == "NOT_REPORTABLE", na.rm = TRUE),
      sum(is.finite(results$ll4_calculated_nt50), na.rm = TRUE),
      sum(is.finite(results$uniroot_calculated_nt50), na.rm = TRUE),
      sum(results$ll4_reportable %in% TRUE, na.rm = TRUE),
      sum(results$uniroot_reportable %in% TRUE, na.rm = TRUE),
      sum(results$both_methods_reportable %in% TRUE, na.rm = TRUE),
      sum(results$normalization_reference_is_plate_matched %in% TRUE, na.rm = TRUE),
      sum(results$normalization_available %in% TRUE &
        !(results$normalization_reference_is_plate_matched %in% TRUE), na.rm = TRUE),
      sum(results$normalization_source == "plate_matched_row_inferred", na.rm = TRUE),
      nrow(plot_manifest),
      n_saved,
      n_failed,
      n_plot_fallback
    )),
    stringsAsFactors = FALSE
  )

  validate_output_table(summary, "run_summary.csv")
  readr::write_csv(summary, file.path(output_dir, "run_summary.csv"), na = "")
}

run_analysis <- function(config = NULL) {
  if (is.null(config)) {
    config <- build_run_config()
  }

  check_r_version()
  check_required_packages()
  config <- validate_config(config)
  config <- prepare_runtime_paths(config)
  generated_at_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)

  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  message("DirectView-FRA analysis v", SCRIPT_VERSION)
  message("Configuration source: ", config$runtime$configuration_source)
  if (nzchar(config$runtime$cli_overrides)) {
    message("CLI overrides: ", config$runtime$cli_overrides)
  }
  message("Working directory: ", config$runtime$working_directory)
  message("Input folder: ", config$paths$input_folder)
  message("Reading raw counts: ", config$input$raw_counts_file)
  raw_data <- standardize_raw_counts(config$input$raw_counts_file)

  lookup_data <- standardize_lookup(config$input$lookup_file)
  raw_data <- merge_lookup_metadata(raw_data, lookup_data)
  raw_data <- classify_controls(raw_data)
  raw_data <- validate_raw_counts(raw_data)
  raw_data <- assign_control_reference_viruses(raw_data, config)
  raw_data <- derive_replicate_index(raw_data)

  gen5_data <- standardize_gen5_fits(config$input$gen5_fits_file)

  normalization_qc <- compute_normalization_qc(raw_data, config)
  normalized_data <- normalize_to_controls(raw_data, normalization_qc)
  dilution_summary <- summarize_dilutions(normalized_data)

  analysis <- analyze_all_curves(
    normalized_data,
    dilution_summary,
    normalization_qc,
    gen5_data,
    config
  )

  analysis_results_output <- select_output_fields(
    analysis$results,
    ANALYSIS_RESULTS_OUTPUT_FIELDS,
    "analysis_results.csv"
  )

  validate_output_table(
    analysis_results_output,
    "analysis_results.csv",
    required_fields = ANALYSIS_RESULT_REQUIRED_FIELDS
  )
  validate_output_table(normalized_data, "normalized_well_data.csv")
  validate_output_table(dilution_summary, "dilution_summary.csv")
  validate_output_table(normalization_qc, "normalization_qc.csv")

  readr::write_csv(
    analysis_results_output,
    file.path(config$output_dir, "analysis_results.csv"),
    na = ""
  )
  readr::write_csv(
    normalized_data,
    file.path(config$output_dir, "normalized_well_data.csv"),
    na = ""
  )
  readr::write_csv(
    dilution_summary,
    file.path(config$output_dir, "dilution_summary.csv"),
    na = ""
  )
  readr::write_csv(
    normalization_qc,
    file.path(config$output_dir, "normalization_qc.csv"),
    na = ""
  )

  plot_manifest <- save_plots(analysis$objects, config$output_dir, config)
  write_plot_manifest(plot_manifest, config$output_dir, config)

  write_run_metadata(
    config,
    config$output_dir,
    list(
      raw_counts = config$input$raw_counts_file,
      lookup = config$input$lookup_file,
      gen5_fits = config$input$gen5_fits_file
    ),
    generated_at_utc
  )
  write_run_summary(
    analysis$results,
    plot_manifest,
    config$output_dir,
    config,
    generated_at_utc
  )

  n_reportable_clean <- sum(
    analysis$results$report_status == "REPORTABLE",
    na.rm = TRUE
  )
  n_review <- sum(
    analysis$results$report_status == "REPORTABLE_REVIEW_REQUIRED",
    na.rm = TRUE
  )
  n_not_reportable <- sum(
    analysis$results$report_status == "NOT_REPORTABLE",
    na.rm = TRUE
  )

  message("Analysis complete: ", nrow(analysis$results), " curve(s).")
  message("  Reportable without review flags: ", n_reportable_clean)
  message("  Reportable, review required: ", n_review)
  message("  Not reportable: ", n_not_reportable)
  message(
    "  Plot files written: ",
    if (nrow(plot_manifest) == 0L) 0L else {
      sum(plot_manifest$save_status == "saved")
    }
  )
  message(
    "Results directory: ",
    normalizePath(config$output_dir, winslash = "/", mustWork = FALSE)
  )

  invisible(list(
    config = config,
    results = analysis_results_output,
    full_internal_results = analysis$results,
    normalized_well_data = normalized_data,
    dilution_summary = dilution_summary,
    normalization_qc = normalization_qc,
    analysis_objects = analysis$objects,
    plot_manifest = plot_manifest
  ))
}

run_directview_fra <- function(args = character()) {
  cli <- parse_cli_args(args)
  config <- build_run_config(cli)
  run_analysis(config)
}

if (sys.nframe() == 0L) {
  run_directview_fra(commandArgs(trailingOnly = TRUE))
}

