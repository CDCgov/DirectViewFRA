# Detailed usage

DirectViewFRA can be run from an interactive R environment or from the command line. Both routes use the same validation, normalization, model fitting, reporting, and plotting functions.

## Interactive mode: hardcoded settings

The default is:

```r
USE_HARDCODED_SETTINGS <- TRUE
```

Edit the single `HARDCODED_SETTINGS` list near the top of `directview_fra_analysis.R`, then source and run the script:

```r
source("directview_fra_analysis.R")
results <- run_directview_fra()
```

The returned object contains the effective configuration, curve-level results, normalized well data, dilution summaries, normalization QC, plot objects, and plot manifest.

### Paths

```r
paths = list(
  set_working_directory_to_script_folder = TRUE,
  input_folder = "example-data"
),
output_dir = "directview_fra_results"
```

When `set_working_directory_to_script_folder = TRUE`, the script attempts to use the folder containing the R file as the working directory. A relative `input_folder` is resolved from that location. A relative `output_dir` is resolved from the input folder.  This was added in testing as most errors derived from users incorrectly setting paths.

Windows examples:

```r
input_folder = "C:/DirectViewFRA/Run_01"
input_folder = "C:\\DirectViewFRA\\Run_01"
```

Linux example:

```r
input_folder = "/home/user/DirectViewFRA/Run_01"
```

Forward slashes are recommended on Windows because they do not require escaping.

### Input files

```r
input = list(
  raw_counts_file = "mixed_quality_example_raw_counts.csv",
  lookup_file = NULL,
  gen5_fits_file = NULL
)
```

- `raw_counts_file` is required.
- `lookup_file` is optional and supplies plate-layout metadata when the raw CSV is not fully annotated.
- `gen5_fits_file` is optional and adds external four-parameter fit comparison fields. It is not required for the R LL4 or observed-value uniroot calculations.

- Detailed information on the expected fields for input can be found in [INPUT_FORMAT](INPUT_FORMAT.md) documentation.


### Normalization settings

```r
normalization = list(
  infer_control_virus_from_well_row = FALSE,
  allow_cross_plate_virus_controls = FALSE,
  cross_plate_vc_plate_numbers = NULL,
  allow_cross_plate_vc_for_reporting = FALSE,
  pool_equal_weight_by_plate = TRUE
)
```

`infer_control_virus_from_well_row` permits row-based control-virus assignment only for a validated row-block plate layout. Explicit metadata are otherwise preferred.

`allow_cross_plate_virus_controls` permits same-virus controls from other plates when no plate-matched control is available. This is disabled by default.

`cross_plate_vc_plate_numbers` optionally limits the eligible source plates. Values must match the `plate_number` strings in the input exactly:

```r
cross_plate_vc_plate_numbers = c("Plate 10", "Plate 11", "Plate 12")
```

`allow_cross_plate_vc_for_reporting` is a separate reporting permission. A cross-plate denominator can be used for diagnostic normalization while remaining non-reportable.

When `pool_equal_weight_by_plate = TRUE`, each contributing plate mean receives equal weight. When `FALSE`, all eligible VC wells are pooled directly.

### Analysis settings

```r
analysis = list(
  threshold_pct = 50,
  reporting_method = "ll4",
  uniroot_interpolation_scale = "linear",
  ll4_ci_method = "delta",
  confidence_level = 0.95,
  min_unique_dilutions = 4L,
  min_total_observations = 8L,
  prediction_grid_points = 400L
)
```

- `reporting_method` is `"ll4"` or `"uniroot"`. Both calculated estimates remain in the output when available.
- `uniroot_interpolation_scale` is `"linear"` or `"log2"` and applies consistently to the full run.
- `ll4_ci_method` is `"none"`, `"delta"`, or `"inv"`.
- `prediction_grid_points` controls fitted-curve plotting resolution.

The observed-value uniroot calculation checks the complete ordered dilution-mean series. It returns an estimate only when there is exactly one threshold crossing and never extrapolates beyond the tested dilutions.

### QC settings

```r
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
)
```

A `NULL` optional threshold means the diagnostic is calculated and exported but is not used to gate reportability. Set a threshold only when a defensible rule has been selected for the analysis plan.

The LL4 model is first fitted without a Hill constraint. A curve-local refit with `fixed_hill` is attempted only when the initial Hill estimate lies outside `hill_range`.

### Plot settings

Each named branch writes to a separate subfolder under `plots/`.

```r
individual_review_png = list(
  enabled = TRUE,
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
)
```

Important choices:

- `point_display`: `"raw"`, `"summary"`, or `"both"`.
- `nt50_display`: `"primary"`, `"ll4"`, `"uniroot"`, `"both"`, or `"none"`.
- `show_comparison_details = FALSE`: removes the top-right method detail box completely.
- `y_axis = "data"`: displays the complete plotted data range.
- `y_axis = "fixed"`: applies `y_limits` with `coord_cartesian()`.
- `fallback_on_plot_error = TRUE`: retries only the failed plot after disabling optional overlays and intervals; it does not change later plots or other branches.

When normalization is unavailable, the individual branch displays raw focus counts rather than an empty percent-neutralization panel.

## Command-line mode

Use `--use-hardcoded-settings false` to ignore the hardcoded settings and start from conservative built-in defaults.

```sh
Rscript directview_fra_analysis.R \
  --use-hardcoded-settings false \
  --input-folder "/home/user/run-01" \
  --raw-counts "Raw Counts.csv" \
  --lookup none \
  --gen5-fits none \
  --output-dir "../results" \
  --reporting-method ll4 \
  --uniroot-scale linear \
  --ll4-ci delta \
  --plot-branches individual_review_png
```

Windows PowerShell can use the same command on one line with forward-slash paths:

```powershell
Rscript directview_fra_analysis.R --use-hardcoded-settings false --input-folder "C:/DirectViewFRA/Run_01" --raw-counts "Raw Counts.csv" --lookup none --gen5-fits none --output-dir "../results" --reporting-method ll4 --uniroot-scale linear --ll4-ci delta --plot-branches individual_review_png
```

### Cross-plate diagnostic example

```sh
Rscript directview_fra_analysis.R \
  --use-hardcoded-settings false \
  --input-folder example-data \
  --raw-counts mixed_quality_example_raw_counts.csv \
  --allow-cross-plate-vc true \
  --cross-plate-vc-plates "Plate 10,Plate 11,Plate 12" \
  --allow-cross-plate-vc-reporting false \
  --output-dir ../cross-plate-review
```

### Optional external fit comparison

```sh
Rscript directview_fra_analysis.R \
  --use-hardcoded-settings false \
  --input-folder example-data \
  --raw-counts mixed_quality_example_raw_counts.csv \
  --gen5-fits mixed_quality_example_ll4_fit_curve.csv \
  --output-dir ../results-with-external-comparison
```

Use `Rscript directview_fra_analysis.R --help` for the current command-line option list.
