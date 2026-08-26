# Output files

Each run writes its outputs to the configured `output_dir`.

## Main result files

### `analysis_results.csv`

One row per unique `plate_number + plate_id + virus + sample_id` curve. The columns are ordered for manual review:

1. Pair identifiers.
2. NT50 values, intervals, and direct method comparisons.
3. Reportability and QC flags.
4. Method and normalization status.
5. Tested-range and observed-curve descriptors.
6. Normalization details.
7. LL4 parameters, fit diagnostics, warnings, and errors.
8. Observed-value uniroot diagnostics.
9. Optional external fit coefficients and diagnostics.
10. QC settings and script version.

Key fields:

| Field | Meaning |
|---|---|
| `plate_number`, `sample_id`, `virus`, `plate_id`, `curve_id` | Curve identifiers. |
| `reported_nt50` | The selected primary result, populated only when the configured primary method is reportable. |
| `reported_ci_lower`, `reported_ci_upper` | Reported LL4 interval limits when available and inside the tested range. The observed-value uniroot method does not borrow an LL4 interval. |
| `ll4_calculated_nt50` | Calculated absolute 50% crossing from the final R LL4 model when available. |
| `uniroot_calculated_nt50` | Calculated crossing from the observed dilution-level means when a unique in-range crossing exists. |
| `gen5_calculated_nt50` | Optional crossing from the supplied external four-parameter fit. |
| `reportable` | Whether the configured primary method passes the hard reporting gates. |
| `report_status` | `REPORTABLE`, `REPORTABLE_REVIEW_REQUIRED`, or `NOT_REPORTABLE`. |
| `hard_qc_flags` | Conditions that prevent reporting the primary result. |
| `review_qc_flags` | Conditions requiring review but not automatically preventing reporting. |
| `primary_reporting_method` | The run-wide configured method, `ll4` or `uniroot`. |
| `normalization_source` | The selected VC source, such as plate-matched explicit controls or an enabled cross-plate fallback. |
| `normalization_reference_plate_numbers` | Plate labels contributing to the denominator. |
| `virus_control_mean` | Denominator used to calculate percent neutralization. |
| `ll4_fit_strategy` | `unconstrained`, `fixed_hill_refit`, or `failed`. |
| `ll4_hill`, `ll4_lower_asymptote`, `ll4_upper_asymptote`, `ll4_midpoint_parameter` | Final LL4 parameters. |
| `uniroot_method`, `uniroot_interpolation_scale`, `uniroot_n_crossings` | Definition and status of the observed-value interpolation. |
| `analysis_script_version` | Script version that produced the row. |

A calculated value may be retained even when it is not reportable. Use the method-specific status fields and the primary hard/review flags when interpreting a row.

### `normalized_well_data.csv`

One row per input well, with standardized metadata, control classification, normalization source, replicate index, and `percent_neutralization` when a valid denominator was available.

This file is the primary trace-back from a curve result to the individual focus counts.

### `dilution_summary.csv`

One row per curve and tested dilution. It includes:

- replicate count;
- mean and SD focus count;
- mean, SD, and SE percent neutralization;
- minimum and maximum percent neutralization.

The observed-value uniroot calculation uses `dilution` and `mean_percent_neutralization` from this table.

### `normalization_qc.csv`

One row per `plate_number + plate_id + virus` normalization group. It records:

- whether a plate-matched VC was available;
- whether row inference was required;
- whether cross-plate controls were allowed and used;
- source plate numbers and IDs;
- the number, mean, SD, CV, minimum, and maximum of the selected controls;
- working-range checks;
- cell-control summaries; and
- normalization hard and review flags.

## Reproducibility and run files

### `run_summary.csv`

A metric/value table describing the effective run configuration and summary counts, including the number of curves, reportable results, available LL4 and uniroot estimates, normalization types, and plot outcomes.

### `analysis_config.txt`

The complete effective R configuration after hardcoded settings and command-line overrides have been combined and validated.

### `session_info.txt`

The R version, platform, and loaded package versions from `sessionInfo()`.

### `input_manifest.csv`

The raw-count, lookup, and external-fit paths used for the run, with existence status, size, modification time, and MD5 hash.

### `plot_manifest.csv`

One row per attempted plot. Fields include branch, layout, format, curve or virus scope, relative output path, saved width and height, status, local fallback use, file existence, and file size.

All successful plot paths should end with the configured `.png` or `.pdf` extension.

## Plot directories

Each enabled branch writes to:

```text
output_dir/
  plots/
    <branch_name>/
```

The default branch is:

```text
plots/individual_review_png/
```

Individual percent-neutralization plots can include raw wells, dilution means with SD, the LL4 fit, NT50 lines, an LL4 confidence interval, the 50% reference line, a top-right comparison box, and wrapped QC flags. Every component is controlled by the branch settings.

If no percent-neutralization value is available for a curve, the individual review branch displays raw focus counts and identifies normalization as unavailable. It does not apply percent-neutralization y-axis limits to that fallback plot.

## Missing values

CSV files are written with blank cells for R `NA` values. A blank calculated NT50 means the method did not return a finite estimate. A blank reported NT50 means the configured primary result was unavailable or failed a hard reporting rule.
