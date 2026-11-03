# Example data

This folder contains an anonymized mixed-quality dataset, an optional lookup table and external fit export, and two archives of generated results for demonstrating the DirectViewFRA workflow.

The input CSVs use `Serum` as the sample-identifier column. The script maps this to `sample_id` internally and in outputs. The legacy `Global.Name` header is also accepted, but only one sample-identifier column should be included in each file.

These files demonstrate software behavior and output review. They are not a substitute for local assay validation or acceptance criteria.

## Input files

### `mixed_quality_example_raw_counts.csv`

Required per-well input for the default example. It contains 1,152 rows across 12 plates, including plate, well, virus, serum, dilution, focus-count, and control metadata.

### `mixed_quality_example_lookup.csv`

Optional plate-layout lookup table containing 48 unique mappings: `CC`, `VC`, `Sera1`, and `Sera2` for each of the 12 plates. It includes `Plate ID`, `Plate Number`, `Well ID`, `Virus`, `Serum`, and `Control Type`.

The script joins the lookup to the raw counts using `Plate ID` and `Well ID`. Here, `Well ID` is the local label, such as `Sera1`, `VC`, or `CC`. Physical well positions, such as `A1`, remain in the raw-count file's `Well` column.

The lookup must cover every raw-count mapping, including controls. It can supply missing metadata; conflicting nonblank values between the raw counts and lookup cause the run to stop for review.

### `mixed_quality_example_ll4_fit_curve.csv`

Optional commercial imaging-system software fit export, validated using Agilent Gen5 software. The default analysis does not read this file because `gen5_fits_file = NULL`. When enabled, it supplies external fit coefficients and comparison results alongside the R-derived estimates.

## Running the examples

Keep the CSVs in the configured `input_folder`. Filenames are resolved relative to that folder; absolute file paths are also accepted. With `USE_HARDCODED_SETTINGS <- TRUE`, edit the existing `input` entry inside `HARDCODED_SETTINGS`, leaving the other settings intact.

For the default raw-count-only example:

```r
set_working_directory_to_script_folder = TRUE,
input_folder = "example-data"

  raw_counts_file = "mixed_quality_example_raw_counts.csv",
  lookup_file = NULL,
  gen5_fits_file = NULL
```

For the example using both optional inputs:

```r
set_working_directory_to_script_folder = TRUE,
input_folder = "example-data"

  raw_counts_file = "mixed_quality_example_raw_counts.csv",
  lookup_file = "mixed_quality_example_lookup.csv",
  gen5_fits_file = "mixed_quality_example_ll4_fit_curve.csv"

```

## Expected output archives

### `ExpectedResultsWithNoChange.zip`

Reference outputs from the default example using the raw-count CSV, without a lookup table or external fit input.

### `ExpectedResultsWithLookupAndLL4FitFiles.zip`

Reference outputs using the same raw-count CSV together with the lookup table and external LL4 fit export. This run populates the Gen5 fit and comparison fields in `analysis_results.csv` and records the lookup-supplied control designations as explicit.

The R-derived LL4 and uniroot estimates and reportability classifications agree between the two archived runs. Both contain 24 curves: 19 classified as `REPORTABLE`, four as `REPORTABLE_REVIEW_REQUIRED`, and one as `NOT_REPORTABLE`.

### Archive contents and comparison

Each archive contains a `directview_fra_results/` folder with:

- Curve-level results in `analysis_results.csv`.
- Per-well and dilution-level data in `normalized_well_data.csv` and `dilution_summary.csv`.
- Control and normalization information in `normalization_qc.csv`.
- 24 generated curve-review PNGs in `plots/individual_review_png/`.
- Run documentation in `analysis_config.txt`, `run_summary.csv`, `session_info.txt`, `input_manifest.csv`, and `plot_manifest.csv`.

Extract the archives into separate folders because both use the same top-level output-folder name. Configure `input_folder` for your own system rather than copying the archived machine paths.

The included PNGs are identical between the two runs. Loading the external fit file does not automatically display its curves: `show_gen5_overlay = FALSE` in both archived configurations. To display those curves, enable `show_gen5_overlay = TRUE` in the desired plot branch and rerun with the external fit input selected.

Use these archives as reference outputs. Exact numerical precision and plot appearance may vary with the R environment; consult `session_info.txt` when comparing runs.

## Plate 09 diagnostic normalization example

Plate 09 includes serum `2025-086` tested against `V3X2_H1N2v Example 3` without a plate-matched virus control. Both archived runs use `allow_cross_plate_virus_controls = FALSE`, so no cross-plate denominator is selected. This curve has no usable virus-control reference, is flagged `virus_control_reference_unavailable`, and is not reportable. Its individual review plot displays raw focus counts because normalization is unavailable.

The lookup preserves the different virus assignments on Plate 09. Adding the lookup or external fit file does not supply the missing virus-control reference.

To demonstrate cross-plate diagnostic normalization, set `allow_cross_plate_virus_controls = TRUE`. Same-virus controls from other plates can then be used for diagnostic normalization and plotting, while `allow_cross_plate_vc_for_reporting = FALSE` keeps that denominator from producing a reportable result. This is a separate configuration from either archived example.
