# Start here

This guide runs the included DirectViewFRA example with the fewest required edits.

## 1. Install R and the required packages

Use R 4.1.0 or newer. In R or RStudio, run:

```r
install.packages(c("drc", "dplyr", "ggplot2", "readr", "scales", "tibble"))
```

## 2. Keep the repository structure intact

The default example expects this layout:

```text
DirectViewFRA/
  directview_fra_analysis.R
  example-data/
    mixed_quality_example_raw_counts.csv
    mixed_quality_example_ll4_fit_curve.csv
```

The external fit CSV is included for optional comparison but is not read by default.

## 3. Run from RStudio or Jupyter

Open or change to the repository folder, then run:

```r
source("directview_fra_analysis.R")
results <- run_directview_fra()
```

The script uses the editable `HARDCODED_SETTINGS` block near the top. The default input folder is `example-data`, the required raw-count file is selected, and the output folder is:

```text
example-data/directview_fra_results/
```

## 4. Run from a terminal

Windows PowerShell:

```powershell
Rscript directview_fra_analysis.R --use-hardcoded-settings false --input-folder "example-data" --raw-counts "mixed_quality_example_raw_counts.csv" --output-dir "../directview_fra_results"
```

Linux or macOS shell:

```sh
Rscript directview_fra_analysis.R \
  --use-hardcoded-settings false \
  --input-folder "example-data" \
  --raw-counts "mixed_quality_example_raw_counts.csv" \
  --output-dir "../directview_fra_results"
```

The command-line route starts from conservative built-in defaults. Optional inputs and cross-plate controls are not enabled unless requested.

## 5. Review these outputs first

Open the output directory and review:

1. `analysis_results.csv` — one row per plate, virus, and sample curve.
2. `normalization_qc.csv` — the control source and denominator used for each plate-virus combination.
3. `normalized_well_data.csv` — every well with its calculated percent neutralization when available.
4. `plots/individual_review_png/` — individual curve review plots.
5. `run_summary.csv` — a compact description of the run configuration and result counts.

A calculated NT50 is not necessarily reportable. Use the `reportable`, `report_status`, `hard_qc_flags`, and `review_qc_flags` fields together.

## 6. Analyze another dataset

Edit these values near the top of `directview_fra_analysis.R`:

```r
USE_HARDCODED_SETTINGS <- TRUE

# Change these entries inside the existing HARDCODED_SETTINGS list:
input_folder = "C:/path/to/your/input-folder"
output_dir = "directview_fra_results"
raw_counts_file = "Raw Counts.csv"
lookup_file = NULL
gen5_fits_file = NULL
```

On Windows, forward slashes are recommended in R strings. A Linux path can be supplied directly, for example `/home/user/run-01`.

See [Detailed usage](docs/USAGE.md), [Input files](docs/INPUT_FORMAT.md), and [Output files](docs/OUTPUTS.md) for the complete reference.
