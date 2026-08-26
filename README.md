# DirectViewFRA

DirectViewFRA is an R workflow for analysis of per-well focus counts from the DirectView focus reduction assay for influenza neutralizing antibodies. The workflow performs virus-control normalization, four-parameter log-logistic fitting, NT50 estimation, quality-control review, and configurable plot and CSV export.

**Software version:** 1.0.0

## Responsible team

- **Main point of contact:** Bin Zhou (Lead and Technical Contact), Influenza Division, National Center for Immunization and Respiratory Diseases, CDC — bzhou@cdc.gov
- **Code/programming contact:** Michael Currier, Influenza Division, National Center for Immunization and Respiratory Diseases, CDC — ovn5@cdc.gov


## What the workflow does

- Accepts a raw per-well focus-count CSV, with an optional plate-layout lookup CSV.
- Supports multiple viruses on a physical plate and matches virus controls by virus-specific metadata.
- Uses plate- and virus-matched virus controls whenever they are available.
- The included hardcoded example enables same-virus cross-plate controls for diagnostic normalization and plotting when a plate-matched control is absent; cross-plate results remain non-reportable unless separate reporting permission is enabled.
- Calculates an LL4 NT50 and an independent observed-value uniroot NT50 when the data support each method.
- Does not promote an extrapolated NT50 when the observed dilution series does not bracket the selected threshold.
- Fits an unconstrained LL4 first, then refits only the affected curve with Hill fixed at 1.0 when the initial Hill estimate is outside 0.5-2.0.
- Writes one comprehensive curve-level result CSV plus per-well, dilution-summary, normalization-QC, run-summary, and plot-manifest files.
- Runs one or more explicitly enabled plot branches, each in its own output folder.

## Quick start

The repository is arranged so the default settings run the included example:

```text
DirectViewFRA/
  directview_fra_analysis.R
  example-data/
    mixed_quality_example_raw_counts.csv
    mixed_quality_example_ll4_fit_curve.csv
```

Install the required R packages:

```r
install.packages(c("drc", "dplyr", "ggplot2", "readr", "scales", "tibble"))
```

### Jupyter or RStudio

```r
source("directview_fra_analysis.R")
results <- run_directview_fra()
```

The settings block near the top of the script is enabled by default:

```r
USE_HARDCODED_SETTINGS <- TRUE
```

Edit `HARDCODED_SETTINGS` to select files, normalization rules, analysis settings, optional QC thresholds, and plot branches.

### Command line

```sh
Rscript directview_fra_analysis.R \
  --use-hardcoded-settings false \
  --input-folder /path/to/input \
  --raw-counts "Raw Counts.csv" \
  --output-dir results
```

Explicit command-line values override the selected starting configuration. Run the following for all recognized command-line options:

```sh
Rscript directview_fra_analysis.R --help
```

## Important analysis rules

The primary reporting method is selected once for the run. Both calculated method values remain available in `analysis_results.csv`, but the workflow does not silently switch the reporting method for individual curves.

Plate- and virus-matched controls are preferred. In the included hardcoded example, `allow_cross_plate_virus_controls = TRUE` keeps curves visible for diagnostic review when a matched control is absent. The separate `allow_cross_plate_vc_for_reporting` setting remains `FALSE`, so those calculations do not pass the reporting gate. The built-in command-line defaults remain plate-only unless cross-plate use is explicitly enabled.

The observed-value uniroot method evaluates the complete ordered dilution-mean series for threshold crossings. It reports a value only when there is one unique in-range crossing and performs no interval extension beyond the two adjacent observed dilution means containing that crossing.

Review the individual well values, curve plots, normalization source, and quality-control flags before reporting an NT50.

## Plot controls

Each plot branch has independent settings. Plot files are written with an explicit `.png` or `.pdf` filename extension and graphics device. The default `y_axis = "data"` displays the full plotted range; set `y_axis = "fixed"` to apply `y_limits` instead.

The top-right LL4/uniroot detail box is controlled by:

```r
show_comparison_details = TRUE
```

Set it to `FALSE` to omit the box entirely while retaining any other requested plot elements, such as fitted curves, NT50 lines, confidence intervals, points, and captions. From the command line, `--show-comparison-details false` applies the setting to all configured branches.

## Documentation

- [Start here](START_HERE.md)
- [Detailed usage](docs/USAGE.md)
- [Input files](docs/INPUT_FORMAT.md)
- [Output files](docs/OUTPUTS.md)
- [Example data](example-data/README.md)
- [Contributing](CONTRIBUTING.md)
- [Disclaimer](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md)
- [Code of conduct](https://github.com/CDCgov/template/blob/main/code-of-conduct.md)

## Citation

When using this workflow, cite the software repository and the DirectView-FRA method publication:

Feng C, Rowe T, Currier M, Dong R, Huang Y, Atteberry G, Wang L, Hatta M, Davis CT, Wentworth DE, Zhou B. *A DirectView focus reduction assay for high-throughput quantification of neutralizing antibodies against influenza A and B viruses.* Cell Reports Methods. 2026;101479. doi:10.1016/j.crmeth.2026.101479.

Citation metadata are also provided in [`CITATION.cff`](CITATION.cff).

## Public domain standard notice

This repository constitutes a work of the United States Government and is not subject to domestic copyright protection under 17 USC § 105. This repository is in the public domain within the United States, and copyright and related rights in the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/). All contributions to this repository will be released under the CC0 dedication. By submitting a pull request, you agree to comply with this waiver of copyright interest.

## License standard notice

The source code is licensed under the Apache Software License, Version 2.0 or later. See [`LICENSE`](LICENSE).

The source code is distributed in the hope that it will be useful, but without warranty of any kind, including any implied warranty of merchantability or fitness for a particular purpose.

## Privacy standard notice

This repository contains only non-sensitive, publicly available data and information. Do not submit or store personally identifiable information, protected health information, sensitive laboratory information, or other non-public data. See the [Disclaimer](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md) and [Code of Conduct](https://github.com/CDCgov/template/blob/main/code-of-conduct.md).

## Contributing standard notice

Anyone is encouraged to contribute by [forking](https://help.github.com/articles/fork-a-repo) and submitting a pull request. By contributing, you grant a world-wide, royalty-free, perpetual, irrevocable, non-exclusive, transferable license to all users under the terms of the [Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or later. All submissions received through CDC GitHub may be subject to applicable federal law, including the Federal Records Act, and may be archived. See [`CONTRIBUTING`](CONTRIBUTING.md).

## Records management standard notice

This repository is not a source of government records, but is a copy intended to increase collaboration and collaborative potential. Government records are maintained through applicable CDC records-management processes.

## Additional standard notices

Please refer to [CDC's Template Repository](https://github.com/CDCgov/template) for more information about contributing, public domain notices and disclaimers, and the code of conduct.
