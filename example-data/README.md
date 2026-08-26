# Example data

This folder contains an anonymized mixed-quality dataset for demonstrating the DirectViewFRA workflow.

## Files

### `mixed_quality_example_raw_counts.csv`

Required per-well input for the default example. It includes plate, well, virus, sample, dilution, focus-count, and control metadata.

### `mixed_quality_example_ll4_fit_curve.csv`

Optional commercial imaging-system software fit export.  This was validated on Angilent's Gen5 software. The default analysis does not read this file because `gen5_fits_file = NULL`. To include the comparison, set:

```r
gen5_fits_file = "mixed_quality_example_ll4_fit_curve.csv"
```

The example data are intended for software demonstration and review of output behavior. They are not a substitute for local assay validation or acceptance criteria.

## Plate 09 diagnostic normalization example

The example includes a Plate 09 curve for V3X2_H1N2v Example 3 without a plate-matched virus control. Under the default settings, `allow_cross_plate_virus_controls = FALSE`, so no cross-plate denominator is selected: this curve has no usable virus-control reference, is flagged virus_control_reference_unavailable, and is not reportable. 

To demonstrate cross-plate diagnostic normalization instead, set `allow_cross_plate_virus_controls = TRUE`; same-virus controls from other plates are then used for diagnostic normalization and plotting, while `allow_cross_plate_vc_for_reporting = FALSE` keeps that denominator from producing a reportable result.
