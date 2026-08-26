# Input files

DirectViewFRA reads comma-separated value files. Column names are normalized to lowercase underscore form, and several common aliases are accepted.

## Required raw-count CSV

The raw-count file contains one row per physical well.

### Required fields

| Canonical field | Accepted examples | Requirement |
|---|---|---|
| `well_id` | `well_id`, `wellid`, `local_name`, `localname` | Well role or local identifier, such as a sample label, `VC`, or `CC`. |
| `well` | `well`, `well_position`, `well_name` | Physical well position, such as `A1`. |
| `dilution` | `conc_dil`, `concentration_dilution`, `dilution`, `dilution_factor`, `reciprocal_dilution` | Positive numeric value for sample wells. |
| `focus_count` | `cell_count`, `focus_count`, `foci_count`, `spot_count` | Finite, nonnegative count. |
| `plate_id` | `plate_id`, `plateid` | Unique plate identifier. |
| `plate_number` | `plate_number`, `plate_no`, `platenumber` | Human-readable plate label. |

### Sample metadata

These fields may be present in the raw-count file or supplied through the optional lookup file:

| Canonical field | Accepted examples |
|---|---|
| `virus` | `virus`, `virus_name` |
| `sample_id` | `global_name`, `sample_id`, `sample`, `sera`, `antibody_id` |

Both values are required for rows classified as samples. Explicitly designated control rows may omit sample or virus text when control assignment can be resolved from the other metadata.

### Optional fields

| Canonical field | Accepted examples | Notes |
|---|---|---|
| `control_type` | `control_type`, `control_designation`, `well_type`, `sample_type` | Accepted values include `sample`, `virus_control`, `VC`, `cell_control`, and `CC`. |
| `reading_datetime` | `reading_date_time`, `reading_datetime`, `read_datetime` | Retained as input metadata when supplied. |

When `control_type` is absent, the script infers controls from `well_id` or the sample name using `VC`, `Virus Control`, `CC`, or `Cell Control` conventions. **An explicit control designation that conflicts with those names causes the run to stop for review.**

### Raw-data validation

The script checks that:

- each physical `plate_id + well` combination is unique;
- each `plate_id` maps to only one `plate_number`;
- focus counts are finite and nonnegative;
- sample dilutions are finite and greater than zero;
- sample rows have virus and sample identifiers after optional lookup merging; and
- normalized column names do not become duplicates.

## Optional plate-layout lookup CSV

Use a lookup file when the raw-count export lacks complete virus, sample, or control metadata.

Required lookup fields:

| Canonical field | Accepted examples |
|---|---|
| `well_id` | `well_id`, `local_name`, `localname` |
| `plate_id` | `plate_id`, `plateid` |
| `sample_id` | `sample_id`, `global_name`, `sera`, `sample`, `antibody_id` |
| `virus` | `virus`, `virus_name` |

Optional lookup fields:

| Canonical field | Accepted examples |
|---|---|
| `plate_number` | `plate_number`, `plate_no`, `platenumber` |
| `control_type` | `control_type`, `control_designation`, `well_type`, `sample_type` |

The lookup must contain one unique row per `plate_id + well_id` mapping and must cover every raw-count mapping. When both files provide the same metadata field, conflicting nonblank values cause the run to stop rather than silently choosing one.

## Optional external four-parameter fit CSV

This file is not required for the R LL4 or observed-value uniroot estimates. It supplies optional comparison coefficients and fitted curves from compatible external software.

Required fields:

| Canonical field | Accepted examples |
|---|---|
| `a` | `a` |
| `b` | `b` |
| `c` | `c` |
| `d` | `d` |
| `plate_number` | `plate_number`, `plate_no`, `platenumber` |
| `plate_id` | `plate_id`, `plateid` |
| `virus` | `virus`, `virus_name` |
| `sample_id` | `global_name`, `sample_id`, `sample`, `sera`, `antibody_id` |

Optional fields:

| Canonical field | Accepted examples |
|---|---|
| `curve_name` | `curve_name` |
| `local_name` | `name`, `local_name` |
| `curve_formula` | `curve_formula`, `y_formula` |
| `r2` | `r2`, `r_squared` |
| `fit_f_prob` | `fit_f_prob`, `fit_f_probability` |

Only one external fit row may exist for each `plate_number + plate_id + virus + sample_id` key.

## Multiple viruses on one plate

Multiple sample viruses are allowed on one physical plate. Virus controls are assigned to a normalization virus using explicit metadata by default. Row-based reassignment is disabled unless `infer_control_virus_from_well_row = TRUE` is explicitly selected for a validated row-block layout.

When a plate lacks a matched control for a sample virus, the default behavior is to leave normalization unavailable. Same-virus controls from other plates are considered only when `allow_cross_plate_virus_controls = TRUE`.

## Example files

See [`example-data/README.md`](../example-data/README.md) for the supplied demonstration files and their intended use.
