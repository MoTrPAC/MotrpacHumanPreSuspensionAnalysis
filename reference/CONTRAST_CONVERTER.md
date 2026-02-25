# Contrast Converter

A `data.table` containing information about contrasts.

## Usage

``` r
CONTRAST_CONVERTER
```

## Format

A `data.table` with 33 rows and 5 columns, where each row is a different
contrast:

- contrast_order:

  integer; the order of the contrast.

- contrast:

  factor; full contrast (33 levels).

- contrast_short:

  factor; shortened version of the contrasts.

- contrast_type:

  factor; one of "exercise_with_controls", "exercise_no_controls",
  "Endur_vs_Resist", "baseline", or "control_only".

- contrast_category:

  factor; one of "EE-CON", "RE-CON", "EE-EE", "RE-RE", "EE-RE", or
  "CON-CON".
