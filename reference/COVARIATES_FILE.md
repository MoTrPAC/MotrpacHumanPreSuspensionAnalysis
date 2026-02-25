# Covariates used for statistical analysis

Table of covariates used in the differential abundance analysis.
`visitcode` is only used for training analysis. This is also used for
any batch correction for regressing out technical covariates in the
structure.

## Usage

``` r
COVARIATES_FILE
```

## Format

A `data.frame` object with the columns:

- ome:

  character; the ome that the covariates are for

- tech_or_design:

  factor; either technical or Design, to designate if individual
  covariates should be regressed out in the norm-qc tables

- data_type:

  factor; which data type a covariate is (e.g. factor)

- covariate:

  character; the name of the covariate

- tissue:

  factor; one of `tissue_available_list`.
