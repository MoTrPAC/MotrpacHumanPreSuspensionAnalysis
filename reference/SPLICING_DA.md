# Differential Alternative Splicing Results

A data object containing results from a differential alternative
splicing analysis. The analytical workflow used to generate these
results—including read alignment, isoform quantification, splicing event
detection, and statistical testing—is described in detail in:
"Characterization of exercise-modulated alternative splicing landscape
in human skeletal muscle, adipose tissue, and blood in the MoTrPAC
Study" To accommodate file size and distribution constraints, this
object includes only statistically significant splicing isoforms,
defined by a false discovery rate (FDR) less than 0.05. As such, the
object represents a filtered subset of all tested isoforms and is
intended for downstream analysis, visualization, and integration with
other molecular datasets rather than for reproducing the full splicing
analysis pipeline.

## Usage

``` r
SPLICING_DA
```

## Format

A data frame containing differential splicing results for significant
isoforms only (FDR \< 0.05).

## Source

Generated as described in the splicing companion paper
