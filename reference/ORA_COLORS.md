# Color ramp for over-representation analysis (ORA) heatmaps

Two-color ramp for ORA significance heatmaps, from white at no
enrichment to purple at the strongest `-log10(p)`. Pass it as the
`colors` argument of
[`TMSig::enrichmap()`](https://rdrr.io/pkg/TMSig/man/enrichmap.html).

## Usage

``` r
ORA_COLORS
```

## Format

A named character vector of length 2: `low` (`"white"`) and `high`
(`"#543483"`).

## Examples

``` r
if (FALSE) { # \dontrun{
TMSig::enrichmap(ora_results, colors = ORA_COLORS)
} # }
```
