# fmrireport

`fmrireport` generates publication-oriented PDF reports for fitted
[`fmrireg::fmri_lm`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.html)
models. It also provides a standalone
[`cluster_table()`](https://bbuchsbaum.github.io/fmrireport/reference/cluster_table.md)
function for building COBIDAS-compliant activation tables from any
statistical brain volume.

## Documentation

Full documentation and vignettes:
<https://bbuchsbaum.github.io/fmrireport/>

## Installation

``` r
# install.packages("remotes")
remotes::install_github("bbuchsbaum/fmrireport")
```

## Quick example

``` r
library(fmrireport)
library(neuroim2)

# Create a synthetic t-stat volume
set.seed(42)
sp <- NeuroSpace(c(64, 64, 40), spacing = c(3, 3, 3))
arr <- array(rnorm(64 * 64 * 40), c(64, 64, 40))
arr[20:28, 25:33, 15:22] <- arr[20:28, 25:33, 15:22] + 5

vol <- NeuroVol(arr, sp)
ct <- cluster_table(vol, threshold = 3.0, stat_type = "t", df = 50)
print(ct)
```

## Full PDF report

``` r
library(fmrireg)

# fit <- fmri_lm(...)
# report(fit, output_file = "my_analysis.pdf")
```

## Synthetic demo reports

``` r
Rscript inst/examples/fmri_lm_report_demo.R /tmp/fmrireport-report-examples
```

This writes three complete example reports (matrix, spatial,
spatial+atlas).
