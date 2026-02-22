# fmrireport

`fmrireport` generates publication-oriented PDF reports for fitted
`fmrireg::fmri_lm` models.

## Example

```r
library(fmrireg)
library(fmrireport)

# fit <- fmri_lm(...)
# report(fit, output_file = "fmri_lm_report.pdf")
```

## Synthetic demo reports

Run:

```r
Rscript inst/examples/fmri_lm_report_demo.R /tmp/fmrireport-report-examples
```

This writes three complete example reports (matrix, spatial, spatial+atlas).
