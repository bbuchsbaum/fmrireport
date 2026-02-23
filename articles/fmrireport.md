# Getting started with fmrireport

## Why fmrireport?

After fitting an fMRI GLM with `fmrireg`, you typically need to inspect
the model specification, visualize contrast maps, and report peak
activation coordinates. Doing this manually — extracting coefficient
summaries, building brain montages, assembling cluster tables — is
tedious and error-prone.

**fmrireport** generates publication-ready PDF reports from a fitted
`fmri_lm` object in a single function call. It also provides a
standalone
[`cluster_table()`](https://bbuchsbaum.github.io/fmrireport/reference/cluster_table.md)
function for building COBIDAS-compliant activation tables from any
statistical brain volume.

## Package overview

The package has two main entry points:

- **[`report()`](https://bbuchsbaum.github.io/fmrireport/reference/report.md)**
  — renders a full PDF report (model info, design matrix, HRF shapes,
  parameter estimates, contrast maps with peak tables, and diagnostics).
- **[`cluster_table()`](https://bbuchsbaum.github.io/fmrireport/reference/cluster_table.md)**
  — standalone function that detects suprathreshold clusters, identifies
  sub-peaks, computes p-values, and optionally labels regions with a
  brain atlas.

## Quick example: cluster_table()

You don’t need a full `fmri_lm` fit to use
[`cluster_table()`](https://bbuchsbaum.github.io/fmrireport/reference/cluster_table.md).
Any `NeuroVol` containing a statistical map will do.

``` r
library(fmrireport)
library(neuroim2)
#> Loading required package: Matrix
#> 
#> Attaching package: 'neuroim2'
#> The following object is masked from 'package:base':
#> 
#>     scale

# Create a synthetic t-stat volume with two activation clusters
set.seed(42)
sp <- NeuroSpace(c(64, 64, 40), spacing = c(3, 3, 3))
arr <- array(rnorm(64 * 64 * 40), c(64, 64, 40))

# Cluster 1: large blob in left hemisphere
arr[20:28, 25:33, 15:22] <- arr[20:28, 25:33, 15:22] + 5

# Cluster 2: smaller blob in right hemisphere
arr[40:46, 30:36, 20:26] <- arr[40:46, 30:36, 20:26] + 4

vol <- NeuroVol(arr, sp)
```

Now detect clusters and print the result:

``` r
ct <- cluster_table(vol, threshold = 3.0, stat_type = "t", df = 50)
print(ct)
#> cluster_table: 2 cluster(s)  |  threshold = 3.00 (t)  |  space = Unknown
#> 
#>   Cluster 1  (k = 638, 17226 mm3)
#>     Peak:   x = 60.0  y = 81.0  z = 51.0  stat = 8.208  p = 7.936e-11
#>       sub-peak 2:  x = 60.0  y = 72.0  z = 45.0  stat = 7.870  p = 2.647e-10
#>       sub-peak 3:  x = 66.0  y = 87.0  z = 63.0  stat = 7.834  p = 3.008e-10
#> 
#>   Cluster 2  (k = 281, 7587 mm3)
#>     Peak:   x = 132.0  y = 96.0  z = 75.0  stat = 6.562  p = 2.905e-08
#>       sub-peak 2:  x = 132.0  y = 87.0  z = 57.0  stat = 6.416  p = 4.923e-08
#>       sub-peak 3:  x = 135.0  y = 105.0  z = 63.0  stat = 6.139  p = 1.329e-07
```

The output shows each cluster with its size, world coordinates, peak
statistic, and p-value. Sub-peaks (local maxima within the same cluster)
are indented beneath their parent.

## Flat table for rendering

Use [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) to
get a flat table suitable for inclusion in documents:

``` r
df <- as.data.frame(ct)
head(df, 10)
#>   Cluster   k Vol_mm3   x   y  z     Stat            p Region Hemi
#> 1       1 638   17226  60  81 51 8.208263 7.936489e-11     --   --
#> 2       1  NA      NA  60  72 45 7.869639 2.647068e-10     --   --
#> 3       1  NA      NA  66  87 63 7.833813 3.008181e-10     --   --
#> 4       2 281    7587 132  96 75 6.562311 2.904745e-08     --   --
#> 5       2  NA      NA 132  87 57 6.415751 4.923016e-08     --   --
#> 6       2  NA      NA 135 105 63 6.139384 1.328658e-07     --   --
```

Sub-peak rows have `NA` in the `k` and `Vol_mm3` columns, making it easy
to distinguish them from cluster-level rows.

## MNI coordinate labels

When you specify a known coordinate space, column names adapt
automatically:

``` r
ct_mni <- cluster_table(vol, threshold = 3.0,
                        stat_type = "t", df = 50,
                        coord_space = "MNI152")
df_mni <- as.data.frame(ct_mni)
names(df_mni)
#>  [1] "Cluster" "k"       "Vol_mm3" "MNI_x"   "MNI_y"   "MNI_z"   "Stat"   
#>  [8] "p"       "Region"  "Hemi"
```

## Tinytable output

For Quarto or R Markdown reports,
[`format_cluster_tt()`](https://bbuchsbaum.github.io/fmrireport/reference/format_cluster_tt.md)
produces a formatted table with a descriptive caption:

``` r
format_cluster_tt(ct)
```

| Cluster | k   | Vol_mm3 | x   | y   | z   | Stat  | p   | Region | Hemi |
|---------|-----|---------|-----|-----|-----|-------|-----|--------|------|
| 1       | 638 | 17226   | 60  | 81  | 51  | 8.208 | 0   | --     | --   |
| 1       | NA  | NA      | 60  | 72  | 45  | 7.870 | 0   | --     | --   |
| 1       | NA  | NA      | 66  | 87  | 63  | 7.834 | 0   | --     | --   |
| 2       | 281 | 7587    | 132 | 96  | 75  | 6.562 | 0   | --     | --   |
| 2       | NA  | NA      | 132 | 87  | 57  | 6.416 | 0   | --     | --   |
| 2       | NA  | NA      | 135 | 105 | 63  | 6.139 | 0   | --     | --   |

Threshold: t = 3 \| df = 50 \| Space: Unknown

## Full PDF report workflow

If you have a fitted `fmri_lm` object from `fmrireg`, generating a full
report is a single call:

``` r
library(fmrireg)

# Fit a model (shown for context; not run here)
fit <- fmri_lm(
  onset ~ hrf(condition, contrasts = my_contrasts),
  block = ~run,
  dataset = my_dataset
)

# Generate the report
report(fit,
       output_file = "my_analysis.pdf",
       title = "Visual Localizer Analysis",
       author = "Jane Doe",
       cluster_thresh = 3.5,
       min_cluster_size = 20L,
       atlas = my_atlas)
```

The report includes:

| Section         | Contents                                                     |
|-----------------|--------------------------------------------------------------|
| **Model**       | Formula, dataset class, voxel count, residual df, conditions |
| **Design**      | Design matrix heatmap, regressor correlation plot            |
| **HRF**         | Fitted hemodynamic response curves per condition             |
| **Estimates**   | Mean/SD/median betas, max t-stat, % significant voxels       |
| **Contrasts**   | Brain map montage, hierarchical activation table             |
| **Diagnostics** | AR parameters, residual sigma summary                        |

You can select which sections to include with the `sections` argument.

## Controlling the activation table

Two new parameters on
[`report()`](https://bbuchsbaum.github.io/fmrireport/reference/report.md)
control the cluster table:

- `local_maxima_dist` — minimum distance (mm) between sub-peaks within a
  cluster (default: 15)
- `max_sub_peaks` — maximum sub-peaks reported per cluster (default: 3)

``` r
report(fit,
       output_file = "report.pdf",
       local_maxima_dist = 12,
       max_sub_peaks = 5)
```

## Next steps

- [`vignette("cluster-tables")`](https://bbuchsbaum.github.io/fmrireport/articles/cluster-tables.md)
  — in-depth guide to
  [`cluster_table()`](https://bbuchsbaum.github.io/fmrireport/reference/cluster_table.md)
  with synthetic Gaussian activations and atlas labeling
- [`?cluster_table`](https://bbuchsbaum.github.io/fmrireport/reference/cluster_table.md)
  — full parameter reference
- [`?report.fmri_lm`](https://bbuchsbaum.github.io/fmrireport/reference/report.md)
  — all report customization options
