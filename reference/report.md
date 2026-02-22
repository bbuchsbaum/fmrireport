# Render Analysis Reports

Generic report interface and `fmri_lm` implementation that renders a
publication-style PDF report through Quarto.

## Usage

``` r
report(x, ...)

# Default S3 method
report(x, ...)

# S3 method for class 'fmri_lm'
report(
  x,
  output_file = "fmri_lm_report.pdf",
  title = "fMRI GLM Analysis Report",
  author = NULL,
  sections = c("model", "design", "hrf", "estimates", "contrasts", "diagnostics"),
  brain_map_stat = c("tstat", "estimate", "prob"),
  slice_axis = 3L,
  n_slices = 9L,
  threshold = NULL,
  bg_vol = NULL,
  atlas = NULL,
  cluster_thresh = 3,
  min_cluster_size = 10L,
  max_peaks = 15L,
  local_maxima_dist = 15,
  max_sub_peaks = 3L,
  open = interactive(),
  quiet = TRUE,
  ...
)
```

## Arguments

- x:

  Object to report on.

- ...:

  Additional arguments passed to methods.

- output_file:

  Output PDF path.

- title:

  Report title.

- author:

  Optional author string.

- sections:

  Sections to include.

- brain_map_stat:

  Statistic for contrast maps: `"tstat"`, `"estimate"`, or `"prob"`.

- slice_axis:

  Axis for montage slicing (1, 2, or 3).

- n_slices:

  Number of slices to render in contrast maps.

- threshold:

  Optional hard threshold for contrast map rendering.

- bg_vol:

  Optional background `NeuroVol` for overlay plots.

- atlas:

  Optional atlas object for peak labeling.

- cluster_thresh:

  Absolute statistic threshold for peak clustering.

- min_cluster_size:

  Minimum cluster size (voxels) in peak table.

- max_peaks:

  Maximum peaks per contrast table.

- local_maxima_dist:

  Minimum distance (mm) between local maxima within a cluster for
  sub-peak detection. Default 15.

- max_sub_peaks:

  Maximum number of sub-peaks per cluster in the activation table.
  Default 3.

- open:

  If `TRUE`, open the generated PDF after rendering.

- quiet:

  If `TRUE`, suppress Quarto render logs.

## Value

Method-specific report output. For `fmri_lm` fits, returns the output
PDF path invisibly.
