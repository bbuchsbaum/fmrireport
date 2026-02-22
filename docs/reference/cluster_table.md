# Build a Publication-Quality Activation Cluster Table

Detects suprathreshold clusters in a statistical brain volume and
returns a structured table of cluster-level and peak-level statistics
suitable for publication (NeuroImage / HBM / COBIDAS compliant).

## Usage

``` r
cluster_table(
  vol,
  threshold,
  atlas = NULL,
  stat_type = c("t", "z", "F", "other"),
  df = NULL,
  coord_space = NULL,
  min_cluster_size = 10L,
  local_maxima_dist = 15,
  max_peaks = 3L,
  connectivity = c("26-connect", "18-connect", "6-connect"),
  sort_by = c("size", "peak_stat")
)
```

## Arguments

- vol:

  A `NeuroVol` object containing the statistical map.

- threshold:

  Numeric threshold for cluster detection (absolute value).

- atlas:

  Optional atlas object (e.g., a `volatlas` from neuroatlas) for
  anatomical labeling.

- stat_type:

  Type of statistic in the volume: `"t"`, `"z"`, `"F"`, or `"other"`.

- df:

  Degrees of freedom for p-value computation. A single value for t/z
  statistics, or a length-2 vector `c(df1, df2)` for F statistics.
  `NULL` omits p-values.

- coord_space:

  Character string naming the coordinate space (e.g., `"MNI152"`,
  `"MNI305"`). If `NULL`, inferred from `atlas` or set to `"Unknown"`.

- min_cluster_size:

  Minimum number of voxels for a cluster to be retained.

- local_maxima_dist:

  Minimum distance (mm) between local maxima within a cluster.

- max_peaks:

  Maximum number of sub-peaks reported per cluster.

- connectivity:

  Voxel connectivity for cluster detection.

- sort_by:

  Sort clusters by `"size"` (descending) or `"peak_stat"` (descending
  absolute value).

## Value

An S3 object of class `"cluster_table"` with elements:

- clusters:

  A data.frame with one row per cluster.

- peaks:

  A data.frame with sub-peaks (local maxima) per cluster.

- threshold:

  The threshold used.

- stat_type:

  The statistic type.

- df:

  Degrees of freedom.

- n_clusters:

  Number of clusters.

- coord_space:

  Coordinate space label.

- atlas_name:

  Name of atlas used, or `NA`.
