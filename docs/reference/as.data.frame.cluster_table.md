# Coerce cluster_table to a Flat data.frame

Produces a flat data.frame suitable for table rendering. Cluster rows
have all columns populated; sub-peak rows have `k` and `volume_mm3` set
to `NA`.

## Usage

``` r
# S3 method for class 'cluster_table'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  max_clusters = 30L,
  max_sub_peaks = 3L,
  ...
)
```

## Arguments

- x:

  A `cluster_table` object.

- row.names:

  Ignored.

- optional:

  Ignored.

- max_clusters:

  Maximum clusters to include.

- max_sub_peaks:

  Maximum sub-peaks per cluster.

- ...:

  Ignored.

## Value

A data.frame.
