# Print a cluster_table Object

CLI-formatted hierarchical output showing clusters with sub-peaks
indented.

## Usage

``` r
# S3 method for class 'cluster_table'
print(x, max_clusters = 20L, max_sub_peaks = 3L, ...)
```

## Arguments

- x:

  A `cluster_table` object.

- max_clusters:

  Maximum number of clusters to display.

- max_sub_peaks:

  Maximum sub-peaks per cluster to display.

- ...:

  Ignored.

## Value

Invisibly returns `x`.
