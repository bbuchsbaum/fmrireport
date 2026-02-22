# Format a cluster_table as a tinytable

Returns a
[`tinytable::tt()`](https://vincentarelbundock.github.io/tinytable/man/tt.html)
object with appropriate caption and formatting for Quarto/Typst
rendering.

## Usage

``` r
format_cluster_tt(x, max_clusters = 30L, max_sub_peaks = 3L, digits = 3L, ...)
```

## Arguments

- x:

  A `cluster_table` object.

- max_clusters:

  Maximum clusters to include.

- max_sub_peaks:

  Maximum sub-peaks per cluster.

- digits:

  Number of digits for rounding numeric columns.

- ...:

  Passed to
  [`tinytable::tt()`](https://vincentarelbundock.github.io/tinytable/man/tt.html).

## Value

A `tinytable` object.
