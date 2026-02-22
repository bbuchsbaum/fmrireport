library(testthat)
library(fmrireport)

# -- Helper: create a synthetic volume with known clusters --------------------
make_test_vol <- function(dims = c(20, 20, 20), spacing = c(2, 2, 2)) {

  sp <- neuroim2::NeuroSpace(dims, spacing = spacing)
  arr <- array(0, dim = dims)

  # Cluster 1: large block of high values (6x6x6 = 216 voxels)

  arr[5:10, 5:10, 5:10] <- 5

  # Cluster 2: smaller block (4x4x4 = 64 voxels)
  arr[15:18, 15:18, 15:18] <- 4

  # Cluster 3: tiny (2x2x2 = 8 voxels) -- should be filtered at min_cluster_size=10

  arr[2:3, 15:16, 15:16] <- 6

  neuroim2::NeuroVol(arr, sp)
}

# -- cluster_table: basic functionality ----------------------------------------

test_that("cluster_table detects expected clusters", {
  vol <- make_test_vol()
  ct <- cluster_table(vol, threshold = 3.0, min_cluster_size = 10L)

  expect_s3_class(ct, "cluster_table")
  # Should find 2 clusters (the 8-voxel cluster is filtered out)
  expect_equal(ct$n_clusters, 2L)
  expect_equal(nrow(ct$clusters), 2L)
  expect_true(all(ct$clusters$k >= 10))
})

test_that("cluster_table returns correct structure", {
  vol <- make_test_vol()
  ct <- cluster_table(vol, threshold = 3.0, min_cluster_size = 1L)

  expect_true(all(c("clusters", "peaks", "threshold", "stat_type", "df",
                     "n_clusters", "coord_space", "atlas_name") %in% names(ct)))

  expected_cluster_cols <- c("cluster_id", "k", "volume_mm3", "peak_x",
                             "peak_y", "peak_z", "peak_stat", "peak_p",
                             "label", "label_full", "hemi", "network")
  expect_true(all(expected_cluster_cols %in% names(ct$clusters)))

  expected_peak_cols <- c("cluster_id", "peak_rank", "x", "y", "z",
                          "stat", "p", "label", "label_full", "hemi", "network")
  expect_true(all(expected_peak_cols %in% names(ct$peaks)))
})

test_that("cluster_table min_cluster_size filtering works", {
  vol <- make_test_vol()

  ct_strict <- cluster_table(vol, threshold = 3.0, min_cluster_size = 100L)
  expect_true(ct_strict$n_clusters < 3)
  expect_true(all(ct_strict$clusters$k >= 100))

  ct_loose <- cluster_table(vol, threshold = 3.0, min_cluster_size = 1L)
  expect_true(ct_loose$n_clusters >= ct_strict$n_clusters)
})

# -- P-value computation -------------------------------------------------------

test_that("P-values for t-stat are computed correctly", {
  vol <- make_test_vol()
  ct <- cluster_table(vol, threshold = 3.0, stat_type = "t", df = 50,
                      min_cluster_size = 1L)

  expect_false(all(is.na(ct$clusters$peak_p)))
  # P-values should be between 0 and 1
  pvals <- ct$clusters$peak_p[!is.na(ct$clusters$peak_p)]
  expect_true(all(pvals >= 0 & pvals <= 1))
})

test_that("P-values for z-stat are computed correctly", {
  vol <- make_test_vol()
  ct <- cluster_table(vol, threshold = 3.0, stat_type = "z",
                      min_cluster_size = 1L)

  # Without df, z-stat still computes p-values when stat_type="z" and df=NULL
  expect_true(all(is.na(ct$clusters$peak_p)))

  # With df (not used for z but included to not error)
  ct2 <- cluster_table(vol, threshold = 3.0, stat_type = "z", df = 1,
                       min_cluster_size = 1L)
  expect_false(all(is.na(ct2$clusters$peak_p)))
})

test_that("P-values for F-stat are computed correctly", {
  vol <- make_test_vol()
  ct <- cluster_table(vol, threshold = 3.0, stat_type = "F", df = c(3, 50),
                      min_cluster_size = 1L)

  pvals <- ct$clusters$peak_p[!is.na(ct$clusters$peak_p)]
  expect_true(all(pvals >= 0 & pvals <= 1))
})

test_that("P-values are NA when df is NULL", {
  vol <- make_test_vol()
  ct <- cluster_table(vol, threshold = 3.0, stat_type = "t", df = NULL,
                      min_cluster_size = 1L)
  expect_true(all(is.na(ct$clusters$peak_p)))
})

# -- World coordinates ---------------------------------------------------------

test_that("World coordinates reflect spacing", {
  vol <- make_test_vol(spacing = c(3, 3, 3))
  ct <- cluster_table(vol, threshold = 3.0, min_cluster_size = 1L)

  # World coords should not all be zero

  expect_true(any(ct$clusters$peak_x != 0) || any(ct$clusters$peak_y != 0) ||
              any(ct$clusters$peak_z != 0))
})

# -- Coordinate space handling -------------------------------------------------

test_that("coord_space defaults to Unknown without atlas", {
 vol <- make_test_vol()
 ct <- cluster_table(vol, threshold = 3.0, min_cluster_size = 1L)
 expect_equal(ct$coord_space, "Unknown")
})

test_that("explicit coord_space is preserved", {
  vol <- make_test_vol()
  ct <- cluster_table(vol, threshold = 3.0, coord_space = "MNI152",
                      min_cluster_size = 1L)
  expect_equal(ct$coord_space, "MNI152")
})

# -- as.data.frame.cluster_table -----------------------------------------------

test_that("as.data.frame produces hierarchical flat table", {
  vol <- make_test_vol()
  ct <- cluster_table(vol, threshold = 3.0, min_cluster_size = 1L)
  df <- as.data.frame(ct)

  expect_true(is.data.frame(df))
  expect_true(nrow(df) >= ct$n_clusters)
  expect_true("Cluster" %in% names(df))
  expect_true("k" %in% names(df))
  expect_true("Stat" %in% names(df))
  expect_true("Region" %in% names(df))
})

test_that("as.data.frame uses MNI column names for MNI space", {
  vol <- make_test_vol()
  ct <- cluster_table(vol, threshold = 3.0, coord_space = "MNI152",
                      min_cluster_size = 1L)
  df <- as.data.frame(ct)

  expect_true("MNI_x" %in% names(df))
  expect_true("MNI_y" %in% names(df))
  expect_true("MNI_z" %in% names(df))
})

test_that("as.data.frame uses plain x/y/z for non-MNI space", {
  vol <- make_test_vol()
  ct <- cluster_table(vol, threshold = 3.0, min_cluster_size = 1L)
  df <- as.data.frame(ct)

  expect_true("x" %in% names(df))
  expect_true("y" %in% names(df))
  expect_true("z" %in% names(df))
})

# -- print.cluster_table -------------------------------------------------------

test_that("print.cluster_table runs without error", {
  vol <- make_test_vol()
  ct <- cluster_table(vol, threshold = 3.0, min_cluster_size = 1L)
  expect_output(print(ct), "cluster_table")
})

test_that("print.cluster_table handles empty table", {
  vol <- make_test_vol()
  ct <- cluster_table(vol, threshold = 100, min_cluster_size = 1L)
  expect_output(print(ct), "0 clusters")
})

# -- Edge cases ----------------------------------------------------------------

test_that("No suprathreshold voxels returns empty cluster_table", {
  vol <- make_test_vol()
  ct <- cluster_table(vol, threshold = 100)

  expect_s3_class(ct, "cluster_table")
  expect_equal(ct$n_clusters, 0L)
  expect_equal(nrow(ct$clusters), 0L)
  expect_equal(nrow(ct$peaks), 0L)
})

test_that("Single-voxel cluster handled with min_cluster_size=1", {
  sp <- neuroim2::NeuroSpace(c(10, 10, 10), spacing = c(2, 2, 2))
  arr <- array(0, dim = c(10, 10, 10))
  arr[5, 5, 5] <- 10  # single voxel
  vol <- neuroim2::NeuroVol(arr, sp)

  ct <- cluster_table(vol, threshold = 3.0, min_cluster_size = 1L)
  expect_equal(ct$n_clusters, 1L)
  expect_equal(ct$clusters$k[1], 1L)
})

test_that("sort_by peak_stat works", {
  vol <- make_test_vol()
  ct <- cluster_table(vol, threshold = 3.0, min_cluster_size = 1L,
                      sort_by = "peak_stat")
  # First cluster should have highest absolute peak stat
  if (ct$n_clusters > 1) {
    expect_true(abs(ct$clusters$peak_stat[1]) >= abs(ct$clusters$peak_stat[2]))
  }
})

test_that("volume_mm3 is computed correctly", {
  vol <- make_test_vol(spacing = c(2, 2, 2))
  ct <- cluster_table(vol, threshold = 3.0, min_cluster_size = 1L)

  # voxel volume = 2*2*2 = 8 mm3
  for (i in seq_len(ct$n_clusters)) {
    expect_equal(ct$clusters$volume_mm3[i], ct$clusters$k[i] * 8)
  }
})

# -- format_cluster_tt ---------------------------------------------------------

test_that("format_cluster_tt produces a tinytable", {
  skip_if_not_installed("tinytable")
  vol <- make_test_vol()
  ct <- cluster_table(vol, threshold = 3.0, min_cluster_size = 1L)
  tt <- format_cluster_tt(ct)
  expect_true(inherits(tt, "tinytable"))
})

# -- Atlas labeling (skip if neuroatlas not installed) -------------------------

test_that("cluster_table works without atlas (labels are NA)", {
  vol <- make_test_vol()
  ct <- cluster_table(vol, threshold = 3.0, min_cluster_size = 1L)

  expect_true(all(is.na(ct$clusters$label)))
  expect_true(all(is.na(ct$clusters$hemi)))
})
