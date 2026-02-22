## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  dpi = 96
)

## ----make-synthetic-----------------------------------------------------------
library(fmrireport)
library(neuroim2)

# Define a brain-like volume in MNI-ish space
# 2mm isotropic, 91x109x91 grid, origin at (-90, -126, -72)
sp <- NeuroSpace(c(91, 109, 91),
                 spacing = c(2, 2, 2),
                 origin = c(-90, -126, -72))

# Start with Gaussian noise (simulating a t-stat map under the null)
set.seed(2024)
arr <- array(rnorm(prod(dim(sp))), dim(sp))

## ----add-blobs----------------------------------------------------------------
# Helper: add a 3D Gaussian blob centered at grid coordinate (cx, cy, cz)
add_gaussian_blob <- function(arr, cx, cy, cz, peak_value, sigma = 4) {
  dims <- dim(arr)
  # Work in a local window for efficiency
  radius <- ceiling(3 * sigma)
  x_range <- max(1, cx - radius):min(dims[1], cx + radius)
  y_range <- max(1, cy - radius):min(dims[2], cy + radius)
  z_range <- max(1, cz - radius):min(dims[3], cz + radius)

  for (xi in x_range) {
    for (yi in y_range) {
      for (zi in z_range) {
        d2 <- (xi - cx)^2 + (yi - cy)^2 + (zi - cz)^2
        arr[xi, yi, zi] <- arr[xi, yi, zi] +
          peak_value * exp(-d2 / (2 * sigma^2))
      }
    }
  }
  arr
}

# Place blobs at grid coordinates that map to known cortical MNI regions.
# With origin=(-90,-126,-72) and spacing=2mm:
#   grid (27,51,65) -> MNI (-38,-26,56) = Left SomMot
#   grid (65,51,65) -> MNI ( 38,-26,56) = Right SomMot
#   grid (46,55,67) -> MNI (  0,-18,62) = Medial motor
#   grid (41,21,39) -> MNI (-10,-86, 4) = Left V1
#   grid (51,21,39) -> MNI ( 10,-86, 4) = Right V1

# Blob 1: Left somatomotor cortex (large, strong activation)
arr <- add_gaussian_blob(arr, cx = 27, cy = 51, cz = 65, peak_value = 8, sigma = 5)

# Blob 2: Right somatomotor cortex (moderate)
arr <- add_gaussian_blob(arr, cx = 65, cy = 51, cz = 65, peak_value = 6, sigma = 4)

# Blob 3: Medial motor / SMA (moderate)
arr <- add_gaussian_blob(arr, cx = 46, cy = 55, cz = 67, peak_value = 5.5, sigma = 3)

# Blob 4: Left primary visual cortex (small, strong)
arr <- add_gaussian_blob(arr, cx = 41, cy = 21, cz = 39, peak_value = 7, sigma = 3)

# Blob 5: Right primary visual cortex (small, moderate)
arr <- add_gaussian_blob(arr, cx = 51, cy = 21, cz = 39, peak_value = 5, sigma = 3)

vol <- NeuroVol(arr, sp)

## ----plot-montage, fig.height = 3.5-------------------------------------------
# Show axial slices at the levels where we placed blobs
slices <- c(39, 55, 60, 65, 67)
plot_montage(vol, zlevels = slices, along = 3L,
             cmap = "inferno", ncol = 5L,
             title = "Synthetic t-statistic map")

## ----basic-cluster------------------------------------------------------------
ct <- cluster_table(vol,
                    threshold = 3.0,
                    stat_type = "t",
                    df = 100,
                    min_cluster_size = 20)
print(ct)

## ----subpeak-dist-------------------------------------------------------------
# Tight distance: many sub-peaks
ct_tight <- cluster_table(vol, threshold = 3.0, stat_type = "t", df = 100,
                          local_maxima_dist = 8, max_peaks = 5)
cat("Sub-peaks with 8mm distance:\n")
nrow(ct_tight$peaks)

# Wide distance: fewer sub-peaks
ct_wide <- cluster_table(vol, threshold = 3.0, stat_type = "t", df = 100,
                         local_maxima_dist = 25, max_peaks = 5)
cat("Sub-peaks with 25mm distance:\n")
nrow(ct_wide$peaks)

## ----sort-options-------------------------------------------------------------
ct_by_stat <- cluster_table(vol, threshold = 3.0, stat_type = "t", df = 100,
                            sort_by = "peak_stat")
# First cluster has the highest absolute peak stat
ct_by_stat$clusters[, c("cluster_id", "k", "peak_stat")]

## ----pvalues------------------------------------------------------------------
# t-statistic
ct_t <- cluster_table(vol, threshold = 3.0, stat_type = "t", df = 100)
head(ct_t$clusters[, c("peak_stat", "peak_p")])

# z-statistic (provide df to enable p-value computation)
ct_z <- cluster_table(vol, threshold = 3.0, stat_type = "z", df = 1)
head(ct_z$clusters[, c("peak_stat", "peak_p")])

# Without df, p-values are NA
ct_nodf <- cluster_table(vol, threshold = 3.0, stat_type = "t", df = NULL)
head(ct_nodf$clusters[, c("peak_stat", "peak_p")])

## ----flat-df------------------------------------------------------------------
df <- as.data.frame(ct)
df

## ----mni-columns--------------------------------------------------------------
ct_mni <- cluster_table(vol, threshold = 3.0, stat_type = "t", df = 100,
                        coord_space = "MNI152")
df_mni <- as.data.frame(ct_mni)
names(df_mni)

## ----tinytable-output, eval = requireNamespace("tinytable", quietly = TRUE)----
tt <- format_cluster_tt(ct_mni, max_sub_peaks = 2, digits = 2)
tt

## ----atlas-label, eval = requireNamespace("neuroatlas", quietly = TRUE)-------
library(neuroatlas)

# Load the Schaefer 200-parcel, 7-network atlas
atlas <- get_schaefer_atlas(parcels = "200", networks = "7")

ct_labeled <- cluster_table(vol, threshold = 3.0,
                            stat_type = "t", df = 100,
                            atlas = atlas,
                            coord_space = "MNI152")
print(ct_labeled)

## ----atlas-df, eval = requireNamespace("neuroatlas", quietly = TRUE)----------
as.data.frame(ct_labeled)

## ----overlay-plot, fig.height = 4---------------------------------------------
# Create a thresholded overlay for visualization
arr_thresh <- as.array(vol)
arr_thresh[abs(arr_thresh) < 3.0] <- NA
vol_thresh <- NeuroVol(arr_thresh, sp)

# Create a background volume (uniform gray for this demo)
bg_arr <- array(1000, dim(sp))
bg <- NeuroVol(bg_arr, sp)

plot_overlay(bg, vol_thresh,
             zlevels = c(39, 55, 60, 65),
             along = 3L,
             ov_cmap = "inferno",
             ov_thresh = 0,
             ncol = 4L,
             title = "Thresholded activation (t > 3.0)")

## ----connectivity-------------------------------------------------------------
ct_6 <- cluster_table(vol, threshold = 3.0, stat_type = "t", df = 100,
                      connectivity = "6-connect")
ct_26 <- cluster_table(vol, threshold = 3.0, stat_type = "t", df = 100,
                       connectivity = "26-connect")
cat("6-connect clusters:", ct_6$n_clusters, "\n")
cat("26-connect clusters:", ct_26$n_clusters, "\n")

## ----edge-cases---------------------------------------------------------------
# No suprathreshold voxels
ct_empty <- cluster_table(vol, threshold = 100)
ct_empty$n_clusters
print(ct_empty)

# Very liberal threshold (many small clusters filtered by min_cluster_size)
ct_liberal <- cluster_table(vol, threshold = 1.5, min_cluster_size = 50)
cat("Clusters with k >= 50:", ct_liberal$n_clusters, "\n")

## ----structure----------------------------------------------------------------
str(ct, max.level = 1)

