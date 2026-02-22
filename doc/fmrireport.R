## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4,
  dpi = 96
)

## ----quick-example------------------------------------------------------------
library(fmrireport)
library(neuroim2)

# Create a synthetic t-stat volume with two activation clusters
set.seed(42)
sp <- NeuroSpace(c(64, 64, 40), spacing = c(3, 3, 3))
arr <- array(rnorm(64 * 64 * 40), c(64, 64, 40))

# Cluster 1: large blob in left hemisphere
arr[20:28, 25:33, 15:22] <- arr[20:28, 25:33, 15:22] + 5

# Cluster 2: smaller blob in right hemisphere
arr[40:46, 30:36, 20:26] <- arr[40:46, 30:36, 20:26] + 4

vol <- NeuroVol(arr, sp)

## ----cluster-detect-----------------------------------------------------------
ct <- cluster_table(vol, threshold = 3.0, stat_type = "t", df = 50)
print(ct)

## ----flat-table---------------------------------------------------------------
df <- as.data.frame(ct)
head(df, 10)

## ----mni-coords---------------------------------------------------------------
ct_mni <- cluster_table(vol, threshold = 3.0,
                        stat_type = "t", df = 50,
                        coord_space = "MNI152")
df_mni <- as.data.frame(ct_mni)
names(df_mni)

## ----tinytable, eval = requireNamespace("tinytable", quietly = TRUE)----------
format_cluster_tt(ct)

## ----report-example, eval = FALSE---------------------------------------------
# library(fmrireg)
# 
# # Fit a model (shown for context; not run here)
# fit <- fmri_lm(
#   onset ~ hrf(condition, contrasts = my_contrasts),
#   block = ~run,
#   dataset = my_dataset
# )
# 
# # Generate the report
# report(fit,
#        output_file = "my_analysis.pdf",
#        title = "Visual Localizer Analysis",
#        author = "Jane Doe",
#        cluster_thresh = 3.5,
#        min_cluster_size = 20L,
#        atlas = my_atlas)

## ----report-params, eval = FALSE----------------------------------------------
# report(fit,
#        output_file = "report.pdf",
#        local_maxima_dist = 12,
#        max_sub_peaks = 5)

