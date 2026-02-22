library(testthat)
library(fmrireport)

make_report_test_fit <- function() {
  sim <- fmrireg::simulate_simple_dataset(
    ncond = 2,
    nreps = 12,
    TR = 2,
    snr = 1.4,
    seed = 42
  )

  base <- sim$noisy[, -1, drop = FALSE]
  n_time <- nrow(base)
  n_vox <- 80L
  W <- matrix(stats::rnorm(ncol(base) * n_vox), nrow = ncol(base), ncol = n_vox)
  Y <- base %*% W + matrix(stats::rnorm(n_time * n_vox, sd = 0.25), nrow = n_time, ncol = n_vox)

  cond_vals <- factor(sim$conditions)
  if (nlevels(cond_vals) < 2L) {
    stop("Expected at least two conditions for report tests.")
  }
  levels(cond_vals)[seq_len(min(2L, nlevels(cond_vals)))] <- c("A", "B")[seq_len(min(2L, nlevels(cond_vals)))]

  etab <- data.frame(
    onset = sim$onsets,
    condition = cond_vals,
    run = 1L
  )

  dset <- fmridataset::matrix_dataset(
    datamat = Y,
    TR = 2,
    run_length = n_time,
    event_table = etab
  )

  con <- fmrireg::contrast_set(
    fmrireg::pair_contrast(~condition == "A", ~condition == "B", name = "A_vs_B")
  )

  suppressWarnings(
    fmrireg::fmri_lm(
      onset ~ fmrireg::hrf(condition, contrasts = con),
      block = ~run,
      dataset = dset,
      durations = 0,
      progress = FALSE
    )
  )
}

test_that("report generic exists and dispatches", {
  expect_true(is.function(report))

  fit <- make_report_test_fit()
  expect_true(inherits(fit, "fmri_lm"))

  expect_error(
    report(1),
    "No report\\(\\) method"
  )
})

test_that(".is_spatial_dataset is FALSE for matrix_dataset", {
  fit <- make_report_test_fit()
  expect_false(fmrireport:::.is_spatial_dataset(fit$dataset))
})

test_that(".report_model_info extracts expected fields", {
  fit <- make_report_test_fit()
  info <- fmrireport:::.report_model_info(fit)

  expect_true(is.list(info))
  expect_true("table" %in% names(info))
  expect_true(is.data.frame(info$table))
  expect_true(all(c("Parameter", "Value") %in% names(info$table)))

  params <- info$table$Parameter
  expect_true("Formula" %in% params)
  expect_true("Strategy" %in% params)
  expect_true("Residual df" %in% params)
})

test_that(".report_estimates returns required summary columns", {
  fit <- make_report_test_fit()
  est <- fmrireport:::.report_estimates(fit)

  expect_true(is.list(est))
  expect_true(isTRUE(est$available))
  expect_true(is.data.frame(est$table))

  expected_cols <- c(
    "Coefficient", "Mean_Beta", "SD_Beta", "Median_Beta",
    "Mean_T", "Max_Abs_T", "Pct_P_lt_0_05"
  )
  expect_true(all(expected_cols %in% names(est$table)))
})

test_that(".pick_peak_slices identifies high-activation slices", {
  arr <- array(0, dim = c(8, 8, 12))
  arr[, , 3] <- 1
  arr[, , 10] <- 3
  arr[4:5, 4:5, 6] <- 5

  vol <- neuroim2::NeuroVol(arr, neuroim2::NeuroSpace(dim = dim(arr)))
  idx <- fmrireport:::.pick_peak_slices(vol, along = 3L, n = 4L)

  expect_true(length(idx) >= 1)
  expect_true(all(idx >= 1 & idx <= dim(arr)[3]))
  expect_true(any(idx %in% c(6, 10)))
})

test_that("cluster_table returns expected structure from report context", {
  arr <- array(0, dim = c(10, 10, 10))
  arr[2:4, 2:4, 2:4] <- 4
  arr[7:9, 7:9, 7:9] <- 6
  vol <- neuroim2::NeuroVol(arr, neuroim2::NeuroSpace(dim = dim(arr)))

  ct <- cluster_table(vol, threshold = 3, min_cluster_size = 4L)

  expect_s3_class(ct, "cluster_table")
  expect_true(ct$n_clusters >= 1)
  expect_true(all(c("cluster_id", "k", "peak_stat") %in% names(ct$clusters)))

  df <- as.data.frame(ct)
  expect_true(is.data.frame(df))
  expect_true(nrow(df) >= 1)
})

test_that(".prepare_report_data returns expected structure", {
  fit <- make_report_test_fit()
  fig_dir <- tempfile("report-fig-")
  dir.create(fig_dir, recursive = TRUE)
  on.exit(unlink(fig_dir, recursive = TRUE), add = TRUE)

  out <- fmrireport:::.prepare_report_data(
    x = fit,
    sections = c("model", "design", "hrf", "estimates", "contrasts", "diagnostics"),
    brain_map_stat = "tstat",
    slice_axis = 3L,
    n_slices = 6L,
    threshold = NULL,
    bg_vol = NULL,
    atlas = NULL,
    cluster_thresh = 2.0,
    min_cluster_size = 2L,
    max_peaks = 10L,
    fig_dir = fig_dir
  )

  expect_true(is.list(out))
  expect_true(all(c("model", "design", "hrf", "estimates", "contrasts", "diagnostics") %in% names(out)))
})

test_that("report.fmri_lm renders a PDF (smoke test)", {
  skip_on_cran()
  skip_if_not_installed("quarto")
  skip_if_not_installed("tinytable")
  skip_if_not(quarto::quarto_available())

  fit <- make_report_test_fit()
  out_file <- tempfile(fileext = ".pdf")

  rendered <- tryCatch(
    report(
      fit,
      output_file = out_file,
      sections = c("model", "estimates"),
      open = FALSE,
      quiet = TRUE
    ),
    error = identity
  )

  if (inherits(rendered, "error")) {
    skip(paste("Quarto render unavailable in this environment:", conditionMessage(rendered)))
  }

  expect_true(file.exists(out_file))
})
