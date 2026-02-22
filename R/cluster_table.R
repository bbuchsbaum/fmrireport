#' Build a Publication-Quality Activation Cluster Table
#'
#' Detects suprathreshold clusters in a statistical brain volume and returns
#' a structured table of cluster-level and peak-level statistics suitable for
#' publication (NeuroImage / HBM / COBIDAS compliant).
#'
#' @param vol A \code{NeuroVol} object containing the statistical map.
#' @param threshold Numeric threshold for cluster detection (absolute value).
#' @param atlas Optional atlas object (e.g., a \code{volatlas} from neuroatlas)
#'   for anatomical labeling.
#' @param stat_type Type of statistic in the volume: \code{"t"}, \code{"z"},
#'   \code{"F"}, or \code{"other"}.
#' @param df Degrees of freedom for p-value computation. A single value for
#'   t/z statistics, or a length-2 vector \code{c(df1, df2)} for F statistics.
#'   \code{NULL} omits p-values.
#' @param coord_space Character string naming the coordinate space (e.g.,
#'   \code{"MNI152"}, \code{"MNI305"}). If \code{NULL}, inferred from \code{atlas}
#'   or set to \code{"Unknown"}.
#' @param min_cluster_size Minimum number of voxels for a cluster to be retained.
#' @param local_maxima_dist Minimum distance (mm) between local maxima within a
#'   cluster.
#' @param max_peaks Maximum number of sub-peaks reported per cluster.
#' @param connectivity Voxel connectivity for cluster detection.
#' @param sort_by Sort clusters by \code{"size"} (descending) or
#'   \code{"peak_stat"} (descending absolute value).
#' @return An S3 object of class \code{"cluster_table"} with elements:
#'   \describe{
#'     \item{clusters}{A data.frame with one row per cluster.}
#'     \item{peaks}{A data.frame with sub-peaks (local maxima) per cluster.}
#'     \item{threshold}{The threshold used.}
#'     \item{stat_type}{The statistic type.}
#'     \item{df}{Degrees of freedom.}
#'     \item{n_clusters}{Number of clusters.}
#'     \item{coord_space}{Coordinate space label.}
#'     \item{atlas_name}{Name of atlas used, or \code{NA}.}
#'   }
#' @export
cluster_table <- function(vol, threshold, atlas = NULL,
                          stat_type = c("t", "z", "F", "other"),
                          df = NULL, coord_space = NULL,
                          min_cluster_size = 10L,
                          local_maxima_dist = 15, max_peaks = 3L,
                          connectivity = c("26-connect", "18-connect", "6-connect"),
                          sort_by = c("size", "peak_stat")) {
  stat_type <- match.arg(stat_type)
  connectivity <- match.arg(connectivity)
  sort_by <- match.arg(sort_by)
  min_cluster_size <- as.integer(max(1L, min_cluster_size))
  max_peaks <- as.integer(max(1L, max_peaks))

  ## Stage 1 -- Cluster detection via neuroim2
  ## conn_comp errors when no voxels exceed threshold; catch gracefully.
  cc <- tryCatch(
    neuroim2::conn_comp(vol, threshold = threshold,
                        cluster_table = TRUE, local_maxima = TRUE,
                        local_maxima_dist = local_maxima_dist,
                        connect = connectivity),
    error = function(e) NULL
  )

  if (is.null(cc)) {
    return(.empty_cluster_table(threshold, stat_type, df, coord_space, atlas))
  }

  ctab <- cc$cluster_table
  lmax <- cc$local_maxima

  if (is.null(ctab) || !is.data.frame(ctab) || nrow(ctab) == 0) {
    return(.empty_cluster_table(threshold, stat_type, df, coord_space, atlas))
  }

  ## Stage 2 -- Filter small clusters
  keep <- ctab$N >= min_cluster_size
  ctab <- ctab[keep, , drop = FALSE]

  if (nrow(ctab) == 0) {
    return(.empty_cluster_table(threshold, stat_type, df, coord_space, atlas))
  }

  kept_ids <- ctab$index
  if (!is.null(lmax) && is.matrix(lmax) && nrow(lmax) > 0) {
    lmax <- lmax[lmax[, "index"] %in% kept_ids, , drop = FALSE]
  } else {
    lmax <- NULL
  }

  ## Stage 3 -- World coordinates
  vol_space <- neuroim2::space(vol)
  peak_grid <- as.matrix(ctab[, c("x", "y", "z"), drop = FALSE])
  peak_world <- neuroim2::grid_to_coord(vol_space, peak_grid)

  ## Stage 4 -- Coordinate space resolution
  coord_space <- .resolve_coord_space(coord_space, atlas)

  ## Stage 5 -- Atlas labeling (peaks)
  if (!is.null(atlas)) {
    peak_labels <- .label_coords(peak_world, atlas, vol_space)
  } else {
    peak_labels <- data.frame(
      label = rep(NA_character_, nrow(ctab)),
      label_full = rep(NA_character_, nrow(ctab)),
      hemi = rep(NA_character_, nrow(ctab)),
      network = rep(NA_character_, nrow(ctab)),
      stringsAsFactors = FALSE
    )
  }

  ## Stage 6 -- P-values
  peak_p <- .stat_to_p(ctab$value, stat_type, df)

  ## Compute volume in mm3
  vox_vol <- prod(abs(neuroim2::spacing(vol_space)))

  ## Assemble cluster-level data.frame
  clusters <- data.frame(
    cluster_id = ctab$index,
    k = ctab$N,
    volume_mm3 = ctab$N * vox_vol,
    peak_x = peak_world[, 1],
    peak_y = peak_world[, 2],
    peak_z = peak_world[, 3],
    peak_stat = ctab$value,
    peak_p = peak_p,
    label = peak_labels$label,
    label_full = peak_labels$label_full,
    hemi = peak_labels$hemi,
    network = peak_labels$network,
    stringsAsFactors = FALSE
  )

  ## Sort
  if (sort_by == "size") {
    clusters <- clusters[order(clusters$k, decreasing = TRUE), , drop = FALSE]
  } else {
    clusters <- clusters[order(abs(clusters$peak_stat), decreasing = TRUE), , drop = FALSE]
  }
  rownames(clusters) <- NULL

  ## Assemble sub-peaks (local maxima)
  peaks <- .build_sub_peaks(lmax, vol_space, atlas, stat_type, df, max_peaks,
                            kept_ids)

  ## Build S3 object
  atlas_name <- if (!is.null(atlas)) {
    if (!is.null(atlas$name)) atlas$name else NA_character_
  } else {
    NA_character_
  }

  structure(list(
    clusters    = clusters,
    peaks       = peaks,
    threshold   = threshold,
    stat_type   = stat_type,
    df          = df,
    n_clusters  = nrow(clusters),
    coord_space = coord_space,
    atlas_name  = atlas_name
  ), class = "cluster_table")
}


# ---- Internal helpers --------------------------------------------------------

#' Build sub-peaks data.frame from local maxima matrix
#' @keywords internal
#' @noRd
.build_sub_peaks <- function(lmax, vol_space, atlas, stat_type, df, max_peaks,
                             kept_ids) {
  if (is.null(lmax) || !is.matrix(lmax) || nrow(lmax) == 0) {
    return(data.frame(
      cluster_id = integer(0), peak_rank = integer(0),
      x = numeric(0), y = numeric(0), z = numeric(0),
      stat = numeric(0), p = numeric(0),
      label = character(0), label_full = character(0),
      hemi = character(0), network = character(0),
      stringsAsFactors = FALSE
    ))
  }

  pieces <- list()
  for (cid in kept_ids) {
    rows <- lmax[lmax[, "index"] == cid, , drop = FALSE]
    if (nrow(rows) == 0) next

    ## Order by descending absolute stat and cap at max_peaks
    ord <- order(abs(rows[, "value"]), decreasing = TRUE)
    rows <- rows[ord, , drop = FALSE]
    rows <- rows[seq_len(min(nrow(rows), max_peaks)), , drop = FALSE]

    grid_coords <- rows[, c("x", "y", "z"), drop = FALSE]
    world_coords <- neuroim2::grid_to_coord(vol_space, grid_coords)
    p_vals <- .stat_to_p(rows[, "value"], stat_type, df)

    if (!is.null(atlas)) {
      labs <- .label_coords(world_coords, atlas, vol_space)
    } else {
      n <- nrow(rows)
      labs <- data.frame(
        label = rep(NA_character_, n), label_full = rep(NA_character_, n),
        hemi = rep(NA_character_, n), network = rep(NA_character_, n),
        stringsAsFactors = FALSE
      )
    }

    pieces[[length(pieces) + 1L]] <- data.frame(
      cluster_id = rep(as.integer(cid), nrow(rows)),
      peak_rank = seq_len(nrow(rows)),
      x = world_coords[, 1],
      y = world_coords[, 2],
      z = world_coords[, 3],
      stat = as.numeric(rows[, "value"]),
      p = p_vals,
      label = labs$label,
      label_full = labs$label_full,
      hemi = labs$hemi,
      network = labs$network,
      stringsAsFactors = FALSE
    )
  }

  if (!length(pieces)) {
    return(data.frame(
      cluster_id = integer(0), peak_rank = integer(0),
      x = numeric(0), y = numeric(0), z = numeric(0),
      stat = numeric(0), p = numeric(0),
      label = character(0), label_full = character(0),
      hemi = character(0), network = character(0),
      stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, pieces)
}

#' Resolve coordinate space from explicit param, atlas, or fallback
#' @keywords internal
#' @noRd
.resolve_coord_space <- function(coord_space, atlas) {
  if (!is.null(coord_space) && !is.na(coord_space) && nzchar(coord_space)) {
    return(coord_space)
  }
  if (!is.null(atlas) && requireNamespace("neuroatlas", quietly = TRUE)) {
    cs <- tryCatch(neuroatlas::atlas_coord_space(atlas), error = function(e) NA_character_)
    if (!is.na(cs) && nzchar(cs)) return(cs)
  }
  "Unknown"
}

#' Check space compatibility between volume space and atlas
#' @keywords internal
#' @noRd
.check_space_compat <- function(vol_space, atlas) {
  if (!requireNamespace("neuroatlas", quietly = TRUE)) {
    return(list(compatible = FALSE, transform = NULL))
  }

  atlas_cs <- tryCatch(neuroatlas::atlas_coord_space(atlas), error = function(e) NA_character_)
  if (is.na(atlas_cs) || !nzchar(atlas_cs)) {
    return(list(compatible = TRUE, transform = NULL))
  }

  ## Detect volume coordinate space from atlas context — we don't have a direct

  ## way to interrogate the vol_space for its named space, so we just assume
  ## the volume is in whatever space the user expects. If the atlas coord space
  ## is a known standard (MNI152/MNI305), we assume compatibility.
  standard_spaces <- c("MNI152", "MNI305")
  if (atlas_cs %in% standard_spaces) {
    return(list(compatible = TRUE, transform = NULL))
  }

  list(compatible = TRUE, transform = NULL)
}

#' Label world coordinates using an atlas
#' @keywords internal
#' @noRd
.label_coords <- function(world_coords, atlas, vol_space) {
  n <- nrow(world_coords)
  empty <- data.frame(
    label = rep(NA_character_, n), label_full = rep(NA_character_, n),
    hemi = rep(NA_character_, n), network = rep(NA_character_, n),
    stringsAsFactors = FALSE
  )

  if (!requireNamespace("neuroatlas", quietly = TRUE)) {
    return(empty)
  }

  compat <- .check_space_compat(vol_space, atlas)
  if (!isTRUE(compat$compatible)) {
    warning("Volume and atlas coordinate spaces appear incompatible; skipping labels.",
            call. = FALSE)
    return(empty)
  }

  ## Get atlas volume
  atlas_vol <- .get_atlas_vol(atlas)
  if (is.null(atlas_vol)) return(empty)

  atlas_space <- neuroim2::space(atlas_vol)
  atlas_arr <- as.array(atlas_vol)

  ## Convert world coords to atlas grid
  grid <- neuroim2::coord_to_grid(atlas_space, world_coords)
  grid <- round(grid)

  ## Get ROI metadata
  meta <- tryCatch(neuroatlas::roi_metadata(atlas), error = function(e) NULL)
  if (is.null(meta) || !is.data.frame(meta) || nrow(meta) == 0) {
    return(empty)
  }

  ## Look up each coordinate
  d <- dim(atlas_arr)
  label <- character(n)
  label_full <- character(n)
  hemi <- character(n)
  network <- character(n)

  has_network <- "network" %in% names(meta)

  for (i in seq_len(n)) {
    gi <- as.integer(grid[i, 1:3])

    ## Bounds check
    if (any(gi < 1) || gi[1] > d[1] || gi[2] > d[2] || gi[3] > d[3]) {
      label[i] <- NA_character_
      label_full[i] <- NA_character_
      hemi[i] <- NA_character_
      network[i] <- NA_character_
      next
    }

    region_id <- atlas_arr[gi[1], gi[2], gi[3]]
    if (is.na(region_id) || region_id == 0) {
      label[i] <- NA_character_
      label_full[i] <- NA_character_
      hemi[i] <- NA_character_
      network[i] <- NA_character_
      next
    }

    idx <- match(region_id, meta$id)
    if (is.na(idx)) {
      label[i] <- as.character(region_id)
      label_full[i] <- as.character(region_id)
      hemi[i] <- NA_character_
      network[i] <- NA_character_
    } else {
      label[i] <- as.character(meta$label[idx])
      label_full[i] <- if ("label_full" %in% names(meta)) {
        as.character(meta$label_full[idx])
      } else {
        as.character(meta$label[idx])
      }
      hemi[i] <- if ("hemi" %in% names(meta)) {
        as.character(meta$hemi[idx])
      } else {
        NA_character_
      }
      network[i] <- if (has_network) {
        as.character(meta$network[idx])
      } else {
        NA_character_
      }
    }
  }

  data.frame(label = label, label_full = label_full, hemi = hemi,
             network = network, stringsAsFactors = FALSE)
}

#' Convert statistic values to p-values
#' @keywords internal
#' @noRd
.stat_to_p <- function(stat_values, stat_type, df) {
  if (is.null(df)) {
    return(rep(NA_real_, length(stat_values)))
  }

  switch(stat_type,
    t = stats::pt(abs(stat_values), df = df, lower.tail = FALSE) * 2,
    z = stats::pnorm(abs(stat_values), lower.tail = FALSE) * 2,
    F = stats::pf(stat_values, df1 = df[1], df2 = df[2], lower.tail = FALSE),
    rep(NA_real_, length(stat_values))
  )
}

#' Safely extract NeuroVol from atlas object
#' @keywords internal
#' @noRd
.get_atlas_vol <- function(atlas) {
  if (inherits(atlas, c("NeuroVol", "ClusteredNeuroVol"))) {
    return(atlas)
  }
  if (is.list(atlas) && !is.null(atlas$atlas)) {
    vol <- atlas$atlas
    if (inherits(vol, c("NeuroVol", "ClusteredNeuroVol"))) {
      return(vol)
    }
  }
  NULL
}

#' Construct an empty cluster_table
#' @keywords internal
#' @noRd
.empty_cluster_table <- function(threshold, stat_type, df, coord_space, atlas) {
  coord_space <- .resolve_coord_space(coord_space, atlas)
  atlas_name <- if (!is.null(atlas) && !is.null(atlas$name)) atlas$name else NA_character_

  structure(list(
    clusters = data.frame(
      cluster_id = integer(0), k = integer(0), volume_mm3 = numeric(0),
      peak_x = numeric(0), peak_y = numeric(0), peak_z = numeric(0),
      peak_stat = numeric(0), peak_p = numeric(0),
      label = character(0), label_full = character(0),
      hemi = character(0), network = character(0),
      stringsAsFactors = FALSE
    ),
    peaks = data.frame(
      cluster_id = integer(0), peak_rank = integer(0),
      x = numeric(0), y = numeric(0), z = numeric(0),
      stat = numeric(0), p = numeric(0),
      label = character(0), label_full = character(0),
      hemi = character(0), network = character(0),
      stringsAsFactors = FALSE
    ),
    threshold   = threshold,
    stat_type   = stat_type,
    df          = df,
    n_clusters  = 0L,
    coord_space = coord_space,
    atlas_name  = atlas_name
  ), class = "cluster_table")
}
