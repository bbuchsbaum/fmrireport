#' Print a cluster_table Object
#'
#' CLI-formatted hierarchical output showing clusters with sub-peaks indented.
#'
#' @param x A \code{cluster_table} object.
#' @param max_clusters Maximum number of clusters to display.
#' @param max_sub_peaks Maximum sub-peaks per cluster to display.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}.
#' @export
print.cluster_table <- function(x, max_clusters = 20L, max_sub_peaks = 3L, ...) {
  if (x$n_clusters == 0L) {
    cat("cluster_table: 0 clusters (threshold =", x$threshold, ")\n")
    return(invisible(x))
  }

  mni <- .is_mni_space(x$coord_space)
  coord_prefix <- if (mni) "MNI " else ""

  cat(sprintf("cluster_table: %d cluster(s)  |  threshold = %.2f (%s)  |  space = %s",
              x$n_clusters, x$threshold, x$stat_type, x$coord_space))
  if (!is.na(x$atlas_name)) cat("  |  atlas =", x$atlas_name)
  cat("\n\n")

  clusters <- x$clusters
  n_show <- min(nrow(clusters), max_clusters)

  for (i in seq_len(n_show)) {
    row <- clusters[i, ]
    label_str <- if (!is.na(row$label)) paste0("  [", row$label, "]") else ""
    hemi_str <- if (!is.na(row$hemi)) paste0(" (", row$hemi, ")") else ""

    cat(sprintf("  Cluster %d  (k = %d, %.0f mm3)\n",
                row$cluster_id, row$k, row$volume_mm3))
    cat(sprintf("    Peak: %s%s  %sx = %.1f  %sy = %.1f  %sz = %.1f  stat = %.3f",
                label_str, hemi_str,
                coord_prefix, row$peak_x,
                coord_prefix, row$peak_y,
                coord_prefix, row$peak_z,
                row$peak_stat))
    if (!is.na(row$peak_p)) cat(sprintf("  p = %.4g", row$peak_p))
    cat("\n")

    ## Sub-peaks
    sub <- x$peaks[x$peaks$cluster_id == row$cluster_id, , drop = FALSE]
    if (nrow(sub) > 1) {
      sub <- sub[sub$peak_rank > 1, , drop = FALSE]
      n_sub <- min(nrow(sub), max_sub_peaks)
      for (j in seq_len(n_sub)) {
        sr <- sub[j, ]
        sl <- if (!is.na(sr$label)) paste0("  [", sr$label, "]") else ""
        cat(sprintf("      sub-peak %d:%s  %sx = %.1f  %sy = %.1f  %sz = %.1f  stat = %.3f",
                    sr$peak_rank, sl,
                    coord_prefix, sr$x, coord_prefix, sr$y, coord_prefix, sr$z,
                    sr$stat))
        if (!is.na(sr$p)) cat(sprintf("  p = %.4g", sr$p))
        cat("\n")
      }
    }
    cat("\n")
  }

  if (n_show < nrow(clusters)) {
    cat(sprintf("  ... and %d more cluster(s)\n", nrow(clusters) - n_show))
  }

  invisible(x)
}

#' Coerce cluster_table to a Flat data.frame
#'
#' Produces a flat data.frame suitable for table rendering. Cluster rows have
#' all columns populated; sub-peak rows have \code{k} and \code{volume_mm3}
#' set to \code{NA}.
#'
#' @param x A \code{cluster_table} object.
#' @param max_clusters Maximum clusters to include.
#' @param max_sub_peaks Maximum sub-peaks per cluster.
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param ... Ignored.
#' @return A data.frame.
#' @export
as.data.frame.cluster_table <- function(x, row.names = NULL, optional = FALSE,
                                        max_clusters = 30L,
                                        max_sub_peaks = 3L, ...) {
  if (x$n_clusters == 0L) {
    return(.empty_cluster_df(x$coord_space))
  }

  mni <- .is_mni_space(x$coord_space)
  xn <- if (mni) "MNI_x" else "x"
  yn <- if (mni) "MNI_y" else "y"
  zn <- if (mni) "MNI_z" else "z"

  clusters <- x$clusters
  n_show <- min(nrow(clusters), max_clusters)
  rows <- list()

  for (i in seq_len(n_show)) {
    cl <- clusters[i, ]
    row <- data.frame(
      Cluster = cl$cluster_id,
      k = cl$k,
      Vol_mm3 = cl$volume_mm3,
      Stat = cl$peak_stat,
      p = cl$peak_p,
      Region = if (is.na(cl$label)) "" else cl$label,
      Hemi = if (is.na(cl$hemi)) "" else cl$hemi,
      stringsAsFactors = FALSE
    )
    row[[xn]] <- cl$peak_x
    row[[yn]] <- cl$peak_y
    row[[zn]] <- cl$peak_z
    rows[[length(rows) + 1L]] <- row

    ## Sub-peaks
    sub <- x$peaks[x$peaks$cluster_id == cl$cluster_id & x$peaks$peak_rank > 1,
                   , drop = FALSE]
    n_sub <- min(nrow(sub), max_sub_peaks)
    if (n_sub > 0) {
      for (j in seq_len(n_sub)) {
        sr <- sub[j, ]
        srow <- data.frame(
          Cluster = cl$cluster_id,
          k = NA_integer_,
          Vol_mm3 = NA_real_,
          Stat = sr$stat,
          p = sr$p,
          Region = if (is.na(sr$label)) "" else sr$label,
          Hemi = if (is.na(sr$hemi)) "" else sr$hemi,
          stringsAsFactors = FALSE
        )
        srow[[xn]] <- sr$x
        srow[[yn]] <- sr$y
        srow[[zn]] <- sr$z
        rows[[length(rows) + 1L]] <- srow
      }
    }
  }

  out <- do.call(rbind, rows)
  ## Reorder columns for readability
  coord_cols <- c(xn, yn, zn)
  front <- c("Cluster", "k", "Vol_mm3")
  back <- c("Stat", "p", "Region", "Hemi")
  col_order <- c(front, coord_cols, back)
  col_order <- col_order[col_order %in% names(out)]
  out <- out[, col_order, drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Format a cluster_table as a tinytable
#'
#' Returns a \code{tinytable::tt()} object with appropriate caption and
#' formatting for Quarto/Typst rendering.
#'
#' @param x A \code{cluster_table} object.
#' @param max_clusters Maximum clusters to include.
#' @param max_sub_peaks Maximum sub-peaks per cluster.
#' @param digits Number of digits for rounding numeric columns.
#' @param ... Passed to \code{tinytable::tt()}.
#' @return A \code{tinytable} object.
#' @export
format_cluster_tt <- function(x, max_clusters = 30L, max_sub_peaks = 3L,
                              digits = 3L, ...) {
  if (!requireNamespace("tinytable", quietly = TRUE)) {
    stop("Package 'tinytable' is required for format_cluster_tt().", call. = FALSE)
  }

  df <- as.data.frame(x, max_clusters = max_clusters,
                      max_sub_peaks = max_sub_peaks)

  ## Round numeric columns
  for (nm in names(df)) {
    if (is.numeric(df[[nm]])) {
      df[[nm]] <- round(df[[nm]], digits = digits)
    }
  }

  ## Build caption
  parts <- character(0)
  parts <- c(parts, paste0("Threshold: ", x$stat_type, " = ", x$threshold))
  if (!is.null(x$df)) {
    parts <- c(parts, paste0("df = ", paste(x$df, collapse = ", ")))
  }
  if (!is.na(x$coord_space) && nzchar(x$coord_space)) {
    parts <- c(parts, paste0("Space: ", x$coord_space))
  }
  if (!is.na(x$atlas_name)) {
    parts <- c(parts, paste0("Atlas: ", x$atlas_name))
  }
  caption <- paste(parts, collapse = "  |  ")

  tinytable::tt(df, caption = caption, ...)
}


# ---- Internal helpers --------------------------------------------------------

#' @keywords internal
#' @noRd
.is_mni_space <- function(coord_space) {
  !is.null(coord_space) && !is.na(coord_space) &&
    grepl("^MNI", coord_space, ignore.case = TRUE)
}

#' @keywords internal
#' @noRd
.empty_cluster_df <- function(coord_space) {
  mni <- .is_mni_space(coord_space)
  xn <- if (mni) "MNI_x" else "x"
  yn <- if (mni) "MNI_y" else "y"
  zn <- if (mni) "MNI_z" else "z"

  df <- data.frame(
    Cluster = integer(0), k = integer(0), Vol_mm3 = numeric(0),
    stringsAsFactors = FALSE
  )
  df[[xn]] <- numeric(0)
  df[[yn]] <- numeric(0)
  df[[zn]] <- numeric(0)
  df$Stat <- numeric(0)
  df$p <- numeric(0)
  df$Region <- character(0)
  df$Hemi <- character(0)
  df
}
