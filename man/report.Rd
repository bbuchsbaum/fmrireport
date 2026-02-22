\name{report}
\alias{report}
\alias{report.default}
\alias{report.fmri_lm}
\title{Render Analysis Reports}
\usage{
report(x, ...)

\method{report}{default}(x, ...)

\method{report}{fmri_lm}(
  x,
  output_file = "fmri_lm_report.pdf",
  title = "fMRI GLM Analysis Report",
  author = NULL,
  sections = c("model", "design", "hrf", "estimates", "contrasts", "diagnostics"),
  brain_map_stat = c("tstat", "estimate", "prob"),
  slice_axis = 3L,
  n_slices = 9L,
  threshold = NULL,
  bg_vol = NULL,
  atlas = NULL,
  cluster_thresh = 3,
  min_cluster_size = 10L,
  max_peaks = 15L,
  local_maxima_dist = 15,
  max_sub_peaks = 3L,
  open = interactive(),
  quiet = TRUE,
  ...
)
}
\arguments{
\item{x}{Object to report on.}

\item{...}{Additional arguments passed to methods.}

\item{output_file}{Output PDF path.}

\item{title}{Report title.}

\item{author}{Optional author string.}

\item{sections}{Sections to include.}

\item{brain_map_stat}{Statistic for contrast maps: \code{"tstat"}, \code{"estimate"}, or \code{"prob"}.}

\item{slice_axis}{Axis for montage slicing (1, 2, or 3).}

\item{n_slices}{Number of slices to render in contrast maps.}

\item{threshold}{Optional hard threshold for contrast map rendering.}

\item{bg_vol}{Optional background \code{NeuroVol} for overlay plots.}

\item{atlas}{Optional atlas object for peak labeling.}

\item{cluster_thresh}{Absolute statistic threshold for peak clustering.}

\item{min_cluster_size}{Minimum cluster size (voxels) in peak table.}

\item{max_peaks}{Maximum peaks per contrast table.}

\item{local_maxima_dist}{Minimum distance (mm) between local maxima within a cluster for sub-peak detection. Default 15.}

\item{max_sub_peaks}{Maximum number of sub-peaks per cluster in the activation table. Default 3.}

\item{open}{If \code{TRUE}, open the generated PDF after rendering.}

\item{quiet}{If \code{TRUE}, suppress Quarto render logs.}
}
\value{
Method-specific report output. For \code{fmri_lm} fits, returns the output PDF
path invisibly.
}
\description{
Generic report interface and \code{fmri_lm} implementation that renders a
publication-style PDF report through Quarto.
}
