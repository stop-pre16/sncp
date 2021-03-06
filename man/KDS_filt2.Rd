% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kds_filt_functions.R
\name{KDS_filt2}
\alias{KDS_filt2}
\title{Function to do KDS-Filt using random draws for number of points}
\usage{
KDS_filt2(
  sim_pattern,
  xwin,
  ywin,
  bw_method_c = "bw.ppl",
  thresh = 0.5,
  n_sim_data = 20,
  edge_correct = T
)
}
\arguments{
\item{sim_pattern}{Data frame with standard long format CT info for just diseased tissue}

\item{xwin}{Range of x values}

\item{ywin}{Range of y values}

\item{bw_method_c}{Either numeric smoothing bandwidth, or character string for methods (bw_ppl, bw_diggle, bw_CvL)}

\item{thresh}{Quantile to use for filtering noise}

\item{n_sim_data}{Number of simulation datasets to use for constructing null distribution}

\item{edge_correct}{Logical variable denoting if kernel intensity estimates should be edge corrected}
}
\value{
Returns a data frame with x and y coordinates, and the feature or noise classification for each point
}
\description{
Wrapper function that does KDS-Filt, then gets cluster assignments
}
