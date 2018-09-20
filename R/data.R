####    Documentation for data sets in the sncp R package    ####

#' Boundry pixels of the lung in an example axial slice from a CT scan.
#'
#' A dataframe that contains the x and y coordinates of the boundary pixels for a single lung in an axial CT slice
#'
#' @format A data frame with 614 rows and 2 variables:
#' \describe{
#'   \item{x}{x coordinate}
#'   \item{y}{y coordinate}
#'   ...
#' }
"boundry_total"


#' All pixels of lung tissue in an example axial slice from a CT scan.
#'
#' A dataframe that contains the x and y coordinates for all of the pixels in a single lung in one axial CT slice
#'
#' @format A data frame with 25,473 rows and 2 variables:
#' \describe{
#'   \item{x}{x coordinate}
#'   \item{y}{y coordinate}
#'   ...
#' }
"lung_data"

#' Simulated point pattern and associated model parameters for fitting SNCP to
#'
#' A named list that contains a simulated 2D point pattern of clusters (i.e. emphysema airholes), along with the true model parameters that went into generating the simulated data set
#'
#' @format A list with 11 components:
#' \describe{
#'   \item{centers}{Data frame with x and y coordinates for true cluster center locations}
#'   \item{logAlphas}{Vector of true values of log(alpha) associated with each center (i.e. log of expected number of pixels in each cluster)}
#'   \item{sigmaList}{List of true covariance matrices sigma associated with each center (i.e. covariance matrix of multivariate normal dispersion kernal centered at each true cluster center)}
#'   \item{nPointVec}{Vector that has the number of simulated points that go with each cluster center.  The last value is the number of noise/random scatter points}
#'   \item{pointPattern}{Data frame with x and y coordinates for simulated point pattern of diseased tissue}
#'   \item{beta}{Intensity value for homogeneous Cox Process that generated the random scatter points}
#'   \item{meanMuAlpha}{Mean for normal distribution from which muAlpha was drawn for this data set}
#'   \item{varMuAlpha}{Variance for normal distribution from which muAlpha was drawn for this data set}
#'   \item{muAlpha}{Actual value for muAlpha used in this simulated data set (i.e. mean of normal distribution from which the log(alpha)'s were drawn from)}
#'   \item{varAlpha}{Variance of normal distribution from which the log(alpha)'s were drawn from}
#'
#' }
"sim_data"
