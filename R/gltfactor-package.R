#' gltfactor: Bayesian Factor Models with GLT Process Priors
#'
#' Fits a Bayesian Gaussian factor model under a Generalized Lower Triangular
#' (GLT) process prior using an MCMC sampler implemented in C++ via Rcpp.
#' The package provides tools for model fitting, posterior summarisation,
#' identification post-processing, and visualisation of loading matrices and
#' variance decompositions.
#'
#' The main entry point is [gltfa()], which accepts a data matrix `y` (T x m)
#' and returns a `"gltfit"` object containing posterior draws and metadata.
#'
#' @section Main functions:
#' \describe{
#'   \item{[gltfa()]}{Fit a GLT factor model via MCMC.}
#'   \item{[get_variance()]}{Compute posterior draws of the factor and total
#'     covariance matrices.}
#'   \item{[get_pivots()]}{Extract pivot row indices from a binary allocation
#'     matrix `Delta`.}
#'   \item{[plot_real_matrix()]}{Visualise a numeric matrix as a colour image
#'     using a diverging palette.}
#' }
#'
#' @useDynLib glt, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom sparvaride sparvaride
"_PACKAGE"