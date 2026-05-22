#' Fit a Bayesian factor model with GLT process prior
#'
#' Fits a Bayesian Gaussian factor model under a Generalized Lower Triangular
#' (GLT) process prior using an MCMC sampler implemented in C++ via Rcpp.
#' This is the main user-facing entry point of the package.
#'
#' The model is:
#' \deqn{y_t = \Lambda \eta_t + \epsilon_t, \quad \epsilon_t \sim N_m(0, \Sigma)}
#' where \eqn{\Lambda} is a sparse lower-triangular loading matrix controlled
#' by a binary allocation matrix \eqn{\Delta}, and
#' \eqn{\Sigma = \mathrm{diag}(\sigma^2_1, \ldots, \sigma^2_m)}.
#'
#' Notation and dimensions:
#' \itemize{
#'   \item \code{y} is a numeric matrix of dimension \eqn{T \times m}.
#'   \item \code{Lambda} is a loading matrix of dimension \eqn{m \times H}.
#'   \item \code{Delta} is a binary allocation matrix of dimension \eqn{m \times H}.
#'   \item \code{Eta} is a factor score matrix of dimension \eqn{H \times T}.
#' }
#'
#' @param y A numeric \eqn{T \times m} data matrix with observations in rows
#'   and variables in columns. Missing values are not supported.
#' @param mcmc A named list of MCMC settings. Recognised fields:
#'   \describe{
#'     \item{`niter`}{Positive integer. Total number of MCMC iterations.
#'       Default \code{5000}.}
#'     \item{`nburn`}{Non-negative integer. Number of burn-in iterations to
#'       discard. Must be strictly less than \code{niter}. Default \code{2500}.}
#'     \item{`thin`}{Positive integer. Thinning interval. Default \code{1}.}
#'   }
#' @param Hmax A positive integer giving a conservative upper bound for the
#'   number of active factors. Default \code{ncol(y)}.
#' @param prior A named list of prior hyperparameters. Recognised fields:
#'   \describe{
#'     \item{`a_nu`}{Positive numeric. First shape of the Beta prior on the
#'       global inclusion probability \eqn{\nu}. Default \code{1}.}
#'     \item{`b_nu`}{Positive numeric. Second shape of the Beta prior on
#'       \eqn{\nu}. Default \code{1}.}
#'     \item{`alpha`}{Positive numeric. Shape of the Gamma prior on pivot
#'       loadings. Default \code{1}.}
#'     \item{`beta`}{Positive numeric. Rate of the Gamma prior on pivot
#'       loadings. Default \code{1}.}
#'     \item{`kappa`}{Positive numeric. Precision of the normal prior on
#'       non-pivot loadings. Default \code{1}.}
#'     \item{`a_sigma`}{Positive numeric. Shape of the inverse-Gamma prior on
#'       idiosyncratic variances \eqn{\sigma^2_i}. Default \code{2}.}
#'     \item{`b_sigma`}{Positive numeric. Scale of the inverse-Gamma prior on
#'       \eqn{\sigma^2_i}. Default \code{2}.}
#'   }
#' @param init A named list of initial values for the sampler. Recognised
#'   fields:
#'   \describe{
#'     \item{`H`}{Single non-negative integer. Initial number of active factors.}
#'     \item{`Lambda`}{Numeric \eqn{m \times H} loading matrix.}
#'     \item{`Delta`}{Integer \eqn{m \times H} binary allocation matrix with
#'       entries in \eqn{\{0,1\}}.}
#'     \item{`Eta`}{Numeric \eqn{H \times T} factor score matrix.}
#'     \item{`ell`}{Integer vector of length \eqn{H} with pivot row indices.}
#'     \item{`tau`}{Numeric vector of length \eqn{H} with initial
#'       column-specific inclusion probabilities.}
#'     \item{`nu`}{Single numeric in \eqn{(0,1)}. Initial global sparsity
#'       parameter.}
#'     \item{`sigma2`}{Positive numeric vector of length \eqn{m}. Initial
#'       idiosyncratic variances.}
#'   }
#' @param control A named list of algorithmic and storage options. Recognised
#'   fields:
#'   \describe{
#'     \item{`store_draws`}{Logical. Whether to store posterior draws of
#'       \code{Lambda}, \code{Delta}, \code{sigma2}, and \code{H}.
#'       Default \code{TRUE}.}
#'     \item{`store_eta`}{Logical. Whether to additionally store draws of the
#'       factor score matrix \code{Eta}. Default \code{FALSE}.}
#'     \item{`print_every`}{Positive integer. Print progress every this many
#'       iterations (only when \code{verbose = TRUE}). Default \code{100}.}
#'     \item{`seed`}{Single integer or \code{NULL}. Random seed for the C++
#'       sampler. Default \code{NULL}.}
#'   }
#' @param fixed A named list of logical flags indicating which MCMC updates to
#'   skip. For example, \code{list(H = TRUE)} keeps the number of factors fixed
#'   at its initial value throughout sampling.
#' @param verbose Logical. Whether to print sampling progress to the console.
#'   Defaults to \code{interactive()}.
#'
#' @return An object of class \code{"gltfit"}, a named list with components:
#' \describe{
#'   \item{`draws`}{A named list of posterior draws (if \code{control$store_draws = TRUE}):
#'     \code{Lambda} (\eqn{m \times H_{\max} \times n_{\text{save}}}),
#'     \code{Delta} (\eqn{m \times H_{\max} \times n_{\text{save}}}),
#'     \code{sigma2} (\eqn{m \times n_{\text{save}}}),
#'     \code{H} (\eqn{n_{\text{save}}}), and optionally
#'     \code{Eta} (\eqn{H_{\max} \times T \times n_{\text{save}}}).}
#'   \item{`meta`}{A named list with \code{nsave}, \code{Hmax}, \code{T}, \code{m}.}
#'   \item{`input`}{A named list echoing the processed \code{mcmc}, \code{prior},
#'     \code{control}, and \code{fixed} arguments.}
#'   \item{`call`}{The matched call.}
#' }
#'
#' @seealso [get_variance()], [get_pivots()], [plot_real_matrix()]
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' Y <- matrix(rnorm(200 * 10), nrow = 200, ncol = 10)
#'
#' fit <- gltfa(
#'   y       = Y,
#'   mcmc    = list(niter = 2000, nburn = 1000, thin = 1),
#'   Hmax    = 5,
#'   prior   = list(a_sigma = 2, b_sigma = 1),
#'   init    = list(H = 2),
#'   control = list(store_draws = TRUE, seed = 1L, random_scan = TRUE),
#'   verbose = FALSE
#' )
#' }
#'
#' @export
gltfa <- function(
  y,
  mcmc = list(),
  Hmax = ncol(y),
  prior = list(),
  init = list(),
  control = list(),
  fixed = list(),
  verbose = interactive()
) {
  y <- as.matrix(y)

  if (!is.numeric(y)) {
    stop("`y` must be a numeric matrix.", call. = FALSE)
  }

  if (anyNA(y)) {
    stop("Missing values in `y` are not supported in the current version.", call. = FALSE)
  }

  T <- nrow(y)
  m <- ncol(y)

  if (is.null(T) || T < 2L) {
    stop("`y` must have at least 2 rows.", call. = FALSE)
  }

  if (is.null(m) || m < 1L) {
    stop("`y` must have at least 1 column.", call. = FALSE)
  }

  if (!is.list(mcmc)) {
    stop("`mcmc` must be a named list.", call. = FALSE)
  }
  if (!is.list(prior)) {
    stop("`prior` must be a named list.", call. = FALSE)
  }
  if (!is.list(init)) {
    stop("`init` must be a named list.", call. = FALSE)
  }
  if (!is.list(control)) {
    stop("`control` must be a named list.", call. = FALSE)
  }
  if (!is.list(fixed)) {
    stop("`fixed` must be a named list.", call. = FALSE)
  }

  # if H is fixed, nu is fixed as well
  if (isTRUE(fixed$H)) {
    fixed$nu <- TRUE
  }

  # if pivots are fixed, H and nu are fixed as well
  if (isTRUE(fixed$pivots)) {
    fixed$H <- TRUE
    fixed$nu <- TRUE
  }

  # if Delta is fixed, pivots, H, tau, and nu are fixed as well
  if (isTRUE(fixed$Delta)) {
    fixed$pivots <- TRUE
    fixed$H <- TRUE
    fixed$nu <- TRUE
    fixed$tau <- TRUE
  }
  
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be a single TRUE or FALSE.", call. = FALSE)
  }

  mcmc <- gltfa_process_mcmc(mcmc)
  prior <- gltfa_process_prior(prior)
  control <- gltfa_process_control(control, verbose = verbose)

  init <- gltfa_process_init(
    init = init,
    T = T,
    m = m,
    Hmax = Hmax,
    y = y
  )

   model <- list(
    T = T,
    m = m,
    Hmax = Hmax
  )

  fit <- gltfa_cpp(
    y = y,
    mcmc = mcmc,
    model = model,
    prior = prior,
    init = init,
    control = control,
    fixed = fixed
  )

  fit$call <- match.call()
  fit$input <- list(
    mcmc = mcmc,
    model = model,
    prior = prior,
    control = control,
    fixed = fixed
  )

  class(fit) <- "gltfit"
  fit

}