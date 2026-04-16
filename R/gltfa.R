#' Fit a Bayesian factor model with GLT process prior
#'
#' Fits a Bayesian Gaussian factor model under a GLT process prior using MCMC.
#' This is the main user-facing entry point of the package and will call
#' an internal Rcpp sampler.
#'
#' Notation and dimensions:
#' \itemize{
#'   \item `y` is a numeric matrix of dimension `T x m`.
#'   \item `Lambda` is a loading matrix of dimension `m x H`.
#'   \item `Delta` is a binary allocation matrix of dimension `m x H`.
#'   \item `Eta` is a factor score matrix of dimension `H x T`.
#' }
#'
#' @param y A numeric data matrix with observations in rows and variables in columns. `NA`` not admitted for now
#' @param mcmc A named list containing MCMC settings:
#'   \describe{
#'     \item{`niter`}{total number of MCMC iterations.}
#'     \item{`nburn`}{number of burn-in iterations.}
#'     \item{`thin`}{thinning interval.}
#'   }
#' @param Hmax A conservative upper bound for the number of factors `H`.
#'   The default is `2 * ncol(y)`.
#' @param prior A named list of prior hyperparameters.
#' @param init A named list of initial values for the sampler.
#' @param control A named list of algorithmic and storage options.
#' @param fixed A named list of logical flags indicating which updates to skip (e.g. list(H=TRUE)).
#' @param verbose Logical; whether to print sampling progress.
#'
#' @return An object of class `"gltfit"` containing processed inputs,
#'   model metadata, and later the output of the internal sampler.
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
    Hmax = Hmax
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