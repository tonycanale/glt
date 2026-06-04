#' Validate and complete initial values for the MCMC sampler
#'
#' Processes the user-supplied `init` list, fills in defaults, checks
#' cross-dimensional consistency, and coerces types as required by
#' [gltfa_cpp()].
#'
#' @param init A named list of initial values. Recognised fields are:
#' \describe{
#'   \item{`H`}{Single nonnegative integer. Initial number of active factors.}
#'   \item{`Lambda`}{Numeric `m x H` loading matrix.}
#'   \item{`Delta`}{Integer `m x H` binary allocation matrix with entries in `{0,1}`.}
#'   \item{`Eta`}{Numeric `H x T` factor score matrix.}
#'   \item{`ell`}{Integer vector of length `H` with pivot row indices.}
#'   \item{`tau`}{Numeric vector of length `H` with initial inclusion probabilities.}
#'   \item{`nu`}{Single numeric in `(0, 1)`. Initial global sparsity parameter.}
#'   \item{`sigma2`}{Positive numeric vector of length `m`. Initial idiosyncratic variances.}
#' }
#' @param T Integer. Number of observations (`nrow(y)`).
#' @param m Integer. Number of variables (`ncol(y)`).
#' @param Hmax Integer. Maximum number of factors.
#'
#' @return A validated and coerced named list ready to be passed to
#'   [gltfa_cpp()].
#'
#' @keywords internal
#' @noRd
gltfa_process_init <- function(init, T, m, Hmax, y = NULL) {
  # ...existing code...
  init_defaults <- list(
    H = NULL,
    nu = NULL,
    ell = NULL,
    tau = NULL,
    Delta = NULL,
    Lambda = NULL,
    Eta = NULL,
    sigma2 = NULL
  )
  init <- utils::modifyList(init_defaults, init)

  H_candidates <- integer(0)

  if (!is.null(init$H)) {
    if (!is.numeric(init$H) || length(init$H) != 1L || is.na(init$H)) {
      stop("`init$H` must be NULL or a single nonnegative integer.", call. = FALSE)
    }
    init$H <- as.integer(init$H)
    if (init$H < 0L) {
      stop("`init$H` must be nonnegative.", call. = FALSE)
    }
    H_candidates <- c(H_candidates, init$H)
  }

  if (!is.null(init$Lambda)) {
    if (!is.matrix(init$Lambda) || !is.numeric(init$Lambda)) {
      stop("`init$Lambda` must be NULL or a numeric matrix.", call. = FALSE)
    }
    if (nrow(init$Lambda) != m) {
      stop("`init$Lambda` must have `m = ncol(y)` rows.", call. = FALSE)
    }
    H_candidates <- c(H_candidates, ncol(init$Lambda))
  }

  if (!is.null(init$Delta)) {
    if (!is.matrix(init$Delta)) {
      stop("`init$Delta` must be NULL or a matrix.", call. = FALSE)
    }
    if (nrow(init$Delta) != m) {
      stop("`init$Delta` must have `m = ncol(y)` rows.", call. = FALSE)
    }
    H_candidates <- c(H_candidates, ncol(init$Delta))
  }

  if (!is.null(init$Eta)) {
    if (!is.matrix(init$Eta) || !is.numeric(init$Eta)) {
      stop("`init$Eta` must be NULL or a numeric matrix.", call. = FALSE)
    }
    if (ncol(init$Eta) != T) {
      stop("`init$Eta` must have `T = nrow(y)` columns.", call. = FALSE)
    }
    H_candidates <- c(H_candidates, nrow(init$Eta))
  }

  if (!is.null(init$ell)) {
    if (!is.numeric(init$ell)) {
      stop("`init$ell` must be NULL or a numeric vector.", call. = FALSE)
    }
    H_candidates <- c(H_candidates, length(init$ell))
  }

  if (!is.null(init$tau)) {
    if (!is.numeric(init$tau)) {
      stop("`init$tau` must be NULL or a numeric vector.", call. = FALSE)
    }
    H_candidates <- c(H_candidates, length(init$tau))
  }

  H_candidates <- unique(as.integer(H_candidates))

  if (length(H_candidates) > 1L) {
    stop(
      paste0(
        "Inconsistent initial values: supplied objects imply different values of `H` (",
        paste(H_candidates, collapse = ", "),
        ")."
      ),
      call. = FALSE
    )
  }

  if (length(H_candidates) == 0L) {
    init$H <- 0L
  } else {
    init$H <- H_candidates
  }

  if (init$H > Hmax) {
    stop("The implied initial value of `H` exceeds `Hmax`.", call. = FALSE)
  }

  if (!is.null(init$Lambda) && ncol(init$Lambda) != init$H) {
    stop("`init$Lambda` must have `H` columns.", call. = FALSE)
  }

  if (!is.null(init$Delta) && ncol(init$Delta) != init$H) {
    stop("`init$Delta` must have `H` columns.", call. = FALSE)
  }

  if (!is.null(init$Eta) && nrow(init$Eta) != init$H) {
    stop("`init$Eta` must have `H` rows.", call. = FALSE)
  }

  if (!is.null(init$ell)) {
    if (length(init$ell) != init$H) {
      stop("`init$ell` must have length `H`.", call. = FALSE)
    }
    if (anyNA(init$ell)) {
      stop("`init$ell` cannot contain missing values.", call. = FALSE)
    }
    init$ell <- as.integer(init$ell)
  }

  if (!is.null(init$tau)) {
    if (length(init$tau) != init$H) {
      stop("`init$tau` must have length `H`.", call. = FALSE)
    }
    if (anyNA(init$tau)) {
      stop("`init$tau` cannot contain missing values.", call. = FALSE)
    }
  }

  if (!is.null(init$nu)) {
    if (!is.numeric(init$nu) || length(init$nu) != 1L || is.na(init$nu)) {
      stop("`init$nu` must be NULL or a single numeric value.", call. = FALSE)
    }
    if (init$nu <= 0 || init$nu >= 1) {
      stop("`init$nu` must lie strictly between 0 and 1.", call. = FALSE)
    }
  }

  if (!is.null(init$sigma2)) {
    if (!is.numeric(init$sigma2) || length(init$sigma2) != m || anyNA(init$sigma2)) {
      stop("`init$sigma2` must be NULL or a numeric vector of length `m`.", call. = FALSE)
    }
    if (any(init$sigma2 <= 0)) {
      stop("`init$sigma2` must contain positive values.", call. = FALSE)
    }
  }

  if (!is.null(init$Delta)) {
    vals <- unique(as.vector(init$Delta))
    vals <- vals[!is.na(vals)]
    if (!all(vals %in% c(0, 1))) {
      stop("`init$Delta` must contain only 0/1 entries.", call. = FALSE)
    }
    storage.mode(init$Delta) <- "integer"
  }

if (!is.null(init$Delta)) {
    if(is.null(init$tau)){
      init$tau <- colSums(init$Delta)/(nrow(init$Delta) - get_pivots(init$Delta)+1)
    }
  if(is.null(init$ell)){
      init$ell <- get_pivots(init$Delta)
    }
}


  # -- initialise Eta from its full conditional if not supplied --
  if (is.null(init$Eta) && !is.null(init$Lambda) &&
      !is.null(init$sigma2) && !is.null(y)) {

    H0      <- ncol(init$Lambda)
    Lam     <- init$Lambda                        # m x H
    sig2    <- init$sigma2                        # m
    Y       <- t(y)                               # m x T  (y is T x m)

    SiLam   <- Lam / sig2                         # m x H  (Sigma^{-1} Lambda, row-wise)
    V_inv   <- diag(H0) + crossprod(Lam, SiLam)  # H x H
    V       <- solve(V_inv)                       # H x H
    mu_eta  <- V %*% crossprod(SiLam, Y)          # H x T  (posterior mean, one col per obs)

    L       <- t(chol(V))                         # lower Cholesky H x H
    init$Eta <- mu_eta + L %*% matrix(rnorm(H0 * T), H0, T)
  }

  init
}

#' Validate and complete MCMC settings
#'
#' Processes the user-supplied `mcmc` list, fills in defaults, validates types
#' and ranges, and computes `nsave`.
#'
#' @param mcmc A named list with any subset of:
#' \describe{
#'   \item{`niter`}{Positive integer. Total number of MCMC iterations.
#'     Default `5000`.}
#'   \item{`nburn`}{Non-negative integer. Number of burn-in iterations to
#'     discard. Must be strictly less than `niter`. Default `2500`.}
#'   \item{`thin`}{Positive integer. Thinning interval. Default `1`.}
#' }
#'
#' @return A validated named list with entries `niter`, `nburn`, and `thin`
#'   (all integers). The derived quantity `nsave` is computed internally but
#'   not returned (it is recomputed inside [gltfa_cpp()]).
#'
#' @keywords internal
#' @noRd
gltfa_process_mcmc <- function(mcmc) {
  # ...existing code...
  if (!is.list(mcmc)) {
    stop("`mcmc` must be a named list.", call. = FALSE)
  }

  defaults <- list(
    niter = 5000L,
    nburn = 2500L,
    thin  = 1L
  )

  mcmc <- utils::modifyList(defaults, mcmc)

  req_names <- names(defaults)
  bad_names <- setdiff(names(mcmc), req_names)
  if (length(bad_names) > 0L) {
    stop(
      paste0(
        "Unknown entries in `mcmc`: ",
        paste(bad_names, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  mcmc$niter <- as.integer(mcmc$niter)
  mcmc$nburn <- as.integer(mcmc$nburn)
  mcmc$thin  <- as.integer(mcmc$thin)

  if (length(mcmc$niter) != 1L || is.na(mcmc$niter) || mcmc$niter <= 0L) {
    stop("`mcmc$niter` must be a positive integer.", call. = FALSE)
  }
  if (length(mcmc$nburn) != 1L || is.na(mcmc$nburn) || mcmc$nburn < 0L) {
    stop("`mcmc$nburn` must be a nonnegative integer.", call. = FALSE)
  }
  if (length(mcmc$thin) != 1L || is.na(mcmc$thin) || mcmc$thin <= 0L) {
    stop("`mcmc$thin` must be a positive integer.", call. = FALSE)
  }
  if (mcmc$nburn >= mcmc$niter) {
    stop("`mcmc$nburn` must be strictly smaller than `mcmc$niter`.", call. = FALSE)
  }

  mcmc$nsave <- as.integer((mcmc$niter - mcmc$nburn) %/% mcmc$thin)

  if (mcmc$nsave <= 0L) {
    stop(
      "The combination of `niter`, `nburn`, and `thin` yields no saved draws.",
      call. = FALSE
    )
  }

  mcmc$nsave <- NULL

  mcmc
}


#' Validate and complete prior hyperparameters
#'
#' Processes the user-supplied `prior` list, fills in defaults, and validates
#' all hyperparameters. Aliases `a_sigma`/`b_sigma` as `a_s`/`b_s` for
#' internal C++ functions.
#'
#' @param prior A named list with any subset of:
#' \describe{
#'   \item{`a_nu`}{Positive numeric. First shape of the Beta prior on `nu`.
#'     Default `1`.}
#'   \item{`b_nu`}{Positive numeric. Second shape of the Beta prior on `nu`.
#'     Default `1`.}
#'   \item{`alpha`}{Positive numeric. Shape parameter of the Gamma prior on
#'     factor loadings. Default `1`.}
#'   \item{`beta`}{Positive numeric. Rate parameter of the Gamma prior on
#'     factor loadings. Default `1`.}
#'   \item{`kappa`}{Positive numeric. Precision of the normal prior on
#'     non-pivot loadings. Default `1`.}
#'   \item{`a_sigma`}{Positive numeric. Shape of the inverse-Gamma prior on
#'     idiosyncratic variances. Default `2`.}
#'   \item{`b_sigma`}{Positive numeric. Scale of the inverse-Gamma prior on
#'     idiosyncratic variances. Default `2`.}
#' }
#'
#' @return A validated named list with all hyperparameters plus aliases
#'   `a_s = a_sigma` and `b_s = b_sigma`.
#'
#' @keywords internal
#' @noRd
gltfa_process_prior <- function(prior) {
  # ...existing code...
  if (!is.list(prior)) {
    stop("`prior` must be a named list.", call. = FALSE)
  }

  defaults <- list(
    a_nu    = 1,
    b_nu    = 1,
    alpha   = 1,
    beta    = 1,
    kappa   = 1,
    a_sigma     = 2,
    b_sigma     = 2,
    mh_sd_alpha = 0.2,    # RW std dev on log(alpha) scale; tune for ~30-40% acceptance
    a_alpha = NULL,
    b_alpha = NULL
    # no a_alpha / b_alpha defaults since the default is to have alpha fixed
  )

  prior <- utils::modifyList(defaults, prior)

  # Aliases: marg_lik / update_pivots_cpp use short names a_s / b_s
  # while gltfa_cpp uses a_sigma / b_sigma. Keep both in sync.
  prior$a_s <- prior$a_sigma
  prior$b_s <- prior$b_sigma

  req_names <- c(setdiff(names(defaults), c("a_alpha", "b_alpha")), "a_s", "b_s")
  opt_names <- c("a_alpha", "b_alpha")   # optional: only present when alpha is updated

  bad_names <- setdiff(names(prior), c(req_names, opt_names))
  if (length(bad_names) > 0L) {
    stop(
      paste0(
        "Unknown entries in `prior`: ",
        paste(bad_names, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  for (nm in req_names) {
    if (!is.numeric(prior[[nm]]) || length(prior[[nm]]) != 1L || is.na(prior[[nm]])) {
      stop(
        paste0("`prior$", nm, "` must be a single numeric value."),
        call. = FALSE
      )
    }
  }

  # Validate optional hyperprior parameters only when supplied
  if (!is.null(prior$a_alpha)) {
    if (!is.numeric(prior$a_alpha) || length(prior$a_alpha) != 1L || is.na(prior$a_alpha) || prior$a_alpha <= 0)
      stop("`prior$a_alpha` must be a single positive numeric value.", call. = FALSE)
  }
  if (!is.null(prior$b_alpha)) {
    if (!is.numeric(prior$b_alpha) || length(prior$b_alpha) != 1L || is.na(prior$b_alpha) || prior$b_alpha <= 0)
      stop("`prior$b_alpha` must be a single positive numeric value.", call. = FALSE)
  }

  if (prior$a_nu <= 0)        stop("`prior$a_nu` must be positive.",        call. = FALSE)
  if (prior$b_nu <= 0)        stop("`prior$b_nu` must be positive.",        call. = FALSE)
  if (prior$alpha <= 0)       stop("`prior$alpha` must be positive.",       call. = FALSE)
  if (prior$beta <= 0)        stop("`prior$beta` must be positive.",        call. = FALSE)
  if (prior$kappa <= 0)       stop("`prior$kappa` must be positive.",       call. = FALSE)
  if (prior$a_sigma <= 0)     stop("`prior$a_sigma` must be positive.",     call. = FALSE)
  if (prior$b_sigma <= 0)     stop("`prior$b_sigma` must be positive.",     call. = FALSE)
  if (prior$mh_sd_alpha <= 0) stop("`prior$mh_sd_alpha` must be positive.", call. = FALSE)

  prior
}


#' Validate and complete algorithmic control options
#'
#' Processes the user-supplied `control` list, fills in defaults, validates
#' types, and appends the `verbose` flag.
#'
#' @param control A named list with any subset of:
#' \describe{
#'   \item{`store_draws`}{Logical. Whether to store posterior draws of
#'     `Lambda`, `Delta`, `sigma2`, and `H`. Default `TRUE`.}
#'   \item{`store_eta`}{Logical. Whether to additionally store draws of the
#'     factor score matrix `Eta`. Default `FALSE`.}
#'   \item{`print_every`}{Positive integer. Print a progress message every
#'     this many iterations (only when `verbose = TRUE`). Default `100`.}
#'   \item{`seed`}{Single integer or `NULL`. Random seed passed to the C++
#'     sampler for reproducibility. Default `NULL`.}
#' }
#' @param verbose Logical. Whether to print sampling progress. Inherited from
#'   the `verbose` argument of [gltfa()].
#'
#' @return A validated named list with all control options plus
#'   `control$verbose`.
#'
#' @keywords internal
#' @noRd
gltfa_process_control <- function(control, verbose) {
  # ...existing code...
  if (!is.list(control)) {
    stop("`control` must be a named list.", call. = FALSE)
  }

  defaults <- list(
    store_draws = TRUE,
    store_eta   = FALSE,
    print_every = 100L,
    seed        = NULL,
    random_scan = TRUE
  )

  control <- utils::modifyList(defaults, control)

  req_names <- names(defaults)
  bad_names <- setdiff(names(control), req_names)
  if (length(bad_names) > 0L) {
    stop(
      paste0(
        "Unknown entries in `control`: ",
        paste(bad_names, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be a single TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(control$store_draws) ||
      length(control$store_draws) != 1L ||
      is.na(control$store_draws)) {
    stop("`control$store_draws` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(control$store_eta) ||
      length(control$store_eta) != 1L ||
      is.na(control$store_eta)) {
    stop("`control$store_eta` must be TRUE or FALSE.", call. = FALSE)
  }

  control$print_every <- as.integer(control$print_every)

  if (length(control$print_every) != 1L ||
      is.na(control$print_every) ||
      control$print_every <= 0L) {
    stop("`control$print_every` must be a positive integer.", call. = FALSE)
  }

  if (!is.null(control$seed)) {
    if (!is.numeric(control$seed) || length(control$seed) != 1L || is.na(control$seed)) {
      stop("`control$seed` must be NULL or a single numeric value.", call. = FALSE)
    }
    control$seed <- as.integer(control$seed)
  }

  control$verbose <- verbose

  control
}


#' Extract pivot row indices from a Delta matrix
#'
#' For each column of `Delta`, returns the row index of the first `1`, which
#' defines the pivot row for that factor.
#'
#' @param Delta An `m x H` integer matrix with entries in `{0, 1}`. Each column
#'   must contain at least one `1`.
#'
#' @return An integer vector of length `H`.
#'
#' @export
get_pivots <- function(Delta) {
  if (!is.matrix(Delta)) {
    stop("`Delta` must be a matrix.", call. = FALSE)
  }
  pivots <- apply(Delta, 2L, match, x = 1L)
  if (length(pivots) != length(unique(pivots))) {
    stop("Pivots need to be on different rows.", call. = FALSE)
  }
  if (anyNA(pivots)) {
    stop("Each column of `Delta` must contain at least one 1 to define a pivot.", call. = FALSE)
  }
  as.integer(pivots)
}