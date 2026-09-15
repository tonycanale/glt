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
    a_nu        = 1,
    b_nu        = 1,
    alpha       = 1,
    beta        = 1,
    kappa       = 1,
    a_sigma     = 2,
    b_sigma     = 2,
    mh_sd_alpha = 0.2,
    mh_sd_beta  = 0.2,
    mh_sd_a_nu  = 0.2,
    a_alpha     = NULL,
    b_alpha     = NULL,
    b_beta      = NULL,
    a_anu       = NULL,
    b_anu       = NULL
  )

  prior <- utils::modifyList(defaults, prior)

  # Aliases: marg_lik / update_pivots_cpp use short names a_s / b_s
  # while gltfa_cpp uses a_sigma / b_sigma. Keep both in sync.
  prior$a_s <- prior$a_sigma
  prior$b_s <- prior$b_sigma

  req_names <- c(setdiff(names(defaults), c("a_alpha", "b_alpha", "a_beta", "b_beta", "a_anu", "b_anu")), "a_s", "b_s")
  opt_names <- c("a_alpha", "b_alpha", "a_beta", "b_beta", "a_anu", "b_anu")

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
  if (!is.null(prior$a_beta)) {
    if (!is.numeric(prior$a_beta) || length(prior$a_beta) != 1L || is.na(prior$a_beta) || prior$a_beta <= 0)
      stop("`prior$a_beta` must be a single positive numeric value.", call. = FALSE)
  }
  if (!is.null(prior$b_beta)) {
    if (!is.numeric(prior$b_beta) || length(prior$b_beta) != 1L || is.na(prior$b_beta) || prior$b_beta <= 0)
      stop("`prior$b_beta` must be a single positive numeric value.", call. = FALSE)
  }
  if (!is.null(prior$a_anu)) {
    if (!is.numeric(prior$a_anu) || length(prior$a_anu) != 1L || is.na(prior$a_anu) || prior$a_anu <= 0)
      stop("`prior$a_anu` must be a single positive numeric value.", call. = FALSE)
  }
  if (!is.null(prior$b_anu)) {
    if (!is.numeric(prior$b_anu) || length(prior$b_anu) != 1L || is.na(prior$b_anu) || prior$b_anu <= 0)
      stop("`prior$b_anu` must be a single positive numeric value.", call. = FALSE)
  }

  if (prior$a_nu <= 0)       stop("`prior$a_nu` must be positive.",       call. = FALSE)
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


#' Build default STAR thresholds for a single column
#'
#' Default threshold sequence is integer-spaced by 1: interior cutpoints
#' `1, 2, ..., max(y_col) + 1` (0-indexed levels `0, ..., max(y_col)`,
#' i.e. the original floor-rounding behaviour).
#'
#' @param y_col Integer vector of 0-indexed observed levels for one column.
#'
#' @return Numeric vector of interior cutpoints.
#'
#' @keywords internal
#' @noRd
star_default_thresholds <- function(y_col) {
  kmax <- max(y_col)
  if (kmax < 0L) {
    stop("Observed levels must be nonnegative integers.", call. = FALSE)
  }
  seq_len(kmax + 1L) # 1, 2, ..., kmax + 1
}

#' Validate and complete `thresholds` / `g_type` for [gltfa_star()]
#'
#' @param y Integer T x m matrix of 0-indexed observed levels.
#' @param thresholds NULL, or a list of length `ncol(y)`; each element is
#'   either NULL (use the default integer-spaced cutpoints for that column)
#'   or a strictly increasing numeric vector of interior cutpoints.
#' @param g_type NULL, or a character vector recycled/validated against
#'   `ncol(y)`, each entry one of `"identity"`, `"exp"`.
#'
#' @return A list with elements `thresholds` (list of numeric vectors) and
#'   `g_type` (integer vector of transform ids, 0 = identity, 1 = exp).
#'
#' @keywords internal
#' @noRd
gltfa_process_star_link <- function(y, thresholds, g_type) {
  m <- ncol(y)

  if (is.null(g_type)) {
    g_type <- rep("exp", m)
  }
  if (!is.character(g_type)) {
    stop("`g_type` must be a character vector.", call. = FALSE)
  }
  if (length(g_type) == 1L) {
    g_type <- rep(g_type, m)
  }
  if (length(g_type) != m) {
    stop("`g_type` must have length 1 or `ncol(y)`.", call. = FALSE)
  }
  if (!all(g_type %in% c("identity", "exp"))) {
    stop("`g_type` entries must be 'identity' or 'exp'.", call. = FALSE)
  }
  g_type_id <- as.integer(match(g_type, c("identity", "exp")) - 1L)

  if (is.null(thresholds)) {
    thresholds <- vector("list", m)
  }
  if (!is.list(thresholds) || length(thresholds) != m) {
    stop("`thresholds` must be NULL or a list of length `ncol(y)`.", call. = FALSE)
  }

  for (j in seq_len(m)) {
    if (is.null(thresholds[[j]])) {
      thresholds[[j]] <- star_default_thresholds(y[, j])
    } else {
      a_j <- as.numeric(thresholds[[j]])
      if (anyNA(a_j) || length(a_j) < 1L) {
        stop(sprintf("`thresholds[[%d]]` must be a non-empty numeric vector.", j), call. = FALSE)
      }
      if (is.unsorted(a_j, strictly = TRUE)) {
        stop(sprintf("`thresholds[[%d]]` must be strictly increasing.", j), call. = FALSE)
      }
      kmax_valid <- length(a_j) # thresholds imply valid levels 0, ..., length(a_j)
      kmax_j <- max(y[, j])
      if (kmax_j > kmax_valid) {
        stop(
          sprintf(
            "`thresholds[[%d]]` has %d interior cutpoint(s), which only supports levels 0..%d, but column %d of `y` contains a level as high as %d. Supply more cutpoints (one per boundary between consecutive levels).",
            j, kmax_valid, kmax_valid, j, kmax_j
          ),
          call. = FALSE
        )
      }
      thresholds[[j]] <- a_j
    }
  }

  list(thresholds = thresholds, g_type = g_type_id)
}

#' Validate and complete initial values for the STAR MCMC sampler
#'
#' Extends [gltfa_process_init()] with STAR-specific initial state: the
#' regression coefficients `Beta` (p x m) and the latent matrix `z` (T x m).
#'
#' @param init A named list, see [gltfa_process_init()] plus `Beta`, `z`.
#' @param T,m,Hmax,y As in [gltfa_process_init()].
#' @param p Integer. Number of covariates (`ncol(X)`, may be 0).
#'
#' @return A validated and coerced named list ready to be passed to
#'   [gltfa_star_cpp()].
#'
#' @keywords internal
#' @noRd
gltfa_process_init_star <- function(init, T, m, Hmax, p, y = NULL) {
  Beta <- init$Beta
  z    <- init$z
  init$Beta <- NULL
  init$z    <- NULL

  init <- gltfa_process_init(init, T = T, m = m, Hmax = Hmax, y = y)

  if (!is.null(Beta)) {
    if (!is.matrix(Beta) || !is.numeric(Beta)) {
      stop("`init$Beta` must be NULL or a numeric matrix.", call. = FALSE)
    }
    if (nrow(Beta) != p || ncol(Beta) != m) {
      stop("`init$Beta` must have dimension p x m.", call. = FALSE)
    }
  } else {
    Beta <- matrix(0, nrow = p, ncol = m)
  }
  init$Beta <- Beta

  if (!is.null(z)) {
    if (!is.matrix(z) || !is.numeric(z)) {
      stop("`init$z` must be NULL or a numeric matrix.", call. = FALSE)
    }
    if (nrow(z) != T || ncol(z) != m) {
      stop("`init$z` must have dimension T x m.", call. = FALSE)
    }
    init$z <- z
  }

  init
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