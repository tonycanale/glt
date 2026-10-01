#' Fit a STAR-extended GLT factor model for discrete/mixed-scale data
#'
#' Fits the Kowal & Canale (2020) simultaneous transformation and rounding
#' (STAR) extension of [gltfa()] to discrete-valued (e.g. counts, ordinal)
#' observed data. A continuous latent matrix `z` is linked to the observed
#' matrix `y` via a column-specific monotone transform `G_j` and a
#' column-specific threshold sequence, and modelled with a regression mean
#' plus the same sparse GLT factor model used by [gltfa()]:
#'
#' \deqn{y_{ij} = \ell \iff G_j(z_{ij}) \in [a_{j,\ell}, a_{j,\ell+1})}
#' \deqn{z_{ij} = x_i'\beta_j + \epsilon_{ij}, \quad
#'   \epsilon_t \sim N_m(0, \Lambda\Lambda' \text{-structured } \Sigma)}
#'
#' Two extra Gibbs steps run before the original 8-step sampler on each
#' iteration: imputing `z` from its truncated-normal full conditional (Step
#' 0a), and updating the regression coefficients `Beta` (Step 0b, conjugate
#' normal). The residual `epsilon = z - X Beta` is then fed into the
#' unchanged Delta/pivots/H/Lambda/sigma2/Eta updates.
#'
#' @param y An integer (or coercible) \eqn{T \times m} matrix of observed,
#'   0-indexed levels (\eqn{0, 1, \ldots, K_j - 1} for column j). Missing
#'   values are not supported.
#' @param X A numeric \eqn{T \times p} covariate matrix, or `NULL` (default)
#'   for a zero-mean latent process (\eqn{p = 0}).
#' @param thresholds A list of length `ncol(y)`. Element `j` is either `NULL`
#'   (use the default integer-spaced cutpoints \eqn{1, 2, \ldots, \max(y_{\cdot
#'   j}) + 1}, reproducing floor rounding) or a strictly increasing numeric
#'   vector of \eqn{\max(y_{\cdot j}) + 1} interior cutpoints
#'   \eqn{(a_{j,1}, \ldots, a_{j,K_j - 1})}. Columns may have different
#'   numbers of levels \eqn{K_j}. Default `NULL` uses the default for every
#'   column.
#' @param g_type Character vector of length 1 or `ncol(y)`, each entry one of
#'   `"identity"` or `"exp"`, giving the column-specific monotone transform
#'   \eqn{G_j}. Default `"exp"` for every column (count-data STAR).
#' @param mcmc,Hmax,control,fixed,verbose As in [gltfa()].
#' @param prior A named list of prior hyperparameters; see [gltfa()] for the
#'   factor-model hyperparameters (`a_nu`, `b_nu`, `alpha`, `beta`, `kappa`,
#'   `a_sigma`, `b_sigma`, ...), plus regression hyperparameters:
#'   \describe{
#'     \item{`b0_beta`}{Numeric vector of length `ncol(X)`. Prior mean for
#'       every column's regression coefficients. Default `0`.}
#'     \item{`V_beta_inv`}{`p x p` prior precision matrix (common across
#'       columns), i.e. `beta_j ~ N(b0_beta, sigma2_j * V_beta_inv^{-1})`.
#'       Default `NULL`, in which case `V_beta_inv = I_p / c0_beta`.}
#'     \item{`c0_beta`}{Positive numeric. Prior variance scale used when
#'       `V_beta_inv` is not supplied. Default `100`.}
#'   }
#' @param init A named list of initial values; see [gltfa()] for the
#'   factor-model fields, plus:
#'   \describe{
#'     \item{`Beta`}{Numeric \eqn{p \times m} matrix. Default all zero.}
#'     \item{`z`}{Numeric \eqn{T \times m} matrix. Default a bound-consistent
#'       midpoint initialisation.}
#'   }
#'
#' @return An object of class `c("gltfit_star", "gltfit")`. In addition to
#'   the components documented in [gltfa()], `draws$Beta` holds `p x m`
#'   draws of the regression coefficients (always stored), and `last$z` /
#'   `last$epsilon` hold the final latent-z / residual draws (not stored
#'   across iterations).
#'
#' @references
#' Canale, A. and Frühwirth-Schnatter, S. (2026).
#' "The generalized lower triangular process prior."
#' \emph{Technical Report}, 
#' 
#' Kowal, D.  and Canale, A. (2020).
#' "Simultaneous transformation and rounding (STAR) models for integer-valued data."
#' \emph{Electronic Journal of Statistics}, 14(1), 1744-1772.
#' 
#' @seealso [gltfa()]
#'
#' @export
gltfa_star <- function(
  y,
  X = NULL,
  thresholds = NULL,
  g_type = "exp",
  mcmc = list(),
  Hmax = ncol(y),
  prior = list(),
  init = list(),
  control = list(),
  fixed = list(),
  verbose = interactive()
) {
  y <- as.matrix(y)
  storage.mode(y) <- "integer"

  if (anyNA(y)) {
    stop("Missing values in `y` are not supported in the current version.", call. = FALSE)
  }
  if (any(y < 0L)) {
    stop("`y` must contain only nonnegative integer levels (0-indexed).", call. = FALSE)
  }

  T <- nrow(y)
  m <- ncol(y)

  if (is.null(T) || T < 2L) {
    stop("`y` must have at least 2 rows.", call. = FALSE)
  }
  if (is.null(m) || m < 1L) {
    stop("`y` must have at least 1 column.", call. = FALSE)
  }

  if (!is.numeric(Hmax) || length(Hmax) != 1L || is.na(Hmax) || Hmax <= 1) {
    stop("`Hmax` must be a single number greater than 1.", call. = FALSE)
  }
  Hmax <- as.integer(Hmax)

  if (is.null(X)) {
    X <- matrix(0, nrow = T, ncol = 0)
  } else {
    X <- as.matrix(X)
    if (!is.numeric(X)) {
      stop("`X` must be NULL or a numeric matrix.", call. = FALSE)
    }
    if (nrow(X) != T) {
      stop("`X` must have `nrow(X) == nrow(y)`.", call. = FALSE)
    }
  }
  p <- ncol(X)

  if (!is.list(mcmc))    stop("`mcmc` must be a named list.", call. = FALSE)
  if (!is.list(prior))   stop("`prior` must be a named list.", call. = FALSE)
  if (!is.list(init))    stop("`init` must be a named list.", call. = FALSE)
  if (!is.list(control)) stop("`control` must be a named list.", call. = FALSE)
  if (!is.list(fixed))   stop("`fixed` must be a named list.", call. = FALSE)

  # Same hyperprior wiring as gltfa(): alpha/beta/a_nu fixed unless both
  # shape hyperparameters are supplied.
  if (is.null(prior$a_alpha) && is.null(prior$b_alpha)) {
    fixed$alpha   <- TRUE
    prior$a_alpha     <- 1.0
    prior$b_alpha <- 1.0
    prior$mh_sd_alpha <- prior$mh_sd_alpha %||% 0.2
  } else if (is.null(prior$a_alpha) || is.null(prior$b_alpha)) {
    stop("`prior$a_alpha` and `prior$b_alpha` must be both provided or both NULL.", call. = FALSE)
  } else {
    if (prior$a_alpha <= 0) stop("`prior$a_alpha` must be positive.", call. = FALSE)
    if (prior$b_alpha <= 0) stop("`prior$b_alpha` must be positive.", call. = FALSE)
    fixed$alpha <- FALSE
  }

  if (is.null(prior$a_beta) && is.null(prior$b_beta)) {
    fixed$beta   <- TRUE
    prior$a_beta <- 1.0
    prior$b_beta <- 1.0
    prior$mh_sd_beta <- prior$mh_sd_beta %||% 0.2
  } else if (is.null(prior$a_beta) || is.null(prior$b_beta)) {
    stop("`prior$a_beta` and `prior$b_beta` must be both provided or both NULL.", call. = FALSE)
  } else {
    if (prior$a_beta <= 0) stop("`prior$a_beta` must be positive.", call. = FALSE)
    if (prior$b_beta <= 0) stop("`prior$b_beta` must be positive.", call. = FALSE)
    fixed$beta <- FALSE
  }

  if (is.null(prior$a_anu) && is.null(prior$b_anu)) {
    fixed$a_nu  <- TRUE
    prior$a_anu <- 1.0
    prior$b_anu <- 1.0
    prior$mh_sd_a_nu <- prior$mh_sd_a_nu %||% 0.2
  } else if (is.null(prior$a_anu) || is.null(prior$b_anu)) {
    stop("`prior$a_anu` and `prior$b_anu` must be both provided or both NULL.", call. = FALSE)
  } else {
    if (prior$a_anu <= 0) stop("`prior$a_anu` must be positive.", call. = FALSE)
    if (prior$b_anu <= 0) stop("`prior$b_anu` must be positive.", call. = FALSE)
    fixed$a_nu <- FALSE
  }

  if (isTRUE(fixed$H)) {
    fixed$nu   <- TRUE
    fixed$a_nu <- TRUE
  }
  if (isTRUE(fixed$pivots)) {
    fixed$H    <- TRUE
    fixed$nu   <- TRUE
    fixed$a_nu <- TRUE
  }
  if (isTRUE(fixed$Delta)) {
    fixed$pivots <- TRUE
    fixed$H      <- TRUE
    fixed$nu     <- TRUE
    fixed$a_nu   <- TRUE
    fixed$tau    <- TRUE
    fixed$alpha  <- TRUE
    fixed$beta   <- TRUE
  }
  if (isTRUE(fixed$tau)) {
    fixed$alpha <- TRUE
    fixed$beta  <- TRUE
  }
  if (p == 0L) {
    fixed$Beta <- TRUE
  }

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be a single TRUE or FALSE.", call. = FALSE)
  }

  link <- gltfa_process_star_link(y, thresholds, g_type)

  mcmc    <- gltfa_process_mcmc(mcmc)
  prior   <- gltfa_process_prior(prior)
  control <- gltfa_process_control(control, verbose = verbose)

  # Regression-prior defaults (b0_beta / V_beta_inv / c0_beta)
  if (is.null(prior$b0_beta)) {
    prior$b0_beta <- rep(0, p)
  } else if (length(prior$b0_beta) != p) {
    stop("`prior$b0_beta` must have length `ncol(X)`.", call. = FALSE)
  }
  if (!is.null(prior$V_beta_inv)) {
    Vb <- as.matrix(prior$V_beta_inv)
    if (nrow(Vb) != p || ncol(Vb) != p) {
      stop("`prior$V_beta_inv` must be a `p x p` matrix.", call. = FALSE)
    }
    prior$V_beta_inv <- Vb
  }
  if (is.null(prior$c0_beta)) {
    prior$c0_beta <- 100
  } else if (prior$c0_beta <= 0) {
    stop("`prior$c0_beta` must be positive.", call. = FALSE)
  }

  init <- gltfa_process_init_star(
    init = init,
    T = T,
    m = m,
    Hmax = Hmax,
    p = p,
    y = NULL # y is discrete; Eta-from-conditional init not applicable here
  )

  model <- list(
    T = T,
    m = m,
    Hmax = Hmax
  )

  fit <- gltfa_star_cpp(
    y_obs = y,
    X = X,
    thresholds = link$thresholds,
    g_type = link$g_type,
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
    fixed = fixed,
    thresholds = link$thresholds,
    g_type = g_type
  )

  class(fit) <- c("gltfit_star", "gltfit")
  fit
}
