#' Posterior draws of variance decomposition
#'
#' For each saved MCMC draw, computes the \eqn{m \times m} matrices
#' \eqn{\Lambda\Lambda^\top} (factor covariance) and
#' \eqn{\Lambda\Lambda^\top + \Sigma} (total marginal covariance), where
#' \eqn{\Sigma = \mathrm{diag}(\sigma^2_1, \ldots, \sigma^2_m)}.
#'
#' @param fit A `gltfit` object returned by [gltfa()].
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{`LLt`}{An \eqn{m \times m \times S} array of
#'     \eqn{\Lambda_s \Lambda_s^\top} draws.}
#'   \item{`total`}{An \eqn{m \times m \times S} array of
#'     \eqn{\Lambda_s \Lambda_s^\top + \Sigma_s} draws.}
#' }
#'
#' @export
get_variance <- function(fit) {
  if (!inherits(fit, "gltfit")) {
    stop("`fit` must be a `gltfit` object.", call. = FALSE)
  }

  draws <- fit$draws

  if (is.null(draws$Lambda) || is.null(draws$sigma2)) {
    stop(
      "Variance draws require `store_draws = TRUE` in `control`.",
      call. = FALSE
    )
  }

  Lambda_draws <- draws$Lambda   # list of length S, each m x H_s
  sigma2_draws <- draws$sigma2   # S x m matrix

  S <- length(Lambda_draws)
  m <- nrow(Lambda_draws[[1L]])

  LLt   <- array(NA_real_, dim = c(m, m, S))
  total <- array(NA_real_, dim = c(m, m, S))

  for (s in seq_len(S)) {
    L  <- Lambda_draws[[s]]                      # m x H_s
    LL <- tcrossprod(L)                          # L %*% t(L), m x m
    LLt[, , s]   <- LL
    total[, , s] <- LL + diag(sigma2_draws[s, ], nrow = m)
  }

  list(LLt = LLt, total = total)
}
