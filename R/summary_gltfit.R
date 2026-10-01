#' Identify admissible posterior draws of Delta
#'
#' Flags each saved MCMC draw of the allocation matrix \eqn{\Delta} as
#' "admissible" for posterior summarisation. A draw is excluded when, after
#' dropping spurious columns (columns with exactly one nonzero entry, i.e.
#' factors loading on a single variable) and any resulting all-zero rows, the
#' reduced allocation matrix fails the counting ("3-5-7-9") identifiability
#' rule of [sparvaride::counting_rule_holds()].
#'
#' @param Delta_draws A list of length `nsave`, each element an `m x H_s`
#'   integer/binary matrix (as stored in `object$draws$Delta` for a
#'   `"gltfit"` object; `H_s` may vary across draws).
#'
#' @return A logical vector of length `length(Delta_draws)`; `TRUE` marks an
#'   admissible draw.
#'
#' @seealso [summary.gltfit()], [sparvaride::counting_rule_holds()]
#' @keywords internal
#' @noRd
admissible_draws <- function(Delta_draws) {
  nsave <- length(Delta_draws)
  # Count of spurious columns (colSum == 1) per draw, via the C++ routine
  # shared with metrics_cpp.cpp (see src/utils_aux.h::compute_spurious).
  spurious_count <- compute_spurious_cpp(Delta_draws)

  admissible <- logical(nsave)
  for (s in seq_len(nsave)) {
    Del <- Delta_draws[[s]]

    # Drop spurious (single-loading) columns, if any.
    D_sub <- if (spurious_count[s] > 0) {
      Del[, colSums(Del) > 1, drop = FALSE]
    } else {
      Del
    }

    # Drop rows with no remaining loadings.
    if (nrow(D_sub) > 0 && any(rowSums(D_sub) == 0)) {
      D_sub <- D_sub[rowSums(D_sub) > 0, , drop = FALSE]
    }

    admissible[s] <- if (nrow(D_sub) == 0 || ncol(D_sub) == 0) {
      FALSE
    } else {
      sparvaride::counting_rule_holds(D_sub)
    }
  }
  admissible
}


#' Pad an m x H matrix with zero columns up to m x Hmax
#'
#' @param mat An `m x H` numeric matrix, `H <= Hmax`.
#' @param Hmax Target number of columns.
#'
#' @return An `m x Hmax` matrix, zero-padded on the right.
#' @keywords internal
#' @noRd
pad_cols <- function(mat, Hmax) {
  H <- ncol(mat)
  if (H == Hmax) return(mat)
  cbind(mat, matrix(0, nrow = nrow(mat), ncol = Hmax - H))
}


#' Summarize a "gltfit" object
#'
#' Computes posterior averages (across admissible saved MCMC draws) of
#' \code{Lambda}, \code{Delta}, \code{sigma2}, and \code{H}, along with the
#' posterior average of \eqn{\Omega = \Lambda \Lambda^\top + \Sigma}, where
#' \eqn{\Sigma = \mathrm{diag}(\sigma^2)}. The posterior-mean \code{Lambda},
#' \code{Delta}, and \code{Omega} matrices are returned with class
#' \code{"facmatrix"} so they inherit [plot.facmatrix()].
#'
#' Before averaging, draws are filtered using `admissible_draws()`: a draw is
#' dropped when, after removing spurious columns (loadings on a single
#' variable) and any resulting all-zero rows, the allocation matrix fails
#' [sparvaride::counting_rule_holds()]. This mirrors the usual post-processing
#' step for GLT-type factor models, where unidentified draws (e.g. those with
#' too few active variables per factor) should not contribute to posterior
#' summaries of \code{Lambda}/\code{Delta}.
#'
#' \code{object$draws$Lambda} and \code{object$draws$Delta} are stored as
#' ragged lists of length \eqn{n_{\text{save}}}, each element an \eqn{m \times
#' H_s} matrix (\eqn{H_s} is the active number of factors for that draw, see
#' \code{object$draws$H}). Kept draws are zero-padded to \eqn{m \times
#' H_{\max}} before elementwise averaging. \code{sigma2} draws are stored as
#' an \eqn{n_{\text{save}} \times m} matrix, averaged over kept rows.
#' \eqn{\Omega} is computed per kept draw as \eqn{\Lambda_s \Lambda_s^\top +
#' \mathrm{diag}(\sigma^2_s)} (using the draw's own, unpadded \eqn{\Lambda_s})
#' and then averaged.
#'
#' @param object An object of class \code{"gltfit"}, as returned by [gltfa()].
#' @param ... Currently unused.
#'
#' @return A named list with class \code{"summary.gltfit"} containing:
#' \describe{
#'   \item{`Lambda`}{Posterior mean loading matrix (\eqn{m \times H_{\max}}),
#'     class \code{"facmatrix"}.}
#'   \item{`Delta`}{Posterior mean allocation matrix (\eqn{m \times H_{\max}}),
#'     class \code{"facmatrix"}.}
#'   \item{`sigma2`}{Posterior mean idiosyncratic variances (length \eqn{m}).}
#'   \item{`H`}{Posterior mean number of active factors (scalar).}
#'   \item{`Omega`}{Posterior mean of \eqn{\Lambda \Lambda^\top + \Sigma}
#'     (\eqn{m \times m}), class \code{"facmatrix"}.}
#'   \item{`n_admissible`}{Number of draws retained after filtering.}
#'   \item{`n_total`}{Total number of saved draws before filtering.}
#' }
#'
#' @seealso [gltfa()], [plot.facmatrix()], [sparvaride::counting_rule_holds()]
#' @export
summary.gltfit <- function(object, ...) {
  if (is.null(object$draws)) {
    stop("`object$draws` is missing; was `control$store_draws` set to FALSE?", call. = FALSE)
  }

  draws <- object$draws

  Lambda_draws <- draws$Lambda
  Delta_draws  <- draws$Delta
  sigma2_draws <- draws$sigma2
  H_draws      <- draws$H

  if (!is.list(Lambda_draws) || !is.list(Delta_draws)) {
    stop("`object$draws$Lambda`/`Delta` are expected to be lists of per-draw matrices.", call. = FALSE)
  }

  nsave <- length(Lambda_draws)
  m     <- nrow(Lambda_draws[[1]])

  keep <- admissible_draws(Delta_draws)
  n_admissible <- sum(keep)

  if (n_admissible == 0) {
    stop(
      "No posterior draws satisfy the counting-rule admissibility criterion; ",
      "cannot compute a posterior summary.",
      call. = FALSE
    )
  }

  idx <- which(keep)

  # Use the largest H actually explored among the admissible draws, not the
  # user-specified ceiling in `meta$Hmax` (which may never be attained).
  Hmax <- max(vapply(Lambda_draws[idx], ncol, integer(1)))

  Lambda_sum <- matrix(0, nrow = m, ncol = Hmax)
  Delta_sum  <- matrix(0, nrow = m, ncol = Hmax)
  Omega_sum  <- matrix(0, nrow = m, ncol = m)

  for (s in idx) {
    Lambda_s <- Lambda_draws[[s]]
    Delta_s  <- Delta_draws[[s]]
    sigma2_s <- sigma2_draws[s, ]

    Lambda_sum <- Lambda_sum + pad_cols(Lambda_s, Hmax)
    Delta_sum  <- Delta_sum  + pad_cols(Delta_s,  Hmax)
    Omega_sum  <- Omega_sum  + tcrossprod(Lambda_s) + diag(sigma2_s, nrow = m)
  }

  Lambda_mean <- Lambda_sum / n_admissible
  Delta_mean  <- Delta_sum  / n_admissible
  Omega_mean  <- Omega_sum  / n_admissible
  sigma2_mean <- colMeans(sigma2_draws[idx, , drop = FALSE])
  H_mean      <- mean(H_draws[idx])

  out <- list(
    Lambda = facmatrix(Lambda_mean),
    Delta  = facmatrix(Delta_mean),
    sigma2 = sigma2_mean,
    H      = H_mean,
    Omega  = facmatrix(Omega_mean),
    n_admissible = n_admissible,
    n_total = nsave
  )

  class(out) <- "summary.gltfit"
  out
}


#' Print a "summary.gltfit" object
#'
#' @param x A \code{"summary.gltfit"} object, as returned by
#'   [summary.gltfit()].
#' @param ... Currently unused.
#'
#' @return Invisibly `x`.
#'
#' @export
print.summary.gltfit <- function(x, ...) {
  cat("Posterior summary of a \"gltfit\" object\n")
  cat("  Draws kept: ", x$n_admissible, "/", x$n_total,
      "(after counting-rule admissibility filtering)\n")
  cat("  Lambda: ", nrow(x$Lambda), "x", ncol(x$Lambda), "posterior mean loading matrix\n")
  cat("  Delta:  ", nrow(x$Delta), "x", ncol(x$Delta), "posterior mean allocation matrix\n")
  cat("  Omega:  ", nrow(x$Omega), "x", ncol(x$Omega), "posterior mean of Lambda %*% t(Lambda) + Sigma\n")
  cat("  sigma2: posterior mean idiosyncratic variances (length", length(x$sigma2), ")\n")
  cat("  H:      posterior mean number of active factors =", round(x$H, 3), "\n")
  invisible(x)
}
