#' @keywords internal
#' @noRd
gltfa_pivots_from_delta <- function(Delta) {
  if (!is.matrix(Delta)) {
    stop("`Delta` must be a matrix.", call. = FALSE)
  }

  pivots <- apply(Delta, 2L, match, x = 1L)

  if (anyNA(pivots)) {
    stop("Each column of `Delta` must contain at least one 1 to define a pivot.", call. = FALSE)
  }

  as.integer(pivots)
}

#' @keywords internal
#' @noRd
update_delta_column_r <- function(y, Delta, Eta, j, tau_j, hyperpar, rows = NULL) {
  y <- as.matrix(y)
  Delta <- as.matrix(Delta)
  Eta <- as.matrix(Eta)

  m <- nrow(Delta)
  H <- ncol(Delta)

  if (!is.numeric(j) || length(j) != 1L || is.na(j)) {
    stop("`j` must be a single column index.", call. = FALSE)
  }
  j <- as.integer(j)

  if (j < 1L || j > H) {
    stop("`j` is out of bounds for `Delta`.", call. = FALSE)
  }

  if (!is.numeric(tau_j) || length(tau_j) != 1L || is.na(tau_j)) {
    stop("`tau_j` must be a single numeric value.", call. = FALSE)
  }
  if (tau_j <= 0 || tau_j >= 1) {
    stop("`tau_j` must lie strictly between 0 and 1.", call. = FALSE)
  }

  if (nrow(y) != ncol(Eta)) {
    stop("`y` and `Eta` have incompatible dimensions.", call. = FALSE)
  }
  if (ncol(y) != m) {
    stop("`y` and `Delta` have incompatible dimensions.", call. = FALSE)
  }
  if (nrow(Eta) != H) {
    stop("`Eta` and `Delta` must refer to the same number of factors.", call. = FALSE)
  }

  ell <- gltfa_pivots_from_delta(Delta)

  if (is.null(rows)) {
    rows <- seq.int(ell[j] + 1L, m)
  } else {
    if (anyNA(rows) || !is.numeric(rows)) {
      stop("`rows` must be NULL or a numeric vector of row indices.", call. = FALSE)
    }
    rows <- unique(as.integer(rows))
  }

  rows <- rows[rows >= 1L & rows <= m & rows > ell[j]]

  if (length(rows) == 0L) {
    return(list(
      Delta = Delta,
      rows = integer(0),
      log_post_odds = numeric(0),
      accepted = logical(0)
    ))
  }

  log_prior_odds <- log(tau_j) - log1p(-tau_j)
  log_post_odds <- numeric(length(rows))
  accepted <- logical(length(rows))

  for (k in seq_along(rows)) {
    i <- rows[k]

    Delta_one <- Delta
    Delta_zero <- Delta
    Delta_one[i, j] <- 1L
    Delta_zero[i, j] <- 0L

    ll_one <- marg_lik(i, y, Delta_one, Eta, hyperpar)
    ll_zero <- marg_lik(i, y, Delta_zero, Eta, hyperpar)
    log_post_odds[k] <- (ll_one - ll_zero) + log_prior_odds

    log_u <- log(stats::runif(1L))

    if (Delta[i, j] == 0L) {
      accepted[k] <- log_u <= log_post_odds[k]
      if (accepted[k]) {
        Delta[i, j] <- 1L
      }
    } else {
      accepted[k] <- log_u <= -log_post_odds[k]
      if (accepted[k]) {
        Delta[i, j] <- 0L
      }
    }
  }

  list(
    Delta = Delta,
    rows = rows,
    log_post_odds = log_post_odds,
    accepted = accepted
  )
}

#' @keywords internal
#' @noRd
update_delta_r <- function(y, Delta, Eta, tau, hyperpar, column_order = NULL) {
  Delta <- as.matrix(Delta)
  H <- ncol(Delta)

  if (H == 0L) {
    return(list(
      Delta = Delta,
      column_order = integer(0),
      rows = list(),
      log_post_odds = list(),
      accepted = list()
    ))
  }

  if (!is.numeric(tau) || length(tau) != H || anyNA(tau)) {
    stop("`tau` must be a numeric vector of length `ncol(Delta)`.", call. = FALSE)
  }
  if (any(tau <= 0 | tau >= 1)) {
    stop("All entries of `tau` must lie strictly between 0 and 1.", call. = FALSE)
  }

  if (is.null(column_order)) {
    column_order <- sample.int(H)
  } else {
    if (anyNA(column_order) || !is.numeric(column_order)) {
      stop("`column_order` must be NULL or a numeric vector of column indices.", call. = FALSE)
    }
    column_order <- as.integer(column_order)
    if (any(column_order < 1L | column_order > H)) {
      stop("`column_order` contains indices outside `1:ncol(Delta)`.", call. = FALSE)
    }
  }

  accepted <- vector("list", length(column_order))
  log_post_odds <- vector("list", length(column_order))
  rows <- vector("list", length(column_order))

  for (k in seq_along(column_order)) {
    j <- column_order[k]
    step <- update_delta_column_r(
      y = y,
      Delta = Delta,
      Eta = Eta,
      j = j,
      tau_j = tau[j],
      hyperpar = hyperpar
    )

    Delta <- step$Delta
    accepted[[k]] <- step$accepted
    log_post_odds[[k]] <- step$log_post_odds
    rows[[k]] <- step$rows
  }

  names(accepted) <- paste0("j", column_order)
  names(log_post_odds) <- names(accepted)
  names(rows) <- names(accepted)

  list(
    Delta = Delta,
    column_order = column_order,
    rows = rows,
    log_post_odds = log_post_odds,
    accepted = accepted
  )
}
