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
