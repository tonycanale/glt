gltfa_process_init <- function(init, T, m, Hmax) {
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

  init
}

gltfa_process_mcmc <- function(mcmc) {
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

  mcmc
}


gltfa_process_prior <- function(prior) {
  if (!is.list(prior)) {
    stop("`prior` must be a named list.", call. = FALSE)
  }

  defaults <- list(
    a_nu    = 1,
    b_nu    = 1,
    alpha   = 1,
    beta    = 1,
    kappa   = 1,
    a_sigma = 2,
    b_sigma = 2
  )

  prior <- utils::modifyList(defaults, prior)

  # Aliases: marg_lik / update_pivots_cpp use short names a_s / b_s
  # while gltfa_cpp uses a_sigma / b_sigma. Keep both in sync.
  prior$a_s <- prior$a_sigma
  prior$b_s <- prior$b_sigma

  req_names <- c(names(defaults), "a_s", "b_s")
  bad_names <- setdiff(names(prior), req_names)
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

  if (prior$a_nu <= 0) {
    stop("`prior$a_nu` must be positive.", call. = FALSE)
  }
  if (prior$b_nu <= 0) {
    stop("`prior$b_nu` must be positive.", call. = FALSE)
  }
  if (prior$alpha <= 0) {
    stop("`prior$alpha` must be positive.", call. = FALSE)
  }
  if (prior$beta <= 0) {
    stop("`prior$beta` must be positive.", call. = FALSE)
  }
  if (prior$kappa <= 0) {
    stop("`prior$kappa` must be positive.", call. = FALSE)
  }
  if (prior$a_sigma <= 0) {
    stop("`prior$a_sigma` must be positive.", call. = FALSE)
  }
  if (prior$b_sigma <= 0) {
    stop("`prior$b_sigma` must be positive.", call. = FALSE)
  }

  prior
}


gltfa_process_control <- function(control, verbose) {
  if (!is.list(control)) {
    stop("`control` must be a named list.", call. = FALSE)
  }

  defaults <- list(
    store_draws = TRUE,
    store_eta = FALSE,
    print_every = 100L,
    seed = NULL
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