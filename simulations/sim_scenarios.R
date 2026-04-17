############################################################
# sim_gltfa.R
#
# Simulation script for Gaussian factor models with
# generalized lower-triangular flavored loading structures.
#
# It:
#   1. generates Gaussian data under multiple scenarios
#   2. runs gltfa(y, mcmc = list(), prior = list())
#   3. stores results in a reproducible format
#
# NOTE:
# - I leave mcmc = list() and prior = list() blank, as requested.
# - Adjust the orientation of y if your gltfa() expects variables in rows
#   rather than columns. Here y is n x p: rows = observations, cols = variables.
############################################################

rm(list = ls())

## ----------------------------- ##
## 0. User-facing configuration ##
## ----------------------------- ##

set.seed(123)

# Output directory
out_dir <- "sim_gltfa_output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Number of replicates per configuration
n_reps <- 20

# Sample sizes
sample_sizes <- c(50, 100, 250)

# Moderate dimensionalities:
# p = observed dimension, k = latent dimension
dim_grid <- list(
  list(p = 20, k = 5),
  list(p = 30, k = 7),
  list(p = 40, k = 10)
)

# Scenario names
scenario_names <- c(
  "scenario1_few_dense_columns",
  "scenario2_sparse_below",
  "scenario3_block_covariance",
  "scenario4_many_sparse_columns"
)

# Idiosyncratic variances: choose either fixed or random
use_random_uniqueness <- TRUE


## -------------------------------- ##
## 1. Helper functions: basic tools ##
## -------------------------------- ##

# Build a p x k indicator matrix Delta with ordered pivots.
# pivots must be strictly increasing integers in 1:p.
make_delta_from_pivots <- function(p, k, pivots, below_probs) {
  stopifnot(length(pivots) == k)
  stopifnot(length(below_probs) == k)
  stopifnot(all(diff(pivots) > 0))
  stopifnot(all(pivots >= 1), all(pivots <= p))

  Delta <- matrix(0, nrow = p, ncol = k)

  for (j in seq_len(k)) {
    # structural zeros above pivot
    # pivot entry forced to 1
    Delta[pivots[j], j] <- 1

    # below pivot, Bernoulli draws with column-specific sparsity
    if (pivots[j] < p) {
      idx <- (pivots[j] + 1):p
      Delta[idx, j] <- rbinom(length(idx), size = 1, prob = below_probs[j])
    }
  }

  Delta
}

# Given Delta, generate Lambda with Gaussian nonzero loadings.
# Signs are random; pivot entries are made reasonably strong.
make_lambda_from_delta <- function(Delta,
                                   pivot_scale = 1.2,
                                   offpivot_scale = 0.7,
                                   min_pivot_abs = 0.6) {
  p <- nrow(Delta)
  k <- ncol(Delta)

  Lambda <- matrix(0, nrow = p, ncol = k)

  for (j in seq_len(k)) {
    nz <- which(Delta[, j] == 1)
    pivot_j <- min(nz)

    # pivot loading
    val <- rnorm(1, mean = 0, sd = pivot_scale)
    while (abs(val) < min_pivot_abs) {
      val <- rnorm(1, mean = 0, sd = pivot_scale)
    }
    Lambda[pivot_j, j] <- val

    # other nonzero loadings below pivot
    other_nz <- setdiff(nz, pivot_j)
    if (length(other_nz) > 0) {
      Lambda[other_nz, j] <- rnorm(length(other_nz), mean = 0, sd = offpivot_scale)
    }
  }

  Lambda
}

# Generate diagonal uniqueness matrix Sigma
make_sigma_diag <- function(p, random = TRUE) {
  if (random) {
    vars <- runif(p, min = 0.3, max = 0.8)
  } else {
    vars <- rep(0.5, p)
  }
  diag(vars, nrow = p, ncol = p)
}

# Simulate Gaussian factor model data:
# y_t = Lambda eta_t + eps_t
# Returns y as n x p
simulate_factor_data <- function(n, Lambda, Sigma) {
  p <- nrow(Lambda)
  k <- ncol(Lambda)

  Eta <- matrix(rnorm(n * k), nrow = n, ncol = k)   # n x k
  Eps <- matrix(rnorm(n * p), nrow = n, ncol = p)   # n x p

  sigma_sd <- sqrt(diag(Sigma))
  Eps <- sweep(Eps, 2, sigma_sd, FUN = "*")

  Y <- Eta %*% t(Lambda) + Eps
  colnames(Y) <- paste0("y", seq_len(p))

  Y
}

# Convenience wrapper to package truth
make_truth_object <- function(Lambda, Sigma, scenario, n, rep_id, p, k, pivots, below_probs) {
  list(
    scenario = scenario,
    n = n,
    rep = rep_id,
    p = p,
    k = k,
    pivots = pivots,
    below_probs = below_probs,
    Lambda = Lambda,
    Sigma = Sigma,
    Omega = Lambda %*% t(Lambda) + Sigma
  )
}


## ---------------------------------------- ##
## 2. Scenario-specific loading generators  ##
## ---------------------------------------- ##

# Scenario 1:
# lower triangular factor loading not sparse below
# different sparsity between columns, with a few dense columns
generate_scenario1 <- function(p, k) {
  # reasonably spread pivots
  pivots <- round(seq(1, p - (k - 1), length.out = k)) + 0:(k - 1)
  pivots <- unique(pmax(1, pmin(p, pivots)))
  if (length(pivots) < k) pivots <- seq_len(k)

  # few dense columns, others moderately sparse
  below_probs <- rep(0.25, k)
  dense_cols <- seq_len(max(1, floor(k / 3)))
  below_probs[dense_cols] <- 0.80
  if (k >= 4) below_probs[k] <- 0.10

  Delta <- make_delta_from_pivots(p, k, pivots, below_probs)
  Lambda <- make_lambda_from_delta(Delta, pivot_scale = 1.3, offpivot_scale = 0.8)
  Sigma <- make_sigma_diag(p, random = use_random_uniqueness)

  list(
    Delta = Delta,
    Lambda = Lambda,
    Sigma = Sigma,
    pivots = pivots,
    below_probs = below_probs
  )
}

# Scenario 2:
# lower triangular factor loading sparse below
generate_scenario2 <- function(p, k) {
  pivots <- round(seq(1, p - (k - 1), length.out = k)) + 0:(k - 1)
  pivots <- unique(pmax(1, pmin(p, pivots)))
  if (length(pivots) < k) pivots <- seq_len(k)

  # all columns sparse below pivot
  below_probs <- seq(0.05, 0.20, length.out = k)

  Delta <- make_delta_from_pivots(p, k, pivots, below_probs)
  Lambda <- make_lambda_from_delta(Delta, pivot_scale = 1.2, offpivot_scale = 0.6)
  Sigma <- make_sigma_diag(p, random = use_random_uniqueness)

  list(
    Delta = Delta,
    Lambda = Lambda,
    Sigma = Sigma,
    pivots = pivots,
    below_probs = below_probs
  )
}

# Scenario 3:
# factors leading to block covariance structure
generate_scenario3 <- function(p, k) {
  # Use contiguous variable blocks, each factor mostly active in one block.
  # Keep GLT-like ordered pivots by placing pivot at first index of each block.
  block_sizes <- rep(floor(p / k), k)
  remainder <- p - sum(block_sizes)
  if (remainder > 0) block_sizes[seq_len(remainder)] <- block_sizes[seq_len(remainder)] + 1

  block_starts <- cumsum(c(1, block_sizes))[1:k]
  block_ends <- cumsum(block_sizes)

  pivots <- block_starts
  Delta <- matrix(0, nrow = p, ncol = k)
  below_probs <- rep(NA_real_, k)

  for (j in seq_len(k)) {
    Delta[pivots[j], j] <- 1

    # Strong within-block activation from pivot downward inside block
    idx_block <- seq.int(block_starts[j], block_ends[j])
    idx_block <- idx_block[idx_block > pivots[j]]
    if (length(idx_block) > 0) {
      Delta[idx_block, j] <- rbinom(length(idx_block), 1, 0.85)
    }

    # Very occasional spillover below block to allow some imperfect structure
    idx_tail <- seq.int(block_ends[j] + 1, p)
    if (length(idx_tail) > 0) {
      Delta[idx_tail, j] <- rbinom(length(idx_tail), 1, 0.05)
    }

    below_probs[j] <- 0.85
  }

  Lambda <- matrix(0, nrow = p, ncol = k)

  for (j in seq_len(k)) {
    nz <- which(Delta[, j] == 1)
    pivot_j <- pivots[j]

    # pivot and within-block loadings slightly stronger to induce blocks
    pivot_val <- rnorm(1, 1.3, 0.25)
    pivot_val <- pivot_val * sample(c(-1, 1), 1)
    Lambda[pivot_j, j] <- pivot_val

    other_nz <- setdiff(nz, pivot_j)
    if (length(other_nz) > 0) {
      # stronger within-block effects
      Lambda[other_nz, j] <- rnorm(length(other_nz), 0, 0.9)
    }
  }

  Sigma <- make_sigma_diag(p, random = use_random_uniqueness)

  list(
    Delta = Delta,
    Lambda = Lambda,
    Sigma = Sigma,
    pivots = pivots,
    below_probs = below_probs
  )
}

# Scenario 4:
# many columns and many zeroes below columns
# different degrees of sparsity, many sparse columns
generate_scenario4 <- function(p, k) {
  pivots <- seq_len(k)

  # many columns, most very sparse, a few moderately sparse
  base_probs <- seq(0.02, 0.15, length.out = k)
  below_probs <- sample(base_probs, size = k, replace = FALSE)

  Delta <- make_delta_from_pivots(p, k, pivots, below_probs)
  Lambda <- make_lambda_from_delta(Delta, pivot_scale = 1.1, offpivot_scale = 0.5)
  Sigma <- make_sigma_diag(p, random = use_random_uniqueness)

  list(
    Delta = Delta,
    Lambda = Lambda,
    Sigma = Sigma,
    pivots = pivots,
    below_probs = below_probs
  )
}


## -------------------------------- ##
## 3. Scenario dispatcher function  ##
## -------------------------------- ##

generate_scenario <- function(scenario, p, k) {
  switch(
    scenario,
    scenario1_few_dense_columns = generate_scenario1(p, k),
    scenario2_sparse_below      = generate_scenario2(p, k),
    scenario3_block_covariance  = generate_scenario3(p, k),
    scenario4_many_sparse_columns = generate_scenario4(p, k),
    stop("Unknown scenario: ", scenario)
  )
}


## ------------------------------- ##
## 4. Run one simulation instance  ##
## ------------------------------- ##

run_one_simulation <- function(scenario, n, p, k, rep_id) {
  scen <- generate_scenario(scenario, p, k)

  Y <- simulate_factor_data(
    n = n,
    Lambda = scen$Lambda,
    Sigma = scen$Sigma
  )

  # Run your function exactly as requested
  fit <- gltfa(
    y = Y,
    mcmc = list(),
    prior = list()
  )

  truth <- make_truth_object(
    Lambda = scen$Lambda,
    Sigma = scen$Sigma,
    scenario = scenario,
    n = n,
    rep_id = rep_id,
    p = p,
    k = k,
    pivots = scen$pivots,
    below_probs = scen$below_probs
  )

  list(
    meta = data.frame(
      scenario = scenario,
      n = n,
      p = p,
      k = k,
      rep = rep_id,
      stringsAsFactors = FALSE
    ),
    y = Y,
    truth = truth,
    fit = fit
  )
}


## --------------------------------------------- ##
## 5. Full simulation loop over all configurations ##
## --------------------------------------------- ##

all_results <- vector("list", length = 0)
counter <- 1

for (scenario in scenario_names) {
  for (dd in dim_grid) {
    p <- dd$p
    k <- dd$k

    for (n in sample_sizes) {
      for (rep_id in seq_len(n_reps)) {
        cat(
          sprintf(
            "Running: scenario=%s | p=%d | k=%d | n=%d | rep=%d\n",
            scenario, p, k, n, rep_id
          )
        )

        res <- run_one_simulation(
          scenario = scenario,
          n = n,
          p = p,
          k = k,
          rep_id = rep_id
        )

        all_results[[counter]] <- res
        counter <- counter + 1

        # save each run immediately for safety
        file_stub <- sprintf(
          "%s_p%d_k%d_n%d_rep%03d",
          scenario, p, k, n, rep_id
        )
        saveRDS(
          res,
          file = file.path(out_dir, paste0(file_stub, ".rds"))
        )
      }
    }
  }
}

# Save all results together as well
saveRDS(all_results, file = file.path(out_dir, "all_results.rds"))


## -------------------------------- ##
## 6. Build a simple summary table  ##
## -------------------------------- ##

meta_table <- do.call(
  rbind,
  lapply(all_results, function(x) x$meta)
)

write.csv(
  meta_table,
  file = file.path(out_dir, "simulation_index.csv"),
  row.names = FALSE
)

cat("\nSimulation study completed.\n")
cat("Results saved in:", normalizePath(out_dir), "\n")