############################################################
# sim_gltfa.R
# Simulation script for Gaussian factor models with
# generalized lower-triangular flavored loading structures.
#
############################################################

rm(list = ls())

## ----------------------------- ##
## 0. User-facing configuration ##
## ----------------------------- ##

set.seed(123)

# Output directory
out_dir <- "sim_gltfa_output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Number of replicates
n_reps <- 20

# Dimensions 
m_dim <- c(20, 50, 100)
T_dim <- c(m_dim, 2 * m_dim)
H_dim <- matrix(
  c(3,8,10,
    5,10,15,
    5,10,15,
    8, 20, 40),
  nrow = 4, ncol = 3, byrow = TRUE
)

# Scenario names
scenario_names <- c(
  "scenario1_few_dense_columns",
  "scenario2_sparse_below", # similar to Sylvia's scenario 
  "scenario3_block_covariance",
  "scenario4_many_sparse_columns". ## H < (m-1)/2
)

## -------------------------------- ##
## 1. Helper functions: basic tools ##
## -------------------------------- ##

# Build a m x H indicator matrix Delta with ordered pivots.
# pivots must be strictly increasing integers in 1:m.
make_delta_from_pivots <- function(m, H, pivots, below_probs) {
  stopifnot(length(pivots) == H)
  stopifnot(length(below_probs) == H)
  stopifnot(all(diff(pivots) > 0))
  stopifnot(all(pivots >= 1), all(pivots <= m))

  Delta <- matrix(0, nrow = m, ncol = H)

  for (j in seq_len(H)) {
    # structural zeros above pivot
    # pivot entry forced to 1
    Delta[pivots[j], j] <- 1

    # below pivot, Bernoulli draws with column-specific sparsity
    if (pivots[j] < m) {
      idx <- (pivots[j] + 1):m
      Delta[idx, j] <- rbinom(length(idx), size = 1, prob = below_probs[j])
    }
  }

  Delta
}

# Given Delta, generate Lambda with Gaussian nonzero loadings.
# Signs are random; pivot entries are made reasonably strong.
make_lambda_from_delta <- function(Delta, negative_prob = 0.1,
  mean =1, sd = 0.1){
  m <- nrow(Delta)
  H <- ncol(Delta)
  Lambda <- matrix( rnorm(m * H, mean = mean, sd = sd),
    nrow = m, ncol = H)
  Lambda_neg <- matrix(sample(c(-1,1),m * H, replace = TRUE, 
    prob = c(negative_prob, 1 - negative_prob)),
    nrow = m, ncol = H)
  Lambda <- Lambda * Lambda_neg
  Lambda <- Lambda * Delta
  Lambda
}

# Simulate Gaussian factor model data:
# y_t = Lambda eta_t + eps_t
# Returns y as n x p
simulate_factor_data <- function(n, Lambda, Sigma) {
  m <- nrow(Lambda)
  H <- ncol(Lambda)

  Eta <- matrix(rnorm(n * H), nrow = n, ncol = H)   # n x H
  Eps <- matrix(rnorm(n * m), nrow = n, ncol = m)   # n x m

  sigma_sd <- sqrt(diag(Sigma))
  Eps <- sweep(Eps, 2, sigma_sd, FUN = "*")

  Y <- Eta %*% t(Lambda) + Eps
  colnames(Y) <- paste0("y", seq_len(m))

  Y
}

# Convenience wrapper to package truth
make_truth_object <- function(Lambda, Sigma, scenario, n, 
  rep_id, m, H, pivots, below_probs) {
  list(
    scenario = scenario,
    n = n,
    rep = rep_id,
    m = m,
    H = H,
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
generate_scenario1 <- function(m, H) {
  # this is exactly PLT
  pivots <- seq_len(H)
  
  # few dense columns, others moderately sparse
  Delta <- make_delta_from_pivots(m, H, pivots, 0.95)
  Lambda <- make_lambda_from_delta(Delta)
  Sigma <- diag(1,m)

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
generate_scenario2 <- function(m, H,m, pivots) {
  
  # all columns sparse below pivot
  below_probs <- rep(0.5, H)

  Delta <- make_delta_from_pivots(m, H, pivots, below_probs)
  Lambda <- make_lambda_from_delta(Delta)
  Sigma <- diag(1,m)

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
generate_delta_scenario3 <- function(m, H, blocks_membership) {
  
  nblocks <- unique(blocks_membership) |> length()
  Delta0 <- model.matrix(~ factor(blocks_membership) - 1) |>
    matrix(nrow = m, ncol = nblocks)
  # Distribute H factors across the nblocks as evenly as possible
  cols_per_block <- rep(floor(H / nblocks), nblocks)
  remainder_H <- H - sum(cols_per_block)
  if (remainder_H > 0) {
    cols_per_block[seq_len(remainder_H)] <- cols_per_block[seq_len(remainder_H)] + 1
  }

    Delta <- Delta0[, rep(seq_len(nblocks), cols_per_block), 
      drop = FALSE]
  Delta
}
Delta1 <- generate_delta_scenario3(m = m_dim[1], H = H_dim[3,1],
  blocks_membership = c(rep(1, 8), rep(2, 6), rep(3, 6)))
Delta1[1,2] <- Delta1[9,4] <- 0
Delta1
Delta2 <- generate_delta_scenario3(m = m_dim[2], H = H_dim[3,2],
  blocks_membership = c(rep(1, 20), rep(2, 15), rep(3, 15)))
Delta2[1,2:4] <- 0
Delta2[2,3:4] <- 0
Delta2[3,4] <- 0
Delta2[21,6:8] <- 0
Delta2[22,7:8] <- 0
Delta2[23,8] <- 0
Delta2[36,9:10] <- 0
Delta2[37,10] <- 0
Delta2
Delta3 <- generate_delta_scenario3(m = m_dim[3], H = H_dim[3,3],
  blocks_membership = c(rep(1, 40), rep(2, 30), rep(3, 30)))
Delta3[1:5,1:5] <- Delta3[1:5,1:5]*lower.tri(Delta3[1:5,1:5],
   diag = TRUE)
Delta3[41:45,6:10] <- Delta3[41:45,6:10]*lower.tri(Delta3[41:45,6:10],
   diag = TRUE)
Delta3[71:75,11:15] <- Delta3[71:75,11:15]*lower.tri(Delta3[71:75,11:15],
   diag = TRUE)


generate_scenario3 <- function(Delta) {
  
  Lambda <- make_lambda_from_delta(Delta, pivot_scale = 1.1, 
    offpivot_scale = 0.5)
  Sigma <- diag(1,nrow(Delta))

  list(
    Delta = Delta,
    Lambda = Lambda,
    Sigma = Sigma,
    pivots = get_pivots(Delta),
    below_probs = 1
  )
}
# Scenario 4:
# many columns and many zeroes below columns
# different degrees of sparsity, many sparse columns
generate_scenario4 <- function(m, H, m, pivots) {
  
  # all columns sparse below pivot
  below_probs <- seq(0.15, 0.40, length.out = H)

  Delta <- make_delta_from_pivots(m, H, pivots, below_probs)
  Lambda <- make_lambda_from_delta(Delta)
  Sigma <- diag(1,m)

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

generate_scenario <- function(scenario, m, H) {
  switch(
    scenario,
    scenario1_few_dense_columns = generate_scenario1(m, H),
    scenario2_sparse_below      = generate_scenario2(m, H),
    scenario3_block_covariance  = generate_scenario3(m, H),
    scenario4_many_sparse_columns = generate_scenario4(m, H),
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