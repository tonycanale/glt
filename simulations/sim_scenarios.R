############################################################
#
# Simulate scenario and data
#
############################################################

rm(list = ls())

## ----------------------------- ##
## 0. User-facing configuration ##
## ----------------------------- ##ì

# Output directory
out_dir <- "simulated_data"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Number of replicates
n_reps <- 20

# Dimensions 
m_dim <- c(20, 50, 100)
T_dim <- cbind(m_dim, 2 * m_dim)
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
  "scenario4_many_sparse_columns" ## H < (m-1)/2
)

## -------------------------------- ##
## 1. Helper functions: basic tools ##
## -------------------------------- ##

# Build a m x H indicator matrix Delta with ordered pivots and given nonzero probabilities below pivots.
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

# Given Delta, generate Lambda with Gaussian nonzero loadings and negative signes with some probability.
make_lambda_from_delta <- function(Delta, negative_prob = 0.1, mean =1, sd = 0.1){
  m <- nrow(Delta)
  H <- ncol(Delta)
  Lambda <- matrix(rnorm(m * H, mean = mean, sd = sd),
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

# Wrapper to package truth
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
# lower triangular dense factor loading
generate_scenario1 <- function(m, H) {
  # this is exactly PLT
  pivots <- seq_len(H)
  
  # few dense columns, others moderately sparse
  Delta <- make_delta_from_pivots(m, H, pivots, rep(0.95,H))
  Lambda <- make_lambda_from_delta(Delta)
  Sigma <- diag(1,m)

  list(
    Delta = Delta,
    Lambda = Lambda,
    Sigma = Sigma,
    pivots = pivots,
    below_probs = rep(0.95,H)
  )
}

# Scenario 2:
# lower triangular factor loading sparse below
generate_scenario2 <- function(m, H) {

  pivots <- c(1, 3, 5, 7, 9)
  if(H==10) pivots <- c(pivots, 10, 11, 13, 14, 17)
  if(H==15) pivots <- c(pivots, 10, 11, 13, 14, 17, 19, 20,21,22, 23)
  
  # all columns sparse below pivot
  below_probs <- rep(0.5, H)

  Delta <- make_delta_from_pivots(m, H, pivots, below_probs)
  Delta[which(rowSums(Delta) == 0),1] <- 1 # check for empty rows, if any, put a 1
  while(!counting_rule_holds(Delta)) {
    Delta <- make_delta_from_pivots(m, H, pivots, below_probs)
    Delta[which(rowSums(Delta) == 0),1] <- 1 # check for empty rows, if any, put a 1
  }

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
  
  Lambda <- make_lambda_from_delta(Delta)
    
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
generate_scenario4 <- function(m, H) {
  
  # this is exactly PLT
  pivots <- seq_len(H)
  

  # all columns sparse below pivot
  below_probs <- seq(0.15, 0.40, length.out = H)


  Delta <- make_delta_from_pivots(m, H, pivots, below_probs)
  Delta[which(rowSums(Delta) == 0),1] <- 1 # check for empty rows, if any, put a 1
  while(!counting_rule_holds(Delta)) {
    Delta <- make_delta_from_pivots(m, H, pivots, below_probs)  
  Delta[which(rowSums(Delta) == 0),1] <- 1 # check for empty rows, if any, put a 1
  }
  
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



## ------------------------------------------------ ##
## 3. Simulate data over all configurations         ##
## ------------------------------------------------ ##

# Scenario 1 -------------------------------------------------------
for(dimens in 1:3){
  for(samplesize in 1:2){
    cat("\nSimulating scenario 1 with m =", m_dim[dimens], "H =", H_dim[1,dimens], "T =", T_dim[dimens,samplesize], "\n")
    for (rep in seq_len(n_reps)) {
      set.seed(dimens + samplesize + rep) # ensure different seed for each scenario and replicate
      sc <- generate_scenario1(m_dim[dimens], H_dim[1,dimens])
      Y  <- simulate_factor_data(T_dim[dimens,samplesize], sc$Lambda, sc$Sigma)
      write.csv(Y, file = file.path(out_dir,
        sprintf("scenario1_Y_m%d_H%d_T%d_rep%02d.csv", m_dim[dimens], H_dim[1,dimens], T_dim[dimens,samplesize], rep)),
        row.names = FALSE)
      write.csv(sc$Delta, file = file.path(out_dir,
        sprintf("scenario1_Delta_m%d_H%d_T%d_rep%02d.csv", m_dim[dimens], H_dim[1,dimens], T_dim[dimens,samplesize], rep)),
        row.names = FALSE)
      write.csv(sc$Lambda, file = file.path(out_dir,
        sprintf("scenario1_Lambda_m%d_H%d_T%d_rep%02d.csv", m_dim[dimens], H_dim[1,dimens], T_dim[dimens,samplesize], rep)),
        row.names = FALSE)
    }
    }
  }

# Scenario 2 -------------------------------------------------------
# H_dim[2,]: H values for scenario2; H must be 5, 10, or 15 (pivots are hardcoded inside)
for(dimens in 1:3){
  for(samplesize in 1:2){
    cat("\nSimulating scenario 2 with m =", m_dim[dimens], "H =", H_dim[2,dimens], "T =", T_dim[dimens,samplesize], "\n")
    for (rep in seq_len(n_reps)) {
      set.seed(dimens + samplesize + rep) # ensure different seed for each scenario and replicate
      sc <- generate_scenario2(m_dim[dimens], H_dim[2,dimens])
      Y  <- simulate_factor_data(T_dim[dimens,samplesize], sc$Lambda, sc$Sigma)
      write.csv(Y, file = file.path(out_dir,
        sprintf("scenario2_Y_m%d_H%d_T%d_rep%02d.csv", m_dim[dimens], H_dim[2,dimens], T_dim[dimens,samplesize], rep)),
        row.names = FALSE)
      write.csv(sc$Delta, file = file.path(out_dir,
        sprintf("scenario2_Delta_m%d_H%d_T%d_rep%02d.csv", m_dim[dimens], H_dim[2 ,dimens], T_dim[dimens,samplesize], rep)),
        row.names = FALSE)
      write.csv(sc$Lambda, file = file.path(out_dir,
        sprintf("scenario2_Lambda_m%d_H%d_T%d_rep%02d.csv", m_dim[dimens], H_dim[2,dimens], T_dim[dimens,samplesize], rep)),
        row.names = FALSE)
    }
    }
  }


# Scenario 3 -------------------------------------------------------
# Uses pre-built Delta1 (m=20), Delta2 (m=50), Delta3 (m=100)
# H_dim[3,] values are implicit in Delta dimensions
for(dimens in 1:3){
  for(samplesize in 1:2){
    cat("\nSimulating scenario 3 with m =", m_dim[dimens], "H =", H_dim[2,dimens], "T =", T_dim[dimens,samplesize], "\n")
    for (rep in seq_len(n_reps)) {
      set.seed(dimens + samplesize + rep) # ensure different seed for each scenario and replicate
      if(dimens==1) sc <- generate_scenario3(Delta1)
      if(dimens==2) sc <- generate_scenario3(Delta2)
      if(dimens==3) sc <- generate_scenario3(Delta3)
      Y  <- simulate_factor_data(T_dim[dimens,samplesize], sc$Lambda, sc$Sigma)
      write.csv(Y, file = file.path(out_dir,
        sprintf("scenario3_Y_m%d_H%d_T%d_rep%02d.csv", m_dim[dimens], H_dim[3,dimens], T_dim[dimens,samplesize], rep)),
        row.names = FALSE)
      write.csv(sc$Delta, file = file.path(out_dir,
        sprintf("scenario3_Delta_m%d_H%d_T%d_rep%02d.csv", m_dim[dimens], H_dim[3,dimens], T_dim[dimens,samplesize], rep)),
        row.names = FALSE)
      write.csv(sc$Lambda, file = file.path(out_dir,
        sprintf("scenario3_Lambda_m%d_H%d_T%d_rep%02d.csv", m_dim[dimens], H_dim[3,dimens], T_dim[dimens,samplesize], rep)),
        row.names = FALSE)
    }
    }
  }

# Scenario 4 -------------------------------------------------------
# H_dim[4,]: H values for scenario4
# Pivots: PLT
for(dimens in 1:3){
  for(samplesize in 1:2){
    cat("\nSimulating scenario 4 with m =", m_dim[dimens], "H =", H_dim[4,dimens], "T =", T_dim[dimens,samplesize], "\n")
    for (rep in seq_len(n_reps)) {
      set.seed(dimens + samplesize + rep) # ensure different seed for each scenario and replicate
      sc <- generate_scenario4(m_dim[dimens], H_dim[4,dimens])
      Y  <- simulate_factor_data(T_dim[dimens,samplesize], sc$Lambda, sc$Sigma)
      write.csv(Y, file = file.path(out_dir,
        sprintf("scenario4_Y_m%d_H%d_T%d_rep%02d.csv", m_dim[dimens], H_dim[4,dimens], T_dim[dimens,samplesize], rep)),
        row.names = FALSE)
      write.csv(sc$Delta, file = file.path(out_dir,
        sprintf("scenario4_Delta_m%d_H%d_T%d_rep%02d.csv", m_dim[dimens], H_dim[4,dimens], T_dim[dimens,samplesize], rep)),
        row.names = FALSE)
      write.csv(sc$Lambda, file = file.path(out_dir,
        sprintf("scenario4_Lambda_m%d_H%d_T%d_rep%02d.csv", m_dim[dimens], H_dim[4,dimens], T_dim[dimens,samplesize], rep)),
        row.names = FALSE)
    }
    }
  }