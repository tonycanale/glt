############################################################
#
# Simulate scenario and data
#
############################################################

rm(list = ls())
library(sparvaride)
## ----------------------------- ##
## 0. User-facing configuration ##
## ----------------------------- ##ì

# Output directory
out_dir <- "simulated_data_4"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Number of replicates
n_reps <- 20

# Dimensions 
m_dim <- c(20, 50, 100)
T_dim <- cbind(rep(100,3), rep(200,3))
H_dim <- matrix(
  c(3,8,10,
    4,10,15,
    4,10,15,
    7, 15, 30),
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
make_lambda_from_delta_old <- function(Delta, negative_prob = 0.1, mean =1, sd = sqrt(0.1)){
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

make_lambda_from_delta_old2 <- function(Delta, mean = 0.5, sd = 1, thresh= 0.02){
  m <- nrow(Delta)
  H <- ncol(Delta)
  Lambda <- matrix(rnorm(m * H, mean = mean, sd = sd),
    nrow = m, ncol = H) * Delta
  totvar <- (rowSums(Lambda^2) + 1)
  smallish <- abs(Lambda / matrix(rep(totvar, H),m,H) )
  small <- (smallish < thresh) & (Lambda != 0)
  if(any(small)) {
  Lambda[small] <- (min(abs(totvar)) + 0.000001) * sign(Lambda[small])
  } 
  Lambda
}

make_lambda_from_delta <- function(Delta, sd = sqrt(0.5), thresh= 0.02){
  m <- nrow(Delta)
  H <- ncol(Delta)
  means <- sample(c(0,1.5,-1.5), size=H, replace=TRUE)
  Lambda <- (matrix(rep(means,each=m),m,H) + matrix(rnorm(m * H, mean = 0, sd = sd),
    nrow = m, ncol = H, byrow = TRUE)) * Delta
  totvar <- (rowSums(Lambda^2) + 1)
  smallish <- abs(Lambda / matrix(rep(totvar, H),m,H) )
  small <- (smallish < thresh) & (Lambda != 0)
  if(any(small)) {
  Lambda[small] <- (min(abs(totvar)) + 0.000001) * sign(Lambda[small])
  } 
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

  pivots <- c(1, 3, 5, 7)
  if(H==10) pivots <- c(pivots, 9, 10, 11, 13, 14, 17)
  if(H==15) pivots <- c(pivots, 9, 10, 11, 13, 14, 17, 19, 20,21,22, 23)
  
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
set.seed(1)
Delta1 <- generate_delta_scenario3(m = m_dim[1], H = H_dim[3,1],
  blocks_membership = c(rep(1, 10), rep(2, 10)))
Delta1[1,2] <- Delta1[11,4] <- 0
Delta1 |> plot_real_matrix()

Delta2 <- generate_delta_scenario3(m = m_dim[2], H = H_dim[3,2],
  blocks_membership = rep(1:2, each=m_dim[2]/2))
Delta2[1,2:5] <- 0
Delta2[2,3:5] <- 0
Delta2[3,4:5] <- 0
Delta2[4,5] <- 0
Delta2[26,7:10] <- 0
Delta2[27,8:10] <- 0
Delta2[28,9:10] <- 0
Delta2[29,10] <- 0
Delta2 |> plot_real_matrix()

Delta3 <- generate_delta_scenario3(m = m_dim[3], H = H_dim[3,3],
  blocks_membership = rep(1:4, each=m_dim[3]/4))
Delta3 |> plot_real_matrix()
Delta3[1:4,1:4] <- Delta3[1:4,1:4]*lower.tri(Delta3[1:4,1:4],
   diag = TRUE)
Delta3[26:28,5:8] <- Delta3[26:28,5:8]*lower.tri(Delta3[26:28,5:8],
   diag = TRUE)
Delta3[51:54,9:12] <- Delta3[51:54,9:12]*lower.tri(Delta3[51:54,9:12],
   diag = TRUE)
Delta3[76:78,13:15] <- Delta3[76:78,13:15]*lower.tri(Delta3[76:78,13:15],
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
generate_scenario4_old <- function(m, H, mintau = 0.25, maxtau = 0.4) {
  
  # this is exactly PLT
  #pivots <- seq_len(H) #OLD

  # Sample H-1 random pivots between 2: H + m/5
  pivots <- c(1,sort(sample(2:(H + floor(m/5)), H-1)))

  # all columns sparse below pivot
  below_probs <- seq(mintau, maxtau, length.out = H)


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


generate_scenario4 <- function(m, H) {
  
  # this is exactly PLT
  #pivots <- seq_len(H) #OLD

  pivots <- seq(1, H + floor(m/5), length=H)
  act_loadings <- max(2, floor(m/5))
  Delta <- make_delta_from_pivots(m, H, pivots, below_probs = rep(0,H))
  #put `blocksize` ones below each pivot
  for (j in seq_len(H)) {
    if(pivots[j] < m) {
      idx <- (pivots[j] + 1):min(m, pivots[j] + act_loadings)
      Delta[idx, j] <- 1
    }
  }
  
  
  Lambda <- make_lambda_from_delta(Delta)
  Sigma <- diag(1,m)

  list(
    Delta = Delta,
    Lambda = Lambda,
    Sigma = Sigma,
    pivots = pivots,
    below_probs = rep(0,H)
  )
  }



s1 <- generate_scenario1(m_dim[2], H_dim[1,2])
s2 <- generate_scenario2(m_dim[2], H_dim[2,2])
s3 <- generate_scenario3(Delta2)
s4 <- generate_scenario4(m_dim[2], H_dim[4,2])

pdf(file.path(out_dir,"example_plots.pdf"))
s1$Delta |> plot_real_matrix()
s1$Lambda |> plot_real_matrix()
s2$Delta |> plot_real_matrix()
s2$Lambda |> plot_real_matrix()
s3$Delta |> plot_real_matrix()
s3$Lambda |> plot_real_matrix()
s4$Delta |> plot_real_matrix()
s4$Lambda |> plot_real_matrix()
dev.off()

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
# Pivots: Randomly sampled for each replicate, but always between 1 and H + m/5
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