############################################################
#
# Run gltfa on each scenario and save results
#
############################################################

devtools::load_all("/Users/antonio/github/glt/", recompile = FALSE)
library(sparvaride)
# ------------------------------------------------------------------ #
# Configuration (must match sim_scenarios.R)                          #
# ------------------------------------------------------------------ #

in_dir  <- "simulated_data_3"
out_dir <- "simulation_results_3"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

n_reps <- 10
m_dim  <- c(20, 50, 100)
T_dim <- cbind(rep(100,3), rep(200,3))
H_dim  <- matrix(
  c( 3,  8, 10,
     5, 10, 15,
     5, 10, 15,
     8, 20, 40),
  nrow = 4, ncol = 3, byrow = TRUE
)


prior_fixed <- list(a_nu = 1, b_nu = 1, alpha = 1, beta = 1,
                    kappa = 1, a_sigma = 1, b_sigma = 0.3)

mcmc_cfg <- list(niter = 7000, nburn = 1000, thin = 2)

# ------------------------------------------------------------------ #
# Metrics                                                              #
# ------------------------------------------------------------------ #

compute_metrics <- function(fit, Lambda_true, Delta_true,
                             Sigma_true = rep(1, nrow(Lambda_true))) {

  m   <- nrow(Lambda_true)
  H0  <- ncol(Lambda_true)
  Omega0       <- Lambda_true %*% t(Lambda_true) + diag(Sigma_true)
  
  Delta_draws  <- fit$draws$Delta    # list of nsave matrices (m x H_s)
  
  # Counting rule: discard draws where any row of Delta is zero, or where the counting rule does not hold
  admissible <- rep(TRUE, fit$meta$nsave)
  for(ite in 1:fit$meta$nsave) {
    if(any(rowSums(Delta_draws[[ite]]) == 0)) {
      admissible[ite] <- FALSE
    } else {
      admissible[ite] <- counting_rule_holds(Delta_draws[[ite]])
    }
  }
  # nsave is now different from fit$meta$nsave, as we discard some draws
  nsave        <- sum(admissible)

  Delta_draws <- Delta_draws[admissible]    # list of nsave matrices (m x H_s)
  Lambda_draws <- fit$draws$Lambda[admissible]   # list of nsave matrices (m x H_s)
  sigma2_draws <- fit$draws$sigma2[admissible, , drop = FALSE]    # nsave x m matrix
  H_draws      <- fit$draws$H[admissible]        # integer vector length nsave

  # -- column alignment: match draw cols to Delta_true cols by Hamming --
  align_cols <- function(D_draw, H_s) {
    D_active <- D_draw[, seq_len(H_s), drop = FALSE]
    n_match  <- min(H_s, H0)
    assigned <- integer(n_match)
    used     <- logical(H_s)
    for (j in seq_len(n_match)) {
      sims       <- sapply(seq_len(H_s), function(k)
        if (used[k]) -Inf else sum(D_active[, k] == Delta_true[, j]))
      best         <- which.max(sims)
      assigned[j]  <- best
      used[best]   <- TRUE
    }
    D_aligned <- D_active[, assigned, drop = FALSE]
    if (H_s < H0)
      D_aligned <- cbind(D_aligned, matrix(0L, m, H0 - H_s))
    D_aligned
  }

  # -- Stein loss: tr(A B^{-1}) - log|det(A B^{-1})| - m --
  stein_loss <- function(Omega_hat) {
    ev <- Re(eigen(Omega_hat %*% solve(Omega0), only.values = TRUE)$values)
    sum(ev) - sum(log(pmax(ev, .Machine$double.eps))) - m
  }

  stein_sum <- 0
  Delta_acc <- matrix(0, m, H0)
  Delta_acc_aligned <- matrix(0, m, H0)

  for (s in seq_len(nsave)) {
    H_s    <- H_draws[s]
    Lam_s  <- Lambda_draws[[s]][, seq_len(H_s), drop = FALSE]
    sig2_s <- sigma2_draws[s, ]
    Omega_s <- Lam_s %*% t(Lam_s) + diag(sig2_s)

    stein_sum <- stein_sum + stein_loss(Omega_s)

    Del_s     <- Delta_draws[[s]]
    Delta_acc <- Delta_acc + Del_s
    Delta_acc_aligned <- Delta_acc_aligned + align_cols(Del_s, H_s)
  }

  avg_stein  <- stein_sum / nsave
  Delta_post <- Delta_acc / nsave
  Delta_post_aligned <- Delta_acc_aligned / nsave

  # -- AUC via Wilcoxon rank-sum --
  labs  <- as.vector(Delta_true)
  scors <- as.vector(Delta_post)
  scors_aligned <- as.vector(Delta_post_aligned)
  n1    <- sum(labs == 1L)
  n0    <- sum(labs == 0L)
  auc   <- tryCatch(
    {
      wilcox.test(scors[labs == 1L], scors[labs == 0L],
                exact = FALSE)$statistic / (n1 * n0)
    }, error = function(e) NA_real_
  )

  auc_aligned   <- tryCatch(
    {
      wilcox.test(scors_aligned[labs == 1L], scors_aligned[labs == 0L],
                exact = FALSE)$statistic / (n1 * n0)
    }, error = function(e) NA_real_
  )
 
  auc <- ifelse(is.null(auc), NA_real_, as.numeric(auc))
  auc_aligned <- ifelse(is.null(auc_aligned), NA_real_, as.numeric(auc_aligned))

  list(H_median = median(H_draws), 
       H_mode = as.integer(names(which.max(table(H_draws)))), 
       H_true_post_prob = mean(H_draws == H0),
       stein_loss = avg_stein, 
       auc = auc,
       auc_aligned = auc_aligned,
       Delta_post_aligned = Delta_post_aligned,
       Delta_post = Delta_post)
}

# ------------------------------------------------------------------ #
# Main loop                                                            #
# ------------------------------------------------------------------ #

n_total <-  3 * 2 * n_reps   # scenarios x dims x sample sizes x reps
results <- matrix(NA_real_, nrow = n_total,
                  ncol = 12,
                  dimnames = list(NULL,
                    c("scenario", "m", "H", "T", "rep",
                      "H_median", "H_mode", "H_true_post_prob", 
                      "stein_loss", "auc", "auc_aligned", "status")))
idx <- 1L