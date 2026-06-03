############################################################
#
# Run gltfa on each scenario and save results
#
############################################################

devtools::load_all("/Users/antonio/github/glt/", recompile = FALSE)
library(sparvaride)
source("metrics.R")

# ------------------------------------------------------------------ #
# Configuration (must match sim_scenarios.R)                          #
# ------------------------------------------------------------------ #

in_dir  <- "simulated_data_5"
out_dir <- "simulation_results_5_allrandom_2"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

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



prior_fixed <- list(a_nu = 1, b_nu = 1, alpha = 1, beta = 1,
                    kappa = 1, a_sigma = 1, b_sigma = 0.3)

mcmc_cfg <- list(niter = 10000, nburn = 6000, thin = 1)

# ------------------------------------------------------------------ #
# Main loop                                                            #
# ------------------------------------------------------------------ #

n_total <-  3 * 2 * n_reps   # scenarios x dims x sample sizes x reps
results <- matrix(NA_real_, nrow = n_total,
                  ncol = 14,
                  dimnames = list(NULL,
                    c("scenario", "m", "H", "T", "rep",
                      "H_median", "H_mode", "H_true_post_prob", 
                      "stein_loss", "stein_loss_sigma", "auc", "auc_aligned", "nsave", "status")))
idx <- 1L