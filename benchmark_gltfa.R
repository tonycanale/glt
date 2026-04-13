## benchmark_gltfa.R
##
## Tests the full gltfa() --> gltfa_cpp() pipeline with a known data-generating
## process. Checks:
##   1. Package loads cleanly
##   2. gltfa() runs without errors
##   3. Output structure is correct
##   4. Posterior means are in the right ballpark (bias check)
##   5. H mixing: chain visits more than one value of H
##   6. sigma2 posterior mean is close to the true value
##   7. Timing benchmark for a longer run
## ---------------------------------------------------------------------------

devtools::load_all(quiet = TRUE)

## ---------------------------------------------------------------------------
## 1. Data-generating process
## ---------------------------------------------------------------------------
set.seed(42)
T_obs  <- 100
m      <- 8
H_true <- 3

# Lower-triangular loading matrix (GLT structure)
Lambda_true <- matrix(0, m, H_true)
Lambda_true[1, 1] <-  1.0
Lambda_true[2, 1] <-  0.8;  Lambda_true[2, 2] <-  1.0
Lambda_true[3, 1] <- -0.6;  Lambda_true[3, 2] <-  0.4;  Lambda_true[3, 3] <-  1.0
Lambda_true[4, 1] <-  0.5;  Lambda_true[4, 2] <- -0.3;  Lambda_true[4, 3] <-  0.7
Lambda_true[5, 2] <-  1.2;  Lambda_true[5, 3] <- -0.5
Lambda_true[6, 2] <- -0.9;  Lambda_true[6, 3] <-  0.3
Lambda_true[7, 3] <-  0.8
Lambda_true[8, 3] <-  1.1

sigma2_true <- rep(0.5, m)
Eta_true    <- matrix(rnorm(H_true * T_obs), H_true, T_obs)
noise       <- matrix(rnorm(T_obs * m, sd = sqrt(sigma2_true)), T_obs, m, byrow = TRUE)
y           <- t(Lambda_true %*% Eta_true) + noise

cat("Data: T =", T_obs, ", m =", m, ", H_true =", H_true, "\n\n")

## ---------------------------------------------------------------------------
## 2. Initial values
## ---------------------------------------------------------------------------
ell_init    <- as.integer(1:H_true)
Delta_init  <- matrix(0L, m, H_true)
for (j in seq_len(H_true)) Delta_init[j:m, j] <- 1L
tau_init    <- rep(0.5, H_true)
sigma2_init <- rep(1.0, m)
Lambda_init <- Lambda_true + matrix(rnorm(m * H_true, sd = 0.1), m, H_true)
Lambda_init[upper.tri(Lambda_init)] <- 0
Eta_init    <- matrix(rnorm(H_true * T_obs), H_true, T_obs)

## ---------------------------------------------------------------------------
## 3. Short run — correctness checks
## ---------------------------------------------------------------------------
cat("=== Short run (500 iter, 250 burn-in) ===\n")

t_short <- system.time({
  fit_short <- gltfa(
    y       = y,
    mcmc    = list(niter = 500L, nburn = 250L, thin = 1L),
    Hmax    = 8L,
    prior   = list(a_nu = 1, b_nu = 1, alpha = 1, beta = 1,
                   kappa = 10, a_sigma = 2, b_sigma = 1),
    init    = list(H      = H_true,
                   nu     = 0.5,
                   ell    = ell_init,
                   tau    = tau_init,
                   Delta  = Delta_init,
                   Lambda = Lambda_init,
                   Eta    = Eta_init,
                   sigma2 = sigma2_init),
    control = list(store_draws = TRUE, store_eta = TRUE,
                   print_every = 100L, seed = 1L)
  )
})

# -- Structure checks ---------------------------------------------------------
cat("\n--- Structure checks ---\n")
stopifnot(inherits(fit_short, "gltfit"))
stopifnot(all(c("draws", "last", "meta", "accept") %in% names(fit_short)))
stopifnot(length(fit_short$draws$H)      == fit_short$meta$nsave)
stopifnot(length(fit_short$draws$Lambda) == fit_short$meta$nsave)
cat("  [OK] output structure\n")

# -- H mixing -----------------------------------------------------------------
H_chain <- fit_short$draws$H
cat("\n--- H chain ---\n")
cat("  min =", min(H_chain), "| max =", max(H_chain),
    "| mean =", round(mean(H_chain), 2), "\n")
cat("  values visited:", sort(unique(H_chain)), "\n")
if (length(unique(H_chain)) > 1L) {
  cat("  [OK] H is mixing\n")
} else {
  cat("  [WARN] H is stuck at", unique(H_chain), "\n")
}

# -- sigma2 posterior mean ----------------------------------------------------
sigma2_pm <- colMeans(fit_short$draws$sigma2)
cat("\n--- sigma2 posterior mean vs truth ---\n")
df_sig <- data.frame(
  var      = paste0("y", seq_len(m)),
  truth    = sigma2_true,
  post_mn  = round(sigma2_pm, 3),
  abs_err  = round(abs(sigma2_pm - sigma2_true), 3)
)
print(df_sig, row.names = FALSE)
rmse_s <- sqrt(mean((sigma2_pm - sigma2_true)^2))
cat("  RMSE(sigma2):", round(rmse_s, 4), "\n")
if (rmse_s < 0.5) cat("  [OK] sigma2 in right ballpark\n") else cat("  [WARN] sigma2 far from truth\n")

# -- Lambda last draw ---------------------------------------------------------
cat("\n--- Lambda last draw (H =", fit_short$last$H, ") ---\n")
print(round(fit_short$last$Lambda, 3))
cat("  Lambda truth:\n")
print(round(Lambda_true, 3))

# -- H acceptance rate --------------------------------------------------------
cat("\n--- H birth/death acceptance rate:",
    round(mean(fit_short$accept[, 1]), 3), "---\n")
cat("  Short run wall time:", round(t_short["elapsed"], 2), "sec\n")

## ---------------------------------------------------------------------------
## 4. Medium benchmark run
## ---------------------------------------------------------------------------
cat("\n=== Medium run (2000 iter, 1000 burn-in, thin = 2) ===\n")

t_med <- system.time({
  fit_med <- gltfa(
    y       = y,
    mcmc    = list(niter = 2000L, nburn = 1000L, thin = 2L),
    Hmax    = 8L,
    prior   = list(a_nu = 1, b_nu = 1, alpha = 1, beta = 1,
                   kappa = 10, a_sigma = 2, b_sigma = 1),
    init    = list(H      = H_true,
                   nu     = 0.5,
                   ell    = ell_init,
                   tau    = tau_init,
                   Delta  = Delta_init,
                   Lambda = Lambda_init,
                   Eta    = Eta_init,
                   sigma2 = sigma2_init),
    control = list(store_draws = TRUE, store_eta = FALSE,
                   print_every = 500L, seed = 42L)
  )
})

H_chain_med   <- fit_med$draws$H
sigma2_pm_med <- colMeans(fit_med$draws$sigma2)

cat("  H chain: mean =", round(mean(H_chain_med), 2),
    "| range [", min(H_chain_med), ",", max(H_chain_med), "]\n")
cat("  sigma2 posterior mean:", round(sigma2_pm_med, 3), "\n")
cat("  sigma2 truth:         ", sigma2_true, "\n")
cat("  RMSE(sigma2):", round(sqrt(mean((sigma2_pm_med - sigma2_true)^2)), 4), "\n")
cat("  Medium run wall time:", round(t_med["elapsed"], 2), "sec\n")

## ---------------------------------------------------------------------------
## 5. Timing summary
## ---------------------------------------------------------------------------
cat("\n=== Timing summary ===\n")
cat("  Short run (500 iter)  :", round(t_short["elapsed"], 2), "sec\n")
cat("  Medium run (2000 iter):", round(t_med["elapsed"],   2), "sec\n")
cat("  Iterations/sec (medium):",
    round(2000 / t_med["elapsed"], 1), "\n")

cat("\n=== Done ===\n")
