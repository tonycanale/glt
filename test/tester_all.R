# =============================================================================
# tester_all.R  –  Full MCMC check for gltfactor
# =============================================================================

devtools::load_all("/Users/antonio/github/gltfactor/", recompile = TRUE)
#devtools::load_all("/Users/antonio/github/gltf/", recompile = TRUE)

# -----------------------------------------------------------------------------
# 1.  Data-generating process
# -----------------------------------------------------------------------------
set.seed(1)
T_obs <- 100; m <- 6; H <- 3

Lambda_true <- matrix(
  c(1.4, 0.0,  0.0,
    0.6, 0.9,  0.0,
    0.8, 1.2,  0.0,
    1.4, 0.0,  1.2,
    0.0,  -1.6, -1.8,
    0.0,   0.8,  1.1),
  nrow = m, ncol = H, byrow = TRUE
)
Eta_true    <- matrix(rnorm(H * T_obs), H, T_obs)
sigma2_true <- rep(0.5, m)

total_true <- tcrossprod(Lambda_true) + diag(sigma2_true)
total_true

y <- t(Lambda_true %*% Eta_true) +
  matrix(rnorm(T_obs * m, sd = sqrt(sigma2_true)), T_obs, m, byrow = TRUE)
Delta <- matrix(as.integer(Lambda_true != 0), m, H)

prior_full <- list(a_nu = 1, b_nu = 2, alpha = 2, beta = 2,
                   kappa = 5, a_sigma = 1, b_sigma = 2)

# -----------------------------------------------------------------------------
# 2.  Run the sampler
# -----------------------------------------------------------------------------
fit <- gltfa(
  y       = y,
  mcmc    = list(niter = 6000, nburn = 1000, thin = 2),
  Hmax    = 6,
  prior   = prior_full,
  init    = list(H = H, nu = 0.5, ell = c(1, 2, 4), tau = 0.75*rep(1,H),
                 Delta = Delta, Lambda = Lambda_true,
                 Eta = Eta_true, sigma2 = sigma2_true),
  control = list(store_draws = TRUE, store_eta = TRUE,
                 print_every = 250L,  seed = 1L),
  fixed   = list(pivots=TRUE, H=TRUE, Eta=TRUE, LambdaSigma=TRUE)
)


par(mfrow=c(1,2))
image(Delta)
apply(array(unlist(fit$draws$Delta),dim=c(6,3,2500)),c(1,2),mean) |> image()
apply(array(unlist(fit$draws$Delta),dim=c(6,3,2500)),c(1,2),mean)
par(mfrow=c(1,1))




plot(fit$draws$H, type = "l", main = "Trace: H", xlab = "iter", ylab = "H", col = "steelblue")

var <- get_variance(fit)

par(mfrow=c(1,2))
image(tcrossprod(Lambda_true))
estLLt <- apply(var$LLt,c(1,2),mean)
image(estLLt)
image(tcrossprod(Lambda_true) + sigma2_true)
estOmega <- apply(var$total,c(1,2),mean)
image(estOmega)
par(mfrow=c(1,1))


nsave <- fit$meta$nsave
cat(sprintf("\n=== Meta: niter=%d  nburn=%d  thin=%d  nsave=%d ===\n\n",
            fit$meta$niter, fit$meta$nburn, fit$meta$thin, nsave))

# -----------------------------------------------------------------------------
# 3.  Acceptance rate (nu MH step)
# -----------------------------------------------------------------------------
accept_rate <- mean(fit$accept)
cat(sprintf("Nu MH acceptance rate: %.3f\n\n", accept_rate))

# -----------------------------------------------------------------------------
# 4.  Traceplots  (base-R, 3x3 grid)
# -----------------------------------------------------------------------------
draws <- fit$draws

# 4a. H trace
op <- par(mfrow = c(3, 3), mar = c(3, 3, 2, 1))

plot(draws$H, type = "l", main = "Trace: H", xlab = "iter", ylab = "H", col = "steelblue")
abline(h = H, col = "red", lty = 2)

# 4b. nu trace
plot(draws$nu, type = "l", main = "Trace: nu", xlab = "iter", ylab = "nu", col = "steelblue")

# 4c–4h. sigma2[j] traces  (all 6 in remaining panels)
for (j in seq_len(m)) {
  plot(draws$sigma2[, j], type = "l",
       main = sprintf("Trace: sigma2[%d]", j),
       xlab = "iter", ylab = "", col = "steelblue")
  abline(h = sigma2_true[j], col = "red", lty = 2)
}
par(op)



# 4c. Lambda traceplots for the first 3 factors, diagonal entries
#     (only meaningful draws where H >= 3)
idx3 <- which(draws$H >= 3)
cat(sprintf("Draws with H >= 3: %d / %d\n\n", length(idx3), nsave))

op2 <- par(mfrow = c(H, H), mar = c(3, 3, 2, 1))
for (h in seq_len(H)) {
  for (i in seq_len(H)) {          # row = variable (use first H vars for brevity)
    vals <- sapply(idx3, function(s) {
      Ls <- draws$Lambda[[s]]
      if (ncol(Ls) >= h) Ls[i, h] else NA_real_
    })
    plot(vals, type = "l",
         main = sprintf("Lambda[%d,%d]", i, h),
         xlab = "iter (H>=3)", ylab = "", col = "steelblue")
    abline(h = Lambda_true[i, h], col = "red", lty = 2)
  }
}
par(op2)

# -----------------------------------------------------------------------------
# 5.  Posterior means
# -----------------------------------------------------------------------------

## 5a.  sigma2
sigma2_pm <- colMeans(draws$sigma2)
cat("=== sigma2: posterior mean vs truth ===\n")
print(round(rbind(posterior_mean = sigma2_pm, truth = sigma2_true), 3))

## 5b.  nu
nu_pm <- mean(draws$nu)
cat(sprintf("\n=== nu: posterior mean = %.3f ===\n\n", nu_pm))

## 5c.  H marginal distribution
H_tab <- table(draws$H)
cat("=== Marginal distribution of H ===\n")
print(H_tab)
cat(sprintf("  True H = %d\n\n", H))

## 5d.  Lambda posterior mean (draws where H >= H)
#  Note: label-switching is *not* corrected here; columns may be permuted.
Lambda_pm_raw <- Reduce("+", lapply(idx3, function(s) {
  Ls <- draws$Lambda[[s]]
  Ls[, seq_len(H), drop = FALSE]     # take first H columns
})) / length(idx3)

cat("=== Lambda: posterior mean (first H cols, no sign/permutation correction) ===\n")
print(round(Lambda_pm_raw, 3))
cat("\n=== Lambda: truth ===\n")
print(Lambda_true)

## 5e.  Absolute error per entry (useful at a glance)
cat("\n=== |Lambda_pm - Lambda_true| ===\n")
print(round(abs(Lambda_pm_raw - Lambda_true), 3))

# -----------------------------------------------------------------------------
# 6.  Posterior inclusion probabilities  (Delta)
# -----------------------------------------------------------------------------
Delta_pm <- Reduce("+", lapply(idx3, function(s) {
  Ds <- draws$Delta[[s]]
  Ds[, seq_len(H), drop = FALSE]
})) / length(idx3)

cat("\n=== Posterior inclusion probabilities (Delta) ===\n")
print(round(Delta_pm, 3))
cat("\n=== True Delta ===\n")
print(Delta)

# -----------------------------------------------------------------------------
# 7.  Quick visual: posterior mean Lambda vs truth (side-by-side image)
# -----------------------------------------------------------------------------
op3 <- par(mfrow = c(1, 2))
image(t(Lambda_pm_raw), main = "Posterior mean Lambda\n(first H cols)",
      col = hcl.colors(20, "RdBu", rev = TRUE),
      zlim = c(-1.5, 1.5), axes = FALSE)
image(t(Lambda_true), main = "True Lambda",
      col = hcl.colors(20, "RdBu", rev = TRUE),
      zlim = c(-1.5, 1.5), axes = FALSE)
par(op3)

# -----------------------------------------------------------------------------
# 8.  Effective sample sizes  (coda if available, else naive)
# -----------------------------------------------------------------------------
if (requireNamespace("coda", quietly = TRUE)) {
  ess_sigma2 <- coda::effectiveSize(coda::mcmc(draws$sigma2))
  ess_nu     <- coda::effectiveSize(coda::mcmc(draws$nu))
  cat("\n=== ESS: sigma2 ===\n");  print(round(ess_sigma2, 0))
  cat(sprintf("ESS: nu = %.0f\n", ess_nu))
} else {
  cat("\n(install 'coda' for effective sample sizes)\n")
}
