devtools::load_all("/Users/antonio/Dropbox/1_Ricerca/SylviaFS/R/gltfactor", quiet = TRUE)

set.seed(1)
T <- 50; m <- 6; H <- 3
Lambda_true <- matrix(c(1,0, 0, 
                      0.6,0.9, 0,
                         0.8,0.2, 0,
                         0.4,0, 1,
                         0, -0.6, -0.8, 
                         0, 0.1, 0.2), nrow=m, ncol=H, byrow=TRUE)
Eta_true    <- matrix(rnorm(H*T), H, T)
sigma2_true <- rep(0.5, m)
y <- t(Lambda_true %*% Eta_true) +
     matrix(rnorm(T*m, sd=sqrt(sigma2_true)), T, m, byrow=TRUE)
Delta <- matrix(as.integer(Lambda_true != 0), m, H)

prior_full <- list(a_nu=1, b_nu=1, alpha=1, beta=1,
                   kappa=10, a_sigma=2, b_sigma=1, a_s=2, b_s=1)

fit <- gltfa_cpp(
  y       = y,
  mcmc    = list(niter=5000, nburn=1000, thin=2),
  model   = list(T=T, m=m, Hmax=5L),
  prior   = prior_full,
  init    = list(H=H, nu=0.5, ell=c(1,2,3), tau=c(0.5,0.5, 0.5),
                 Delta=Delta, Lambda=Lambda_true,
                 Eta=Eta_true, sigma2=sigma2_true),
  control = list(store_draws=TRUE, store_eta=TRUE,
                 print_every=25L, verbose=FALSE, seed=1L)
)

fit <- gltfa_cpp(
  y       = y,
  mcmc    = list(niter=5000, nburn=1000, thin=2),
  model   = list(T=T, m=m, Hmax=5L),
  prior   = prior_full,
  init    = list(H=H, nu=0.5, ell=c(1,2,3), tau=c(0.5,0.5, 0.5),
                 Delta=Delta, Lambda=Lambda_true,
                 Eta=Eta_true, sigma2=sigma2_true),
  control = list(store_draws=TRUE, store_eta=TRUE,
                 print_every=25L, verbose=FALSE, seed=1L)
)
fit <- gltfa_cpp(
  y       = y,
  mcmc    = list(niter=5000, nburn=1000, thin=2),
  model   = list(T=T, m=m, Hmax=5L),
  prior   = prior_full,
  init    = list(H=H, nu=0.5, ell=c(1,2,3), tau=c(0.5,0.5, 0.5),
                 Delta=Delta, Lambda=Lambda_true,
                 Eta=Eta_true, sigma2=sigma2_true),
  control = list(store_draws=TRUE, store_eta=TRUE,
                 print_every=25L, verbose=FALSE, seed=1L)
)

cat("=== Lambda (last draw) ===\n")
print(round(fit$draws$Lambda[[25]], 3))
cat("\n=== True Lambda ===\n")
print(Lambda_true)
cat("\n=== sigma2 posterior mean ===\n")
print(round(colMeans(fit$draws$sigma2), 3))
cat("\n=== Eta[,1:5] (last draw) ===\n")
print(round(fit$draws$Eta[[25]][, 1:5], 3))
cat("\n=== True Eta[,1:5] ===\n")
print(round(Eta_true[, 1:5], 3))
