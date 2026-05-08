set.seed(123)
sc <- generate_scenario4(50, 20, mintau = 0.25)
plot_real_matrix(sc$Delta)
plot_real_matrix(sc$Lambda)
Y  <- simulate_factor_data(100, sc$Lambda, sc$Sigma)


prior_fixed <- list(a_nu = 1, b_nu = 1, alpha = 1, beta = 1,
                    kappa = 1, a_sigma = 1, b_sigma = 0.3)

mcmc_cfg <- list(niter = 5000, nburn = 1000, thin = 2)

fit <- gltfa(
            y       = Y,
            mcmc    = mcmc_cfg,
            Hmax    = 30,
            prior   = prior_fixed,
            init    = list(H = 20, Delta=sc$Delta, Lambda = sc$Lambda, Sigma = sc$Sigma),
            control = list(store_draws = TRUE, store_eta = FALSE,
                           print_every = 1000L,
                           seed        = 12345L),
            fixed = list(H = TRUE, pivots = TRUE)
          )

met <- compute_metrics(fit, sc$Lambda, sc$Delta)
plot_real_matrix(met$Delta_post_aligned)
plot_real_matrix(met$Delta_post)
met$auc_aligned




set.seed(123)
sc <- generate_scenario4(50, 20, mintau = 0.15)
Y  <- simulate_factor_data(100, sc$Lambda, sc$Sigma)


prior_fixed <- list(a_nu = 1, b_nu = 1, alpha = 1, beta = 1,
                    kappa = 1, a_sigma = 1, b_sigma = 0.3)

mcmc_cfg <- list(niter = 5000, nburn = 1000, thin = 2)

fit2 <- gltfa(
            y       = Y,
            mcmc    = mcmc_cfg,
            Hmax    = 30,
            prior   = prior_fixed,
            init    = list(H = 20, Delta=sc$Delta, Lambda = sc$Lambda, Sigma = sc$Sigma),
            control = list(store_draws = TRUE, store_eta = FALSE,
                           print_every = 1000L,
                           seed        = 12345L),
            fixed = list(H = TRUE, pivots = TRUE)
          )

met2 <- compute_metrics(fit2, sc$Lambda, sc$Delta)
plot_real_matrix(met2$Delta_post_aligned)
plot_real_matrix(met2$Delta_post)
chieilpivot <- function(x) which(x==1) 
apply(met2$Delta_post,1, chieilpivot)
met2$auc_aligned

