############################################################
#
# Run gltfa fixing different things, to understand where the errors are coming from
#
############################################################

source("run_simulations_1234.R")

mcmc_cfg <- list(niter = 15000, nburn = 0, thin = 1)

scen <- 2
dimens <- 3
rep <- 1
samplezie <- 1

m <- m_dim[dimens]
H <- H_dim[scen, dimens]
Hmax_run <- min(2L * H, m - 1L)
TT <- T_dim[dimens, samplesize]

        tag <- sprintf("scenario%d_m%d_H%d_T%d_rep%02d", scen, m, H, TT, rep)
        cat("\n[", format(Sys.time(), "%H:%M:%S"), "]  Fitting", tag, "\n")

        # -- load data --
        Y           <- as.matrix(read.csv(file.path(in_dir,
          sprintf("scenario%d_Y_m%d_H%d_T%d_rep%02d.csv",      scen, m, H, TT, rep))))
        Delta_true  <- as.matrix(read.csv(file.path(in_dir,
          sprintf("scenario%d_Delta_m%d_H%d_T%d_rep%02d.csv",  scen, m, H, TT, rep))))
        Lambda_true <- as.matrix(read.csv(file.path(in_dir,
          sprintf("scenario%d_Lambda_m%d_H%d_T%d_rep%02d.csv", scen, m, H, TT, rep))))
        prior_fixed$alpha <- prior_fixed$alpha / H
        # -- fit --
        mcmc_cfg <- list(niter = 10000, nburn = 2000, thin = 1)
#fattori a zero, sequential scan       
fit1 <- tryCatch(
          gltfa(
            y       = Y,
            mcmc    = mcmc_cfg,
            Hmax    = Hmax_run,
            prior   = prior_fixed,
            init    = list(H = ncol(Delta_true), 
                      Delta=Delta_true, 
                      ell=get_pivots(Delta_true),
                      Lambda = Lambda_true,
                      Eta = matrix(0, nrow = H, ncol = TT),
                      sigma2 = rep(1, m)
                    ),
            control = list(store_draws = TRUE, store_eta = FALSE,
                           print_every = 1000L,
                           random_scan = FALSE,
                           seed        = scen * 1000L + dimens * 100L +
                                         samplesize * 10L + rep),
            fixed = list(H = TRUE, pivots = TRUE)
          ),
          error = function(e) { message("ERROR: ", e$message); NULL }
        )

#fattori a zero, random scan       
fit2 <- tryCatch(
          gltfa(
            y       = Y,
            mcmc    = mcmc_cfg,
            Hmax    = Hmax_run,
            prior   = prior_fixed,
            init    = list(H = ncol(Delta_true), 
                      Delta=Delta_true, 
                      ell=get_pivots(Delta_true),
                      Lambda = Lambda_true,
                      Eta = matrix(0, nrow = H, ncol = TT),
                      sigma2 = rep(1, m)
                    ),
            control = list(store_draws = TRUE, store_eta = FALSE,
                           print_every = 1000L,
                           random_scan = TRUE,
                           seed        = scen * 1000L + dimens * 100L +
                                         samplesize * 10L + rep),
            fixed = list(H = TRUE, pivots = TRUE)
          ),
          error = function(e) { message("ERROR: ", e$message); NULL }
        )


#fattori random from full condit, sequential scan       
fit3 <- tryCatch(
          gltfa(
            y       = Y,
            mcmc    = mcmc_cfg,
            Hmax    = Hmax_run,
            prior   = prior_fixed,
            init    = list(H = ncol(Delta_true), 
                      Delta=Delta_true, 
                      ell=get_pivots(Delta_true),
                      Lambda = Lambda_true,
                      sigma2 = rep(1, m)
                    ),
            control = list(store_draws = TRUE, store_eta = FALSE,
                           print_every = 1000L,
                           random_scan = FALSE,
                           seed        = scen * 1000L + dimens * 100L +
                                         samplesize * 10L + rep),
            fixed = list(H = TRUE, pivots = TRUE)
          ),
          error = function(e) { message("ERROR: ", e$message); NULL }
        )

#fattori random from full condit, random scan       
fit4 <- tryCatch(
          gltfa(
            y       = Y,
            mcmc    = mcmc_cfg,
            Hmax    = Hmax_run,
            prior   = prior_fixed,
            init    = list(H = ncol(Delta_true), 
                      Delta=Delta_true, 
                      ell=get_pivots(Delta_true),
                      Lambda = Lambda_true,
                      sigma2 = rep(1, m)
                    ),
            control = list(store_draws = TRUE, store_eta = FALSE,
                           print_every = 1000L,
                           random_scan = TRUE,
                           seed        = scen * 1000L + dimens * 100L +
                                         samplesize * 10L + rep),
            fixed = list(H = TRUE, pivots = TRUE)
          ),
          error = function(e) { message("ERROR: ", e$message); NULL }
        )


save.image(file = "/Users/antonio/Dropbox/check22Maggio_a.RData")


stein1 <- stein_mcmc(fit1, Lambda_true, Delta_true)
stein2 <- stein_mcmc(fit2, Lambda_true, Delta_true)
stein3 <- stein_mcmc(fit3, Lambda_true, Delta_true) 
stein4 <- stein_mcmc(fit4, Lambda_true, Delta_true)

plot(stein1, type = "l", col = "red", ylim=c(0, max(c(stein1, stein2, stein3, stein4))), xlab = "MCMC iteration", ylab = "Stein loss")
lines(stein2, col = "blue")
lines(stein3, col = "green")
lines(stein4, col = "orange")
legend("topright", legend = c("Eta 0 - Sequential", "Eta 0 - Random", "Full Cond - Sequential", "Full Cond - Random"), 
col = c("red", "blue", "green", "orange"), lty = 1)

save.image(file = "/Users/antonio/Dropbox/check22Maggio_a.RData")

