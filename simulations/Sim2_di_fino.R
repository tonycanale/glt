############################################################
#
# Run gltfa fixing different things, to understand where the errors are coming from
#
############################################################

source("run_simulations_1234.R")

scen <- 2
dimens <- 3

m <- m_dim[dimens]
H <- H_dim[scen, dimens]
Hmax_run <- min(2L * H, m - 1L)

for(samplesize in 1:2) {

TT <- T_dim[dimens, samplesize]

      for (rep in 1:5) {

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
        fit <- tryCatch(
          gltfa(
            y       = Y,
            mcmc    = mcmc_cfg,
            Hmax    = Hmax_run,
            prior   = prior_fixed,
            init    = list(H = ncol(Delta_true), 
                      Delta=Delta_true, 
                      ell=get_pivots(Delta_true),
                      Lambda = Lambda_true,
                      Eta = matrix(0, nrow = H, ncol = TT)),
            control = list(store_draws = TRUE, store_eta = FALSE,
                           print_every = 1000L,
                           seed        = scen * 1000L + dimens * 100L +
                                         samplesize * 10L + rep),
            fixed = list(H = TRUE, pivots = TRUE)
          ),
          error = function(e) { message("ERROR: ", e$message); NULL }
        )

if (is.null(fit)) {
          results[idx, ] <- c(scen, m, H, TT, rep, NA, NA, NA, NA, NA, NA, NA, 0)
        } else {
          met            <- compute_metrics(fit, Lambda_true, Delta_true)
            if( length(met$H_mode)==0) met$H_mode <- NA

          results[idx, ] <- c(scen, m, H, TT, rep,
                               met$H_median, 
                               met$H_mode, 
                               met$H_true_post_prob,
                               met$stein_loss, 
                               met$stein_loss_sigma, 
                               met$auc, 
                               met$auc_aligned, 
                               1)
}
print(idx)
        idx <- idx + 1L
    }
  }  
results

save.image(results, file = file.path("/Users/antonio/Dropbox/che21Maggio.RData"))
#results_newiniteta <- as.data.frame(results[1:10,])
#results_deltafixed <- results
#results_lambdafixed <- results


res_stloss <- ggplot(results, aes(x = factor(m), y = stein_loss, fill=factor(T))) +
  geom_boxplot() +
  theme_bw() +
  labs(x ="m", fill = "T", y="") +
  ggtitle("Stein loss full")
  res_stloss

res_stloss_sigma <- ggplot(results, aes(x = factor(m), y = stein_loss_sigma, fill=factor(T))) +
  geom_boxplot() +
  theme_bw() +
  labs(x ="m", fill = "T", y="") +
  ggtitle("Stein loss sigma")
  res_stloss_sigma

