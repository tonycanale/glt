############################################################
#
# Run gltfa on each scenario and save results
#
############################################################

source("run_simulations_1234.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Devi fornire scen da terminale, ad esempio: --args 1")
}

scen <- as.integer(args[1])

if (is.na(scen)) {
  stop("scen deve essere un numero intero")
}

cat("Scenario:", scen, "\n")

for (dimens in 3:1) {
    m <- m_dim[dimens]
    H <- H_dim[scen, dimens]
    Hmax_run <- m/2

    for (samplesize in 1:2) {
      TT <- T_dim[dimens, samplesize]

      for (rep in seq_len(n_reps)) {

        tag <- sprintf("scenario%d_m%d_H%d_T%d_rep%02d", scen, m, H, TT, rep)
        cat("\n[", format(Sys.time(), "%H:%M:%S"), "]  Fitting", tag, "\n")

        # -- load data --
        Y           <- as.matrix(read.csv(file.path(in_dir,
          sprintf("scenario%d_Y_m%d_H%d_T%d_rep%02d.csv",      scen, m, H, TT, rep))))
        Delta_true  <- as.matrix(read.csv(file.path(in_dir,
          sprintf("scenario%d_Delta_m%d_H%d_T%d_rep%02d.csv",  scen, m, H, TT, rep))))
        Lambda_true <- as.matrix(read.csv(file.path(in_dir,
          sprintf("scenario%d_Lambda_m%d_H%d_T%d_rep%02d.csv", scen, m, H, TT, rep))))

        # -- fit --
        fit <- tryCatch(
           gltfa(
            y       = Y,
            mcmc    = mcmc_cfg,
            Hmax    = Hmax_run,
            prior   = prior_fixed,
            init    = list(H = H, 
              Delta=Delta_true, 
              Lambda=Lambda_true,
              ell=get_pivots(Delta_true),
              sigma2 = rep(1, m)),
            control = list(store_draws = TRUE, store_eta = FALSE,
                           print_every = 1000L,
                           seed        = scen * 1000L + dimens * 100L +
                                         samplesize * 10L + rep)#,
            #fixed = list(H = TRUE, pivots = TRUE)
          ),
          error = function(e) { message("ERROR: ", e$message); NULL }
        )


if (is.null(fit)) {
          results[idx, ] <- c(scen, m, H, TT, rep, 
            NA, NA,NA, NA, NA, NA, NA, NA, 0)
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
                               met$nsave, 
                               1 )
}
          saveRDS(fit, file = file.path(out_dir, sprintf("fit_scen%d_m%d_H%d_T%d_rep%02d.rds", scen, m, H, TT, rep)))


        idx <- idx + 1L
      } # rep

      # -- save block after each (scen, dimens, samplesize) --
      block <- results[results[, "scenario"] == scen &
                       results[, "m"]        == m    &
                       results[, "T"]        == TT, , drop = FALSE]
      write.csv(block,
        file      = file.path(out_dir,
          sprintf("scenario%d_m%d_H%d_T%d_metrics.csv", scen, m, H, TT)),
        row.names = FALSE
      )

    } # samplesize
  } # dimens
