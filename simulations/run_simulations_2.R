############################################################
#
# Run gltfa on each scenario and save results
#
############################################################

devtools::load_all("/Users/antonio/github/glt/", recompile = FALSE)
source("run_simulations_1234.R")

scen <- 2
  for (dimens in 1:3) {

    m <- m_dim[dimens]
    H <- H_dim[scen, dimens]
    Hmax_run <- min(2L * H, m - 1L)

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
            init    = list(H = H, Delta=Delta_true),
            control = list(store_draws = TRUE, store_eta = FALSE,
                           print_every = 1000L,
                           seed        = scen * 1000L + dimens * 100L +
                                         samplesize * 10L + rep),
            fixed = list(H = TRUE, pivots = TRUE)
          ),
          error = function(e) { message("ERROR: ", e$message); NULL }
        )

if (is.null(fit)) {
          results[idx, ] <- c(scen, m, H, TT, rep, NA, NA, NA, NA, NA, 0)
        } else {
          met            <- compute_metrics(fit, Lambda_true, Delta_true)
          results[idx, ] <- c(scen, m, H, TT, rep,
                               met$H_median, 
                               met$H_mode, 
                               met$H_true_post_prob,
                               met$stein_loss, met$auc, 1)
        }

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
