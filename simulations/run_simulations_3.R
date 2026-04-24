############################################################
#
# Run gltfa on each scenario and save results
#
############################################################

devtools::load_all("/Users/antonio/github/glt/", recompile = FALSE)

# ------------------------------------------------------------------ #
# Configuration (must match sim_scenarios.R)                          #
# ------------------------------------------------------------------ #

in_dir  <- "../simulated_data"
out_dir <- "simulation_results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

n_reps <- 20
m_dim  <- c(20, 50, 100)
T_dim  <- cbind(m_dim, 2 * m_dim)
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
  Lambda_draws <- fit$draws$Lambda   # list of nsave matrices (m x H_s)
  Delta_draws  <- fit$draws$Delta    # list of nsave matrices (m x H_s)
  sigma2_draws <- fit$draws$sigma2   # nsave x m matrix
  H_draws      <- fit$draws$H        # integer vector length nsave
  nsave        <- fit$meta$nsave

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

  for (s in seq_len(nsave)) {
    H_s    <- H_draws[s]
    Lam_s  <- Lambda_draws[[s]][, seq_len(H_s), drop = FALSE]
    sig2_s <- sigma2_draws[s, ]
    Omega_s <- Lam_s %*% t(Lam_s) + diag(sig2_s)

    stein_sum <- stein_sum + stein_loss(Omega_s)

    Del_s     <- Delta_draws[[s]]
    Delta_acc <- Delta_acc + align_cols(Del_s, H_s)
  }

  avg_stein  <- stein_sum / nsave
  Delta_post <- Delta_acc / nsave

  # -- AUC via Wilcoxon rank-sum --
  labs  <- as.vector(Delta_true)
  scors <- as.vector(Delta_post)
  n1    <- sum(labs == 1L)
  n0    <- sum(labs == 0L)
  auc   <- as.numeric(
    wilcox.test(scors[labs == 1L], scors[labs == 0L],
                exact = FALSE)$statistic / (n1 * n0)
  )

  list(H_median = median(H_draws), stein_loss = avg_stein, auc = auc)
}

# ------------------------------------------------------------------ #
# Main loop                                                            #
# ------------------------------------------------------------------ #

n_total <-  3 * 2 * n_reps   # scenarios x dims x sample sizes x reps
results <- matrix(NA_real_, nrow = n_total,
                  ncol = 9,
                  dimnames = list(NULL,
                    c("scenario", "m", "H", "T", "rep",
                      "H_median", "stein_loss", "auc", "status")))
idx <- 1L

scen <- 3
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
            init    = list(H = H),
            control = list(store_draws = TRUE, store_eta = FALSE,
                           print_every = 1000L,
                           seed        = scen * 1000L + dimens * 100L +
                                         samplesize * 10L + rep)
          ),
          error = function(e) { message("ERROR: ", e$message); NULL }
        )

        if (is.null(fit)) {
          results[idx, ] <- c(scen, m, H, TT, rep, NA, NA, NA, 0)
        } else {
          met            <- compute_metrics(fit, Lambda_true, Delta_true)
          results[idx, ] <- c(scen, m, H, TT, rep,
                               met$H_median, met$stein_loss, met$auc, 1)
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
  