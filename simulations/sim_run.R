############################################################
#
# Run gltfac on each scenario and save results
#
############################################################
compute_metrics <- function(fit, Lambda_true, Delta_true, Sigma_true = rep(1, nrow(Lambda_true))) {

  m      <- nrow(Lambda_true)
  H0     <- ncol(Lambda_true)
  nsave  <- fit$meta$nsave
  Hmax   <- fit$meta$Hmax

  Omega0 <- Lambda_true %*% t(Lambda_true) + diag(Sigma_true)

  # Extract draws
  Lambda_draws <- fit$draws$Lambda   # m x Hmax x nsave  (columns beyond H_s are zeros)
  Delta_draws  <- fit$draws$Delta    # m x Hmax x nsave
  sigma2_draws <- fit$draws$sigma2   # m x nsave
  H_draws      <- fit$draws$H        # nsave

  # ------------------------------------------------------------------
  # Helper: greedily match columns of D_draw (m x H_draw) to Delta_true
  # (m x H0) by maximum Hamming overlap; pad with 0 cols if H_draw < H0,
  # truncate after matching if H_draw > H0.
  # ------------------------------------------------------------------
  align_cols <- function(D_draw, H_draw) {
    D_active <- D_draw[, seq_len(H_draw), drop = FALSE]

    n_match  <- min(H_draw, H0)
    assigned <- integer(n_match)
    used     <- logical(H_draw)

    for (j in seq_len(n_match)) {
      sims       <- sapply(seq_len(H_draw), function(k) {
        if (used[k]) -Inf else sum(D_active[, k] == Delta_true[, j])
      })
      best         <- which.max(sims)
      assigned[j]  <- best
      used[best]   <- TRUE
    }

    D_aligned <- D_active[, assigned, drop = FALSE]

    # Pad with zeros if fewer active columns than H0
    if (H_draw < H0) {
      D_aligned <- cbind(D_aligned, matrix(0L, m, H0 - H_draw))
    }
    D_aligned
  }

  # ------------------------------------------------------------------
  # Stein loss: tr(A B^{-1}) - log|det(A B^{-1})| - m
  # ------------------------------------------------------------------
  stein_loss <- function(Omega_hat) {
    ev <- Re(eigen(Omega_hat %*% solve(Omega0), only.values = TRUE)$values)
    sum(ev) - sum(log(pmax(ev, .Machine$double.eps))) - m
  }

  # ------------------------------------------------------------------
  # Accumulate over MCMC draws
  # ------------------------------------------------------------------
  stein_sum <- 0
  Delta_acc <- matrix(0, m, H0)

  for (s in seq_len(nsave)) {
    H_s     <- H_draws[s]
    Lam_s   <- Lambda_draws[[s]]
    sig2_s  <- sigma2_draws[s,]
    Omega_s <- Lam_s %*% t(Lam_s) + diag(sig2_s)

    stein_sum <- stein_sum + stein_loss(Omega_s)

    Del_s     <- Delta_draws[[s]]
    Delta_acc <- Delta_acc + align_cols(Del_s, H_s)
  }

  avg_stein  <- stein_sum / nsave
  Delta_post <- Delta_acc / nsave   # posterior inclusion probability, column-aligned

  # ------------------------------------------------------------------
  # AUC: scores = posterior inclusion probs; labels = Delta_true entries
  # Computed via Wilcoxon rank-sum (equivalent to empirical AUC)
  # ------------------------------------------------------------------
  labels <- as.vector(Delta_true)
  scores <- as.vector(Delta_post)
  n1     <- sum(labels == 1L)
  n0     <- sum(labels == 0L)
  auc    <- wilcox.test(scores[labels == 1L], scores[labels == 0L],
                        exact = FALSE)$statistic / (n1 * n0)

  list(
    H_median   = median(H_draws),
    stein_loss = avg_stein,
    auc        = as.numeric(auc)
  )
}

