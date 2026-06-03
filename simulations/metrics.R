# ------------------------------------------------------------------ #
# Metrics  (Rcpp-accelerated)                                          #
# ------------------------------------------------------------------ #
# Compile once per session with:
#   Rcpp::sourceCpp("simulations/metrics_cpp.cpp")
# ------------------------------------------------------------------ #

compute_metrics <- function(fit, Lambda_true, Delta_true,
                             Sigma_true = rep(1, nrow(Lambda_true))) {

  m   <- nrow(Lambda_true)
  H0  <- ncol(Lambda_true)
  Omega0 <- Lambda_true %*% t(Lambda_true) + diag(Sigma_true)

  Delta_draws <- fit$draws$Delta   # list of nsave matrices (m x H_s)

  # -- spurious count (C++) --
  spurious   <- compute_spurious_cpp(Delta_draws)

  # -- admissibility filter (R; calls user-supplied counting_rule_holds) --
  admissible <- rep(TRUE, fit$meta$nsave)
  for (ite in seq_len(fit$meta$nsave)) {
    Del <- Delta_draws[[ite]]
    if (!spurious[ite]) {
      D_sub <- if (any(rowSums(Del) == 0))
        Del[rowSums(Del) > 0, , drop = FALSE] else Del
      admissible[ite] <- counting_rule_holds(D_sub)
    } else {
      D_ns  <- Del[, colSums(Del) > 1, drop = FALSE]
      D_sub <- if (any(rowSums(D_ns) == 0))
        D_ns[rowSums(D_ns) > 0, , drop = FALSE] else D_ns
      admissible[ite] <- counting_rule_holds(D_sub)
    }
  }

  # nsave is now different from fit$meta$nsave, as we discard some draws
  nsave        <- sum(admissible)
  Delta_draws  <- Delta_draws[admissible]
  Lambda_draws <- fit$draws$Lambda[admissible]
  sigma2_draws <- fit$draws$sigma2[admissible, , drop = FALSE]
  H_draws      <- as.integer((fit$draws$H - spurious)[admissible])

  # -- inner C++ loop: Stein losses + Delta accumulation --
  inner <- compute_metrics_inner_cpp(
    Lambda_draws  = Lambda_draws,
    Delta_draws   = Delta_draws,
    sigma2_draws  = sigma2_draws,
    H_draws       = H_draws,
    Delta_true    = Delta_true,
    Sigma_true    = Sigma_true,
    Omega0_r      = Omega0
  )

  avg_stein          <- inner$stein_sum         / nsave
  avg_stein_sigma    <- inner$stein_sum_sigma   / nsave
  Delta_post         <- inner$Delta_acc         / nsave
  Delta_post_aligned <- inner$Delta_acc_aligned / nsave

  # -- AUC via Wilcoxon rank-sum --
  labs          <- as.vector(Delta_true)
  scors         <- as.vector(Delta_post)
  scors_aligned <- as.vector(Delta_post_aligned)
  n1 <- sum(labs == 1L);  n0 <- sum(labs == 0L)

  auc <- tryCatch(
    wilcox.test(scors[labs == 1L], scors[labs == 0L],
                exact = FALSE)$statistic / (n1 * n0),
    error = function(e) NA_real_
  )
  auc_aligned <- tryCatch(
    wilcox.test(scors_aligned[labs == 1L], scors_aligned[labs == 0L],
                exact = FALSE)$statistic / (n1 * n0),
    error = function(e) NA_real_
  )
  auc         <- ifelse(is.null(auc),         NA_real_, as.numeric(auc))
  auc_aligned <- ifelse(is.null(auc_aligned), NA_real_, as.numeric(auc_aligned))
  

  list(
    H_median          = median(H_draws),
    H_mode            = as.integer(names(which.max(table(H_draws)))),
    H_true_post_prob  = mean(H_draws == H0),
    stein_loss        = avg_stein,
    stein_loss_sigma  = avg_stein_sigma,
    auc               = auc,
    auc_aligned       = auc_aligned,
    Delta_post_aligned = Delta_post_aligned,
    Delta_post        = Delta_post,
    nsave             = nsave,
    spurious          = spurious,
    admissible        = admissible
  )
}