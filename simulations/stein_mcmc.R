
stein_mcmc <- function(fit, Lambda_true, Delta_true,
                             Sigma_true = rep(1, nrow(Lambda_true))) {

  m   <- nrow(Lambda_true)
  H0  <- ncol(Lambda_true)
  Omega0       <-   Lambda_true %*% t(Lambda_true) + diag(Sigma_true)
  
  Delta_draws  <- fit$draws$Delta    # list of nsave matrices (m x H_s)
  
  # Counting rule: discard draws where any row of Delta is zero, or where the counting rule does not hold
  admissible <- rep(TRUE, fit$meta$nsave)
  for(ite in 1:fit$meta$nsave) {
    if(any(rowSums(Delta_draws[[ite]]) == 0)) {
      admissible[ite] <- counting_rule_holds(Delta_draws[[ite]][rowSums(Delta_draws[[ite]]) > 0, , drop = FALSE])
    } else {
      admissible[ite] <- counting_rule_holds(Delta_draws[[ite]])
    }
  }
  # nsave is now different from fit$meta$nsave, as we discard some draws
  nsave        <- sum(admissible)

  Delta_draws <- Delta_draws[admissible]    # list of nsave matrices (m x H_s)
  Lambda_draws <- fit$draws$Lambda[admissible]   # list of nsave matrices (m x H_s)
  sigma2_draws <- fit$draws$sigma2[admissible, , drop = FALSE]    # nsave x m matrix
  H_draws      <- fit$draws$H[admissible]        # integer vector length nsave


  # -- Stein loss: tr(A B^{-1}) - log|det(A B^{-1})| - m --
  stein_loss <- function(Omega_hat, B=Omega0) {
    ev <- Re(eigen(Omega_hat %*% solve(B), only.values = TRUE)$values)
    sum(ev) - sum(log(pmax(ev, .Machine$double.eps))) - nrow(Omega_hat)
  }

  stein_all <- rep(NA, nsave)
  
  for (s in seq_len(nsave)) {
    H_s    <- H_draws[s]
    Lam_s  <- Lambda_draws[[s]][, seq_len(H_s), drop = FALSE]
    sig2_s <- sigma2_draws[s, ]
    Omega_s <- Lam_s %*% t(Lam_s) + diag(sig2_s)
    stein_all[s] <-  stein_loss(Omega_s, B=Omega0)
 }

  stein_all
}
