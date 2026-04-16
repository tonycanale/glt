# language: R
# Diagnostics for H and variance estimates (run with your fit in workspace)

diag_glt <- function(fit, Lambda_true = NULL, sigma2_true = NULL, true_H = NULL, burn = 0) {
  draws <- fit$draws
  Hvec <- draws$H
  if (burn > 0) {
    keep <- (seq_along(Hvec) > burn)
    Hvec <- Hvec[keep]
  }
  res <- list()

  # basic H summaries
  res$H_table <- table(Hvec)
  res$prop_trueH <- if (!is.null(true_H)) mean(Hvec == true_H) else NA_real_
  res$H_mean <- mean(Hvec)
  res$H_median <- median(Hvec)

  # run-lengths (stickiness)
  rleH <- rle(Hvec)
  res$runlength_summary <- summary(rleH$lengths)
  res$mean_runlength <- mean(rleH$lengths)

  # change statistics (inferred accepted moves)
  diffs <- diff(Hvec)
  res$n_changes <- sum(diffs != 0)
  res$n_increases <- sum(diffs > 0)
  res$n_decreases <- sum(diffs < 0)
  res$frac_increase <- res$n_increases / (res$n_increases + res$n_decreases)

  # boundaries
  res$prop_at_1 <- mean(Hvec == 1)
  m <- nrow(draws$Lambda[[1]])
  res$prop_at_m <- mean(Hvec == m)

  # acceptance matrix if present
  if (!is.null(fit$accept)) {
    acc <- fit$accept
    if (is.matrix(acc) || is.data.frame(acc)) res$accept_colmeans <- colMeans(acc)
    else res$accept_colmeans <- NA
  } else res$accept_colmeans <- NA

  # ESS / autocorr for H
  if (requireNamespace("coda", quietly = TRUE)) {
    res$ess_H <- as.numeric(coda::effectiveSize(coda::mcmc(Hvec)))
    res$acf_H <- stats::acf(Hvec, plot = FALSE)$acf
  } else res$ess_H <- NA

  # inferred LL' and total covariances if Lambda_true provided
  if (!is.null(Lambda_true) && !is.null(draws$Lambda)) {
    # compute posterior mean of L L' and L L' + Sigma
    nsave <- length(draws$Lambda)
    m <- nrow(draws$Lambda[[1]])
    LLt_sum <- matrix(0, m, m)
    Total_sum <- matrix(0, m, m)
    count <- 0L
    for (s in seq_len(nsave)) {
      Ls <- draws$Lambda[[s]]
      ss <- draws$sigma2[s, ]
      LLt_sum <- LLt_sum + tcrossprod(Ls)
      Total_sum <- Total_sum + (tcrossprod(Ls) + diag(as.vector(ss), m))
    }
    res$LLt_pm <- LLt_sum / nsave
    res$Total_pm <- Total_sum / nsave
    if (!is.null(Lambda_true) && !is.null(sigma2_true)) {
      trueTotal <- tcrossprod(Lambda_true) + diag(as.vector(sigma2_true), nrow(Lambda_true))
      res$LLt_frob_err <- sqrt(sum((res$LLt_pm - tcrossprod(Lambda_true))^2))
      res$Total_frob_err <- sqrt(sum((res$Total_pm - trueTotal)^2))
    }
  }

  return(res)
}

# Example run (adjust burn if needed)
d <- diag_glt(fit, Lambda_true = Lambda_true, sigma2_true = sigma2_true, true_H = 3, burn = 1000L)
print(d$H_table)
cat("Prop equal to true H:", d$prop_trueH, "\n")
cat("Mean run-length:", d$mean_runlength, "\n")
cat("n increases / decreases:", d$n_increases, d$n_decreases, "\n")
print(d$runlength_summary)
print(d$accept_colmeans)
cat("ESS(H):", d$ess_H, "\n")
if (!is.null(d$LLt_frob_err)) {
  cat("Frobenius error LLt:", d$LLt_frob_err, " Total:", d$Total_frob_err, "\n")
}

# Quick plots
par(mfrow = c(2,2))
plot(fit$draws$H, type = "l", main = "Trace: H"); abline(h = 3, col = "red", lty = 2)
hist(fit$draws$H, main = "Histogram H", xlab = "H")
if (exists("d$LLt_pm")) {
  image(t(d$LLt_pm[nrow(d$LLt_pm):1, ]), main = "Posterior mean LLt")
  image(t(tcrossprod(Lambda_true)[nrow(Lambda_true):1, ]), main = "True LLt")
}