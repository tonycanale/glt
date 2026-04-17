############################################################
# summarize_gltfa_sim.R
#
# Post-processing script for simulation results produced by
# sim_gltfa.R
#
# It computes:
#   1. latent dimension recovery metrics
#   2. covariance matrix error (Omega)
#   3. Delta recovery metrics using ROC / AUC style summaries
#
# IMPORTANT:
# You will likely need to edit only the extractor functions
# in Section 1 so they match the structure of your gltfa() output.
############################################################

rm(list = ls())

## ----------------------------- ##
## 0. User-facing configuration ##
## ----------------------------- ##

results_dir <- "sim_gltfa_output"
summary_dir <- file.path(results_dir, "summaries")
if (!dir.exists(summary_dir)) dir.create(summary_dir, recursive = TRUE)

# Thresholds for turning posterior probabilities into hard 0/1 calls
delta_thresholds <- c(0.1, 0.25, 0.5, 0.75, 0.9)

# Tiny value to avoid division-by-zero in edge cases
eps_num <- 1e-12


## ---------------------------------------------------- ##
## 1. EDIT HERE IF NEEDED: extract pieces from gltfa()  ##
## ---------------------------------------------------- ##

# This function should return the estimated latent dimension.
# Edit according to your fitted object.
extract_k_hat <- function(fit) {
  # Examples of possible structures:
  # if (!is.null(fit$H)) return(fit$H)
  # if (!is.null(fit$k)) return(fit$k)
  # if (!is.null(fit$n_factors)) return(fit$n_factors)

  if (!is.null(fit$H)) return(as.integer(fit$H))
  if (!is.null(fit$k)) return(as.integer(fit$k))
  if (!is.null(fit$n_factors)) return(as.integer(fit$n_factors))

  stop("Could not extract estimated latent dimension from fit. Edit extract_k_hat().")
}

# This function should return an estimated covariance matrix Omega_hat, p x p.
# Edit according to your fitted object.
extract_Omega_hat <- function(fit) {
  # Examples:
  # if (!is.null(fit$Omega)) return(fit$Omega)
  # if (!is.null(fit$Omega_hat)) return(fit$Omega_hat)
  # if (!is.null(fit$Lambda) && !is.null(fit$Sigma)) {
  #   return(fit$Lambda %*% t(fit$Lambda) + diag(fit$Sigma))
  # }

  if (!is.null(fit$Omega)) return(fit$Omega)
  if (!is.null(fit$Omega_hat)) return(fit$Omega_hat)

  if (!is.null(fit$Lambda) && !is.null(fit$Sigma)) {
    Sig <- fit$Sigma
    if (is.vector(Sig)) Sig <- diag(Sig)
    return(fit$Lambda %*% t(fit$Lambda) + Sig)
  }

  stop("Could not extract Omega_hat from fit. Edit extract_Omega_hat().")
}

# This function should return a matrix of posterior scores for Delta, p x k_est_or_ref.
# BEST choice: posterior inclusion probabilities in [0,1].
#
# Examples:
#   - posterior mean of Delta samples
#   - inclusion probability matrix
#   - frequency of nonzero loadings across MCMC draws
#
# If you only have a single hard estimate 0/1 matrix, ROC/AUC is much less informative,
# but you can still return a binary matrix and AUC becomes based on ties.
extract_delta_scores <- function(fit) {
  # Examples:
  # if (!is.null(fit$Delta_prob)) return(fit$Delta_prob)
  # if (!is.null(fit$post_inclusion_prob)) return(fit$post_inclusion_prob)
  # if (!is.null(fit$Delta_mean)) return(fit$Delta_mean)
  # if (!is.null(fit$Delta)) return(fit$Delta)

  if (!is.null(fit$Delta_prob)) return(fit$Delta_prob)
  if (!is.null(fit$post_inclusion_prob)) return(fit$post_inclusion_prob)
  if (!is.null(fit$Delta_mean)) return(fit$Delta_mean)
  if (!is.null(fit$Delta)) return(fit$Delta)

  stop("Could not extract Delta score matrix from fit. Edit extract_delta_scores().")
}


## -------------------------------- ##
## 2. General utility/helper functions ##
## -------------------------------- ##

mse_matrix <- function(A, B) {
  mean((A - B)^2)
}

safe_div <- function(num, den) {
  ifelse(abs(den) < eps_num, NA_real_, num / den)
}

# Pad or crop a matrix to p x k_target with zeros on the right if needed.
align_matrix_ncol <- function(M, p, k_target) {
  stopifnot(nrow(M) == p)

  k_now <- ncol(M)

  if (k_now == k_target) return(M)

  if (k_now > k_target) {
    return(M[, seq_len(k_target), drop = FALSE])
  }

  out <- matrix(0, nrow = p, ncol = k_target)
  out[, seq_len(k_now)] <- M
  out
}

# True Delta from true Lambda
lambda_to_delta <- function(Lambda) {
  1L * (Lambda != 0)
}

# Hard-call metrics at a given threshold
delta_binary_metrics <- function(delta_true, delta_scores, threshold = 0.5) {
  pred <- 1L * (delta_scores >= threshold)

  tp <- sum(pred == 1 & delta_true == 1)
  fp <- sum(pred == 1 & delta_true == 0)
  tn <- sum(pred == 0 & delta_true == 0)
  fn <- sum(pred == 0 & delta_true == 1)

  tpr <- safe_div(tp, tp + fn)   # sensitivity / recall
  fpr <- safe_div(fp, fp + tn)
  tnr <- safe_div(tn, tn + fp)   # specificity
  ppv <- safe_div(tp, tp + fp)   # precision
  npv <- safe_div(tn, tn + fn)
  acc <- safe_div(tp + tn, tp + tn + fp + fn)
  f1  <- safe_div(2 * tp, 2 * tp + fp + fn)

  data.frame(
    threshold = threshold,
    tp = tp, fp = fp, tn = tn, fn = fn,
    tpr = tpr, fpr = fpr, tnr = tnr,
    ppv = ppv, npv = npv, acc = acc, f1 = f1
  )
}

# ROC curve by threshold sweep
compute_roc_curve <- function(delta_true, delta_scores) {
  truth <- as.integer(as.vector(delta_true))
  score <- as.numeric(as.vector(delta_scores))

  ord <- order(score, decreasing = TRUE)
  truth <- truth[ord]
  score <- score[ord]

  P <- sum(truth == 1)
  N <- sum(truth == 0)

  # Edge case: cannot define ROC if all 0 or all 1
  if (P == 0 || N == 0) {
    return(list(
      roc = data.frame(
        threshold = c(Inf, -Inf),
        tpr = c(0, 1),
        fpr = c(0, 1)
      ),
      auc = NA_real_
    ))
  }

  # Unique thresholds, plus sentinels
  thresholds <- sort(unique(score), decreasing = TRUE)
  thresholds <- c(Inf, thresholds, -Inf)

  roc_list <- lapply(thresholds, function(th) {
    pred <- as.integer(score >= th)
    tp <- sum(pred == 1 & truth == 1)
    fp <- sum(pred == 1 & truth == 0)
    fn <- sum(pred == 0 & truth == 1)
    tn <- sum(pred == 0 & truth == 0)

    tpr <- tp / (tp + fn)
    fpr <- fp / (fp + tn)

    data.frame(threshold = th, tpr = tpr, fpr = fpr)
  })

  roc_df <- do.call(rbind, roc_list)

  # Remove duplicated (fpr,tpr) rows caused by ties
  roc_df <- roc_df[!duplicated(roc_df[, c("fpr", "tpr")]), , drop = FALSE]
  roc_df <- roc_df[order(roc_df$fpr, roc_df$tpr), , drop = FALSE]

  # Trapezoidal AUC
  auc <- 0
  if (nrow(roc_df) >= 2) {
    for (i in 2:nrow(roc_df)) {
      x1 <- roc_df$fpr[i - 1]
      x2 <- roc_df$fpr[i]
      y1 <- roc_df$tpr[i - 1]
      y2 <- roc_df$tpr[i]
      auc <- auc + (x2 - x1) * (y1 + y2) / 2
    }
  }

  list(roc = roc_df, auc = auc)
}

# Optional precision-recall style summary if useful later
compute_pr_curve <- function(delta_true, delta_scores) {
  truth <- as.integer(as.vector(delta_true))
  score <- as.numeric(as.vector(delta_scores))

  thresholds <- sort(unique(score), decreasing = TRUE)
  thresholds <- c(Inf, thresholds, -Inf)

  out <- lapply(thresholds, function(th) {
    pred <- as.integer(score >= th)

    tp <- sum(pred == 1 & truth == 1)
    fp <- sum(pred == 1 & truth == 0)
    fn <- sum(pred == 0 & truth == 1)

    precision <- safe_div(tp, tp + fp)
    recall <- safe_div(tp, tp + fn)

    data.frame(threshold = th, precision = precision, recall = recall)
  })

  do.call(rbind, out)
}


## ------------------------------------------------ ##
## 3. One-run summary: k, Omega, Delta ROC/AUC etc ##
## ------------------------------------------------ ##

summarize_one_result <- function(res, run_id = NA_integer_) {
  fit <- res$fit
  truth <- res$truth

  p <- truth$p
  k_true <- truth$k
  scenario <- truth$scenario
  n <- truth$n
  rep_id <- truth$rep

  # Truth
  Lambda_true <- truth$Lambda
  Omega_true <- truth$Omega
  Delta_true <- lambda_to_delta(Lambda_true)

  # Estimates
  k_hat <- extract_k_hat(fit)
  Omega_hat <- extract_Omega_hat(fit)
  Delta_scores_raw <- extract_delta_scores(fit)

  # Align score matrix to truth dimension p x k_true
  # This is the simplest baseline choice.
  # If later you implement matching/permutation alignment, replace here.
  if (nrow(Delta_scores_raw) != p) {
    stop("Delta score matrix has wrong number of rows.")
  }
  Delta_scores <- align_matrix_ncol(Delta_scores_raw, p = p, k_target = k_true)

  # If scores are not in [0,1], try converting hard 0/1 style
  if (min(Delta_scores) < 0 || max(Delta_scores) > 1) {
    warning("Delta scores not in [0,1]. Converting to binary via != 0.")
    Delta_scores <- 1 * (Delta_scores != 0)
  }

  # Basic fit summaries
  k_err <- k_hat - k_true
  abs_k_err <- abs(k_err)
  omega_mse <- mse_matrix(Omega_true, Omega_hat)

  # ROC/AUC
  roc_obj <- compute_roc_curve(Delta_true, Delta_scores)
  auc <- roc_obj$auc
  roc_df <- roc_obj$roc
  roc_df$run_id <- run_id
  roc_df$scenario <- scenario
  roc_df$n <- n
  roc_df$p <- p
  roc_df$k_true <- k_true
  roc_df$rep <- rep_id

  # Thresholded binary summaries
  thresh_df <- do.call(
    rbind,
    lapply(delta_thresholds, function(th) {
      out <- delta_binary_metrics(Delta_true, Delta_scores, threshold = th)
      out$run_id <- run_id
      out$scenario <- scenario
      out$n <- n
      out$p <- p
      out$k_true <- k_true
      out$rep <- rep_id
      out
    })
  )

  # Main per-run summary row
  main_df <- data.frame(
    run_id = run_id,
    scenario = scenario,
    n = n,
    p = p,
    k_true = k_true,
    k_hat = k_hat,
    k_err = k_err,
    abs_k_err = abs_k_err,
    omega_mse = omega_mse,
    delta_auc = auc,
    rep = rep_id,
    stringsAsFactors = FALSE
  )

  list(
    main = main_df,
    roc = roc_df,
    thresholds = thresh_df
  )
}


## ------------------------------------------- ##
## 4. Read all run files and summarize them    ##
## ------------------------------------------- ##

rds_files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
rds_files <- rds_files[basename(rds_files) != "all_results.rds"]

if (length(rds_files) == 0) {
  stop("No per-run .rds files found in results_dir.")
}

all_main <- list()
all_roc <- list()
all_thresholds <- list()

for (i in seq_along(rds_files)) {
  cat(sprintf("Summarizing %d / %d: %s\n", i, length(rds_files), basename(rds_files[i])))

  res <- readRDS(rds_files[i])
  sm <- summarize_one_result(res, run_id = i)

  all_main[[i]] <- sm$main
  all_roc[[i]] <- sm$roc
  all_thresholds[[i]] <- sm$thresholds
}

main_df <- do.call(rbind, all_main)
roc_df <- do.call(rbind, all_roc)
threshold_df <- do.call(rbind, all_thresholds)


## ---------------------------------------------------- ##
## 5. Aggregate summaries across reps / configurations ##
## ---------------------------------------------------- ##

agg_main <- aggregate(
  cbind(k_hat, k_err, abs_k_err, omega_mse, delta_auc) ~ scenario + n + p + k_true,
  data = main_df,
  FUN = mean
)

agg_main_sd <- aggregate(
  cbind(k_hat, k_err, abs_k_err, omega_mse, delta_auc) ~ scenario + n + p + k_true,
  data = main_df,
  FUN = sd
)

# threshold-wise averages
agg_threshold <- aggregate(
  cbind(tpr, fpr, tnr, ppv, npv, acc, f1) ~ scenario + n + p + k_true + threshold,
  data = threshold_df,
  FUN = mean
)

# merge means and sds for main summaries
names(agg_main)[5:ncol(agg_main)] <- paste0(names(agg_main)[5:ncol(agg_main)], "_mean")
names(agg_main_sd)[5:ncol(agg_main_sd)] <- paste0(names(agg_main_sd)[5:ncol(agg_main_sd)], "_sd")

summary_main <- merge(
  agg_main,
  agg_main_sd,
  by = c("scenario", "n", "p", "k_true"),
  all = TRUE
)


## ----------------------------- ##
## 6. Save outputs to disk       ##
## ----------------------------- ##

write.csv(main_df, file.path(summary_dir, "per_run_summary.csv"), row.names = FALSE)
write.csv(roc_df, file.path(summary_dir, "per_run_roc_points.csv"), row.names = FALSE)
write.csv(threshold_df, file.path(summary_dir, "per_run_threshold_metrics.csv"), row.names = FALSE)
write.csv(summary_main, file.path(summary_dir, "summary_main_mean_sd.csv"), row.names = FALSE)
write.csv(agg_threshold, file.path(summary_dir, "summary_threshold_mean.csv"), row.names = FALSE)

saveRDS(main_df, file.path(summary_dir, "per_run_summary.rds"))
saveRDS(roc_df, file.path(summary_dir, "per_run_roc_points.rds"))
saveRDS(threshold_df, file.path(summary_dir, "per_run_threshold_metrics.rds"))
saveRDS(summary_main, file.path(summary_dir, "summary_main_mean_sd.rds"))
saveRDS(agg_threshold, file.path(summary_dir, "summary_threshold_mean.rds"))


## ----------------------------- ##
## 7. Optional plotting section  ##
## ----------------------------- ##

# Average ROC by scenario/n/p/k_true:
# This averages TPR over a common FPR grid after interpolation.

interp_roc <- function(df_one, fpr_grid) {
  # df_one must contain fpr and tpr
  o <- order(df_one$fpr, df_one$tpr)
  x <- df_one$fpr[o]
  y <- df_one$tpr[o]

  # remove duplicates in x for approx()
  keep <- !duplicated(x)
  x <- x[keep]
  y <- y[keep]

  # need at least two points
  if (length(x) < 2) {
    return(rep(NA_real_, length(fpr_grid)))
  }

  approx(x = x, y = y, xout = fpr_grid, ties = mean, rule = 2)$y
}

fpr_grid <- seq(0, 1, length.out = 201)

split_keys <- unique(roc_df[, c("scenario", "n", "p", "k_true", "rep", "run_id")])

roc_interp_list <- vector("list", nrow(split_keys))

for (i in seq_len(nrow(split_keys))) {
  key <- split_keys[i, ]
  sub <- roc_df[
    roc_df$run_id == key$run_id &
      roc_df$rep == key$rep &
      roc_df$scenario == key$scenario &
      roc_df$n == key$n &
      roc_df$p == key$p &
      roc_df$k_true == key$k_true,
    ,
    drop = FALSE
  ]

  tpr_interp <- interp_roc(sub, fpr_grid)

  roc_interp_list[[i]] <- data.frame(
    scenario = key$scenario,
    n = key$n,
    p = key$p,
    k_true = key$k_true,
    rep = key$rep,
    run_id = key$run_id,
    fpr = fpr_grid,
    tpr = tpr_interp
  )
}

roc_interp_df <- do.call(rbind, roc_interp_list)

avg_roc_df <- aggregate(
  tpr ~ scenario + n + p + k_true + fpr,
  data = roc_interp_df,
  FUN = function(x) mean(x, na.rm = TRUE)
)

write.csv(avg_roc_df, file.path(summary_dir, "average_roc_curves.csv"), row.names = FALSE)

# Base R plot example: one pdf per scenario
pdf(file.path(summary_dir, "average_roc_curves.pdf"), width = 9, height = 7)

scenario_levels <- unique(avg_roc_df$scenario)
for (sc in scenario_levels) {
  sub_sc <- avg_roc_df[avg_roc_df$scenario == sc, , drop = FALSE]

  plot(
    c(0, 1), c(0, 1),
    type = "n",
    xlab = "False Positive Rate",
    ylab = "True Positive Rate",
    main = paste("Average ROC -", sc)
  )
  abline(0, 1, lty = 2)

  combos <- unique(sub_sc[, c("n", "p", "k_true")])

  for (i in seq_len(nrow(combos))) {
    ss <- sub_sc[
      sub_sc$n == combos$n[i] &
      sub_sc$p == combos$p[i] &
      sub_sc$k_true == combos$k_true[i],
      ,
      drop = FALSE
    ]

    lines(ss$fpr, ss$tpr)
  }

  legend(
    "bottomright",
    legend = apply(combos, 1, function(z) paste0("n=", z[1], ", p=", z[2], ", k=", z[3])),
    lty = 1,
    cex = 0.7,
    bty = "n"
  )
}
dev.off()

cat("\nPost-processing completed.\n")
cat("Summaries saved in:", normalizePath(summary_dir), "\n")