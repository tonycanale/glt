############################################################
#
# Analize results of simulations
#
############################################################


res_dir <- "simulations/simulation_results"


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

out <- NULL
for(scen in 1:4) {
for (dimens in 1:3) {
    m <- m_dim[dimens]
    H <- H_dim[scen, dimens]
    Hmax_run <- min(2L * H, m - 1L)

    for (samplesize in 1:2) {
      TT <- T_dim[dimens, samplesize]

        out_temp <- tryCatch(
          read.csv(file.path(res_dir,
          sprintf("scenario%d_m%d_H%d_T%d_metrics.csv", scen, m, H, TT)))[1:n_reps, ],
          error = function(e) NULL
        )
        if (!is.null(out_temp)) {
          out <- rbind(out, out_temp)
        }

    } # samplesize
  } # dimens
}

out$scenario <- factor(out$scenario, levels = 1:4,
                       labels = paste("Scenario", 1:4))
out$dimension <- factor(paste("(", out$m, ",", out$T, ")", sep = ""), levels = c("(20,20)", "(20,40)", "(50,50)", "(50,100)", "(100,100)", "(100,200)"))

library(ggplot2)
ggplot(out, aes(x = factor(dimension), y = stein_loss, fill = factor(H))) +
  geom_boxplot() +
  facet_wrap(~ scenario) +
  theme_bw() +
  labs(x = "(m,T)", fill = "H", y="") +
  ggtitle("Stein loss")

ggplot(out, aes(x = factor(dimension), y = auc , fill = factor(H))) +
  geom_boxplot() +
  facet_wrap(~ scenario) +
  theme_bw() +
    labs(x = "(m,T)", fill = "H", y="") +
  ggtitle("Area Under the ROC Curve (AUC)") +
   geom_hline(yintercept = 0.5, linetype = "dashed", color = "blue")

ggplot(out, aes(x = factor(dimension), y = H-H_median, fill = factor(H))) +
  geom_boxplot() +
  facet_wrap(~ scenario) +
  theme_bw() +
    labs(x = "(m,T)", fill = "H", y="") +
  ggtitle("Difference H-H_median") + geom_hline(yintercept = 0, linetype = "dashed", color = "red")