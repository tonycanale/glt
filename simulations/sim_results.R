############################################################
#
# Analize results of simulations
#
############################################################


res_dir <- "simulation_results_3"


n_reps <- 10
m_dim  <- c(20, 50, 100)
T_dim  <- cbind(rep(100,3), rep(200,3))#cbind(m_dim, 2 * m_dim)
H_dim  <- matrix(
  c( 3,  8, 10,
     5, 10, 15,
     5, 10, 15,
     8, 20, 40),
  nrow = 4, ncol = 3, byrow = TRUE
)

out <- NULL
for(scen in 1:4) {
  if(scen<3) names <- c("scenario", "m", "H", "T", "rep", 
          "H_median", "H_mode", "H_true_post_prob", "stein_loss", "auc", "status")
  if(scen>=3) names <- c("scenario", "m", "H", "T", "rep", 
          "H_median", "H_mode", "H_true_post_prob", "stein_loss", "auc_aligned", "status")
for (dimens in 1:3) {
    m <- m_dim[dimens]
    H <- H_dim[scen, dimens]
    Hmax_run <- min(2L * H, m - 1L)

    for (samplesize in 1:2) {
      TT <- T_dim[dimens, samplesize]

        out_temp <- tryCatch(
          read.csv(file.path(res_dir,
          sprintf("scenario%d_m%d_H%d_T%d_metrics.csv", scen, m, H, TT)))[1:n_reps,names ],
          error = function(e) NULL
        )
      
        if (!is.null(out_temp)) {
          names(out_temp) <- c("scenario", "m", "H", "T", "rep", 
                              "H_median", "H_mode", "H_true_post_prob", 
                              "stein_loss", "auc", "status")
          out <- rbind(out, out_temp)
        }

    } # samplesize
  } # dimens
}

out$scenario <- factor(out$scenario, levels = 1:4,
                       labels = paste("Scenario", 1:4))
out$dimension <- factor(paste("(", out$m, ",", out$T, ")", sep = ""), levels = c("(20,20)", "(20,40)", "(50,50)", "(50,100)", "(100,100)", "(100,200)"))

library(ggplot2)
res_stloss <- ggplot(out, aes(x = factor(m), y = stein_loss, fill=interaction(factor(T),factor(H), sep=", "))) +
  geom_boxplot() +
  facet_wrap(~ scenario) +
  theme_bw() +
  labs(x ="m", fill = "T, H", y="", col="H") +
  ggtitle("Stein loss")

res_auc <-ggplot(out, aes(x = factor(m), y = auc , fill=interaction(factor(T),factor(H), sep=", "))) +
  geom_boxplot() +
  facet_wrap(~ scenario) +
  theme_bw() +
    labs(x ="m", fill = "T, H", y="", col="H") +
  ggtitle("Area Under the ROC Curve (AUC)") +
   geom_hline(yintercept = 0.5, linetype = "dashed", color = "blue")

res_diff_H <- ggplot(out, aes(x = factor(m), y = H-H_median, fill=interaction(factor(T),factor(H), sep=", "))) +
  geom_boxplot() +
  facet_wrap(~ scenario) +
  theme_bw() +
  labs(x ="m", fill = "T, H", y="", col="H") +      
  ggtitle("Difference H-H_median") + geom_hline(yintercept = 0, linetype = "dashed", color = "red")

# export to pdf the three plots
#pdf(file = "simulation_results.pdf", width = 12, height = 8)
res_stloss
res_auc
res_diff_H
dev.off()