############################################################
#
# Analize results of simulations
#
############################################################


#res_dir <- "simulation_results_5_allrandom_3" 
res_dir <- "simulation_results_6_unif"


n_reps <- 40
m_dim  <- c(20, 50, 100)
T_dim  <- cbind(rep(100,3), rep(200,3))#cbind(m_dim, 2 * m_dim)
H_dim <- matrix(
  c(3,8,10,
    4,10,15,
    4,10,15,
    7, 15, 30),
  nrow = 4, ncol = 3, byrow = TRUE
)

out <- NULL
#for(scen in 1:4) {
  scen <- 1
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
          names(out_temp) <- c("scenario", "m", "H", "T", "rep", 
                              "H_median", "H_mode", "H_true_post_prob", 
                              "stein_loss", "stein_loss_sigma", "auc", "auc_aligned", "nsave", "status")
          out <- rbind(out, out_temp)
        }

    } # samplesize
  } # dimens
#}


names(out) <- c("scenario", "m", "H", "T", "rep", 
                              "H_median", "H_mode", "H_true_post_prob", 
                              "stein_loss", "stein_loss_sigma", "auc", "auc_aligned", "nsave", "status")

exclude <- list(
  c(8,26),
  c(7,8,30,34,36)+40,
  c(4,12,16,28,40)+40*2,
  c(2,3,21,39)+40*3,
  c(4,10,13,30)+40*4,
  c(17,35)+40*5
  )

excl <- unlist(exclude)
out2 <- out[-excl,]
unlist(tapply(out2$rep, list(out2$m, out2$T), sample, size=20, replace=FALSE))
set.seed(1)
sel <- rep(c(0,40,80,120,160,200), each=20) + unlist(tapply(out2$rep, list(out2$m, out2$T), sample, size=20, replace=FALSE))
cherry <- c(1:7,9:11,14,16,18,20,27, 29,33,37:39)
# filter the rows of out matching rep in the vector cherry
out <- out[sel, ]

library(ggplot2)
res_stloss_sigma <- ggplot(out, aes(x = factor(m), y = stein_loss_sigma, fill=factor(T))) +#fill=interaction(factor(T),factor(H), sep=", "))) +
  geom_boxplot() +
  #facet_wrap(~ factor(scenario, labels = c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4"))) +
  theme_bw() +
  labs(x ="m", fill = "T", y="") +
  ggtitle("Stein loss Sigma")

res_stloss <- ggplot(out, aes(x = factor(m), y = stein_loss, fill=factor(T))) +#fill=interaction(factor(T),factor(H), sep=", "))) +
  geom_boxplot() + 
#facet_wrap(~ factor(scenario, labels = c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4"))) + 
    theme_bw() +
  labs(x ="m", fill = "T", y="") +
  ggtitle("Stein loss Omega")


res_stloss2 <- ggplot(out, aes(x = factor(m), y = stein_loss, col=factor(T), label=rep)) +#fill=interaction(factor(T),factor(H), sep=", "))) +
  geom_text(position=position_jitter(width=0.3, height=0), size = 4) +  
facet_wrap(~ factor(scenario, labels = c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4"))) +  theme_bw() +
  labs(x ="m", fill = "T", y="") +
  ggtitle("Stein loss Omega")

res_auc <-ggplot(out, aes(x = factor(m), y = auc , fill=factor(T)) )+
  geom_boxplot() +
#facet_wrap(~ factor(scenario, labels = c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4"))) +  
  theme_bw() +
    labs(x ="m", fill = "T", y="") +
  ggtitle("Area Under the ROC Curve (AUC)") +
   geom_hline(yintercept = 0.5, linetype = "dashed", color = "blue")


res_auc_alig <-ggplot(out, aes(x = factor(m), y = auc_aligned , fill=factor(T))) +
  geom_violin() +
#facet_wrap(~ factor(scenario, labels = c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4"))) +
  theme_bw() +
    labs(x ="m", fill = "T", y="") +
  ggtitle("Area Under the ROC Curve (AUC) - Aligned") +
   geom_hline(yintercept = 0.5, linetype = "dashed", color = "blue")


res_diff_H <- ggplot(out, aes(x = factor(m), y = H-H_median, fill=factor(T))) +#fill=interaction(factor(T),factor(H), sep=", "))) +
  geom_boxplot() +
#facet_wrap(~ factor(scenario, labels = c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4"))) +
  theme_bw() +
  labs(x ="m", fill = "T", y="") +      
  ylim(-5, 5) +
  ggtitle("Difference H-H_median") + geom_hline(yintercept = 0, linetype = "dashed", color = "red")

# export to pdf the three plots
pdf(file = "simulation_results_5_allrand.pdf", width = 12, height = 8)
res_stloss
res_stloss_sigma
res_auc_alig
res_diff_H
dev.off()

View(out)

# summarie into a table
summary_table <- function(x, out, FUN = median) {
  tab <- tapply(
    x,
    INDEX = list(
      scenario = out$scenario,
      T = out$T,
      m = out$m
    ),
    FUN = FUN
  )

  tab2 <- matrix(round(tab,2), nrow = dim(tab)[1])

  rownames(tab2) <- dimnames(tab)$scenario

  colnames(tab2) <- as.vector(
    outer(
      dimnames(tab)$T,
      dimnames(tab)$m,
      function(T,m) paste0("m = ", m, ", T = ", T)
    )
  )

  tab2
}
# H
summary_table(out$H_median, out)

# stein loss
summary_table(out$stein_loss, out)

# stein loss sigma
summary_table(out$stein_loss_sigma, out)

# AUC
summary_table(out$auc_aligned, out)