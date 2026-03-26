hyperpar <- list(alpha=2, beta=1, kappa=1, sigma=1, a_s = 1, b_s= 1, a_nu = 2, b_nu = 2)
T <- 20
H <- 5
m <- 10
Eta <- matrix(rnorm(T*H),H,T)
set.seed(1)
y <- matrix(rnorm(T*m),T,m)
Lambda <- matrix(rnorm(m*H)*rbinom(m*H,1,prob=0.7),m,H)
Lambda[1:6,4] <- 0
Lambda[1:8,5] <- 0
Lambda[10,5] <- rnorm(1)
Delta <- matrix(as.integer(Lambda != 0),m,H)
Delta
y <- t(Lambda%*% Eta)  + matrix(rnorm(T*m,sd = 0.2),T,m)
nu <- .5
source("aux.R")
source("update.H.R")

updH.step <- update.H(y, Delta, Eta, hyperpar, nu, q = 0.5)
H <- ncol(Delta)
nu <- rbeta(1, hyperpar$a_nu+H, hyperpar$b_nu+m-H)

source("update.pivots.R")
Bad_Delta <- rbind(matrix(0,3,5),diag(1,5),matrix(1,2,5))
Delta <-  Bad_Delta
for(ite in 1:20)
{
cat(ite)
Delta <- update.pivots(y, Delta, Eta, hyperpar)$Delta
}
check_proposal_plot(Bad_Delta,Delta)