FSHL.prior <- function(m, H, gamma, alpha)
{
  tau <- matrix(0,m,H)
  tau_j <- rbeta(H, gamma*alpha/H, gamma)
  l <- sample(1:m, H, replace = F)
  for(j in 1:H)
  {
    tau[l[j],j] <- 1
    tau[(l[j]+1):m,j] <- tau_j[j]
  }
  Delta <- matrix(rbinom(m*H,1,as.double(tau)),m,H)
  return(Delta)
}

set.seed(123)
Delta <- FSHL.prior(10, 5, 5, 4)
Delta