GLT.prior <- function(m, a_nu, b_nu, alpha, beta)
{
nu <- rbeta(1,a_nu, b_nu)
k_i <- c(rbinom(m/2,1,nu), rep(0,m/2))
K_i <- cumsum(k_i)
l_j <- which(k_i==1)
H <- max(K_i)
if(H==0)return(NULL)
tau <- matrix(0,m,H)
tau_j <- rbeta(H, alpha, beta)
for(i in 1:m)
{
  if(k_i[i]==1) tau[i,K_i[i]] <- 1
}
for(j in 1:H){
  tau[(l_j[j]+1):m,j] <- tau_j[j]
}
Delta <- matrix(rbinom(m*H,1,as.double(tau)),m,H)
return(Delta)
}

set.seed(123)
Delta <- GLT.prior(10, 1, 1, 4, 1)
Delta