#@ marginal likelihood function for the single row of Delta (column of y)
marg.lik.i <- function(i, y, Delta, Eta, hyperpar, logarithm=TRUE)
{
  T <- nrow(y)
  kappa <- hyperpar$kappa
  a_s <- hyperpar$a_s
  b_s <- hyperpar$b_s
  yi <- matrix(y[,i],ncol = 1)
  Deltai <- Delta[i,]
  qi <- sum(Deltai)
  if(qi>0){
  eta <- matrix(t(Eta[which(Deltai==1),]),T,qi)
  inner.mat <- diag(1/kappa, qi) + t(eta) %*% eta
  inner.mat.inv <-  solve(inner.mat)
  Qi <- 0.5 * as.double(t(yi) %*% yi - t(yi) %*% eta %*% (inner.mat.inv) %*% t(eta) %*% yi )
  res <- -T/2 * log(2*pi) + lgamma(a_s + T/2) - lgamma(a_s) + 
    a_s * log(b_s) - (a_s + T/2) * log(b_s + Qi) - (qi/2)*log(kappa) - 
    0.5*determinant( inner.mat , log=TRUE )$modulus 
  }
  if(qi==0){
    res <- -T/2 * log(2*pi) + lgamma(a_s + T/2) - lgamma(a_s) + 
      a_s * log(b_s) - (a_s + T/2) * log(b_s+ 0.5*sum(yi^2)) 
  }
  ifelse(logarithm, as.double(res), as.double(exp(res)))
}

#@ marginal likelihood function for a number of rows of Delta (columns of y)
marg.lik <- function(index, y, Delta, Eta, hyperpar, logarithm =TRUE){
  all <- sapply(index, marg.lik.i, y=y, Delta=Delta, Eta=Eta, hyperpar=hyperpar, logarithm=logarithm)
  if(logarithm)
  {
    return((all))
  }
  else
  {
    return((all))
  }
}

#@ log posterior odds with marginal likelihood 
log.post.odds.ij <- function(i, j, y, Delta, Eta, tau, hyperpar, logarithm=TRUE)
{
  Delta_zero <- Delta
  Delta_one <- Delta
  Delta_zero[i,j] <- 0
  Delta_one[i,j] <- 1
  log_prior_odds <- log(tau[j]) - log(1-tau[j])
  ll_one <- marg.lik(i, y, Delta_one, Eta, hyperpar)
  ll_zero <- marg.lik(i, y, Delta_zero, Eta, hyperpar)
  res <- (ll_one - ll_zero) + log_prior_odds
  ifelse(logarithm, as.double(res), as.double(exp(res)))
}

# get pivot function
get_pivots <- function(Delta) {
  if (!is.matrix(Delta)) {
    stop("`Delta` must be a matrix.", call. = FALSE)
  }
  pivots <- apply(Delta, 2L, match, x = 1L)
  if (length(pivots) != length(unique(pivots))) {
    stop("Pivots need to be on different rows.", call. = FALSE)
  }
  if (anyNA(pivots)) {
    stop("Each column of `Delta` must contain at least one 1 to define a pivot.", call. = FALSE)
  }
  as.integer(pivots)
}

rbetabinom <- function(a,b,size){
  probs <- exp(lgamma(size+1)-lgamma(0:size+1)-lgamma(size-0:size+1) + 
    lbeta(a+0:size,b+size:0) - lbeta(a,b))
  sample(0:size, 1, prob=probs)
}

check_proposal_plot <- function(Delta, Deltastar){
  H <- ncol(Delta)
  m <- nrow(Delta)
  Hstar <- ncol(Deltastar)
  ell <- apply(Delta, 2, match, x=1)
  ellstar <- apply(Deltastar, 2, match, x=1)
  par(mfrow=c(1,2))
  Deltaplot <- Delta
  diag(Deltaplot[ell,]) <- 2
  image(1:H,1:m, t(Deltaplot)[,m:1])
  
  Deltastplot <- Deltastar
  diag(Deltastplot[ellstar,]) <- 2
  image(1:Hstar,1:m, t(Deltastplot)[,m:1])
  par(mfrow=c(1,1))
}