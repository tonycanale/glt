update.H <- function(y, Delta, Eta, hyperpar, nu, q = 0.5){
  paralpha <- hyperpar$alpha
  parbeta <- hyperpar$beta
  parkappa <- hyperpar$kappa
  parsigma <- hyperpar$sigma
  m <- ncol(y)
  H <- ncol(Delta)
  ell <- apply(Delta, 2, match, x=1)
  if((H>1) & (H < m)) increase <- round(runif(1))
  if(H==1){
    increase <- 1
    q <- 1
    }
  if(H==m){
    increase <- 0
    q <- 1
  }
  if(increase){
    Hstar <- H + 1
    Deltastar <- cbind(Delta, rep(0,m))
    ellstar <- ifelse(length((ell[H]+1):m)==1, m, sample((ell[H]+1):m,1))
    Deltastar[ellstar,Hstar] <- 1 
    Etastar <- rbind(Eta, rnorm(T,0,1))
    spourious <- ellstar == m
    if(!spourious){
      Sim1 <- Fim1 <- 0 
      for(i in (ellstar+1):m){
        Deltastar[i,Hstar] <- as.numeric(runif(1) < exp(
          lbeta(paralpha*parbeta + Sim1 + 1, parbeta + Fim1) - 
          lbeta(paralpha*parbeta + Sim1, parbeta + Fim1)))
        Sim1 <- Sim1 + Deltastar[i,Hstar]
        Fim1 <- Fim1 + (1-Deltastar[i,Hstar])
      }
    }
    index <- which(Delta[,H]==1)
    log.lik.R <- sum( marg.lik(index, y=y, Delta=Deltastar, Eta=Etastar, hyperpar=hyperpar) -
      marg.lik(index, y=y, Delta=Delta, Eta=Eta, hyperpar=hyperpar) )
    # to speed it up we probably just need to compute the marginal likelihood for the i in which the 
    # additional column is equal to one, i.e. for those i where  \delta_{iH^*}=1
    R <- exp(log.lik.R + log(nu) - log(1-nu) + log(q)) 
    accept <- runif(1) < min(1,R)
  }
  else{
    Hstar <- H - 1
    Deltastar <- Delta[,-H]
    Etastar <- Eta[-H,]
    # do we actually need Eta?
    index <- which(Delta[,H]==1)
    log.lik.R <- sum( marg.lik(index, y=y, Delta=Deltastar, Eta=Etastar, hyperpar=hyperpar) -
                        marg.lik(index, y=y, Delta=Delta, Eta=Eta, hyperpar=hyperpar) )
      R <- exp(log.lik.R + log(1-nu) - log(nu) + log(1-q))
    accept <- runif(1) < min(1,R)
    # to speed it up we probably just need to compute the marginal likelihood for the i in which the 
    # additional column is equal to one, i.e. for those i where  \delta_{iH^*}=1
      }
  if(accept)
    {
    return(list(Delta = Deltastar,Eta = Etastar, 
                increased = increase, accepted = accept))
    }
  else
    {
    return(list(Delta = Delta, Eta = Eta, 
                increased = increase, accepted = accept))
    }
}