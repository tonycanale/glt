update.pivots <- function(y, Delta, Eta, hyperpar){
  paralpha <- hyperpar$alpha
  parbeta <- hyperpar$beta
  parkappa <- hyperpar$kappa
  parsigma <- hyperpar$sigma
  m <- ncol(y)
  H <- ncol(Delta)
  ell <- apply(Delta, 2, match, x=1)
  # propose new pivots (this may be too drastic, 
  # maybe consider only a perturbation?)
  # maybe in terms of k???
  ### OLD PROPOSAL: ellstar <- sort(sample(1:m,H))
  ### NEW 25/4 proposal --> also add a ranom scan of each column
  Deltaold <- Delta
  Etaold <- Eta
  jstar <- sample(1:H,H)
  for(j in jstar)
  {
  range <- max(c(1,ell[j-1]+1), na.rm=TRUE):min(c(ell[j+1]-1,m), na.rm=TRUE)
  if(length(range)==1){
    next
  } 
  elljstar <- sample(range,1) 
  ellstar <- ell 
  ellstar[j] <- elljstar
  if(prod(ell==ellstar)){
    next
  } 
  ### END OF NEW PROPOSAL
  Deltastar <- Delta
  diag(Deltastar[ellstar,]) <- 1
  #changed pivots
  diffs <- ell - ellstar
  cicle <- which(diffs!=0)
  for(j in cicle)
  {
        if(diffs[j]<0){
          # Case 1: l* > l --> remove structural zeroes above l*
          Deltastar[ell[j]:(ellstar[j]-1),j] <- 0
          }
        else{
          if((ellstar[j]+1)==m) next
          # Case 2: l* < l --> sample potential ones between l*+1 and l
          ## Clear all the potential entries 
          Deltastar[(ellstar[j]+1):(ell[j]),j] = 0
          ## Sample the increase in model size
          dj <- sum(Delta[min((ell[j]+1),m):m,j])
          d_a <- rbetabinom(paralpha*parbeta+dj, 
                            parbeta+m-ell[j]-dj,
                            ell[j]-ellstar[j])
          ## add randomly a d_a number of ones in the candidate entries
          if((d_a>0) & (d_a< (ell[j]-ellstar[j]))) {
            new_ones <- sample((ellstar[j]+1):(ell[j]),size = d_a, 
                               replace=FALSE)
            Deltastar[new_ones,j] = 1
          }
          if(d_a == (ell[j]-ellstar[j])) Deltastar[ell[j],j] = 1
        }
  }
  index <- which(rowSums(abs(Delta- Deltastar))!=0)
  log.lik.R <- sum( marg.lik(index, y=y, Delta=Deltastar, Eta=Eta, hyperpar=hyperpar) -
                      marg.lik(index, y=y, Delta=Delta, Eta=Eta, hyperpar=hyperpar) )
  accept <- runif(1) < min(1,exp(log.lik.R))
  if(accept)
  {
    Delta = Deltastar
    Eta = Eta
    ell <- apply(Delta, 2, match, x=1)
  }
  
  }
  check_proposal_plot(Deltaold, Delta)
  return(list(Delta = Delta, Eta = Eta))
}