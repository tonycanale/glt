update.delta <- function(Delta, Eta_star, tau, y,  hyperpar)
{
  ell <- get_pivots(Delta)
  Delta_out <- Delta
  jind <- sample(1:ncol(Delta),ncol(Delta),replace=FALSE)
  # ciclo sulle colonne random
  for(j in jind){
    # ciclo sulle righe sotto il pivot
    accepted <- rep(FALSE, nrow(Delta)-ell[j])
    for(i in (ell[j]+1):nrow(Delta)){
      logu <- log(runif(1,0,1))
      Oij <- log.post.odds.ij(i, j, y, Delta, Eta_star, tau,
        hyperpar, logarithm=TRUE)
      if (Delta[i, j] == 0L) {
        accepted[i-ell[j]] <- logu <= Oij        
      } 
      else {
        accepted[i-ell[j]] <- logu <= -Oij
      }
    }
    #dopo aver deciso di accettare/rifiutare implemento i cambiamenti
    for(i in (ell[j]+1):nrow(Delta)){
      if (Delta[i, j] == 0L & accepted[i-ell[j]]) {
          Delta_out[i, j] <- 1L
      } 
      if (Delta[i, j] == 1L & accepted[i-ell[j]]) {
          Delta_out[i, j] <- 0L
      } 
    }
    print(accepted)
    #passo alla colonna successiva dopo aver aggiornato la precedente
    Delta <- Delta_out
  }
  return(Delta_out)
  }