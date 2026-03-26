
#ifndef GLTFACTOR_UPDATE_H_H
#define GLTFACTOR_UPDATE_H_H

#include <RcppArmadillo.h>
#include "marglik.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Update H via birth/death MCMC step
// Returns list with: Delta, Eta, increased, accepted
inline Rcpp::List update_H_cpp(
    const arma::mat& y,
    const arma::imat& Delta_in,
    const arma::mat& Eta_in,
    const Rcpp::List& hyperpar,
    double nu,
    double q = 0.5
) {
  Rcpp::RNGScope rng;
  
  const int Tt = static_cast<int>(y.n_rows);
  const int m = static_cast<int>(y.n_cols);
  const int H = static_cast<int>(Delta_in.n_cols);
  
  arma::imat Delta = Delta_in;
  arma::mat Eta = Eta_in;
  
  double paralpha = Rcpp::as<double>(hyperpar["alpha"]);
  double parbeta = Rcpp::as<double>(hyperpar["beta"]);
  
  // Compute pivots (1-based)
  std::vector<int> ell(H);
  for (int j = 0; j < H; ++j) {
    int found = m;
    for (int i = 0; i < m; ++i) {
      if (Delta(i, j) == 1) { found = i + 1; break; }
    }
    ell[j] = found;
  }
  
  int increase = 0;
  if ((H > 1) && (H < m)) {
    increase = (int)std::round(Rcpp::runif(1)[0]);
  } else if (H == 1) {
    increase = 1;
    q = 1.0;
  } else if (H == m) {
    increase = 0;
    q = 1.0;
  }
  
  bool accepted = false;
  arma::imat Deltastar;
  arma::mat Etastar;
  int Hstar = H;
  
  if (increase) {
    Hstar = H + 1;
    Deltastar = arma::imat(m, Hstar, arma::fill::zeros);
    for (int j = 0; j < H; ++j) Deltastar.col(j) = Delta.col(j);
    
    int last_ell = ell[H - 1];
    int ellstar = m;
    if (last_ell < m) {
      int start = last_ell + 1;
      int ncand = m - start + 1;
      if (ncand > 1) {
        int pick = (int)std::floor(R::unif_rand() * ncand);
        ellstar = start + pick;
      }
    }
    Deltastar(ellstar - 1, Hstar - 1) = 1;
    
    Etastar = arma::mat(Hstar, Tt, arma::fill::zeros);
    for (int j = 0; j < H; ++j) Etastar.row(j) = Eta.row(j);
    arma::rowvec neweta(Tt);
    for (int t = 0; t < Tt; ++t) neweta[t] = R::norm_rand();
    Etastar.row(Hstar - 1) = neweta;
    
    bool spurious = (ellstar == m);
    if (!spurious) {
      int Sim1 = 0, Fim1 = 0;
      for (int i = ellstar + 1; i <= m; ++i) {
        double a1 = paralpha * parbeta + Sim1 + 1.0;
        double b1 = parbeta + Fim1;
        double a0 = paralpha * parbeta + Sim1;
        double b0 = parbeta + Fim1;
        double log_num = R::lbeta(a1, b1);
        double log_den = R::lbeta(a0, b0);
        double p = std::exp(log_num - log_den);
        double u = R::unif_rand();
        int draw = (u < p) ? 1 : 0;
        Deltastar(i - 1, Hstar - 1) = draw;
        Sim1 += draw;
        Fim1 += (1 - draw);
      }
    }
    
    std::vector<int> idx;
    for (int i = 0; i < m; ++i) if (Delta(i, H - 1) == 1) idx.push_back(i + 1);
    Rcpp::IntegerVector indexR = Rcpp::wrap(idx);
    
    Rcpp::NumericVector new_ll = marg_lik(indexR, y, Deltastar, Etastar, hyperpar);
    Rcpp::NumericVector old_ll = marg_lik(indexR, y, Delta, Eta, hyperpar);
    
    double loglikR = 0.0;
    for (int i = 0; i < (int)new_ll.size(); ++i) loglikR += (new_ll[i] - old_ll[i]);
    
    double Rval = std::exp(loglikR + std::log(nu) - std::log(1.0 - nu) + std::log(q));
    double u = R::unif_rand();
    accepted = (u < std::min(1.0, Rval));
  } else {
    if (H <= 1) {
      return Rcpp::List::create(
        Rcpp::Named("Delta") = Delta,
        Rcpp::Named("Eta") = Eta,
        Rcpp::Named("increased") = increase,
        Rcpp::Named("accepted") = false
      );
    }
    
    Hstar = H - 1;
    Deltastar = arma::imat(m, Hstar, arma::fill::zeros);
    for (int j = 0; j < Hstar; ++j) Deltastar.col(j) = Delta.col(j);
    Etastar = arma::mat(Hstar, Tt, arma::fill::zeros);
    for (int j = 0; j < Hstar; ++j) Etastar.row(j) = Eta.row(j);
    
    std::vector<int> idx;
    for (int i = 0; i < m; ++i) if (Delta(i, H - 1) == 1) idx.push_back(i + 1);
    Rcpp::IntegerVector indexR = Rcpp::wrap(idx);
    
    Rcpp::NumericVector new_ll = marg_lik(indexR, y, Deltastar, Etastar, hyperpar);
    Rcpp::NumericVector old_ll = marg_lik(indexR, y, Delta, Eta, hyperpar);
    
    double loglikR = 0.0;
    for (int i = 0; i < (int)new_ll.size(); ++i) loglikR += (new_ll[i] - old_ll[i]);
    
    double Rval = std::exp(loglikR + std::log(1.0 - nu) - std::log(nu) + std::log(1.0 - q));
    double u = R::unif_rand();
    accepted = (u < std::min(1.0, Rval));
  }
  
  if (accepted) {
    return Rcpp::List::create(
      Rcpp::Named("Delta") = Deltastar,
      Rcpp::Named("Eta") = Etastar,
      Rcpp::Named("increased") = increase,
      Rcpp::Named("accepted") = true
    );
  } else {
    return Rcpp::List::create(
      Rcpp::Named("Delta") = Delta,
      Rcpp::Named("Eta") = Eta,
      Rcpp::Named("increased") = increase,
      Rcpp::Named("accepted") = false
    );
  }
}

#endif

