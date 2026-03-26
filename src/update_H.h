
#ifndef GLTFACTOR_UPDATE_H_H
#define GLTFACTOR_UPDATE_H_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Update H via birth/death MCMC step
// Returns list with: Delta, Eta, increased, accepted
Rcpp::List update_H_cpp(
    const arma::mat& y,
    const arma::imat& Delta_in,
    const arma::mat& Eta_in,
    const Rcpp::List& hyperpar,
    double nu,
    double q = 0.5
);

#endif

