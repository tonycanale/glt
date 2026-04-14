#ifndef GLTFACTOR_UPDATE_PIVOTS_H
#define GLTFACTOR_UPDATE_PIVOTS_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::List update_pivots_cpp(
    const arma::mat& y,
    const arma::imat& Delta_in,
    const arma::mat& Eta,
    const Rcpp::List& hyperpar
);

#endif
