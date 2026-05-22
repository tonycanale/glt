#ifndef GLTFACTOR_UPDATE_DELTA_H
#define GLTFACTOR_UPDATE_DELTA_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::imat update_delta_cpp(
    const arma::mat&       y,
    const arma::imat&      Delta_in,
    const arma::mat&       Eta,
    const arma::vec&       tau,
    const Rcpp::List&      hyperpar,
    const bool             random_scan = true
);

#endif
