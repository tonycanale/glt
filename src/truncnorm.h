#ifndef GLTFACTOR_TRUNCNORM_H
#define GLTFACTOR_TRUNCNORM_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat truncnorm_lg(const arma::mat& y_lower, const arma::mat& y_upper,
                       const arma::mat& mu, const arma::vec& sigma,
                       const arma::mat& u_rand);

#endif // GLTFACTOR_TRUNCNORM_H
