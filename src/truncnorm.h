#ifndef GLTFACTOR_TRUNCNORM_H
#define GLTFACTOR_TRUNCNORM_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Sample from a truncated normal distribution. Samples are drawn
//' componentwise, so each component of the vector is allowed its own
//' mean, standard deviation, and upper and lower limits. The components
//' are assumed to be independent.
//'
//' @param y_lower \code{n x p} matrix of lower endpoints
//' @param y_upper \code{n x p} matrix of upper endpoints
//' @param mu \code{n x p} matrix of conditional expectations
//' @param sigma \code{p x 1} vector of conditional standard deviations
//' @param u_rand \code{n x p} matrix of uniform random variables
//'
//' @return z_star \code{n x p} draw from the truncated normal distribution
//'
//' @note This function uses \code{Rcpp} for computational efficiency.
//' Bounds may be \code{-Inf}/\code{Inf}; these are handled correctly by
//' \code{R::pnorm}/\code{R::qnorm}.
//'
arma::mat truncnorm_lg(const arma::mat& y_lower, const arma::mat& y_upper,
                       const arma::mat& mu, const arma::vec& sigma,
                       const arma::mat& u_rand);

#endif // GLTFACTOR_TRUNCNORM_H
