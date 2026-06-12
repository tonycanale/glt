#ifndef GLTFACTOR_UPDATE_ALPHA_BETA_H
#define GLTFACTOR_UPDATE_ALPHA_BETA_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Update alpha via log-normal MH (beta fixed, tau marginalised out)
double update_alpha_cpp(
    double                        alpha,
    const arma::imat&             Delta,
    const Rcpp::IntegerVector&    ell,
    int                           m,
    double                        beta,
    double                        a_alpha,
    double                        b_alpha,
    double                        mh_sd_alpha
);

// Update beta via log-normal MH (alpha fixed, tau marginalised out)
double update_beta_cpp(
    double                        beta,
    const arma::imat&             Delta,
    const Rcpp::IntegerVector&    ell,
    int                           m,
    double                        alpha,
    double                        a_beta,
    double                        b_beta,
    double                        mh_sd_beta
);

#endif // GLTFACTOR_UPDATE_ALPHA_BETA_H
