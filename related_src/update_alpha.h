#ifndef GLTFACTOR_UPDATE_ALPHA_H
#define GLTFACTOR_UPDATE_ALPHA_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

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

#endif // GLTFACTOR_UPDATE_ALPHA_H
