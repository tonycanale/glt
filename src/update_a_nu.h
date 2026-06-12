#ifndef GLTFACTOR_UPDATE_A_NU_H
#define GLTFACTOR_UPDATE_A_NU_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Update a_nu via log-normal MH (nu marginalised out)
double update_a_nu_cpp(
    double a_nu,
    int    H,
    int    m,
    double b_nu,
    double a_anu,
    double b_anu,
    double mh_sd_a_nu
);

#endif // GLTFACTOR_UPDATE_A_NU_H
