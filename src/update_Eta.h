#ifndef GLTFACTOR_UPDATE_ETA_H
#define GLTFACTOR_UPDATE_ETA_H

#include <RcppArmadillo.h>

// ---------------------------------------------------------------------------
// update_Eta
//
// Step 8: draw eta_t | y_t, Lambda, sigma2  for t = 1,...,T.
//
// Model:  y_t = Lambda * eta_t + eps_t,   eps_t ~ N(0, Sigma)
//         eta_t ~ N(0, I_H)
// where Sigma = diag(sigma2).
//
// Posterior (standard Gaussian linear model, prior N(0, I_H)):
//   eta_t | y_t ~ N( V * Lambda' * Sigma^{-1} * y_t,  V )
//   V = ( I_H + Lambda' * Sigma^{-1} * Lambda )^{-1}
//
// Implementation:
//   Omega = I_H + Lambda' * Sigma^{-1} * Lambda   is computed once.
//   Its Cholesky L (Omega = L L') is reused for all T draws:
//     - posterior mean : solve two triangular systems  L L' mu = rhs
//     - posterior draw : eta_t = mu_t + L'^{-1} z,  z ~ N(0, I_H)
//
// Arguments:
//   y      (T x m)  observed data, rows are y_t'
//   Lambda (m x H)  factor loadings (zeros already applied via Delta)
//   sigma2 (m)      idiosyncratic variances
//   Eta    (H x T)  factors, updated in-place, columns are eta_t
// ---------------------------------------------------------------------------
void update_Eta(
    const arma::mat& y,
    const arma::mat& Lambda,
    const arma::vec& sigma2,
    arma::mat&       Eta
);

#endif // GLTFACTOR_UPDATE_ETA_H
