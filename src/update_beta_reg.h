#ifndef GLTFACTOR_UPDATE_BETA_REG_H
#define GLTFACTOR_UPDATE_BETA_REG_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// ---------------------------------------------------------------------------
// update_Beta_cpp
//
// Step 0b of the STAR sampler: update the regression coefficients Beta
// (p x m) that govern the mean of the latent z, given the current residual
// factor-model draw.
//
// Model (column-wise):
//   z_j = X * beta_j + epsilon_j,     epsilon_j ~ N(0, sigma2_j I_T)
//   beta_j ~ N(b0, sigma2_j * V0)     (conjugate to sigma2_j, as with Lambda)
//
// Posterior:
//   V_post   = ( V0^{-1} + X'X )^{-1}
//   m_post   = V_post * ( V0^{-1} b0 + X' r_j )     where r_j = z_j (residual
//                                                     AFTER removing the
//                                                     factor part, i.e.
//                                                     r_j = z_j - (Lambda Eta)_j )
//   beta_j | sigma2_j, r_j ~ N( m_post, sigma2_j * V_post )
//
// Arguments
// ---------
//   r        : T x m matrix, residual z - (Lambda * Eta)^T (regression target)
//   X        : T x p design matrix (p may be 0, see below)
//   XtX      : p x p precomputed X'X (pass empty matrix if p == 0)
//   sigma2   : m-vector of idiosyncratic variances
//   V0_inv   : p x p prior precision (V0^{-1}) common to all columns
//   b0       : p-vector prior mean common to all columns
//   Beta     : p x m coefficient matrix, updated in place
//
// If p == 0 (no covariates), this function is a no-op (Beta has 0 rows).
// ---------------------------------------------------------------------------
inline void update_Beta_cpp(
    const arma::mat&  r,        // T x m
    const arma::mat&  X,        // T x p
    const arma::mat&  XtX,      // p x p
    const arma::vec&  sigma2,   // m
    const arma::mat&  V0_inv,   // p x p
    const arma::vec&  b0,       // p
    arma::mat&        Beta      // p x m (in/out)
) {
  const int p = X.n_cols;
  const int m = r.n_cols;

  if (p == 0) {
    return; // nothing to update
  }

  const arma::mat V_post_inv = V0_inv + XtX;      // p x p
  const arma::mat L          = arma::chol(V_post_inv, "lower"); // Cholesky: V_post_inv = L L'
  const arma::vec V0_inv_b0  = V0_inv * b0;

  for (int j = 0; j < m; ++j) {
    const arma::vec rhs = V0_inv_b0 + X.t() * r.col(j);   // p
    // Solve V_post_inv * m_post = rhs via the Cholesky factor L
    arma::vec w      = arma::solve(arma::trimatl(L), rhs);
    arma::vec m_post  = arma::solve(arma::trimatu(L.t()), w);

    // Draw beta_j ~ N(m_post, sigma2_j * V_post), V_post = V_post_inv^{-1}
    arma::vec z(p, arma::fill::randn);
    arma::vec draw_noise = arma::solve(arma::trimatu(L.t()), z); // L'^{-1} z has cov V_post
    Beta.col(j) = m_post + std::sqrt(sigma2(j)) * draw_noise;
  }
}

#endif // GLTFACTOR_UPDATE_BETA_REG_H
