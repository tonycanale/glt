#ifndef GLTFACTOR_UPDATE_LAMBDA_SIGMA2_H
#define GLTFACTOR_UPDATE_LAMBDA_SIGMA2_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// ---------------------------------------------------------------------------
// update_sigma2_Lambda
//
// Joint row-wise sampler for idiosyncratic variances (Step 6) and factor
// loadings (Step 7) of the GLT-FA Gibbs sampler.
//
// Model:   y_t = Lambda * eta_t + eps_t,   eps_t ~ N(0, diag(sigma2))
//
// Lambda is m x H, sparse lower-triangular governed by Delta (m x H, 0/1).
//
// For each row i = 0, ..., m-1:
//
//   q_i = sum_j Delta(i,j)   (number of active loadings in row i)
//
//   -----------------------------------------------------------------------
//   CASE q_i = 0  (zero row, analogous to Step P-a in reference):
//     No loading to sample.  Draw sigma2_i from its closed-form marginal:
//
//       sigma2_i | y_i  ~  InvGamma( a_sigma + T/2,
//                                    b_sigma + 0.5 * sum_t y_{ti}^2 )
//   -----------------------------------------------------------------------
//   CASE q_i > 0  (active row, analogous to Step P-b/c):
//     Let F_i  = Eta.rows( active_cols_of_row_i )  -- q_i x T sub-matrix
//         y_i  = y.col(i)                           -- T-vector
//
//     Normal–inverse-Gamma conjugate update (reference: eq 2948-2952 with
//     the notation change beta -> lambda):
//
//       Prior:   lambda_i^delta | sigma2_i ~ N(0, kappa * sigma2_i * I_{q_i})
//                sigma2_i                  ~ InvGamma(a_sigma, b_sigma)
//
//     Posterior sufficient statistics:
//       FtF   = F_i * F_i'                   (q_i x q_i)
//       Fty   = F_i * y_i                    (q_i)
//       B0inv = I_{q_i} / kappa              (prior precision / sigma2_i)
//
//       V_post^{-1} = FtF + B0inv            (q_i x q_i, proportional to sigma2_i)
//       m_post      = V_post * Fty            (q_i)
//
//     Step 6 — draw sigma2_i | y_i  (marginalising lambda_i^delta):
//       c_T = a_sigma + T/2
//       C_T = b_sigma + 0.5*(y_i'y_i - m_post' * V_post^{-1} * m_post)
//             [equivalently b_sigma + 0.5*SSR where SSR = y'y - Fty'*V_post*Fty]
//       sigma2_i ~ InvGamma(c_T, C_T)
//
//     Step 7 — draw lambda_i^delta | sigma2_i, y_i:
//       lambda_i^delta ~ N( m_post,  sigma2_i * V_post )
//       (sample via Cholesky of  V_post^{-1} / sigma2_i)
//
//     Lambda(i, inactive cols) = 0  (structural/indicator zeros enforced).
//
// Arguments
// ---------
//   y       : T x m data matrix
//   Eta     : H x T factor scores
//   Delta   : m x H binary inclusion indicators (0/1)
//   sigma2  : m-vector of idiosyncratic variances  (updated in-place)
//   Lambda  : m x H loading matrix                 (updated in-place)
//   a_sigma : IG shape hyperparameter for sigma2_i
//   b_sigma : IG rate hyperparameter for sigma2_i
//   kappa   : prior variance scalar for active loadings
// ---------------------------------------------------------------------------

void update_sigma2_Lambda(
    const arma::mat&  y,         // T x m
    const arma::mat&  Eta,       // H x T
    const arma::imat& Delta,     // m x H
    arma::vec&        sigma2,    // m      (in/out)
    arma::mat&        Lambda,    // m x H  (in/out)
    double            a_sigma,
    double            b_sigma,
    double            kappa
);

#endif  // GLTFACTOR_UPDATE_LAMBDA_SIGMA2_H
