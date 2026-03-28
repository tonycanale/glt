// src/update_Lambda_sigma2.cpp
#include <RcppArmadillo.h>
#include "update_Lambda_sigma2.h"

// [[Rcpp::depends(RcppArmadillo)]]

static inline double rinvgamma(double shape, double rate) {
  return 1.0 / R::rgamma(shape, 1.0 / rate);
}

void update_sigma2_Lambda(
    const arma::mat&  y,        // T x m
    const arma::mat&  Eta,      // H x T
    const arma::imat& Delta,    // m x H
    arma::vec&        sigma2,   // m      (updated in-place)
    arma::mat&        Lambda,   // m x H  (updated in-place)
    double            a_sigma,
    double            b_sigma,
    double            kappa)
{
  const int m = static_cast<int>(y.n_cols);
  const int T = static_cast<int>(y.n_rows);

  // ------------------------------------------------------------------
  // Precompute the full H x H Gram matrix ONCE.
  // Each row i only needs the qi x qi submatrix EtaEta'(active, active),
  // which is extracted cheaply rather than recomputed as F * F'.
  // ------------------------------------------------------------------
  const arma::mat EtEt = Eta * Eta.t();     // H x H

  Lambda %= arma::conv_to<arma::mat>::from(Delta);

  for (int i = 0; i < m; ++i) {

    arma::uvec active = arma::find(
        arma::conv_to<arma::uvec>::from(Delta.row(i)));
    const int q = static_cast<int>(active.n_elem);

    const arma::vec yi = y.col(i);

    // ------------------------------------------------------------------
    // Case 1 — zero row
    // ------------------------------------------------------------------
    if (q == 0) {
      sigma2(i) = rinvgamma(a_sigma + 0.5 * T,
                            b_sigma + 0.5 * arma::dot(yi, yi));
      Lambda.row(i).zeros();
      continue;
    }

    // ------------------------------------------------------------------
    // Case 2 — active row
    //
    // Omega = EtEt(active,active) + (1/kappa) I_q   (q x q)
    //
    // Single Cholesky L s.t. Omega = L L' is used for:
    //   (a) m_post  = Omega^{-1} F yi   via two triangular solves
    //   (b) SSR     = yi'yi - Fyi' m_post
    //   (c) lambda draw:  lambda = m_post + sqrt(sigma2_i) L^{-T} z
    //       (same L, no second factorisation needed)
    // ------------------------------------------------------------------
    const arma::mat Omega =
        EtEt(active, active) + arma::eye(q, q) / kappa;    // q x q

    arma::mat L;
    if (!arma::chol(L, arma::symmatu(Omega), "lower")) {
      Rcpp::warning("Row %d: Cholesky of Omega failed; skipping.", i);
      continue;
    }

    // F yi  =  Eta(active, :) * yi   (q x 1)
    const arma::vec Fyi    = Eta.rows(active) * yi;
    const arma::vec w      = arma::solve(arma::trimatl(L), Fyi,
                                         arma::solve_opts::fast);
    const arma::vec m_post = arma::solve(arma::trimatu(L.t()), w,
                                         arma::solve_opts::fast);

    // Step 6 — sigma2_i
    const double SSR = arma::dot(yi, yi) - arma::dot(Fyi, m_post);
    const double C_T = b_sigma + 0.5 * SSR;
    if (C_T <= 0.0) {
      Rcpp::warning("Row %d: C_T = %g <= 0; skipping.", i, C_T);
      continue;
    }
    sigma2(i) = rinvgamma(a_sigma + 0.5 * T, C_T);

    // Step 7 — lambda_i | sigma2_i  using the SAME L
    const arma::vec z     = arma::randn<arma::vec>(q);
    const arma::vec noise = arma::solve(arma::trimatu(L.t()), z,
                                        arma::solve_opts::fast);

    Lambda.row(i).zeros();
    const arma::vec lam_i = m_post + std::sqrt(sigma2(i)) * noise;
    for (int k = 0; k < q; ++k)
      Lambda(i, active(k)) = lam_i(k);
  }
}
