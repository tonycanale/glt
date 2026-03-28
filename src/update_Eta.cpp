// src/update_Eta.cpp
#include "update_Eta.h"

void update_Eta(
    const arma::mat& y,
    const arma::mat& Lambda,
    const arma::vec& sigma2,
    arma::mat&       Eta)
{
  const int T = y.n_rows;
  const int H = Lambda.n_cols;

  // ------------------------------------------------------------------
  // LtSi = Lambda' Sigma^{-1}              (H x m, computed once)
  // Omega = I_H + LtSi * Lambda            (H x H, computed once)
  // L     = chol(Omega, "lower")           (H x H, computed once)
  //
  // For all t simultaneously:
  //   RHS  = LtSi * y'           (H x T)  — one BLAS-3 call
  //   W    = L^{-1}  RHS         (H x T)  — one batched forward  solve
  //   MU   = L'^{-1} W           (H x T)  — one batched backward solve
  //   Z    ~ N(0, I_{H x T})
  //   NS   = L'^{-1} Z           (H x T)  — one batched backward solve
  //   Eta  = MU + NS                       — element-wise add
  //
  // This replaces the previous t-loop with three level-3 BLAS operations.
  // ------------------------------------------------------------------

  const arma::vec s_inv = 1.0 / sigma2;
  const arma::mat LtSi  = Lambda.t() * arma::diagmat(s_inv); // H x m
  const arma::mat Omega = arma::eye(H, H) + LtSi * Lambda;   // H x H

  const arma::mat L = arma::chol(arma::symmatu(Omega), "lower");

  // RHS: H x T
  const arma::mat RHS = LtSi * y.t();

  // Mean: solve L L' MU = RHS  →  W = L^{-1} RHS,  MU = L'^{-1} W
  const arma::mat W  = arma::solve(arma::trimatl(L),   RHS,
                                   arma::solve_opts::fast);
  const arma::mat MU = arma::solve(arma::trimatu(L.t()), W,
                                   arma::solve_opts::fast);

  // Noise: solve L' NS = Z,  Z ~ N(0, I_{H x T})
  const arma::mat Z  = arma::randn<arma::mat>(H, T);
  const arma::mat NS = arma::solve(arma::trimatu(L.t()), Z,
                                   arma::solve_opts::fast);

  Eta = MU + NS;
}
