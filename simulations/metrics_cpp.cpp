// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
using namespace Rcpp;
using namespace arma;

// ------------------------------------------------------------------ //
// Stein loss: tr(A B^{-1}) - log|det(A B^{-1})| - m
//
// Uses Cholesky of B for a symmetric eigenvalue decomposition, avoiding
// the costly solve(B) + eigen of a non-symmetric matrix in plain R.
// ------------------------------------------------------------------ //
static double stein_loss_arma(const mat& A, const mat& B) {
  mat L;
  vec ev;
  if (chol(L, B, "lower")) {
    // C = L^{-1} A L^{-T}  (symmetric, same eigenvalues as A B^{-1})
    mat X = solve(trimatl(L), A);
    mat C = solve(trimatl(L), X.t());
    ev = eig_sym(symmatu(C));
  } else {
    ev = real(eig_gen(A * inv(B)));
  }
  ev.clamp(datum::eps, datum::inf);
  return sum(ev) - sum(log(ev)) - static_cast<double>(A.n_rows);
}

// ------------------------------------------------------------------ //
// Greedy column alignment by Hamming similarity
// ------------------------------------------------------------------ //
static imat align_cols_arma(const imat& D, const imat& Dt, int H_s, int H0) {
  int m       = static_cast<int>(D.n_rows);
  int n_match = std::min(H_s, H0);
  imat out(m, H0, arma::fill::zeros);
  std::vector<bool> used(H_s, false);

  for (int j = 0; j < n_match; j++) {
    int best = -1, best_sim = -1;
    for (int k = 0; k < H_s; k++) {
      if (used[k]) continue;
      int sim = 0;
      for (int i = 0; i < m; i++) sim += (D(i, k) == Dt(i, j));
      if (sim > best_sim) { best_sim = sim; best = k; }
    }
    out.col(j) = D.col(best);
    used[best] = true;
  }
  return out;
}

// ------------------------------------------------------------------ //
// Exported: count spurious columns (colSum == 1) for every draw
// ------------------------------------------------------------------ //
// [[Rcpp::export]]
IntegerVector compute_spurious_cpp(List Delta_draws) {
  int n = Delta_draws.size();
  IntegerVector out(n);
  for (int i = 0; i < n; i++) {
    imat D   = as<imat>(Delta_draws[i]);
    int  cnt = 0;
    for (uword k = 0; k < D.n_cols; k++) cnt += (static_cast<int>(accu(D.col(k))) == 1);
    out[i] = cnt;
  }
  return out;
}

// ------------------------------------------------------------------ //
// Exported: inner accumulation loop
//
// Arguments:
//   Lambda_draws   list[nsave] of m x H_draw numeric matrices
//   Delta_draws    list[nsave] of m x H_draw integer matrices
//   sigma2_draws   nsave x m numeric matrix
//   H_draws        integer vector length nsave (already spurious-adjusted)
//   Delta_true     m x H0 integer matrix
//   Sigma_true     numeric vector length m
//   Omega0_r       m x m matrix (Lambda_true %*% t(Lambda_true) + diag(Sigma_true))
// ------------------------------------------------------------------ //
// [[Rcpp::export]]
List compute_metrics_inner_cpp(
    List           Lambda_draws,
    List           Delta_draws,
    NumericMatrix  sigma2_draws,
    IntegerVector  H_draws,
    IntegerMatrix  Delta_true,
    NumericVector  Sigma_true,
    NumericMatrix  Omega0_r
) {
  int nsave = Lambda_draws.size();
  int m     = Delta_true.nrow();
  int H0    = Delta_true.ncol();

  imat Dt     = as<imat>(Delta_true);
  mat  Om0    = as<mat>(Omega0_r);
  mat  diagSt = diagmat(as<vec>(Sigma_true));

  double ss = 0.0, ss_sig = 0.0;
  mat Dacc   (m, H0, fill::zeros);
  mat Dacc_al(m, H0, fill::zeros);

  for (int s = 0; s < nsave; s++) {
    int H_s = H_draws[s];

    // --- Stein losses ---
    mat Lam = as<mat>(Lambda_draws[s]);
    Lam     = Lam.cols(0, H_s - 1);
    vec sig2(m);
    for (int i = 0; i < m; i++) sig2(i) = sigma2_draws(s, i);
    mat Omega_s = Lam * Lam.t() + diagmat(sig2);

    ss     += stein_loss_arma(Omega_s,       Om0);
    ss_sig += stein_loss_arma(diagmat(sig2), diagSt);

    // --- Delta accumulation ---
    imat Del = as<imat>(Delta_draws[s]);
    int  Hd  = static_cast<int>(Del.n_cols);

    // strip spurious columns (colSum == 1) when Hd != H_s
    if (Hd != H_s) {
      uvec keep;
      for (int k = 0; k < Hd; k++)
        if (static_cast<int>(accu(Del.col(k))) > 1)
          keep = join_cols(keep, uvec{static_cast<uword>(k)});
      Del = keep.is_empty() ? imat(m, 0) : Del.cols(keep);
      H_s = static_cast<int>(Del.n_cols);
    }

    if (H_s == 0) continue;

    // aligned accumulation
    Dacc_al += conv_to<mat>::from(align_cols_arma(Del, Dt, H_s, H0));

    // unaligned accumulation (first H0 cols, zero-pad if H_s < H0)
    if (H_s >= H0) {
      Dacc += conv_to<mat>::from(Del.cols(0, H0 - 1));
    } else {
      mat pad(m, H0, fill::zeros);
      pad.cols(0, H_s - 1) = conv_to<mat>::from(Del);
      Dacc += pad;
    }
  }

  return List::create(
    Named("stein_sum")         = ss,
    Named("stein_sum_sigma")   = ss_sig,
    Named("Delta_acc")         = Dacc,
    Named("Delta_acc_aligned") = Dacc_al
  );
}
