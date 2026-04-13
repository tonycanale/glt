#ifndef GLTFACTOR_HELPERS_H
#define GLTFACTOR_HELPERS_H

#include <RcppArmadillo.h>
#include "marglik.h"

// ---------------------------------------------------------------------------
// Scalar extractors from named lists
// ---------------------------------------------------------------------------

inline int get_int(const Rcpp::List& x, const char* name) {
  SEXP obj = x[name];
  if (obj == R_NilValue) Rcpp::stop("Missing entry `%s`.", name);
  return Rcpp::as<int>(obj);
}

inline double get_double(const Rcpp::List& x, const char* name) {
  SEXP obj = x[name];
  if (obj == R_NilValue) Rcpp::stop("Missing entry `%s`.", name);
  return Rcpp::as<double>(obj);
}

inline bool get_bool(const Rcpp::List& x, const char* name) {
  SEXP obj = x[name];
  if (obj == R_NilValue) Rcpp::stop("Missing entry `%s`.", name);
  return Rcpp::as<bool>(obj);
}

// ---------------------------------------------------------------------------
// Helper: extract 1-based pivots from Delta.
// For each column j, ell[j] is the first row in {1,...,m} where Delta == 1.
// Used by update_pivots_cpp and update_delta_cpp.
// ---------------------------------------------------------------------------
inline std::vector<int> pivots_from_delta_1based(const arma::imat& Delta) {
  const int m = static_cast<int>(Delta.n_rows);
  const int H = static_cast<int>(Delta.n_cols);
  std::vector<int> ell(H, m);
  for (int j = 0; j < H; ++j) {
    for (int i = 0; i < m; ++i) {
      if (Delta(i, j) == 1) {
        // remember ell_j are {1, ..., m} but C indexing starts at 0
        ell[j] = i + 1;
        break;
      }
    }
  }
  return ell;
}

// ---------------------------------------------------------------------------
// Helper: set pivot diagonal entries in Delta from a 1-based ell vector.
// For each column j, sets Delta(ell[j]-1, j) = 1.
// Used by update_pivots_cpp.
// ---------------------------------------------------------------------------
inline void set_pivot_diagonal_from_ell_1based(
    arma::imat&             Delta,
    const std::vector<int>& ell
) {
  const int H = static_cast<int>(Delta.n_cols);
  for (int j = 0; j < H; ++j) {
    // remember ell_j are {1, ..., m} but C indexing starts at 0
    Delta(ell[j] - 1, j) = 1;
  }
}

// ---------------------------------------------------------------------------
// Helper: return 1-based row indices where matrices A and B differ.
// Used by update_pivots_cpp to identify which rows changed for marg_lik().
// ---------------------------------------------------------------------------
inline Rcpp::IntegerVector changed_rows_1based(
    const arma::imat& A,
    const arma::imat& B
) {
  const int m = static_cast<int>(A.n_rows);
  std::vector<int> index;
  index.reserve(m);
  for (int i = 0; i < m; ++i) {
    bool different = false;
    for (int j = 0; j < static_cast<int>(A.n_cols); ++j) {
      if (A(i, j) != B(i, j)) { different = true; break; }
    }
    if (different) {
      // marg_lik expects 1-based row indices
      index.push_back(i + 1);
    }
  }
  return Rcpp::wrap(index);
}

// ---------------------------------------------------------------------------
// Helper: log posterior odds for entry (i_1based, j_0based) of Delta.
//
//   log p(delta_ij=1|rest) - log p(delta_ij=0|rest)
//     = log(tau_j / (1-tau_j))
//       + marg_lik(i | Delta_one) - marg_lik(i | Delta_zero)
//
// i_1based : 1-based row index (R convention)
// j_0based : 0-based column index (C++ convention)
// Used by update_delta_cpp.
// ---------------------------------------------------------------------------
inline double log_post_odds_ij(
    int               i_1based,
    int               j_0based,
    const arma::mat&  y,
    const arma::imat& Delta,
    const arma::mat&  Eta,
    const arma::vec&  tau,
    const Rcpp::List& hyperpar
) {
  arma::imat Delta_zero = Delta;
  arma::imat Delta_one  = Delta;
  // remember ell_j are {1, ..., m} but C indexing starts at 0
  Delta_zero(i_1based - 1, j_0based) = 0;
  Delta_one (i_1based - 1, j_0based) = 1;

  Rcpp::IntegerVector idx = Rcpp::IntegerVector::create(i_1based);
  const double ll_one  = Rcpp::as<double>(
      marg_lik(idx, y, Delta_one,  Eta, hyperpar)[0]);
  const double ll_zero = Rcpp::as<double>(
      marg_lik(idx, y, Delta_zero, Eta, hyperpar)[0]);

  const double log_prior_odds =
      std::log(tau[j_0based]) - std::log(1.0 - tau[j_0based]);
  return (ll_one - ll_zero) + log_prior_odds;
}

// ---------------------------------------------------------------------------
// Build a simple default state when H = 0 or when init pieces are NULL
// ---------------------------------------------------------------------------
inline Rcpp::List make_last_state(
    int H, int T, int m,
    SEXP nu_, SEXP ell_, SEXP tau_,
    SEXP Delta_, SEXP Lambda_, SEXP Eta_, SEXP sigma2_
) {
  Rcpp::IntegerVector ell;
  Rcpp::NumericVector tau;
  arma::imat Delta;
  arma::mat  Lambda;
  arma::mat  Eta;
  Rcpp::NumericVector sigma2;

  ell    = (ell_    == R_NilValue) ? Rcpp::IntegerVector(H)              : Rcpp::as<Rcpp::IntegerVector>(ell_);
  tau    = (tau_    == R_NilValue) ? Rcpp::NumericVector(H)              : Rcpp::as<Rcpp::NumericVector>(tau_);
  Delta  = (Delta_  == R_NilValue) ? arma::imat(m, H, arma::fill::zeros) : Rcpp::as<arma::imat>(Delta_);
  Lambda = (Lambda_ == R_NilValue) ? arma::mat(m, H, arma::fill::zeros)  : Rcpp::as<arma::mat>(Lambda_);
  Eta    = (Eta_    == R_NilValue) ? arma::mat(H, T, arma::fill::zeros)  : Rcpp::as<arma::mat>(Eta_);
  sigma2 = (sigma2_ == R_NilValue) ? Rcpp::NumericVector(m, 1.0)         : Rcpp::as<Rcpp::NumericVector>(sigma2_);

  double nu = (nu_ == R_NilValue) ? NA_REAL : Rcpp::as<double>(nu_);

  return Rcpp::List::create(
    Rcpp::Named("H")      = H,
    Rcpp::Named("nu")     = nu,
    Rcpp::Named("ell")    = ell,
    Rcpp::Named("tau")    = tau,
    Rcpp::Named("Delta")  = Delta,
    Rcpp::Named("Lambda") = Lambda,
    Rcpp::Named("Eta")    = Eta,
    Rcpp::Named("sigma2") = sigma2
  );
}

#endif
