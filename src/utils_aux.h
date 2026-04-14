#ifndef GLTFACTOR_AUX_H
#define GLTFACTOR_AUX_H

#include <RcppArmadillo.h>

namespace {

// ============================================================================
// MARGINAL LIKELIHOOD FUNCTIONS (from utils.cpp/marglik.h)
// ============================================================================

//' Compute marginal likelihood for a single row i of Delta
//' 
//' @keywords internal
inline double marg_lik_i(
    int i,
    const arma::mat& y,              // T x m
    const arma::imat& Delta,         // m x H
    const arma::mat& Eta,            // H x T
    double kappa,
    double a_s,
    double b_s,
    bool logarithm = true
) {
  const int T = y.n_rows;
  const arma::vec yi = y.col(i);
  const arma::irowvec Deltai = Delta.row(i);
  
  int qi = 0;
  for (int j = 0; j < (int)Deltai.n_elem; ++j) {
    if (Deltai[j] == 1) qi++;
  }
  
  double res;
  
  if (qi > 0) {
    // Extract rows of Eta where Delta[i, j] == 1
    arma::mat eta(T, qi);
    int col_idx = 0;
    for (int j = 0; j < (int)Deltai.n_elem; ++j) {
      if (Deltai[j] == 1) {
        eta.col(col_idx) = Eta.row(j).t();
        col_idx++;
      }
    }
    
    // Compute inner.mat = diag(1/kappa, qi) + eta^T eta
    arma::mat inner_mat = eta.t() * eta;
    for (int j = 0; j < qi; ++j) {
      inner_mat(j, j) += 1.0 / kappa;
    }
    
    // Compute determinant and inverse
    double det_sign;
    double det_val;
    arma::log_det(det_val, det_sign, inner_mat);
    arma::mat inner_mat_inv = arma::inv(inner_mat);
    
    // Compute Q_i = 0.5 * (yi^T yi - yi^T eta (inner_mat)^{-1} eta^T yi)
    double Qi = 0.5 * as_scalar(
      yi.t() * yi - 
      yi.t() * eta * inner_mat_inv * eta.t() * yi
    );
    
    res = -T / 2.0 * std::log(2.0 * M_PI) + 
      R::lgammafn(a_s + T / 2.0) - R::lgammafn(a_s) + 
      a_s * std::log(b_s) - 
      (a_s + T / 2.0) * std::log(b_s + Qi) - 
      (qi / 2.0) * std::log(kappa) - 
      0.5 * det_val;
  } else {
    // qi == 0: marginal likelihood when no factors active for row i
    double sum_yi2 = 0.0;
    for (int t = 0; t < T; ++t) {
      sum_yi2 += yi[t] * yi[t];
    }
    
    res = -T / 2.0 * std::log(2.0 * M_PI) + 
      R::lgammafn(a_s + T / 2.0) - R::lgammafn(a_s) + 
      a_s * std::log(b_s) - 
      (a_s + T / 2.0) * std::log(b_s + 0.5 * sum_yi2);
  }
  
  if (logarithm) {
    return res;
  } else {
    return std::exp(res);
  }
}

} // anonymous namespace

//' Marginal likelihood for multiple rows
//' 
//' Compute marginal likelihood for specified rows of Delta.
//' Matches OLDSTUFF/aux.R::marg.lik
//' 
//' @param index 1-based indices of rows to compute marginal likelihood for
//' @param y Response matrix (T x m)
//' @param Delta Binary indicator matrix (m x H)
//' @param Eta Factor matrix (H x T)
//' @param hyperpar List with elements: kappa, a_s, b_s
//' @param logarithm If TRUE, return log marginal likelihood; else exponentiate
//' 
//' @return Numeric vector of marginal likelihoods
//' 
//' @keywords internal
//' @noRd
inline Rcpp::NumericVector marg_lik(
    const Rcpp::IntegerVector& index,
    const arma::mat& y,
    const arma::imat& Delta,
    const arma::mat& Eta,
    const Rcpp::List& hyperpar,
    bool logarithm = true
) {
  double kappa = Rcpp::as<double>(hyperpar["kappa"]);
  double a_s = Rcpp::as<double>(hyperpar["a_s"]);
  double b_s = Rcpp::as<double>(hyperpar["b_s"]);
  
  Rcpp::NumericVector result(index.size());
  
  for (int k = 0; k < index.size(); ++k) {
    int i = index[k] - 1;  // Convert to 0-based
    result[k] = marg_lik_i(i, y, Delta, Eta, kappa, a_s, b_s, logarithm);
  }
  
  return result;
}

//' Sample from beta-binomial distribution
//' 
//' @param a Shape parameter 1
//' @param b Shape parameter 2
//' @param size Number of trials
//' 
//' @return Integer sample from beta-binomial(size, a, b)
//' 
//' @keywords internal
//' @noRd
inline int rbetabinom_cpp(double a, double b, int size) {
  // Compute pmf: P(X = k) = C(size, k) * B(a + k, b + size - k) / B(a, b)
  // Use log scale to avoid overflow
  
  arma::vec log_pmf(size + 1);
  
  for (int k = 0; k <= size; ++k) {
    // log C(size, k)
    double log_choose = R::lchoose(size, k);
    
    // log B(a + k, b + size - k)
    double log_beta_num = R::lbeta(a + k, b + size - k);
    
    // log B(a, b)
    double log_beta_den = R::lbeta(a, b);
    
    log_pmf[k] = log_choose + log_beta_num - log_beta_den;
  }
  
  // Normalize to get probabilities
  double max_log = log_pmf.max();
  arma::vec pmf = arma::exp(log_pmf - max_log);
  pmf /= arma::sum(pmf);
  
  // Sample from discrete distribution
  double u = R::unif_rand();
  double cumsum = 0.0;
  for (int k = 0; k <= size; ++k) {
    cumsum += pmf[k];
    if (u < cumsum) return k;
  }
  
  return size;  // Fallback (should not reach here)
}

// ============================================================================
// PARAMETER EXTRACTORS (from helpers.h)
// ============================================================================

//' Extract integer from named list
inline int get_int(const Rcpp::List& x, const char* name) {
  SEXP obj = x[name];
  if (obj == R_NilValue) Rcpp::stop("Missing entry `%s`.", name);
  return Rcpp::as<int>(obj);
}

//' Extract double from named list
inline double get_double(const Rcpp::List& x, const char* name) {
  SEXP obj = x[name];
  if (obj == R_NilValue) Rcpp::stop("Missing entry `%s`.", name);
  return Rcpp::as<double>(obj);
}

//' Extract boolean from named list
inline bool get_bool(const Rcpp::List& x, const char* name) {
  SEXP obj = x[name];
  if (obj == R_NilValue) Rcpp::stop("Missing entry `%s`.", name);
  return Rcpp::as<bool>(obj);
}

// ============================================================================
// DELTA MATRIX UTILITIES 
// ============================================================================

//' Extract 1-based pivots from Delta
//' For each column j, ell[j] is the first row in {1,...,m} where Delta == 1.
//' Used by update_pivots_cpp and update_delta_cpp.
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

//' Set pivot diagonal entries in Delta from a 1-based ell vector
//' For each column j, sets Delta(ell[j]-1, j) = 1.
//' Used by update_pivots_cpp.
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

//' Return 1-based row indices where matrices A and B differ
//' Used by update_pivots_cpp to identify which rows changed for marg_lik().
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

// ============================================================================
// MODEL EVALUATION FUNCTIONS 
// ============================================================================

//' Log posterior odds for entry (i_1based, j_0based) of Delta
//'
//' log p(delta_ij=1|rest) - log p(delta_ij=0|rest)
//'   = log(tau_j / (1-tau_j))
//'     + marg_lik(i | Delta_one) - marg_lik(i | Delta_zero)
//'
//' i_1based : 1-based row index (R convention)
//' j_0based : 0-based column index (C++ convention)
//' Used by update_delta_cpp.
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
  // marg_lik returns NumericVector; [0] gives a double directly — no Rcpp::as<> needed
  const double ll_one  = marg_lik(idx, y, Delta_one,  Eta, hyperpar)[0];
  const double ll_zero = marg_lik(idx, y, Delta_zero, Eta, hyperpar)[0];

  const double log_prior_odds =
      std::log(tau[j_0based]) - std::log(1.0 - tau[j_0based]);
  return (ll_one - ll_zero) + log_prior_odds;
}

// ============================================================================
// STATE INITIALIZATION HELPERS 
// ============================================================================

//' Build a simple default state when H = 0 or when init pieces are NULL
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
