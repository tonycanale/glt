#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

// ============================================================================
// HELPER UTILITIES for parsing and state management
// ============================================================================

//' Parse an integer scalar from a named list entry
int get_int(const Rcpp::List& x, const char* name) {
  SEXP obj = x[name];
  if (obj == R_NilValue) {
    Rcpp::stop("Missing entry `%s`.", name);
  }
  return Rcpp::as<int>(obj);
}

//' Parse a double scalar from a named list entry
double get_double(const Rcpp::List& x, const char* name) {
  SEXP obj = x[name];
  if (obj == R_NilValue) {
    Rcpp::stop("Missing entry `%s`.", name);
  }
  return Rcpp::as<double>(obj);
}

//' Parse a bool scalar from a named list entry
bool get_bool(const Rcpp::List& x, const char* name) {
  SEXP obj = x[name];
  if (obj == R_NilValue) {
    Rcpp::stop("Missing entry `%s`.", name);
  }
  return Rcpp::as<bool>(obj);
}

//' Compute marginal likelihood for a single row i of Delta
//' 
//' Matches OLDSTUFF/aux.R::marg.lik.i
//' 
//' @keywords internal
double marg_lik_i(
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

// ============================================================================
// EXPORTED FUNCTIONS (called from gltfa_cpp)
// ============================================================================

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
// [[Rcpp::export]]
Rcpp::NumericVector marg_lik(
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
//' Matches OLDSTUFF/aux.R::rbetabinom
//' 
//' @param a Shape parameter 1
//' @param b Shape parameter 2
//' @param size Number of trials
//' 
//' @return Integer sample from beta-binomial(size, a, b)
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
int rbetabinom_cpp(double a, double b, int size) {
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