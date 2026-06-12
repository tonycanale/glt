#include <RcppArmadillo.h>
#include "update_alpha_beta.h"
#include "utils_aux.h"

// [[Rcpp::depends(RcppArmadillo)]]

// ---------------------------------------------------------------------------
// Shared helper: compute sufficient statistics d_j and n_j from Delta / ell
// ---------------------------------------------------------------------------
static void compute_dn(
    const arma::imat&          Delta,
    const Rcpp::IntegerVector& ell,
    int                        m,
    arma::ivec&                d,
    arma::ivec&                n
) {
  const int H = static_cast<int>(ell.size());
  d.set_size(H);
  n.set_size(H);
  for (int j = 0; j < H; ++j) {
    const int ell_j = ell[j];   // 1-based
    n[j] = m - ell_j;
    int dj = 0;
    for (int i = ell_j; i < m; ++i) dj += Delta(i, j);
    d[j] = dj;
  }
}

// ---------------------------------------------------------------------------
// update_alpha_cpp  –  MH step for alpha, beta fixed
// ---------------------------------------------------------------------------

//' @keywords internal
//' @noRd
double update_alpha_cpp(
    double                        alpha,
    const arma::imat&             Delta,
    const Rcpp::IntegerVector&    ell,
    int                           m,
    double                        beta,
    double                        a_alpha,
    double                        b_alpha,
    double                        mh_sd_alpha
) {
  arma::ivec d, n;
  compute_dn(Delta, ell, m, d, n);
  const int H = static_cast<int>(ell.size());

  const double log_alpha_cur  = std::log(alpha);
  const double log_alpha_prop = log_alpha_cur + R::rnorm(0.0, mh_sd_alpha);
  const double alpha_prop     = std::exp(log_alpha_prop);

  // log_post_alpha from utils_aux.h; log-Jacobian accounts for log-normal proposal
  const double log_acc =
    log_post_alpha(alpha_prop, d, n, H, beta, a_alpha, b_alpha) -
    log_post_alpha(alpha,      d, n, H, beta, a_alpha, b_alpha) +
    log_alpha_prop - log_alpha_cur;

  if (std::log(R::unif_rand()) < log_acc)
    return alpha_prop;

  return alpha;
}

// ---------------------------------------------------------------------------
// update_beta_cpp  –  MH step for beta, alpha fixed
// ---------------------------------------------------------------------------

//' @keywords internal
//' @noRd
double update_beta_cpp(
    double                        beta,
    const arma::imat&             Delta,
    const Rcpp::IntegerVector&    ell,
    int                           m,
    double                        alpha,
    double                        a_beta,
    double                        b_beta,
    double                        mh_sd_beta
) {
  arma::ivec d, n;
  compute_dn(Delta, ell, m, d, n);
  const int H = static_cast<int>(ell.size());

  const double log_beta_cur  = std::log(beta);
  const double log_beta_prop = log_beta_cur + R::rnorm(0.0, mh_sd_beta);
  const double beta_prop     = std::exp(log_beta_prop);

  // log_post_beta from utils_aux.h; log-Jacobian accounts for log-normal proposal
  const double log_acc =
    log_post_beta(beta_prop, d, n, H, alpha, a_beta, b_beta) -
    log_post_beta(beta,      d, n, H, alpha, a_beta, b_beta) +
    log_beta_prop - log_beta_cur;

  if (std::log(R::unif_rand()) < log_acc)
    return beta_prop;

  return beta;
}
