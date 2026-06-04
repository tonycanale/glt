#include <RcppArmadillo.h>
#include "update_alpha.h"
#include "utils_aux.h"

// [[Rcpp::depends(RcppArmadillo)]]

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
  const int H = static_cast<int>(ell.size());

  // Sufficient statistics: d_j = ones below pivot, n_j = free positions
  arma::ivec d(H);
  arma::ivec n(H);

  for (int j = 0; j < H; ++j) {
    const int ell_j = ell[j];   // 1-based
    n[j] = m - ell_j;
    int dj = 0;
    for (int i = ell_j; i < m; ++i) dj += Delta(i, j);
    d[j] = dj;
  }

  // Log-normal random walk proposal
  const double log_alpha_cur  = std::log(alpha);
  const double log_alpha_prop = log_alpha_cur + R::rnorm(0.0, mh_sd_alpha);
  const double alpha_prop     = std::exp(log_alpha_prop);

  // MH acceptance ratio + log-Jacobian correction; log_post_alpha from utils_aux.h
  const double log_acc =
    log_post_alpha(alpha_prop, d, n, H, beta, a_alpha, b_alpha) -
    log_post_alpha(alpha,      d, n, H, beta, a_alpha, b_alpha) +
    log_alpha_prop - log_alpha_cur;

  if (std::log(R::unif_rand()) < log_acc)
    return alpha_prop;

  return alpha;
}
