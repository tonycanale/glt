#include "truncnorm.h"

// [[Rcpp::depends(RcppArmadillo)]]

//' Sample from a truncated normal distribution
//'
//' Samples are drawn componentwise, so each component of the vector is
//' allowed its own mean, standard deviation, and upper and lower limits. The
//' components are assumed to be independent.
//'
//' @param y_lower `n x p` matrix of lower endpoints
//' @param y_upper `n x p` matrix of upper endpoints
//' @param mu `n x p` matrix of conditional expectations
//' @param sigma `p x 1` vector of conditional standard deviations
//' @param u_rand `n x p` matrix of uniform random variables
//'
//' @return An `n x p` matrix of draws from the truncated normal distribution.
//'
//' @note This function uses Rcpp for computational efficiency. Bounds may be
//'   `-Inf`/`Inf`; these are handled correctly by `R::pnorm`/`R::qnorm`.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat truncnorm_lg(const arma::mat& y_lower, const arma::mat& y_upper,
                       const arma::mat& mu, const arma::vec& sigma,
                       const arma::mat& u_rand) {
  // Dim of matrix:
  int n = y_lower.n_rows;
  int p = y_lower.n_cols;
  // Storage:
  double val = 0;
  arma::mat z_star(n,p);

  for(int t = 0; t < n; ++t) {
    for(int j = 0; j < p; ++j) {
      // Control
      double uptail1 = (y_lower(t,j) - mu(t,j)) * 1 / sigma(j) > 8;
      double uptail2 = (y_upper(t,j) - mu(t,j)) * 1 / sigma(j) > 8;
      // pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
      // true + false = 1, true + true = 2, false + false = 0
      if((uptail1 + uptail2) == 0){
        // Lower and upper limits, transformed via pnorm:
        double F_lower = R::pnorm(y_lower(t,j), mu(t,j), sigma(j), 1, 0);
        double F_upper = R::pnorm(y_upper(t,j), mu(t,j), sigma(j), 1, 0);
        // replace 0 with 0.000001 and 1 with 0.999999
        if (F_lower == 0) {
          F_lower = 0.000001;
        }
        if (F_upper == 0) {
          F_upper = 0.000001;
        } else if (F_upper == 1) {
          F_upper = 0.999999;
        }
        // Corresponding sampled value:
        val = R::qnorm(F_lower + u_rand(t,j) * (F_upper - F_lower), mu(t,j), sigma(j), 1, 0);
      }
      else {
        double F_lower = R::pnorm(y_lower(t,j), mu(t,j), sigma(j), 0, 0);
        double F_upper = R::pnorm(y_upper(t,j), mu(t,j), sigma(j), 0, 0);
        // replace 0 with 0.000001 and 1 with 0.999999
        if (F_lower == 0) {
          F_lower = 0.000001;
        }
        if (F_upper == 0) {
          F_upper = 0.000001;
        } else if (F_upper == 1) {
          F_upper = 0.999999;
        }
        // Corresponding sampled value:
        val = R::qnorm(F_lower + u_rand(t,j) * (F_upper - F_lower), mu(t,j), sigma(j), 0, 0);
      }
      z_star(t,j) = std::min(std::max(y_lower(t,j), val), y_upper(t,j));
    }
  }
  return z_star;
}
