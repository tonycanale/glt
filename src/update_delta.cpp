#include <RcppArmadillo.h>
#include "helpers.h"      // pivots_from_delta_1based, log_post_odds_ij
#include "update_delta.h"

// [[Rcpp::depends(RcppArmadillo)]]

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::imat update_delta_cpp(
    const arma::mat&   y,
    const arma::imat&  Delta_in,
    const arma::mat&   Eta,
    const arma::vec&   tau,
    const Rcpp::List&  hyperpar
) {
  Rcpp::RNGScope rng_scope;

  const int m = static_cast<int>(Delta_in.n_rows);
  const int H = static_cast<int>(Delta_in.n_cols);

  arma::imat Delta = Delta_in;

  // Random column order; Rcpp::sample(H, H, false) returns values in {1,...,H}
  Rcpp::IntegerVector jind = Rcpp::sample(H, H, false);

  for (int k = 0; k < H; ++k) {
    // convert to 0-based column index
    const int j = jind[k] - 1;

    // extract current pivots (1-based) from the current Delta
    std::vector<int> ell = pivots_from_delta_1based(Delta);

    // ell[j] is 1-based; rows strictly below pivot are ell[j]+1,...,m (1-based)
    // if pivot is already at the last row there is nothing to update
    if (ell[j] >= m) continue;

    const int n_free = m - ell[j];   // number of free rows below pivot
    std::vector<bool> accepted(n_free, false);

    // -----------------------------------------------------------------------
    // Pass 1: decide acceptance for each free row
    // -----------------------------------------------------------------------
    for (int h = 0; h < n_free; ++h) {
      // i_1based runs from ell[j]+1 to m  (1-based)
      const int i_1based = ell[j] + 1 + h;

      const double log_u = std::log(R::unif_rand());
      const double Oij   = log_post_odds_ij(
          i_1based, j, y, Delta, Eta, tau, hyperpar);

      // remember ell_j are {1, ..., m} but C indexing starts at 0
      const int cur = Delta(i_1based - 1, j);

      if (cur == 0) {
        accepted[h] = (log_u <= Oij);
      } else {
        accepted[h] = (log_u <= -Oij);
      }
    }

    // -----------------------------------------------------------------------
    // Pass 2: apply accepted flips
    // -----------------------------------------------------------------------
    for (int h = 0; h < n_free; ++h) {
      if (!accepted[h]) continue;

      const int i_1based = ell[j] + 1 + h;

      // remember ell_j are {1, ..., m} but C indexing starts at 0
      const int cur = Delta(i_1based - 1, j);
      Delta(i_1based - 1, j) = 1 - cur;
    }
    // Delta is updated for column j before moving to the next column
  }

  return Delta;
}
