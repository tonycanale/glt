#include <RcppArmadillo.h>
#include "helpers.h"          // pivots_from_delta_1based, set_pivot_diagonal_from_ell_1based, changed_rows_1based
#include "marglik.h"
#include "update.pivots.h"

// [[Rcpp::depends(RcppArmadillo)]]

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List update_pivots_cpp(
    const arma::mat& y,
    const arma::imat& Delta_in,
    const arma::mat& Eta,
    const Rcpp::List& hyperpar
) {
  Rcpp::RNGScope rng;

  const int m = static_cast<int>(Delta_in.n_rows);
  const int H = static_cast<int>(Delta_in.n_cols);

  const double paralpha = Rcpp::as<double>(hyperpar["alpha"]);
  const double parbeta  = Rcpp::as<double>(hyperpar["beta"]);

  arma::imat Delta = Delta_in;
  std::vector<int> ell = pivots_from_delta_1based(Delta);

  // Rcpp::sample(H, H, false) returns values in {1,...,H}
  Rcpp::IntegerVector jstar = Rcpp::sample(H, H, false);

  for (int k = 0; k < H; ++k) {
    // for column 
    // convert to 0-based column index
    const int j_prop = jstar[k] - 1;

    const int lower = (j_prop == 0)     ? 1 : (ell[j_prop - 1] + 1);
    const int upper = (j_prop == H - 1) ? m : (ell[j_prop + 1] - 1);

    if (lower >= upper) continue;

    Rcpp::IntegerVector range        = Rcpp::seq(lower, upper);
    Rcpp::IntegerVector elljstar_draw = Rcpp::sample(range, 1, false);
    const int elljstar = elljstar_draw[0];

    if (elljstar == ell[j_prop]) continue;

    std::vector<int> ellstar = ell;
    ellstar[j_prop] = elljstar;

    arma::imat Deltastar = Delta;
    set_pivot_diagonal_from_ell_1based(Deltastar, ellstar);

    std::vector<int> diffs(H, 0);
    for (int j = 0; j < H; ++j) diffs[j] = ell[j] - ellstar[j];

    for (int j_move = 0; j_move < H; ++j_move) {
      if (diffs[j_move] == 0) continue;

      if (diffs[j_move] < 0) {
        // ell*_j > ell_j: entries from ell_j to ell*_j - 1 become structural zeros
        // remember ell_j are {1, ..., m} but C indexing starts at 0
        for (int i = ell[j_move] - 1; i <= ellstar[j_move] - 2; ++i)
          Deltastar(i, j_move) = 0;

      } else {
        const int interval_size = ell[j_move] - ellstar[j_move];
        if (interval_size <= 0) continue;

        // Clear candidate entries between ell*_j + 1 and ell_j
        // remember ell_j are {1, ..., m} but C indexing starts at 0
        for (int i = ellstar[j_move]; i <= ell[j_move] - 1; ++i)
          Deltastar(i, j_move) = 0;

        int dj = 0;
        const int start_dj = std::min(ell[j_move] + 1, m);
        // sum Delta[(min(ell_j+1,m)):m, j] — 1-based interval, hence -1 offset
        // remember ell_j are {1, ..., m} but C indexing starts at 0
        for (int i = start_dj - 1; i < m; ++i)
          dj += Delta(i, j_move);

        const int d_a = rbetabinom_cpp(
            paralpha * parbeta + dj,
            parbeta + m - ell[j_move] - dj,
            interval_size);

        if (d_a > 0) {
          Rcpp::IntegerVector candidate_rows = Rcpp::seq(ellstar[j_move] + 1, ell[j_move]);
          Rcpp::IntegerVector new_ones = Rcpp::sample(candidate_rows, d_a, false);
          for (int h = 0; h < new_ones.size(); ++h) {
            // candidate rows are 1-based
            // remember ell_j are {1, ..., m} but C indexing starts at 0
            Deltastar(new_ones[h] - 1, j_move) = 1;
          }
        }
      }
    }

    Rcpp::IntegerVector index = changed_rows_1based(Delta, Deltastar);
    if (index.size() == 0) continue;

    Rcpp::NumericVector new_ll = marg_lik(index, y, Deltastar, Eta, hyperpar);
    Rcpp::NumericVector old_ll = marg_lik(index, y, Delta,     Eta, hyperpar);

    double log_lik_R = 0.0;
    for (int h = 0; h < index.size(); ++h)
      log_lik_R += new_ll[h] - old_ll[h];

    const bool accept = (R::unif_rand() < std::min(1.0, std::exp(log_lik_R)));
    if (accept) {
      Delta = Deltastar;
      ell   = ellstar;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("Delta") = Delta,
    Rcpp::Named("Eta")   = Eta,
    Rcpp::Named("ell")   = Rcpp::wrap(ell)
  );
}
