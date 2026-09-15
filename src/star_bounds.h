#ifndef GLTFACTOR_STAR_BOUNDS_H
#define GLTFACTOR_STAR_BOUNDS_H

#include <RcppArmadillo.h>
#include "star_transform.h"
// [[Rcpp::depends(RcppArmadillo)]]

// ---------------------------------------------------------------------------
// compute_star_bounds
//
// Precompute, once (before the MCMC loop starts), the T x m matrices of
// lower/upper bounds on the latent z implied by:
//
//   y_ij = ell   iff   G_j(z_ij) in A_{j,ell} = [a_{j,ell}, a_{j,ell+1})
//
// where {a_{j,0} = -Inf, a_{j,1}, ..., a_{j,K_j-1}, a_{j,K_j} = Inf} is the
// (column-specific) threshold sequence for column j, and G_j is a monotone
// increasing transform (see star_transform.h).
//
// y_obs is 0-indexed: level `ell` for column j corresponds to interior
// cutpoints thresholds[j][ell-1] (lower, when ell >= 1) and
// thresholds[j][ell] (upper, when ell < K_j - 1); the outermost bounds are
// -Inf / Inf.
//
// Arguments
// ---------
//   y_obs      : T x m integer matrix of observed levels, 0-indexed
//                (0, 1, ..., K_j - 1 for column j).
//   thresholds : list of length m; element j is a numeric vector of the
//                K_j - 1 interior cutpoints (a_{j,1}, ..., a_{j,K_j-1}),
//                strictly increasing. (-Inf/Inf bookends are implicit.)
//   g_type     : length-m integer vector of transform ids (see
//                star_transform.h), one per column.
//
// Returns
// -------
//   Rcpp::List with elements `z_lower` and `z_upper`, each a T x m
//   arma::mat (entries may be -Inf/Inf).
// ---------------------------------------------------------------------------
inline Rcpp::List compute_star_bounds(
    const arma::imat&        y_obs,
    const Rcpp::List&        thresholds,
    const Rcpp::IntegerVector& g_type
) {
  const int T = y_obs.n_rows;
  const int m = y_obs.n_cols;

  if (thresholds.size() != m) {
    Rcpp::stop("`thresholds` must be a list of length m = ncol(y_obs).");
  }
  if (g_type.size() != m) {
    Rcpp::stop("`g_type` must have length m = ncol(y_obs).");
  }

  arma::mat z_lower(T, m);
  arma::mat z_upper(T, m);

  for (int j = 0; j < m; ++j) {
    const arma::vec a_j = Rcpp::as<arma::vec>(thresholds[j]); // interior cutpoints, length K_j - 1
    const int Kj_minus_1 = a_j.n_elem;                        // number of interior cutpoints
    const int transform_id = g_type[j];

    for (int t = 0; t < T; ++t) {
      const int ell = y_obs(t, j);
      if (ell < 0 || ell > Kj_minus_1) {
        Rcpp::stop("`y_obs(%d, %d)` = %d is out of range for the supplied thresholds "
                   "(valid levels are 0..%d).", t + 1, j + 1, ell, Kj_minus_1);
      }

      const double a_lower = (ell == 0)            ? -arma::datum::inf : a_j(ell - 1);
      const double a_upper = (ell == Kj_minus_1)    ?  arma::datum::inf : a_j(ell);

      z_lower(t, j) = star_invert_G(transform_id, a_lower);
      z_upper(t, j) = star_invert_G(transform_id, a_upper);
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("z_lower") = z_lower,
    Rcpp::Named("z_upper") = z_upper
  );
}

#endif // GLTFACTOR_STAR_BOUNDS_H
