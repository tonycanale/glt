#ifndef GLTFACTOR_STAR_TRANSFORM_H
#define GLTFACTOR_STAR_TRANSFORM_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// ---------------------------------------------------------------------------
// Column-specific monotone transformation G_j used by the STAR link
//   y_ij = ell   iff   G_j(z_ij) in [a_{j,ell}, a_{j,ell+1})
//
// Currently supported transforms:
//   IDENTITY : G(t) = t
//   EXP      : G(t) = exp(t)
//
// The dispatch is written generically (integer transform id per column) so
// that additional monotone transforms (log1p, Box-Cox, ...) can be added by
// extending the enum and the two switch statements below.
// ---------------------------------------------------------------------------

enum star_transform_id {
  STAR_TRANSFORM_IDENTITY = 0,
  STAR_TRANSFORM_EXP      = 1
};

// Map a user-facing string to the internal transform id.
inline int star_transform_id_from_string(const std::string& g_type) {
  if (g_type == "identity") return STAR_TRANSFORM_IDENTITY;
  if (g_type == "exp")      return STAR_TRANSFORM_EXP;
  Rcpp::stop("Unknown `g_type` value: '%s'. Supported values are 'identity' and 'exp'.",
             g_type.c_str());
}

// G_j(t) forward transform.
inline double star_apply_G(int transform_id, double t) {
  switch (transform_id) {
    case STAR_TRANSFORM_IDENTITY:
      return t;
    case STAR_TRANSFORM_EXP:
      return std::exp(t);
    default:
      Rcpp::stop("Unsupported transform id in star_apply_G().");
  }
}

// G_j^{-1}(a) inverse transform. `a` may be +/-Inf; the inverse must map
// these consistently (e.g. log(0) = -Inf, log(Inf) = Inf) so that thresholds
// at the boundary of the support translate into infinite bounds on z.
inline double star_invert_G(int transform_id, double a) {
  switch (transform_id) {
    case STAR_TRANSFORM_IDENTITY:
      return a;
    case STAR_TRANSFORM_EXP:
      if (a <= 0.0) return -arma::datum::inf;
      return std::log(a);
    default:
      Rcpp::stop("Unsupported transform id in star_invert_G().");
  }
}

#endif // GLTFACTOR_STAR_TRANSFORM_H
