// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "utils_aux.h"
using namespace Rcpp;

// ------------------------------------------------------------------ //
// Exported: count spurious columns (colSum == 1) for every draw
//
// Thin wrapper around the inline compute_spurious() helper in
// utils_aux.h, which is also reused directly (no .Call() round trip) by
// other C++ translation units.
// ------------------------------------------------------------------ //
// [[Rcpp::export]]
IntegerVector compute_spurious_cpp(List Delta_draws) {
  return compute_spurious(Delta_draws);
}
