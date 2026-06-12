#include <RcppArmadillo.h>
#include "update_a_nu.h"
#include "utils_aux.h"

// [[Rcpp::depends(RcppArmadillo)]]

//' @keywords internal
//' @noRd
double update_a_nu_cpp(
    double a_nu,
    int    H,
    int    m,
    double b_nu,
    double a_anu,
    double b_anu,
    double mh_sd_a_nu
) {
  const double log_a_nu_cur  = std::log(a_nu);
  const double log_a_nu_prop = log_a_nu_cur + R::rnorm(0.0, mh_sd_a_nu);
  const double a_nu_prop     = std::exp(log_a_nu_prop);

  // log_post_a_nu from utils_aux.h; log-Jacobian for log-normal proposal
  const double log_acc =
    log_post_a_nu(a_nu_prop, H, m, b_nu, a_anu, b_anu) -
    log_post_a_nu(a_nu,      H, m, b_nu, a_anu, b_anu) +
    log_a_nu_prop - log_a_nu_cur;

  if (std::log(R::unif_rand()) < log_acc)
    return a_nu_prop;

  return a_nu;
}
