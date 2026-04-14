
#ifndef GLTFACTOR_MARGLIK_H
#define GLTFACTOR_MARGLIK_H

#include <RcppArmadillo.h>

namespace {

double marg_lik_i(
    int i,
    const arma::mat& y,
    const arma::imat& Delta,
    const arma::mat& Eta,
    double kappa,
    double a_s,
    double b_s,
    bool logarithm = true
) {
  const int T  = y.n_rows;
  const arma::vec yi = y.col(i);

  // Active factor indices for row i
  arma::uvec active = arma::find(
      arma::conv_to<arma::uvec>::from(Delta.row(i)));
  const int qi = static_cast<int>(active.n_elem);

  double res;

  if (qi == 0) {
    // ------------------------------------------------------------------
    // Marginal: sigma2_i ~ InvGamma  (no factor contribution)
    // ------------------------------------------------------------------
    const double Qi = 0.5 * arma::dot(yi, yi);
    res = -T / 2.0 * std::log(2.0 * M_PI)
        + R::lgammafn(a_s + T / 2.0) - R::lgammafn(a_s)
        + a_s * std::log(b_s)
        - (a_s + T / 2.0) * std::log(b_s + Qi);
  } else {
    // ------------------------------------------------------------------
    // qi > 0:  Omega = F F' + (1/kappa) I_qi   (qi x qi)
    //
    // Use a SINGLE Cholesky to obtain:
    //   (a) log|Omega|  = 2 * sum(log diag(L))
    //   (b) m_post      = Omega^{-1} F y_i   via two triangular solves
    //   (c) Qi          = 0.5 * (yi'yi - m_post' Omega m_post)
    //                   = 0.5 * (yi'yi - Fyi' m_post)   since Omega m_post = Fyi
    // ------------------------------------------------------------------
    const arma::mat F   = Eta.rows(active);          // qi x T
    const arma::mat Omega = F * F.t() +
        arma::eye(qi, qi) / kappa;                   // qi x qi

    // Cholesky  Omega = L L'
    const arma::mat L   = arma::chol(arma::symmatu(Omega), "lower");

    // log|Omega| = 2 * sum(log diag(L))
    const double log_det_Omega = 2.0 * arma::sum(arma::log(L.diag()));

    // m_post = Omega^{-1} (F yi)  via L L' m = F yi
    const arma::vec Fyi    = F * yi;
    const arma::vec w      = arma::solve(arma::trimatl(L), Fyi,
                                         arma::solve_opts::fast);
    const arma::vec m_post = arma::solve(arma::trimatu(L.t()), w,
                                         arma::solve_opts::fast);

    const double Qi = 0.5 * (arma::dot(yi, yi) - arma::dot(Fyi, m_post));

    res = -T / 2.0 * std::log(2.0 * M_PI)
        + R::lgammafn(a_s + T / 2.0) - R::lgammafn(a_s)
        + a_s * std::log(b_s)
        - (a_s + T / 2.0) * std::log(b_s + Qi)
        - (qi / 2.0) * std::log(kappa)
        - 0.5 * log_det_Omega;
  }
  return logarithm ? res : std::exp(res);
}

} // anonymous namespace

inline Rcpp::NumericVector marg_lik(
    const Rcpp::IntegerVector& index,
    const arma::mat& y,
    const arma::imat& Delta,
    const arma::mat& Eta,
    const Rcpp::List& hyperpar,
    bool logarithm = true
) {
  double kappa = Rcpp::as<double>(hyperpar["kappa"]);
  double a_s   = Rcpp::as<double>(hyperpar["a_s"]);
  double b_s   = Rcpp::as<double>(hyperpar["b_s"]);
  Rcpp::NumericVector result(index.size());
  for (int k = 0; k < index.size(); ++k) {
    result[k] = marg_lik_i(index[k] - 1, y, Delta, Eta,
                            kappa, a_s, b_s, logarithm);
  }
  return result;
}

inline int rbetabinom_cpp(double a, double b, int size) {
  arma::vec log_pmf(size + 1);
  for (int k = 0; k <= size; ++k) {
    double log_choose   = R::lchoose(size, k);
    double log_beta_num = R::lbeta(a + k, b + size - k);
    double log_beta_den = R::lbeta(a, b);
    log_pmf[k] = log_choose + log_beta_num - log_beta_den;
  }
  double max_log = log_pmf.max();
  arma::vec pmf  = arma::exp(log_pmf - max_log);
  pmf /= arma::sum(pmf);
  double u = R::unif_rand();
  double cumsum = 0.0;
  for (int k = 0; k <= size; ++k) {
    cumsum += pmf[k];
    if (u < cumsum) return k;
  }
  return size;
}

#endif

