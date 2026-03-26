
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
  const int T = y.n_rows;
  const arma::vec yi = y.col(i);
  const arma::irowvec Deltai = Delta.row(i);
  
  int qi = 0;
  for (int j = 0; j < (int)Deltai.n_elem; ++j) {
    if (Deltai[j] == 1) qi++;
  }
  
  double res;
  if (qi > 0) {
    arma::mat eta(T, qi);
    int col_idx = 0;
    for (int j = 0; j < (int)Deltai.n_elem; ++j) {
      if (Deltai[j] == 1) {
        eta.col(col_idx) = Eta.row(j).t();
        col_idx++;
      }
    }
    
    arma::mat inner_mat = eta.t() * eta;
    for (int j = 0; j < qi; ++j) {
      inner_mat(j, j) += 1.0 / kappa;
    }
    
    double det_val;
    arma::log_det(det_val, arma::inv(inner_mat));
    arma::mat inner_mat_inv = arma::inv(inner_mat);
    
    double Qi = 0.5 * as_scalar(yi.t() * yi - yi.t() * eta * inner_mat_inv * eta.t() * yi);
    
    res = -T / 2.0 * std::log(2.0 * M_PI) + 
      R::lgammafn(a_s + T / 2.0) - R::lgammafn(a_s) + 
      a_s * std::log(b_s) - 
      (a_s + T / 2.0) * std::log(b_s + Qi) - 
      (qi / 2.0) * std::log(kappa) - 0.5 * det_val;
  } else {
    double sum_yi2 = 0.0;
    for (int t = 0; t < T; ++t) {
      sum_yi2 += yi[t] * yi[t];
    }
    res = -T / 2.0 * std::log(2.0 * M_PI) + 
      R::lgammafn(a_s + T / 2.0) - R::lgammafn(a_s) + 
      a_s * std::log(b_s) - 
      (a_s + T / 2.0) * std::log(b_s + 0.5 * sum_yi2);
  }
  return logarithm ? res : std::exp(res);
}

}

inline Rcpp::NumericVector marg_lik(
    const Rcpp::IntegerVector& index,
    const arma::mat& y,
    const arma::imat& Delta,
    const arma::mat& Eta,
    const Rcpp::List& hyperpar,
    bool logarithm = true
) {
  double kappa = Rcpp::as<double>(hyperpar["kappa"]);
  double a_s = Rcpp::as<double>(hyperpar["a_s"]);
  double b_s = Rcpp::as<double>(hyperpar["b_s"]);
  Rcpp::NumericVector result(index.size());
  for (int k = 0; k < index.size(); ++k) {
    int i = index[k] - 1;
    result[k] = marg_lik_i(i, y, Delta, Eta, kappa, a_s, b_s, logarithm);
  }
  return result;
}

inline int rbetabinom_cpp(double a, double b, int size) {
  arma::vec log_pmf(size + 1);
  for (int k = 0; k <= size; ++k) {
    double log_choose = R::lchoose(size, k);
    double log_beta_num = R::lbeta(a + k, b + size - k);
    double log_beta_den = R::lbeta(a, b);
    log_pmf[k] = log_choose + log_beta_num - log_beta_den;
  }
  double max_log = log_pmf.max();
  arma::vec pmf = arma::exp(log_pmf - max_log);
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

