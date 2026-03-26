#ifndef GLTFACTOR_HELPERS_H
#define GLTFACTOR_HELPERS_H

#include <RcppArmadillo.h>

// Parse an integer scalar from a named list entry
inline int get_int(const Rcpp::List& x, const char* name) {
  SEXP obj = x[name];
  if (obj == R_NilValue) {
    Rcpp::stop("Missing entry `%s`.", name);
  }
  return Rcpp::as<int>(obj);
}

// Parse a double scalar from a named list entry
inline double get_double(const Rcpp::List& x, const char* name) {
  SEXP obj = x[name];
  if (obj == R_NilValue) {
    Rcpp::stop("Missing entry `%s`.", name);
  }
  return Rcpp::as<double>(obj);
}

// Parse a bool scalar from a named list entry
inline bool get_bool(const Rcpp::List& x, const char* name) {
  SEXP obj = x[name];
  if (obj == R_NilValue) {
    Rcpp::stop("Missing entry `%s`.", name);
  }
  return Rcpp::as<bool>(obj);
}

// Build a simple default state when H = 0 or when init pieces are NULL
inline Rcpp::List make_last_state(
    int H,
    int T,
    int m,
    SEXP nu_,
    SEXP ell_,
    SEXP tau_,
    SEXP Delta_,
    SEXP Lambda_,
    SEXP Eta_,
    SEXP sigma2_) {

  Rcpp::IntegerVector ell;
  Rcpp::NumericVector tau;
  arma::imat Delta;
  arma::mat Lambda;
  arma::mat Eta;
  Rcpp::NumericVector sigma2;

  if (ell_ == R_NilValue) {
    ell = Rcpp::IntegerVector(H);
  } else {
    ell = Rcpp::as<Rcpp::IntegerVector>(ell_);
  }

  if (tau_ == R_NilValue) {
    tau = Rcpp::NumericVector(H);
  } else {
    tau = Rcpp::as<Rcpp::NumericVector>(tau_);
  }

  if (Delta_ == R_NilValue) {
    Delta = arma::imat(m, H, arma::fill::zeros);
  } else {
    Delta = Rcpp::as<arma::imat>(Delta_);
  }

  if (Lambda_ == R_NilValue) {
    Lambda = arma::mat(m, H, arma::fill::zeros);
  } else {
    Lambda = Rcpp::as<arma::mat>(Lambda_);
  }

  if (Eta_ == R_NilValue) {
    Eta = arma::mat(H, T, arma::fill::zeros);
  } else {
    Eta = Rcpp::as<arma::mat>(Eta_);
  }

  if (sigma2_ == R_NilValue) {
    sigma2 = Rcpp::NumericVector(m, 1.0);
  } else {
    sigma2 = Rcpp::as<Rcpp::NumericVector>(sigma2_);
  }

  double nu = (nu_ == R_NilValue) ? NA_REAL : Rcpp::as<double>(nu_);

  return Rcpp::List::create(
    Rcpp::Named("H") = H,
    Rcpp::Named("nu") = nu,
    Rcpp::Named("ell") = ell,
    Rcpp::Named("tau") = tau,
    Rcpp::Named("Delta") = Delta,
    Rcpp::Named("Lambda") = Lambda,
    Rcpp::Named("Eta") = Eta,
    Rcpp::Named("sigma2") = sigma2
  );
}

#endif
