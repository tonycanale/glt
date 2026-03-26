#include <RcppArmadillo.h>
#include "helpers.h"
#include "marglik.h"
#include "update_H.h"

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

// Parse an integer scalar from a named list entry
int get_int(const Rcpp::List& x, const char* name) {
  SEXP obj = x[name];
  if (obj == R_NilValue) {
    Rcpp::stop("Missing entry `%s`.", name);
  }
  return Rcpp::as<int>(obj);
}

// Parse a double scalar from a named list entry
double get_double(const Rcpp::List& x, const char* name) {
  SEXP obj = x[name];
  if (obj == R_NilValue) {
    Rcpp::stop("Missing entry `%s`.", name);
  }
  return Rcpp::as<double>(obj);
}

// Parse a bool scalar from a named list entry
bool get_bool(const Rcpp::List& x, const char* name) {
  SEXP obj = x[name];
  if (obj == R_NilValue) {
    Rcpp::stop("Missing entry `%s`.", name);
  }
  return Rcpp::as<bool>(obj);
}

// Build a simple default state when H = 0 or when init pieces are NULL
Rcpp::List make_last_state(
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

} // anonymous namespace


//' Internal Rcpp sampler stub for gltfa
//' This is a placeholder implementation used to validate the R-to-C++
//' interface and package compilation. The full sampler will iterate the
//' Section 3 updates in the paper.
//' @noRd
// [[Rcpp::export]]
Rcpp::List gltfa_cpp(
    const arma::mat& y,
    const Rcpp::List& mcmc,
    const Rcpp::List& model,
    const Rcpp::List& prior,
    const Rcpp::List& init,
    const Rcpp::List& control
) {
  // ---------------------------------------------------------------------------
  // Dimensions and basic checks
  // ---------------------------------------------------------------------------
  const int T = get_int(model, "T");
  const int m = get_int(model, "m");
  const int Hmax = get_int(model, "Hmax");

  if (static_cast<int>(y.n_rows) != T || static_cast<int>(y.n_cols) != m) {
    Rcpp::stop("`y` has dimensions %d x %d, but `model` declares T = %d and m = %d.",
               y.n_rows, y.n_cols, T, m);
  }

  // ---------------------------------------------------------------------------
  // MCMC settings
  // ---------------------------------------------------------------------------
  const int niter = get_int(mcmc, "niter");
  const int nburn = get_int(mcmc, "nburn");
  const int thin  = get_int(mcmc, "thin");
  const int nsave = get_int(mcmc, "nsave");

  // ---------------------------------------------------------------------------
  // Prior hyperparameters
  // ---------------------------------------------------------------------------
  const double a_nu    = get_double(prior, "a_nu");
  const double b_nu    = get_double(prior, "b_nu");
  const double alpha   = get_double(prior, "alpha");
  const double beta    = get_double(prior, "beta");
  const double kappa   = get_double(prior, "kappa");
  const double a_sigma = get_double(prior, "a_sigma");
  const double b_sigma = get_double(prior, "b_sigma");

  // Avoid unused variable warnings in the stub.
  (void) a_nu;
  (void) b_nu;
  (void) alpha;
  (void) beta;
  (void) kappa;
  (void) a_sigma;
  (void) b_sigma;

  // ---------------------------------------------------------------------------
  // Control options
  // ---------------------------------------------------------------------------
  const bool store_draws = get_bool(control, "store_draws");
  const bool store_eta   = get_bool(control, "store_eta");
  const int print_every  = get_int(control, "print_every");
  const bool verbose     = get_bool(control, "verbose");

  (void) print_every;

  // Optional seed
  if (control.containsElementNamed("seed")) {
    SEXP seed_ = control["seed"];
    if (seed_ != R_NilValue) {
      const int seed = Rcpp::as<int>(seed_);
      Rcpp::Function set_seed("set.seed");
      set_seed(seed);
    }
  }

  // ---------------------------------------------------------------------------
  // Initial state
  // Conventions:
  //   Lambda : m x H
  //   Delta  : m x H
  //   Eta    : H x T
  //   y      : T x m
  // ---------------------------------------------------------------------------
  const int H = get_int(init, "H");

  if (H < 0) {
    Rcpp::stop("`init$H` must be nonnegative.");
  }
  if (H > Hmax) {
    Rcpp::stop("`init$H` cannot exceed `Hmax`.");
  }

  SEXP nu_     = init["nu"];
  SEXP ell_    = init["ell"];
  SEXP tau_    = init["tau"];
  SEXP Delta_  = init["Delta"];
  SEXP Lambda_ = init["Lambda"];
  SEXP Eta_    = init["Eta"];
  SEXP sigma2_ = init["sigma2"];

  Rcpp::List last = make_last_state(H, T, m, nu_, ell_, tau_, Delta_, Lambda_, Eta_, sigma2_);

  // Extra consistency checks on initialized matrices
  arma::imat Delta = Rcpp::as<arma::imat>(last["Delta"]);
  arma::mat Lambda = Rcpp::as<arma::mat>(last["Lambda"]);
  arma::mat Eta    = Rcpp::as<arma::mat>(last["Eta"]);

  if (static_cast<int>(Delta.n_rows) != m || static_cast<int>(Delta.n_cols) != H) {
    Rcpp::stop("`Delta` must have dimension m x H.");
  }
  if (static_cast<int>(Lambda.n_rows) != m || static_cast<int>(Lambda.n_cols) != H) {
    Rcpp::stop("`Lambda` must have dimension m x H.");
  }
  if (static_cast<int>(Eta.n_rows) != H || static_cast<int>(Eta.n_cols) != T) {
    Rcpp::stop("`Eta` must have dimension H x T.");
  }

  // ---------------------------------------------------------------------------
  // Allocate storage for draws
  // Because H is moving, most draws are stored as lists.
  // ---------------------------------------------------------------------------
  Rcpp::List draws;

  if (store_draws) {
    draws = Rcpp::List::create(
      Rcpp::Named("H")      = Rcpp::IntegerVector(nsave, NA_INTEGER),
      Rcpp::Named("nu")     = Rcpp::NumericVector(nsave, NA_REAL),
      Rcpp::Named("ell")    = Rcpp::List(nsave),
      Rcpp::Named("tau")    = Rcpp::List(nsave),
      Rcpp::Named("Delta")  = Rcpp::List(nsave),
      Rcpp::Named("Lambda") = Rcpp::List(nsave),
      Rcpp::Named("sigma2") = Rcpp::NumericMatrix(nsave, m)
    );

    if (store_eta) {
      draws["Eta"] = Rcpp::List(nsave);
    } else {
      draws["Eta"] = R_NilValue;
    }
  } else {
    draws = Rcpp::List::create();
  }

  // ---------------------------------------------------------------------------
  // Meta information
  // ---------------------------------------------------------------------------
  Rcpp::List meta = Rcpp::List::create(
    Rcpp::Named("T") = T,
    Rcpp::Named("m") = m,
    Rcpp::Named("Hmax") = Hmax,
    Rcpp::Named("niter") = niter,
    Rcpp::Named("nburn") = nburn,
    Rcpp::Named("thin") = thin,
    Rcpp::Named("nsave") = nsave
  );

  // ---------------------------------------------------------------------------
  // Acceptance statistics placeholders
  // ---------------------------------------------------------------------------
  Rcpp::List accept = Rcpp::List::create(
    Rcpp::Named("H_birth") = NA_REAL,
    Rcpp::Named("H_death") = NA_REAL,
    Rcpp::Named("pivots")  = NA_REAL
  );

  // ---------------------------------------------------------------------------
  // Main MCMC loop
  //
  // Section 3.1 sampler order:
  //
  //   1. Update H.
  //   2. Update nu | H ~ Beta(a_nu + H, b_nu + m - H).
  //   3. Update pivot locations given H.
  //   4. Update tau_j | ell_j, Delta for j = 1, ..., H.
  //   5. Update Delta below pivots, column-wise, marginalizing Lambda and Sigma.
  //   6. Update sigma_i^2, row-wise.
  //   7. Update Lambda row-wise given Delta, Eta, sigma_i^2.
  //   8. Update Eta conditionally on Lambda.
  //
  // Section 3.2:
  //   - H is updated with a random-walk Metropolis-Hastings birth/death move.
  //   - Birth adds one column and a new pivot.
  //   - Death removes the last active column.
  //
  // Section 3.3:
  //   - Pivot locations are updated by a Metropolis-Hastings move.
  //   - Conditionally on new pivots, some entries become structural zeros,
  //     while others become newly available and require proposed indicators.
  //
  // This stub does not yet perform the updates. It only validates the interface,
  // allocates output, and can optionally save the unchanged initial state.
  // ---------------------------------------------------------------------------
  if (verbose) {
    Rcpp::Rcout << "Starting gltfa_cpp stub with "
                << "T = " << T << ", m = " << m << ", H = " << H
                << ", Hmax = " << Hmax << ", niter = " << niter << "."
                << std::endl;
  }

  if (store_draws) {
    int save_idx = 0;

    for (int iter = 0; iter < niter; ++iter) {
      // ---------------------------------------------------------------------
      // STEP 1: Update H
      // Placeholder: no update yet.
      // ---------------------------------------------------------------------

      // ---------------------------------------------------------------------
      // STEP 2: Update nu
      // Placeholder: no update yet.
      // ---------------------------------------------------------------------

      // ---------------------------------------------------------------------
      // STEP 3: Update pivot locations ell
      // Placeholder: no update yet.
      // ---------------------------------------------------------------------

      // ---------------------------------------------------------------------
      // STEP 4: Update tau_j, j = 1, ..., H
      // Placeholder: no update yet.
      // ---------------------------------------------------------------------

      // ---------------------------------------------------------------------
      // STEP 5: Update Delta below pivots
      // Placeholder: no update yet.
      // ---------------------------------------------------------------------

      // ---------------------------------------------------------------------
      // STEP 6: Update sigma_i^2, i = 1, ..., m
      // Placeholder: no update yet.
      // ---------------------------------------------------------------------

      // ---------------------------------------------------------------------
      // STEP 7: Update Lambda row-wise
      // Placeholder: no update yet.
      // ---------------------------------------------------------------------

      // ---------------------------------------------------------------------
      // STEP 8: Update Eta
      // Placeholder: no update yet.
      // ---------------------------------------------------------------------

      const bool should_save =
        (iter + 1 > nburn) &&
        (((iter + 1 - nburn) % thin) == 0);

      if (should_save) {
        Rcpp::IntegerVector H_draw = draws["H"];
        Rcpp::NumericVector nu_draw = draws["nu"];
        Rcpp::List ell_draw = draws["ell"];
        Rcpp::List tau_draw = draws["tau"];
        Rcpp::List Delta_draw = draws["Delta"];
        Rcpp::List Lambda_draw = draws["Lambda"];
        Rcpp::NumericMatrix sigma2_draw = draws["sigma2"];

        H_draw[save_idx] = Rcpp::as<int>(last["H"]);
        nu_draw[save_idx] = Rcpp::as<double>(last["nu"]);
        ell_draw[save_idx] = last["ell"];
        tau_draw[save_idx] = last["tau"];
        Delta_draw[save_idx] = last["Delta"];
        Lambda_draw[save_idx] = last["Lambda"];

        if (store_eta) {
          Rcpp::List Eta_draw = draws["Eta"];
          Eta_draw[save_idx] = last["Eta"];
        }

        Rcpp::NumericVector sigma2_last = last["sigma2"];
        for (int i = 0; i < m; ++i) {
          sigma2_draw(save_idx, i) = sigma2_last[i];
        }

        save_idx++;
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("draws") = draws,
    Rcpp::Named("last") = last,
    Rcpp::Named("meta") = meta,
    Rcpp::Named("accept") = accept
  );
}