#include <RcppArmadillo.h>
#include "helpers.h"
#include "marglik.h"
#include "update_H.h"

// [[Rcpp::depends(RcppArmadillo)]]


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
      // Birth/death MCMC for dimension of latent factors
      // ---------------------------------------------------------------------
      Rcpp::List H_step = update_H_cpp(y, Delta, Eta, prior, 
                                       Rcpp::as<double>(last["nu"]), 0.5);
      
      Delta = Rcpp::as<arma::imat>(H_step["Delta"]);
      Eta = Rcpp::as<arma::mat>(H_step["Eta"]);
      int H_old = H;
      H = Delta.n_cols;
      bool H_accepted = Rcpp::as<bool>(H_step["accepted"]);
      int H_increased = Rcpp::as<int>(H_step["increased"]);
      
      // Store accept indicator
      accept(iter, 0) = H_accepted ? 1 : 0;
      
      // Update state to reflect new H, Delta, Eta
      last["H"] = H;
      last["Delta"] = Delta;
      last["Eta"] = Eta;
      
      // Allocate/deallocate vectors that depend on H
      if (H > H_old && H_accepted) {
        // Birth: dimension increased, allocate new storage
        Rcpp::IntegerVector ell = Rcpp::IntegerVector(H);
        Rcpp::NumericVector tau = Rcpp::NumericVector(H);
        arma::mat Lambda_old = Rcpp::as<arma::mat>(last["Lambda"]);
        arma::mat Lambda(m, H, arma::fill::zeros);
        for (int j = 0; j < H_old; ++j) {
          Lambda.col(j) = Lambda_old.col(j);
        }
        last["ell"] = ell;
        last["tau"] = tau;
        last["Lambda"] = Lambda;
      } else if (H < H_old && H_accepted) {
        // Death: dimension decreased, shrink storage
        Rcpp::IntegerVector ell_old = Rcpp::as<Rcpp::IntegerVector>(last["ell"]);
        Rcpp::NumericVector tau_old = Rcpp::as<Rcpp::NumericVector>(last["tau"]);
        arma::mat Lambda_old = Rcpp::as<arma::mat>(last["Lambda"]);
        
        Rcpp::IntegerVector ell = Rcpp::IntegerVector(H);
        Rcpp::NumericVector tau = Rcpp::NumericVector(H);
        arma::mat Lambda(m, H, arma::fill::zeros);
        for (int j = 0; j < H; ++j) {
          ell[j] = ell_old[j];
          tau[j] = tau_old[j];
          Lambda.col(j) = Lambda_old.col(j);
        }
        last["ell"] = ell;
        last["tau"] = tau;
        last["Lambda"] = Lambda;
      }

      // ...existing code...

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