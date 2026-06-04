#include <RcppArmadillo.h>
#include "utils_aux.h"
#include "update_H.h"
#include "update_pivots.h"
#include "update_delta.h"
#include "update_Lambda_sigma2.h"
#include "update_Eta.h"
#include "update_alpha.h"   // <<< NEW

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
Rcpp::List gltfa_cpp(
  const arma::mat& y,
  const Rcpp::List& mcmc,
  const Rcpp::List& model,
  const Rcpp::List& prior,
  const Rcpp::List& init,
  const Rcpp::List& control,
  const Rcpp::List& fixed // NEW: list of logical flags to freeze updates
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
      // nsave is derived, not supplied: (niter - nburn) / thin.
      // Integer division is intentional.
      const int nsave = (niter - nburn) / thin;
      
      // ---------------------------------------------------------------------------
        // Prior hyperparameters
      // ---------------------------------------------------------------------------
        const double a_nu    = get_double(prior, "a_nu");
      const double b_nu    = get_double(prior, "b_nu");
      double       alpha   = get_double(prior, "alpha");  // mutable – updated by MH
      const double beta    = get_double(prior, "beta");
      const double kappa   = get_double(prior, "kappa");
      const double a_sigma = get_double(prior, "a_sigma");
      const double b_sigma = get_double(prior, "b_sigma");
      // Hyperprior alpha ~ Ga(a_alpha, b_alpha); mh_sd_alpha: RW std on log scale
      const double a_alpha     = prior.containsElementNamed("a_alpha")     ? get_double(prior, "a_alpha")     : 1.0;
      const double b_alpha     = prior.containsElementNamed("b_alpha")     ? get_double(prior, "b_alpha")     : 1.0;
      const double mh_sd_alpha = prior.containsElementNamed("mh_sd_alpha") ? get_double(prior, "mh_sd_alpha") : 0.2;
      
      // Avoid unused variable warnings in the stub.
      (void) a_nu;
      (void) b_nu;
      (void) kappa;
      (void) a_sigma;
      (void) b_sigma;
      
      // ---------------------------------------------------------------------------
        // Control options
      // ---------------------------------------------------------------------------
        const bool store_draws  = get_bool(control, "store_draws");
      const bool store_eta    = get_bool(control, "store_eta");
      const int  print_every  = get_int(control, "print_every");
      const bool verbose      = get_bool(control, "verbose");
      const bool random_scan  = get_bool(control, "random_scan");
      
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
        // Fixed flags parsing
      // ---------------------------------------------------------------------------
        bool FIX_all = false;
      bool FIX_H = false;
      bool FIX_nu = false;
      bool FIX_pivots = false;
      bool FIX_tau = false;
      bool FIX_Delta = false;
      bool FIX_LambdaSigma = false; // joint flag for sigma2+Lambda update
      bool FIX_Eta = false;
      bool FIX_alpha = false;   // <<< NEW
      
      if (!fixed.isNULL()) {
        if (fixed.containsElementNamed("all")) FIX_all = Rcpp::as<bool>(fixed["all"]);
        if (fixed.containsElementNamed("H")) FIX_H = Rcpp::as<bool>(fixed["H"]);
        if (fixed.containsElementNamed("nu")) FIX_nu = Rcpp::as<bool>(fixed["nu"]);
        if (fixed.containsElementNamed("pivots")) FIX_pivots = Rcpp::as<bool>(fixed["pivots"]);
        if (fixed.containsElementNamed("tau")) FIX_tau = Rcpp::as<bool>(fixed["tau"]);
        if (fixed.containsElementNamed("alpha")) FIX_alpha = Rcpp::as<bool>(fixed["alpha"]);  // <<< NEW
        if (fixed.containsElementNamed("Delta")) FIX_Delta = Rcpp::as<bool>(fixed["Delta"]);
        if (fixed.containsElementNamed("LambdaSigma")) FIX_LambdaSigma = Rcpp::as<bool>(fixed["LambdaSigma"]);
        if (fixed.containsElementNamed("Eta")) FIX_Eta = Rcpp::as<bool>(fixed["Eta"]);
      }
      
      if (FIX_all) {
        FIX_H = FIX_nu = FIX_pivots = FIX_tau = FIX_Delta = FIX_LambdaSigma = FIX_Eta = true;
      }
      
      // ---------------------------------------------------------------------------
        // Initial state
      // Conventions:
        //   Lambda : m x H
      //   Delta  : m x H
      //   Eta    : H x T
      //   y      : T x m
      // ---------------------------------------------------------------------------
        int H = get_int(init, "H");
      
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
          Rcpp::Named("alpha")  = Rcpp::NumericVector(nsave, NA_REAL),  // <<< NEW
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
        arma::imat accept(niter, 1, arma::fill::zeros);  // column 0: H birth/death
      
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
        // Main MCMC loop starts here -------------------------------------------
          for (int iter = 0; iter < niter; ++iter) {
            // ---------------------------------------------------------------------
              // STEP 1: Update H
            // ---------------------------------------------------------------------
              int H_old = H;
              if (!FIX_H) {
                Rcpp::List H_step = update_H_cpp(y, Delta, Eta, prior,
                                                 Rcpp::as<double>(last["nu"]), 0.5);
                
                Delta = Rcpp::as<arma::imat>(H_step["Delta"]);
                Eta   = Rcpp::as<arma::mat>(H_step["Eta"]);
                H     = Delta.n_cols;
                bool H_accepted = Rcpp::as<bool>(H_step["accepted"]);
                accept(iter, 0) = H_accepted ? 1 : 0;
                
                last["H"]     = H;
                last["Delta"] = Delta;
                last["Eta"]   = Eta;
                
                if (H > H_old && H_accepted) {
                  Rcpp::IntegerVector ell(H);
                  Rcpp::NumericVector tau(H);
                  arma::mat Lambda_old = Rcpp::as<arma::mat>(last["Lambda"]);
                  arma::mat Lambda(m, H, arma::fill::zeros);
                  for (int j = 0; j < H_old; ++j) Lambda.col(j) = Lambda_old.col(j);
                  last["ell"] = ell; last["tau"] = tau; last["Lambda"] = Lambda;
                } else if (H < H_old && H_accepted) {
                  Rcpp::IntegerVector ell_old = Rcpp::as<Rcpp::IntegerVector>(last["ell"]);
                  Rcpp::NumericVector tau_old = Rcpp::as<Rcpp::NumericVector>(last["tau"]);
                  arma::mat Lambda_old = Rcpp::as<arma::mat>(last["Lambda"]);
                  Rcpp::IntegerVector ell(H); Rcpp::NumericVector tau(H);
                  arma::mat Lambda(m, H, arma::fill::zeros);
                  for (int j = 0; j < H; ++j) {
                    ell[j] = ell_old[j]; tau[j] = tau_old[j];
                    Lambda.col(j) = Lambda_old.col(j);
                  }
                  last["ell"] = ell; last["tau"] = tau; last["Lambda"] = Lambda;
                }
              } // end !FIX_H
              
              // State views for this iteration
              Rcpp::IntegerVector ell = Rcpp::as<Rcpp::IntegerVector>(last["ell"]);
              Rcpp::NumericVector tau = Rcpp::as<Rcpp::NumericVector>(last["tau"]);
              arma::mat Lambda = Rcpp::as<arma::mat>(last["Lambda"]);
              Rcpp::NumericVector sigma2 = Rcpp::as<Rcpp::NumericVector>(last["sigma2"]);
              // ---------------------------------------------------------------------
                // STEP 2: Update nu
              // ---------------------------------------------------------------------
                if (!FIX_nu) {
                  last["nu"] = R::rbeta(a_nu + H, b_nu + m - H);
                }
              
              // ---------------------------------------------------------------------
                // STEP 3: Update pivot locations ell
              // ---------------------------------------------------------------------
                if (!FIX_pivots) {
                  Rcpp::List pivots_step = update_pivots_cpp(y, Delta, Eta, prior);
                  Delta = Rcpp::as<arma::imat>(pivots_step["Delta"]);
                  ell   = Rcpp::as<Rcpp::IntegerVector>(pivots_step["ell"]);
                  last["Delta"] = Delta;
                  last["ell"]   = ell;
                }
              
              // ---------------------------------------------------------------------
                // STEP 4: Update tau_j
              // ---------------------------------------------------------------------
                if (!FIX_tau) {
                  for (int j = 0; j < H; ++j) {
                    const int ell_j = ell[j];
                    int d_j = 0;
                    for (int i = ell_j; i < m; ++i) d_j += Delta(i, j);
                    tau[j] = R::rbeta(alpha * beta + d_j,
                                      beta + m - ell_j - d_j);
                  }
                  last["tau"] = tau;
                }
              
              // ---------------------------------------------------------------------
                // STEP 4b: Update alpha via MH (tau_j marginalised out)
              // ---------------------------------------------------------------------
                if (!FIX_alpha) {
                  alpha = update_alpha_cpp(alpha, Delta, ell, m,
                                           beta, a_alpha, b_alpha, mh_sd_alpha);
                  last["alpha"] = alpha;
                }
              
              // ---------------------------------------------------------------------
                // STEP 5: Update Delta entries below pivots
              // ---------------------------------------------------------------------
                if (!FIX_Delta) {
                  const arma::vec tau_vec = Rcpp::as<arma::vec>(
                    Rcpp::NumericVector(last["tau"]));
                  Delta = update_delta_cpp(y, Delta, Eta, tau_vec, prior, random_scan);
                  last["Delta"] = Delta;
                }
              
              // ---------------------------------------------------------------------
                // STEP 6+7: Update sigma2 and Lambda (joint NIG update)
              // ---------------------------------------------------------------------
                if (!FIX_LambdaSigma) {
                  arma::vec sigma2_vec = Rcpp::as<arma::vec>(
                    Rcpp::NumericVector(last["sigma2"]));
                  Lambda = Rcpp::as<arma::mat>(last["Lambda"]);
                  
                  update_sigma2_Lambda(y, Eta, Delta, sigma2_vec, Lambda,
                                       a_sigma, b_sigma, kappa);
                  
                  last["sigma2"] = Rcpp::wrap(sigma2_vec);
                  last["Lambda"] = Lambda;
                  sigma2 = Rcpp::as<Rcpp::NumericVector>(last["sigma2"]);
                }
              
              // ---------------------------------------------------------------------
                // STEP 8: Update Eta
              // ---------------------------------------------------------------------
                if (!FIX_Eta) {
                  arma::vec sigma2_vec = Rcpp::as<arma::vec>(
                    Rcpp::NumericVector(last["sigma2"]));
                  Lambda = Rcpp::as<arma::mat>(last["Lambda"]);
                  update_Eta(y, Lambda, sigma2_vec, Eta);
                  last["Eta"] = Eta;
                }
              
              
              // Save draws if past burn-in and at the correct thinning interval -----
                const bool should_save =
                (iter + 1 > nburn) &&
                (((iter + 1 - nburn) % thin) == 0);
              
              if (should_save) {
                Rcpp::IntegerVector H_draw = draws["H"];
                Rcpp::NumericVector nu_draw = draws["nu"];
                Rcpp::NumericVector alpha_draw = draws["alpha"];  // <<< NEW
                Rcpp::List ell_draw = draws["ell"];
                Rcpp::List tau_draw = draws["tau"];
                Rcpp::List Delta_draw = draws["Delta"];
                Rcpp::List Lambda_draw = draws["Lambda"];
                Rcpp::NumericMatrix sigma2_draw = draws["sigma2"];
                
                H_draw[save_idx]      = Rcpp::as<int>(last["H"]);
                nu_draw[save_idx]    = Rcpp::as<double>(last["nu"]);
                alpha_draw[save_idx] = alpha;                      // <<< NEW
                draws["alpha"]       = alpha_draw;                 // <<< NEW
                ell_draw[save_idx]    = Rcpp::clone(Rcpp::as<Rcpp::IntegerVector>(last["ell"]));
                tau_draw[save_idx]    = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(last["tau"]));
                Delta_draw[save_idx]  = Rcpp::clone(Rcpp::wrap(Rcpp::as<arma::imat>(last["Delta"])));
                Lambda_draw[save_idx] = Rcpp::clone(Rcpp::wrap(Rcpp::as<arma::mat>(last["Lambda"])));
                
                // Write local copies back — Rcpp::List operator[] returns a copy,
                // not a reference, so assignments above do not modify draws in-place.
                draws["ell"]    = ell_draw;
                draws["tau"]    = tau_draw;
                draws["Delta"]  = Delta_draw;
                draws["Lambda"] = Lambda_draw;
                
                if (store_eta) {
                  Rcpp::List Eta_draw = draws["Eta"];
                  Eta_draw[save_idx] = Rcpp::clone(Rcpp::wrap(Rcpp::as<arma::mat>(last["Eta"])));
                  draws["Eta"] = Eta_draw;
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