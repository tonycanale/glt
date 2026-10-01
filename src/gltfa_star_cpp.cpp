#include <RcppArmadillo.h>
#include "utils_aux.h"
#include "update_H.h"
#include "update_pivots.h"
#include "update_delta.h"
#include "update_Lambda_sigma2.h"
#include "update_Eta.h"
#include "update_alpha_beta.h"
#include "update_a_nu.h"
#include "truncnorm.h"
#include "star_transform.h"
#include "star_bounds.h"
#include "update_beta_reg.h"

// [[Rcpp::depends(RcppArmadillo)]]

//' Fit the STAR-extended GLT factor model to discrete/mixed-scale data
//'
//' Implements the Kowal & Canale (2020) simultaneous transformation and
//' rounding (STAR) approach on top of the existing Gaussian GLT factor
//' sampler. The observed T x m matrix `y_obs` (0-indexed integer levels,
//' possibly with a different number of levels per column) is linked to a
//' continuous latent matrix `z` via
//'
//' \deqn{y_{ij} = \ell \iff G_j(z_{ij}) \in [a_{j,\ell}, a_{j,\ell+1})}
//' \deqn{z_{ij} = x_i'\beta_j + \epsilon_{ij}}
//'
//' where `G_j` is a column-specific monotone transform (`identity` or `exp`),
//' the threshold sequence is column-specific, and `epsilon`
//' follows the same sparse factor model as in [gltfa_cpp()]:
//'
//'   epsilon_t = Lambda * eta_t + noise_t,   noise_t ~ N(0, diag(sigma2))
//'
//' Each MCMC iteration performs two extra steps before the original 8-step
//' Gibbs sweep:
//'
//'   Step 0a: impute z | y_obs, Beta, Lambda, Eta, sigma2  (truncated normal)
//'   Step 0b: update Beta | z, Lambda, Eta, sigma2          (conjugate normal)
//'
//' and then runs Steps 1-8 exactly as in [gltfa_cpp()], operating on the
//' residual `epsilon = z - X * Beta` in place of the raw data.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List gltfa_star_cpp(
  const arma::imat& y_obs,        // T x m, 0-indexed levels
  const arma::mat& X,             // T x p design matrix (p may be 0)
  const Rcpp::List& thresholds,   // length-m list of interior cutpoints
  const Rcpp::IntegerVector& g_type, // length-m transform ids
  const Rcpp::List& mcmc,
  const Rcpp::List& model,
  const Rcpp::List& prior,
  const Rcpp::List& init,
  const Rcpp::List& control,
  const Rcpp::List& fixed
) {
  // ---------------------------------------------------------------------------
    // Dimensions and basic checks
  // ---------------------------------------------------------------------------
    const int T = get_int(model, "T");
    const int m = get_int(model, "m");
    const int Hmax = get_int(model, "Hmax");
    const int p = X.n_cols;

    if (static_cast<int>(y_obs.n_rows) != T || static_cast<int>(y_obs.n_cols) != m) {
      Rcpp::stop("`y_obs` has dimensions %d x %d, but `model` declares T = %d and m = %d.",
                 y_obs.n_rows, y_obs.n_cols, T, m);
    }
    if (static_cast<int>(X.n_rows) != T && p > 0) {
      Rcpp::stop("`X` must have T = %d rows.", T);
    }

    // ---------------------------------------------------------------------------
      // Precompute STAR bounds (do not depend on model parameters)
    // ---------------------------------------------------------------------------
      Rcpp::List bounds = compute_star_bounds(y_obs, thresholds, g_type);
      const arma::mat z_lower = Rcpp::as<arma::mat>(bounds["z_lower"]);
      const arma::mat z_upper = Rcpp::as<arma::mat>(bounds["z_upper"]);

      const arma::mat XtX = (p > 0) ? (X.t() * X) : arma::mat(0, 0);

    // ---------------------------------------------------------------------------
      // MCMC settings
    // ---------------------------------------------------------------------------
      const int niter = get_int(mcmc, "niter");
      const int nburn = get_int(mcmc, "nburn");
      const int thin  = get_int(mcmc, "thin");
      const int nsave = (niter - nburn) / thin;

      // ---------------------------------------------------------------------------
        // Prior hyperparameters
      // ---------------------------------------------------------------------------
        double       a_nu    = get_double(prior, "a_nu");
      const double b_nu    = get_double(prior, "b_nu");
      double       alpha   = get_double(prior, "alpha");
      double       beta    = get_double(prior, "beta");
      const double kappa   = get_double(prior, "kappa");
      const double a_sigma = get_double(prior, "a_sigma");
      const double b_sigma = get_double(prior, "b_sigma");
      const double a_alpha     = prior.containsElementNamed("a_alpha")     ? get_double(prior, "a_alpha")     : 1.0;
      const double b_alpha     = prior.containsElementNamed("b_alpha")     ? get_double(prior, "b_alpha")     : 1.0;
      const double mh_sd_alpha = prior.containsElementNamed("mh_sd_alpha") ? get_double(prior, "mh_sd_alpha") : 0.2;
      const double a_beta      = prior.containsElementNamed("a_beta")      ? get_double(prior, "a_beta")      : 1.0;
      const double b_beta      = prior.containsElementNamed("b_beta")      ? get_double(prior, "b_beta")      : 1.0;
      const double mh_sd_beta  = prior.containsElementNamed("mh_sd_beta")  ? get_double(prior, "mh_sd_beta")  : 0.2;
      const double a_anu       = prior.containsElementNamed("a_anu")       ? get_double(prior, "a_anu")       : 1.0;
      const double b_anu       = prior.containsElementNamed("b_anu")       ? get_double(prior, "b_anu")       : 1.0;
      const double mh_sd_a_nu  = prior.containsElementNamed("mh_sd_a_nu")  ? get_double(prior, "mh_sd_a_nu")  : 0.2;

      // Regression prior: beta_j ~ N(b0, sigma2_j * V0), V0 = c0 * I_p by default
      arma::vec b0(p, arma::fill::zeros);
      arma::mat V0_inv(p, p, arma::fill::zeros);
      if (p > 0) {
        if (prior.containsElementNamed("b0_beta")) {
          b0 = Rcpp::as<arma::vec>(prior["b0_beta"]);
        }
        if (prior.containsElementNamed("V_beta_inv")) {
          V0_inv = Rcpp::as<arma::mat>(prior["V_beta_inv"]);
        } else {
          const double c0 = prior.containsElementNamed("c0_beta") ? get_double(prior, "c0_beta") : 100.0;
          V0_inv = arma::eye(p, p) / c0;
        }
      }

      // ---------------------------------------------------------------------------
        // Control options
      // ---------------------------------------------------------------------------
        const bool store_draws  = get_bool(control, "store_draws");
      const bool store_eta    = get_bool(control, "store_eta");
      const int  print_every  = get_int(control, "print_every");
      const bool verbose      = get_bool(control, "verbose");
      const bool random_scan  = get_bool(control, "random_scan");

      (void) print_every;

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
      bool FIX_LambdaSigma = false;
      bool FIX_Eta = false;
      bool FIX_alpha = false;
      bool FIX_beta  = false;
      bool FIX_a_nu  = false;
      bool FIX_z     = false;
      bool FIX_Beta  = false;

      if (!fixed.isNULL()) {
        if (fixed.containsElementNamed("all"))          FIX_all          = Rcpp::as<bool>(fixed["all"]);
        if (fixed.containsElementNamed("H"))            FIX_H            = Rcpp::as<bool>(fixed["H"]);
        if (fixed.containsElementNamed("nu"))           FIX_nu           = Rcpp::as<bool>(fixed["nu"]);
        if (fixed.containsElementNamed("pivots"))       FIX_pivots       = Rcpp::as<bool>(fixed["pivots"]);
        if (fixed.containsElementNamed("tau"))          FIX_tau          = Rcpp::as<bool>(fixed["tau"]);
        if (fixed.containsElementNamed("alpha"))        FIX_alpha        = Rcpp::as<bool>(fixed["alpha"]);
        if (fixed.containsElementNamed("beta"))         FIX_beta         = Rcpp::as<bool>(fixed["beta"]);
        if (fixed.containsElementNamed("a_nu"))         FIX_a_nu         = Rcpp::as<bool>(fixed["a_nu"]);
        if (fixed.containsElementNamed("Delta"))        FIX_Delta        = Rcpp::as<bool>(fixed["Delta"]);
        if (fixed.containsElementNamed("LambdaSigma"))  FIX_LambdaSigma  = Rcpp::as<bool>(fixed["LambdaSigma"]);
        if (fixed.containsElementNamed("Eta"))          FIX_Eta          = Rcpp::as<bool>(fixed["Eta"]);
        if (fixed.containsElementNamed("z"))            FIX_z            = Rcpp::as<bool>(fixed["z"]);
        if (fixed.containsElementNamed("Beta"))         FIX_Beta         = Rcpp::as<bool>(fixed["Beta"]);
      }

      if (FIX_all) {
        FIX_H = FIX_nu = FIX_pivots = FIX_tau = FIX_Delta = FIX_LambdaSigma = FIX_Eta = FIX_z = FIX_Beta = true;
      }
      if (p == 0) {
        FIX_Beta = true; // nothing to update
      }

      // ---------------------------------------------------------------------------
        // Initial state
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

      // Beta initial state (p x m)
      arma::mat Beta(p, m, arma::fill::zeros);
      if (init.containsElementNamed("Beta") && init["Beta"] != R_NilValue) {
        Beta = Rcpp::as<arma::mat>(init["Beta"]);
        if (static_cast<int>(Beta.n_rows) != p || static_cast<int>(Beta.n_cols) != m) {
          Rcpp::stop("`init$Beta` must have dimension p x m.");
        }
      }
      last["Beta"] = Beta;

      // z initial state (T x m): default to a naive transform-consistent guess
      arma::mat z(T, m, arma::fill::zeros);
      if (init.containsElementNamed("z") && init["z"] != R_NilValue) {
        z = Rcpp::as<arma::mat>(init["z"]);
        if (static_cast<int>(z.n_rows) != T || static_cast<int>(z.n_cols) != m) {
          Rcpp::stop("`init$z` must have dimension T x m.");
        }
      } else {
        // Midpoint-ish initialisation, clipped into (z_lower, z_upper)
        for (int t = 0; t < T; ++t) {
          for (int j = 0; j < m; ++j) {
            double lo = z_lower(t, j);
            double hi = z_upper(t, j);
            double val;
            if (std::isfinite(lo) && std::isfinite(hi)) {
              val = 0.5 * (lo + hi);
            } else if (std::isfinite(lo)) {
              val = lo + 1.0;
            } else if (std::isfinite(hi)) {
              val = hi - 1.0;
            } else {
              val = 0.0;
            }
            z(t, j) = val;
          }
        }
      }
      last["z"] = z;

      // ---------------------------------------------------------------------------
        // Allocate storage for draws
      // ---------------------------------------------------------------------------
        Rcpp::List draws;

      if (store_draws) {
        draws = Rcpp::List::create(
          Rcpp::Named("H")      = Rcpp::IntegerVector(nsave, NA_INTEGER),
          Rcpp::Named("nu")     = Rcpp::NumericVector(nsave, NA_REAL),
          Rcpp::Named("a_nu")   = Rcpp::NumericVector(nsave, NA_REAL),
          Rcpp::Named("alpha")  = Rcpp::NumericVector(nsave, NA_REAL),
          Rcpp::Named("beta")   = Rcpp::NumericVector(nsave, NA_REAL),
          Rcpp::Named("ell")    = Rcpp::List(nsave),
          Rcpp::Named("tau")    = Rcpp::List(nsave),
          Rcpp::Named("Delta")  = Rcpp::List(nsave),
          Rcpp::Named("Lambda") = Rcpp::List(nsave),
          Rcpp::Named("Beta")   = Rcpp::List(nsave),
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
          Rcpp::Named("p") = p,
          Rcpp::Named("Hmax") = Hmax,
          Rcpp::Named("niter") = niter,
          Rcpp::Named("nburn") = nburn,
          Rcpp::Named("thin") = thin,
          Rcpp::Named("nsave") = nsave
        );

      // ---------------------------------------------------------------------------
        // Acceptance statistics placeholders
      // ---------------------------------------------------------------------------
        arma::imat accept(niter, 1, arma::fill::zeros);

      if (verbose) {
        Rcpp::Rcout << "Starting gltfa_star_cpp with "
        << "T = " << T << ", m = " << m << ", p = " << p << ", H = " << H
        << ", Hmax = " << Hmax << ", niter = " << niter << "."
        << std::endl;
      }

      if (store_draws) {
        int save_idx = 0;
        for (int iter = 0; iter < niter; ++iter) {

          // ---------------------------------------------------------------------
            // STEP 0a: impute latent z | y_obs, Beta, Lambda, Eta, sigma2
          // ---------------------------------------------------------------------
            Rcpp::NumericVector sigma2_v0 = Rcpp::as<Rcpp::NumericVector>(last["sigma2"]);
            arma::vec sigma2_vec0 = Rcpp::as<arma::vec>(sigma2_v0);
            arma::vec sd0 = arma::sqrt(sigma2_vec0);

            if (!FIX_z) {
              arma::mat mu = (Lambda * Eta).t();      // T x m, factor part
              if (p > 0) {
                mu += X * Beta;                       // add regression part
              }
              arma::mat u_rand = arma::randu(T, m);
              z = truncnorm_lg(z_lower, z_upper, mu, sd0, u_rand);
              last["z"] = z;
            } else {
              z = Rcpp::as<arma::mat>(last["z"]);
            }

          // ---------------------------------------------------------------------
            // STEP 0b: update Beta | z, Lambda, Eta, sigma2
          // ---------------------------------------------------------------------
            if (!FIX_Beta && p > 0) {
              arma::mat r = z - (Lambda * Eta).t();   // T x m, regression target
              update_Beta_cpp(r, X, XtX, sigma2_vec0, V0_inv, b0, Beta);
              last["Beta"] = Beta;
            }

          // Working "data" matrix for the unchanged Gaussian sampler: the
          // residual after removing the regression mean.
          arma::mat epsilon = z;
          if (p > 0) {
            epsilon = z - X * Beta;
          }

          // ---------------------------------------------------------------------
            // STEP 1: Update H
          // ---------------------------------------------------------------------
            int H_old = H;
            if (!FIX_H) {
              Rcpp::List H_step = update_H_cpp(epsilon, Delta, Eta, prior,
                                               Rcpp::as<double>(last["nu"]), Hmax, 0.5);

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
                arma::mat Lambda_new(m, H, arma::fill::zeros);
                for (int j = 0; j < H_old; ++j) Lambda_new.col(j) = Lambda_old.col(j);
                last["ell"] = ell; last["tau"] = tau; last["Lambda"] = Lambda_new;
              } else if (H < H_old && H_accepted) {
                Rcpp::IntegerVector ell_old = Rcpp::as<Rcpp::IntegerVector>(last["ell"]);
                Rcpp::NumericVector tau_old = Rcpp::as<Rcpp::NumericVector>(last["tau"]);
                arma::mat Lambda_old = Rcpp::as<arma::mat>(last["Lambda"]);
                Rcpp::IntegerVector ell(H); Rcpp::NumericVector tau(H);
                arma::mat Lambda_new(m, H, arma::fill::zeros);
                for (int j = 0; j < H; ++j) {
                  ell[j] = ell_old[j]; tau[j] = tau_old[j];
                  Lambda_new.col(j) = Lambda_old.col(j);
                }
                last["ell"] = ell; last["tau"] = tau; last["Lambda"] = Lambda_new;
              }
            } // end !FIX_H

            Rcpp::IntegerVector ell = Rcpp::as<Rcpp::IntegerVector>(last["ell"]);
            Rcpp::NumericVector tau = Rcpp::as<Rcpp::NumericVector>(last["tau"]);
            Lambda = Rcpp::as<arma::mat>(last["Lambda"]);
            Rcpp::NumericVector sigma2 = Rcpp::as<Rcpp::NumericVector>(last["sigma2"]);

          // ---------------------------------------------------------------------
            // STEP 2: Update nu
          // ---------------------------------------------------------------------
            if (!FIX_nu) {
              last["nu"] = R::rbeta(a_nu + H, b_nu + m - H);
            }

          // ---------------------------------------------------------------------
            // STEP 2b: Update a_nu via MH (nu marginalised out)
          // ---------------------------------------------------------------------
            if (!FIX_a_nu) {
              a_nu = update_a_nu_cpp(a_nu, H, m, b_nu, a_anu, b_anu, mh_sd_a_nu);
              last["a_nu"] = a_nu;
            }

          // ---------------------------------------------------------------------
            // STEP 3: Update pivot locations ell
          // ---------------------------------------------------------------------
            if (!FIX_pivots) {
              Rcpp::List pivots_step = update_pivots_cpp(epsilon, Delta, Eta, prior);
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
            // STEP 4c: Update beta via MH (tau_j marginalised out)
          // ---------------------------------------------------------------------
            if (!FIX_beta) {
              beta = update_beta_cpp(beta, Delta, ell, m,
                                     alpha, a_beta, b_beta, mh_sd_beta);
              last["beta"] = beta;
            }

          // ---------------------------------------------------------------------
            // STEP 5: Update Delta entries below pivots
          // ---------------------------------------------------------------------
            if (!FIX_Delta) {
              const arma::vec tau_vec = Rcpp::as<arma::vec>(
                Rcpp::NumericVector(last["tau"]));
              Delta = update_delta_cpp(epsilon, Delta, Eta, tau_vec, prior, random_scan);
              last["Delta"] = Delta;
            }

          // ---------------------------------------------------------------------
            // STEP 6+7: Update sigma2 and Lambda (joint NIG update)
          // ---------------------------------------------------------------------
            if (!FIX_LambdaSigma) {
              arma::vec sigma2_vec = Rcpp::as<arma::vec>(
                Rcpp::NumericVector(last["sigma2"]));
              Lambda = Rcpp::as<arma::mat>(last["Lambda"]);

              update_sigma2_Lambda(epsilon, Eta, Delta, sigma2_vec, Lambda,
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
              update_Eta(epsilon, Lambda, sigma2_vec, Eta);
              last["Eta"] = Eta;
            }

          // Save draws if past burn-in and at the correct thinning interval -----
            const bool should_save =
            (iter + 1 > nburn) &&
            (((iter + 1 - nburn) % thin) == 0);

          if (should_save) {
            Rcpp::IntegerVector H_draw = draws["H"];
            Rcpp::NumericVector nu_draw    = draws["nu"];
            Rcpp::NumericVector a_nu_draw  = draws["a_nu"];
            Rcpp::NumericVector alpha_draw = draws["alpha"];
            Rcpp::NumericVector beta_draw  = draws["beta"];
            Rcpp::List ell_draw = draws["ell"];
            Rcpp::List tau_draw = draws["tau"];
            Rcpp::List Delta_draw = draws["Delta"];
            Rcpp::List Lambda_draw = draws["Lambda"];
            Rcpp::List Beta_draw = draws["Beta"];
            Rcpp::NumericMatrix sigma2_draw = draws["sigma2"];
            H_draw[save_idx]      = Rcpp::as<int>(last["H"]);
            nu_draw[save_idx]    = Rcpp::as<double>(last["nu"]);
            a_nu_draw[save_idx]  = a_nu;
            alpha_draw[save_idx] = alpha;
            beta_draw[save_idx]  = beta;
            ell_draw[save_idx]    = Rcpp::clone(Rcpp::as<Rcpp::IntegerVector>(last["ell"]));
            tau_draw[save_idx]    = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(last["tau"]));
            Delta_draw[save_idx]  = Rcpp::clone(Rcpp::wrap(Rcpp::as<arma::imat>(last["Delta"])));
            Lambda_draw[save_idx] = Rcpp::clone(Rcpp::wrap(Rcpp::as<arma::mat>(last["Lambda"])));
            Beta_draw[save_idx]   = Rcpp::clone(Rcpp::wrap(Rcpp::as<arma::mat>(last["Beta"])));

            draws["ell"]    = ell_draw;
            draws["tau"]    = tau_draw;
            draws["Delta"]  = Delta_draw;
            draws["Lambda"] = Lambda_draw;
            draws["Beta"]   = Beta_draw;
            draws["alpha"]       = alpha_draw;
            draws["beta"]        = beta_draw;

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

      last["epsilon"] = z - (p > 0 ? arma::mat(X * Beta) : arma::mat(T, m, arma::fill::zeros));

      return Rcpp::List::create(
        Rcpp::Named("draws") = draws,
        Rcpp::Named("last") = last,
        Rcpp::Named("meta") = meta,
        Rcpp::Named("accept") = accept
      );
}
