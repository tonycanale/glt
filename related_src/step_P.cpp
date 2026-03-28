#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include "Eigen/src/Core/util/Constants.h"
#include "Eigen/src/OrderingMethods/Ordering.h"
#include "steps.h"
#include "structs.h"
#include "util.h"
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
List sample_facload_one_r(
    const List& data_in,
    const List& model_in,
    const List& prior_in,
    const bool simpost) {
  Data data = list2data(data_in);
  Model model = list2model(model_in);
  Prior prior = list2prior(prior_in);
  Facload_one_out facload_out;
  if (simpost) {
    facload_out = sample_facload_one(data, model, prior, true);
  } else {
    facload_out = sample_facload_one(data, model, prior);
  }
  return List::create(
      _["model"] = model2list(model),
      _["post"] = List::create(
        _["beta"] = List::create(
          _["b"] = wrap(facload_out.beta.b),
          _["B"] = wrap(facload_out.beta.B)),
        _["sigma"] = List::create(
          _["c"] = wrap(facload_out.sigma.c),
          _["C"] = wrap(facload_out.sigma.C))));
}

// [[Rcpp::export]]
List sample_facload_joint_r(
    const List& data_in,
    const List& model_in,
    const List& prior_in,
    const int flag) {
  Data data = list2data(data_in);
  Model model = list2model(model_in);
  Prior prior = list2prior(prior_in);
  Facload_joint_out facload_out;
  double SSBeta = 0;
  if (flag == 0) {
    SSBeta = sample_facload_joint(data, model, prior);
  } else {
    facload_out = posterior_facload_joint(data, model, prior);
  }
  return List::create(
      _["model"] = model2list(model),
      _["post"] = flag == 0 ? R_NilValue : as<SEXP>(List::create(
        _["betadetl"] = wrap(facload_out.betadetl),
        _["sigmalogcr"] = wrap(facload_out.sigmalogcr),
        _["cT"] = wrap(facload_out.cT))),
      _["SSBeta"] = wrap(SSBeta));
}

// Sample all factor loadings and idiosyncratic variances simultaneously for a one-factor model
// Implements simplified version of Step (P) for dedicated rows (Step P-b, eq 2982-3001)
// Used for zero rows in Algorithm P and computes marginal likelihood for all rows
Facload_one_out sample_facload_one(
    const Data& data,
    const Model& model,
    const Prior& prior) {
  // some quick checks
  if (model.delta.n_rows != 1 or
      model.d != 1 or
      data.X.n_rows > 1) {
    stop("Invalid input in 'facload_one'");
  }
  // to facilitate understandability and stronger type checks, we define "aliases"
  const arma::urowvec delta = model.delta.row(0);
  const arma::rowvec X = data.X.as_row();

  Facload_one_out post;
  post.sigma.c = prior.par.sigma.c + 0.5 * data.N;
  post.sigma.C.zeros(data.r);

  const arma::uvec delta_zero = arma::find(delta == 0);
  const arma::uvec delta_one = arma::find(delta == 1);
  if (delta_zero.n_elem > 0) {
    post.sigma.C(delta_zero) = prior.par.sigma.C(delta_zero) + 0.5 * arma::sum(arma::square(data.y.rows(delta_zero)), 1);
  }

  if (delta_one.n_elem > 0) {
    const double XtX = arma::dot(X, X);
    const arma::vec Xty = static_cast<arma::vec>(arma::sum(data.y.each_row() % X, 1))(delta_one);
    switch (prior.type) {
      case Prior::Type::FRAC: {
        post.beta.B = 1 / XtX;
        post.beta.b = post.beta.B * Xty;
        const arma::mat yhat = post.beta.b * X;
        post.sigma.c(delta_one) = post.sigma.c(delta_one)
          - 0.5 * data.N * prior.bfrac;
        post.sigma.C(delta_one) = prior.par.sigma.C(delta_one)
          + 0.5 * (1 - prior.bfrac) * arma::sum(arma::square(data.y.rows(delta_one) - yhat), 1);
        break;
      }
      case Prior::Type::UNIT: {
        const double prior_par_beta_b = 0;
        post.beta.B = 1 / (prior.par.beta.B0inv + XtX);
        post.beta.b = (prior.par.beta.B0inv * prior_par_beta_b + Xty) * post.beta.B;
        const arma::mat yhat = post.beta.b * X;
        post.sigma.C(delta_one) = prior.par.sigma.C(delta_one)
          + 0.5 * arma::square(prior_par_beta_b - post.beta.b) * prior.par.beta.B0inv
          + 0.5 * arma::sum(arma::square(data.y.rows(delta_one) - yhat), 1);
        break;
      }
      default:
        stop("Unknown type of prior");
    }
  }
  return post;
}

Facload_one_out sample_facload_one(
    const Data& data,
    Model& model,
    const Prior& prior,
    const bool simpost) {
  const Facload_one_out post = sample_facload_one(data, model, prior);

  const arma::urowvec delta = model.delta.row(0);
  const unsigned int pd = arma::accu(delta);
  const arma::uvec delta_one = arma::find(delta == 1);

  if (simpost) {
    model.par.sigma = 1.0 / gamrnd<arma::rowvec>(post.sigma.c, 1 / post.sigma.C).t();
    model.par.beta.zeros(1, data.r);
    if (pd > 0) {
      model.par.beta(arma::uvec({0}), delta_one) = arma::trans(post.beta.b + arma::sqrt(model.par.sigma(delta_one) * post.beta.B) % randn(pd));
    }
    if (model.par.beta.has_nan()) {
      Rcpp::stop("NaN in beta");
    }
  }
  return post;
}

struct Facload_joint_help {
  double frac;
  arma::vec cT, CT, mc;
  arma::uvec line_nz, index_col, deltared_one;
  Eigen::SparseMatrix<double> L;
  Eigen::SparseMatrix<double> L_inv;
};

// Helper function for joint sampling of factor loadings and idiosyncratic variances
// Implements core computations for Algorithm P (bfa_version_archive.tex, line 3017)
// Specifically handles steps P-c1 and P-c2 for nonzero rows
Facload_joint_help help_facload_joint(
    const Data& data,
    const Model& model,
    const Prior& prior,
    const arma::uvec& sp) {
  // Identify nonzero rows and columns in the sparse loading matrix
  const arma::uvec line_nz = arma::find(sp);
  const arma::uvec colj = arma::any(model.delta.head_rows(model.d) != 0, 1);
  const arma::uvec index_col = arma::find(colj);
  
  // Construct the factor matrix X^δ and compute F'F (corresponds to X'X in regression)
  // This is part of building the information matrix in Step P-c1
  const arma::mat data1X = data.X.rows(index_col);
  arma::mat FtF = data1X * data1X.t();
  if (prior.type == Prior::Type::UNIT) {
    FtF.diag() += prior.par.beta.B0inv;  // Add prior precision for standard prior (eq 2948)
  }
  //const unsigned int pda = arma::accu(colj) * arma::accu(sp);
  const arma::uvec deltared_one = arma::find(model.delta(index_col, line_nz).as_col());  // TODO is this rowwise or columnwise in matlab?
  const unsigned int FtF_dim = FtF.n_cols;

  // Step P-c1: Construct the information matrix Ω (Omega) - eq (3028-3039)
  // Ω is block diagonal with blocks corresponding to rows i_1, ..., i_n
  // Each block contains (V^δ_{i,T})^{-1} from the posterior (eq 2949 or 2958)
  // The matrix is sparse with bandwidth = max(q_i)
  
  // Count elements per block to optimize sparse matrix construction
  std::vector<unsigned int> block_counts;
  unsigned int sum_squares = 0;
  {
    block_counts.reserve(100);  // should be enough
    unsigned int i = 0;
    while (i < deltared_one.n_elem) {
      unsigned int count = 1;
      const unsigned int i_block = deltared_one(i) / FtF_dim;
      for (i++; i < deltared_one.n_elem and deltared_one(i) / FtF_dim == i_block; i++) {
        count++;
      }
      block_counts.push_back(count);
      sum_squares += count * count;
    }
  }
  
  // Build the sparse information matrix Ω using triplets for efficiency
  Eigen::SparseMatrix<double> Omega(deltared_one.n_elem, deltared_one.n_elem);
  {
    std::vector<Eigen::Triplet<double, unsigned int> > triplets;
    triplets.reserve(sum_squares);
    for (unsigned int i = 0, j = 0; i < deltared_one.n_elem; i += block_counts[j], j++) {
      for (unsigned int ii = i; ii < i + block_counts[j] and ii < deltared_one.n_elem; ii++) {
          {  // diagonal
            const unsigned int irow = deltared_one(ii) % FtF_dim;
            const unsigned int icol = irow;
            triplets.push_back({ii, ii, FtF(irow, icol)});
          }
          for (unsigned int jj = i; jj < ii; jj++) {  // off-diagonal
            const unsigned int irow = deltared_one(ii) % FtF_dim;
            const unsigned int icol = deltared_one(jj) % FtF_dim;
            triplets.push_back({ii, jj, FtF(irow, icol)});
            triplets.push_back({jj, ii, FtF(icol, irow)});
          }
      }
    }
    Omega.setFromTriplets(triplets.cbegin(), triplets.cend());
  }

  // Step P-c1: Construct the covector c - eq (3040-3046)
  // c is stacked from c^δ_{i_l,T} = X'^δ_i * y_i for each nonzero row (eq 2950 or 2959)
  //data1X.print("\ndata1X");
  //data.y.rows(line_nz).print("y");
  Eigen::VectorXd cv(deltared_one.n_elem);
  //deltared_one.t().print("deltared_one");
  for (unsigned int i = 0; i < cv.size(); i++) {  // fill cv
    const unsigned int irow = deltared_one(i) % index_col.n_elem;
    const unsigned int icol = deltared_one(i) / index_col.n_elem;
    cv(i) = arma::dot(data1X.row(irow), data.y.row(line_nz(icol)));
    //Rcout << i << ' ' << irow << ' ' << icol << std::endl;
  }
  //Rcout << "cv\n" << cv.transpose() << std::endl;

  // Step P-c2: Compute Cholesky decomposition Ω = L·L^T (eq 3052)
  // Uses specialized algorithm for band matrices
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int> > Omega_solver;
  Omega_solver.compute(Omega);
  if (Omega_solver.info() != Eigen::Success) {
    stop("Cholesky failed");
  }
  const Eigen::SparseMatrix<double> L {Omega_solver.matrixL()};
#ifndef NDEBUG
  const Eigen::MatrixXd L_print {L};
#endif
  
  // Compute L^{-1} for later use in sampling (Step P-c4)
  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > L_inverter;
  L_inverter.compute(L);
  Eigen::SparseMatrix<double> I(Omega.rows(), Omega.cols());
  I.setIdentity();
  const Eigen::SparseMatrix<double> L_inv = L_inverter.solve(I);

  // Step P-c2: Solve L·m = c for m using triangular solver (eq 3054)
  Eigen::VectorXd mc_eigen = L_inv * cv;
  const arma::vec mc(mc_eigen.data(), mc_eigen.size(), false);
  //Rcout << "\nmc_eigen\n" << mc_eigen.transpose() << std::endl;
  //mc.t().print("mc");

  // Step P-c3: Compute posterior moments for idiosyncratic variances
  // Use m^T·m to efficiently compute C^δ_{i,T} via eq (3060): m^T_{i_l}·m_{i_l} = c^δ'_{i_l,T}·V^δ_{i_l,T}·c^δ_{i_l,T}
  arma::mat xi(arma::accu(colj), arma::accu(sp), arma::fill::zeros);
  xi(deltared_one) = mc;
  const arma::vec mtBm = arma::trans(arma::sum(arma::square(xi)));
  arma::vec cT = prior.par.sigma.c(line_nz) + 0.5 * data.N;
  double frac;

  // Compute posterior parameters c_T and C^δ_{i,T} for σ_i (eq 2951-2952 or 2960-2961)
  switch (prior.type) {
    case Prior::Type::FRAC: {
      cT -= 0.5 * data.N * prior.bfrac;  // Fractional prior adjustment (eq 2960)
      frac = 1 - prior.bfrac;
      break;
    }
    case Prior::Type::UNIT: {
      frac = 1;  // Standard prior (eq 2951)
      break;
    }
    default:
      stop("Unknown type of prior in step (P)");
      break;
  }
  // Compute C_{i,T} using SSR_i = y'y - c'·V·c (eq 2952 or 2961)
  const arma::vec CT = prior.par.sigma.C(line_nz) + 0.5 * frac * (arma::sum(arma::square(data.y.rows(line_nz)), 1) - mtBm);
  if (arma::any(CT <= 0)) {
    Rcpp::stop("negative value in CT");
  }
  return {frac, cT, CT, mc, line_nz, index_col, deltared_one, L, L_inv};
}

// Step (P): Sample factor loadings and idiosyncratic variances jointly
// Implements Algorithm P (bfa_version_archive.tex, line 3017) for sparse Bayesian factor models
// Handles three types of rows: zero rows (P-a), dedicated rows (P-b), and general nonzero rows (P-c)
double sample_facload_joint(
    const Data& data,
    Model& model,
    const Prior& prior) {
  const arma::uvec sd = arma::trans(arma::all(model.delta.head_cols(data.r) == 0));
  double SSBeta = 0;

  const bool TODO_SET_TRUE = true;

  // Step P-a: Sample idiosyncratic variances for zero rows (eq 2936)
  // For rows with q_i = 0, only σ_i needs to be sampled from the "null" model
  if (arma::any(sd)) {
    const arma::uvec line_z = arma::find(sd);
    const Data data1 {false, data.N, line_z.n_elem, data.y.rows(line_z), {}, {}, {}, {}, {}};
    Model model1;
    model1.d = 1;
    model1.delta.zeros(model1.d, data1.r);
    const Prior prior1 {prior.GLT_RJMCMC, prior.lifix, prior.varsel, prior.triangular,
      prior.a0, prior.b0, prior.bfrac, prior.itype, prior.typetau, prior.type,
      {{prior.par.sigma.c(line_z).t(), prior.par.sigma.C(line_z).t()}, {prior.par.beta.B0inv}}};
    sample_facload_one(data1, model1, prior1, true);
    if (TODO_SET_TRUE) {
      model.par.sigma(line_z) = model1.par.sigma;
    }
  }

  // Step P-c: Handle nonzero rows (q_i > 0) using joint block sampling
  // This includes dedicated rows (P-b, q_i = 1) and rows with multiple loadings (q_i > 1)
  const arma::uvec sp = 1u - sd;
  if (arma::any(sp)) {
    const Facload_joint_help help_out = help_facload_joint(data, model, prior, sp);
    const arma::vec &cT = help_out.cT,
                    &CT = help_out.CT,
                    &mc = help_out.mc;
    const arma::uvec &line_nz = help_out.line_nz,
                     &index_col = help_out.index_col,
                     &deltared_one = help_out.deltared_one;
    const Eigen::SparseMatrix<double> &L_inv = help_out.L_inv;
    // if (simpost)
    //const arma::vec sigmasim = 1 / gamrnd<arma::vec>(cT, 1 / CT);  # TODO uncomment
    
    // Step P-c3: Sample idiosyncratic variances σ_{i_1}, ..., σ_{i_n} from inverted Gamma (eq 2942-2943)
    arma::vec sigmasim;
    if (TODO_SET_TRUE) {
      sigmasim = 1 / gamrnd<arma::vec>(cT, 1 / CT);
    } else {
      sigmasim = model.par.sigma(line_nz);
    }
    model.par.sigma(line_nz) = sigmasim;
    // Step P-c4: Sample factor loadings Λ^δ from the joint posterior (eq 3064-3069)
    // Draw z ~ N(0, D) where D = diag(σ_{i_1}·1_{1×q_{i_1}}, ..., σ_{i_n}·1_{1×q_{i_n}})
    const arma::vec z0 = randn(deltared_one.n_elem);

    // Construct m + z with appropriate scaling by idiosyncratic variances
    arma::vec scale(deltared_one.n_elem);
    Eigen::VectorXd draws_help(deltared_one.n_elem);
    for (unsigned int i = 0; i < draws_help.size(); i++) {
      scale(i) = std::sqrt(sigmasim(deltared_one(i) / index_col.n_elem));
      draws_help(i) = mc(i) + scale(i) * z0(i);
    }
    // Solve L^T·Λ^δ = m + z for Λ^δ (eq 3067)
    // This gives a draw from Λ^δ | σ_{i_1}, ..., σ_{i_n}, y, F
    const Eigen::VectorXd betasim = L_inv.transpose() * draws_help;
    //cT.t().print("\ncT");
    //CT.t().print("\nCT");
    //sigmasim.t().print("\nsigmasim");
    //scale.t().print("\nscale");
    //z0.t().print("\nz0");
    //Rcout << "\ndraws_help\n" << draws_help.transpose() << std::endl;
    //Rcout << "\nbetasim\n" << betasim.transpose() << std::endl;
    //model.par.beta.print("\nbeta");

    // Unpack the stacked vector Λ^δ back into the loading matrix
    model.par.beta.zeros();
    for (unsigned int i = 0; i < betasim.size(); i++) {  // fill cv
      const unsigned int irow = deltared_one(i) % index_col.n_elem;
      const unsigned int icol = deltared_one(i) / index_col.n_elem;
      model.par.beta(index_col(irow), line_nz(icol)) = betasim(i);
    }
    if (model.par.beta.has_nan()) {
      Rcpp::stop("NaN in beta");
    }
    SSBeta = arma::dot(z0, z0);
    //model.par.beta.print("\nbeta");
  }
    //stop("hehehe");
  return SSBeta;
}

// Compute posterior distributions for marginal likelihood calculations
// Used in Steps (L) and (D) for model selection with respect to δ (Section 3.8.4, line 3108)
// Returns determinants and log-ratios needed for marginal likelihood computation
Facload_joint_out posterior_facload_joint(
    const Data& data,
    const Model& model,
    const Prior& prior) {
  Facload_joint_out post;
  const arma::uvec sp = arma::trans(arma::any(model.delta.head_cols(data.r) != 0));
  if (arma::any(sp)) {
    const Facload_joint_help help_out = help_facload_joint(data, model, prior, sp);
    const double &frac = help_out.frac;
    const arma::vec &cT = help_out.cT,
                    &CT = help_out.CT,
                    &mc = help_out.mc;
    const arma::uvec &index_col = help_out.index_col,
                     &line_nz = help_out.line_nz,
                     &deltared_one = help_out.deltared_one;
    const Eigen::SparseMatrix<double> &L = help_out.L;
    post.cT = cT.t();
    post.betadetl.set_size(line_nz.n_elem);
    post.betadetl.fill(arma::datum::inf);
    post.sigmalogcr.zeros(CT.n_elem);
    // Extract diagonal elements of L for computing determinants
    // For marginal likelihood computation (eq 3136 or 3190), we need |V^δ_{i,T}|
    // Since L_{i_l} is the Cholesky of (V^δ_{i_l,T})^{-1}, we have |V^δ_{i_l,T}| = 1/|L_{i_l}|^2
    for (unsigned int i = 0; i < deltared_one.n_elem; i++) {
      const unsigned int e = deltared_one(i);
      if (e % index_col.n_elem == index_col.n_elem - 1) {
        post.betadetl(e / index_col.n_elem) = -log(L.coeff(i, i));
        post.sigmalogcr(e / index_col.n_elem) = log(CT(e / index_col.n_elem) + 0.5 * frac * mc(i) * mc(i)) - log(CT(e / index_col.n_elem));
      }
    }
  } else {
    post.betadetl.zeros(data.r);
    post.sigmalogcr.zeros(data.r);
  }
  return post;
}

