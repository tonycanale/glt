#include <RcppArmadillo.h>
#include "Rcpp/Rmath.h"
#include "Rinternals.h"
#include "util.h"
#include "delta.h"
#include "havard.h"
#include "leading_elements.h"
#include "structs.h"
#include "steps.h"

namespace PBFA {
  extern Mcmcout mcmcout;
}

#define NOT_YET_IMPLEMENTED ::Rf_error("nyi")

// TODO remove '_delta' and most 'sample' words from function names

using namespace Rcpp;

// [[Rcpp::export]]
List delta_li_sort_r(
    const arma::umat& delta,
    const CharacterVector& ide_in) {
  const Delta_li_sort_out li_out = delta_li_sort(delta, Charv2ide(ide_in));
  return List::create(
      _["li_glt"] = wrap(li_out.li_glt),
      _["li"] = wrap(li_Cpp2R(li_out.li)),
      _["lisort"] = wrap(li_Cpp2R(li_out.lisort)),
      _["bisort"] = wrap(li_Cpp2R(li_out.bisort)),
      _["is"] = wrap(li_Cpp2R(li_out.is)));
}

struct Delta_odds_likelihood {
  arma::vec odd, odd_lik;
};

Delta_odds_likelihood compute_delta_odds_likelihood(
    const Data& data,
    const Model& model,
    const Prior& prior,
    const unsigned int j,
    const arma::uvec& ichange);

// [[Rcpp::export]]
List compute_delta_odds_likelihood_r(
    const List& data,
    const List& model,
    const List& prior,
    const unsigned int j,
    const arma::uvec& ichange_in) {
  const auto out = compute_delta_odds_likelihood(list2data(data), list2model(model), list2prior(prior), j - 1u, ichange_in - 1u);
  return List::create(
      _["odd"] = wrap(out.odd),
      _["odd_lik"] = wrap(out.odd_lik));
}

double compute_delta_odds_prior(
    const double tau);

// [[Rcpp::export]]
NumericVector compute_delta_odds_prior_r(
    const NumericVector& taus) {  // vectorized variant
  NumericVector prior_odds(taus.length());
  std::transform(taus.cbegin(), taus.cend(), prior_odds.begin(), compute_delta_odds_prior);
  return prior_odds;
}

double compute_delta_move_lead_updown_log_accept_rate(
    const Data& data,
    const Model& model,
    const Prior& prior,
    const unsigned int j,
    const arma::uvec& li,
    const arma::uvec& li_add_in);

// [[Rcpp::export]]
double compute_delta_move_lead_updown_log_accept_rate_r(
    const List& data,
    const List& model,
    const List& prior,
    const unsigned int j,
    const IntegerVector& li,
    const IntegerVector& li_new) {
  return compute_delta_move_lead_updown_log_accept_rate(list2data(data), list2model(model), list2prior(prior), j - 1u, li_R2Cpp(li), li_R2Cpp(li_new));
}

arma::vec compute_delta_flip_log_accept_rate(  // 'flip' as in change from 0 to 1 or vice versa
    const arma::urowvec& delta_j_ichange,
    const arma::vec& odd_ichange,
    const double qratio);

// [[Rcpp::export]]
NumericVector compute_delta_flip_log_accept_rate_r(
    const IntegerVector& delta_j_ichange_in,
    const NumericVector& odd_ichange_in,
    const double qratio) {
  arma::urowvec delta_j_ichange (delta_j_ichange_in.length());
  std::transform(delta_j_ichange_in.cbegin(), delta_j_ichange_in.cend(), delta_j_ichange.begin(),
      [](const int in) -> unsigned int { return in; });
  NumericVector result = wrap(compute_delta_flip_log_accept_rate(delta_j_ichange, Rcpp::as<arma::vec>(odd_ichange_in), qratio));
  result.attr("dim") = R_NilValue;  // remove dim attribute
  return result;
}

double compute_delta_switch_log_accept_rate(  // 'switch' means to switch two columns (between and including their leading indices)
    const Data& data,
    const Model& model,
    const Prior& prior,
    const unsigned int j1,
    const unsigned int j2,
    const arma::uvec& ichange,
    const arma::uvec& li);

// [[Rcpp::export]]
double compute_delta_switch_log_accept_rate_r(
    const List& data,
    const List& model,
    const List& prior,
    const unsigned int j1,
    const unsigned int j2,
    const IntegerVector& ichange,
    const IntegerVector& li) {
  return compute_delta_switch_log_accept_rate(list2data(data), list2model(model), list2prior(prior), j1 - 1u, j2 - 1u, as<arma::uvec>(ichange) - 1u, li_R2Cpp(li));
}

void sample_delta_move_lead_updown(
    const Data& data,
    Model& model,
    const Prior& prior,
    const unsigned int j,
    const unsigned int S_extend,
    const arma::uvec& li);

// [[Rcpp::export]]
List sample_delta_move_lead_updown_r(
    const List& data,
    const List& model_in,
    const List& prior,
    const unsigned int j,
    const unsigned int S_extend,
    const IntegerVector& li) {
  Model model = list2model(model_in);
  sample_delta_move_lead_updown(list2data(data), model, list2prior(prior), j - 1u, S_extend, li_R2Cpp(li));
  return model2list(model);
}

void sample_delta_switch(
    const Data& data,
    Model& model,
    const Prior& prior,
    const unsigned int j1,
    const unsigned int j2,
    const arma::uvec& ichange,
    const arma::uvec& li);

// [[Rcpp::export]]
List sample_delta_switch_r(
    const List& data,
    const List& model_in,
    const List& prior,
    const unsigned int j1,
    const unsigned int j2,
    const IntegerVector& ichange,
    const IntegerVector& li) {
  Model model = list2model(model_in);
  sample_delta_switch(list2data(data), model, list2prior(prior), j1 - 1u, j2 - 1u, as<arma::uvec>(ichange) - 1u, li_R2Cpp(li));
  return model2list(model);
}

// [[Rcpp::export]]
List sample_delta_update_row_r(
    const List& data,
    const List& model_in,
    const List& prior,
    const unsigned int j,
    const IntegerVector& ichange,
    const double qratio) {
  Model model = list2model(model_in);
  sample_delta_update_row(list2data(data), model, list2prior(prior), j - 1u, as<arma::uvec>(ichange) - 1u, qratio);
  return model2list(model);
}

void sample_delta_add_delete_leading_element(
    const Data& data,
    Model& model,
    const Prior& prior,
    const unsigned int j,
    const unsigned int S_extend,
    const arma::uvec& li,
    const double prob_add_no_spurious,
    const bool extended_uniqueness);

// [[Rcpp::export]]
List sample_delta_add_delete_leading_element_r(
    const List& data,
    const List& model_in,
    const List& prior,
    const unsigned int j,
    const unsigned int S_extend,
    const IntegerVector& li,
    const double prob_add_no_spurious,
    const bool extended_uniqueness) {
  Model model = list2model(model_in);
  sample_delta_add_delete_leading_element(list2data(data), model, list2prior(prior), j - 1u, S_extend, li_R2Cpp(li), prob_add_no_spurious, extended_uniqueness);
  return model2list(model);
}

// [[Rcpp::export]]
List sample_delta_move_lead_r(
    const List& data,
    const List& model_in,
    const List& prior,
    const unsigned int j,
    const unsigned int S_extend,
    const IntegerVector& li,  // TODO remove these `li` arguments
    const double prob_switch,
    const double prob_move,
    const double prob_add,
    const bool extended_uniqueness) {
  Model model = list2model(model_in);
  sample_delta_move_lead(list2data(data), model, list2prior(prior), j, S_extend, li_R2Cpp(li), prob_switch, prob_move, prob_add, extended_uniqueness);
  return model2list(model);
}

// C++ implementations

double compute_delta_odds_prior(
    const double tau) {
  return std::log(tau / (1 - tau));
}

arma::vec compute_delta_flip_log_accept_rate(
    const arma::urowvec& delta_j_ichange,
    const arma::vec& odd_ichange,
    const double qratio) {
  // the two cases in the Matlab code are the same; I implement
  // the general case when `size(ichange,1)` can be anything
  arma::vec result (arma::size(odd_ichange));
  std::transform(delta_j_ichange.cbegin(), delta_j_ichange.cend(), result.begin(),
      [](const arma::uword is_one) -> double { return 1. - 2. * is_one; });
  result %= odd_ichange;
  result += qratio;
  return result;
}

Delta_li_sort_out delta_li_sort(
    const arma::umat& delta,
    const Model::Ide ide) {
  arma::uvec li;
  arma::uvec lisort;
  arma::uvec bisort;
  arma::uvec is;
  bool li_glt;

  const unsigned int nm = delta.n_cols;  // number of observations
  const unsigned int kmax = delta.n_rows;  // number of active factors
  const arma::uvec nonzero_columns = arma::sum(delta, 1) != 0;
  const unsigned int nnz = arma::accu(nonzero_columns);  // number of non-zero columns

  if (nnz > 0) {
    li.set_size(kmax);
    bisort.set_size(arma::size(li));
    for (unsigned int i = 0; i < li.n_elem; i++) {
      if (nonzero_columns(i)) {
        li(i) = arma::find(delta.row(i), 1).min();
        bisort(i) = arma::find(delta.row(i)).max();
      } else {
        li(i) = bisort(i) = UWORD_NAN;
      }
    }
    is = arma::sort_index(li);  // UWORD_NAN is larger than other values
    lisort = li(is.head(nnz));
    bisort = bisort(is.head(nnz));
    li_glt = arma::all(arma::diff(lisort) > 0);
    if (li_glt and ide != Model::Ide::NO_BOUND) {
      li_glt = li_check(lisort, nm, ide);
    }
  } else {
    li.set_size(kmax);
    li.fill(UWORD_NAN);
    li_glt = false;
  }

  return {li_glt, li, lisort, bisort, is};
}

double compute_delta_move_lead_updown_log_accept_rate(
    const Data& data,
    const Model& model,
    const Prior& prior,
    const unsigned int j,
    const arma::uvec& li,  // TODO eliminate; li and model.li are always the same
    const arma::uvec& li_new) {
  const arma::vec odd_lik = compute_delta_odds_likelihood(data, model, prior, j, {li_new(j), li(j)}).odd_lik;
  return odd_lik(0) - odd_lik(1);
}

double compute_delta_switch_log_accept_rate(
    const Data& data,
    const Model& model,
    const Prior& prior,
    const unsigned int j1,
    const unsigned int j2,
    const arma::uvec& ichange,
    const arma::uvec& li) {
  // compute odds with a 'hypothetical' delta
  Model model_new {model};
  model_new.delta.elem(arma::uvec{j1}, ichange).zeros();
  model_new.delta.elem(arma::uvec{j2}, ichange).zeros();
  const arma::vec odd_j1 = compute_delta_odds_likelihood(data, model_new, prior, j1, ichange).odd(ichange) +
    compute_delta_odds_prior(model.pitau(j1));
  const arma::vec odd_j2 = compute_delta_odds_likelihood(data, model_new, prior, j2, ichange).odd(ichange) +
    compute_delta_odds_prior(model.pitau(j2));

  // return to the original delta
  // TODO just accumulate instead of extra memory allocation
  arma::vec odd_mh(ichange.n_elem, arma::fill::zeros);
  odd_mh = odd_j2 - odd_j1;
  odd_mh(arma::find(model.delta(arma::uvec{j2}, ichange))) *= -1;

  return arma::accu(odd_mh);
}

// TODO find better name
// Runs li_add but excludes the jth index,
// so it determines the set of all feasible
// leading elements if the jth column is removed.
// Then, it finds feasible leading elements that
// are smaller (so above) than li_boundary.
arma::uvec li_add_exclude(
    const arma::uvec& li,
    const unsigned int j,
    const unsigned int nnz,
    const Model& model) {
  const arma::uvec result = li_add(
      li(arma::find(arma::regspace<arma::uvec>(0, model.delta.n_rows - 1u) != j)),
      nnz, model.ide, true).liadd;
  return result(arma::find(result != UWORD_NAN));
}

void sample_delta_move_lead_updown(
    const Data& data,
    Model& model,
    const Prior& prior,
    const unsigned int j,
    const unsigned int S_extend,
    const arma::uvec& li) {
  // not tested explicitly

  // beginning
  PBFA::mcmcout.mh_lead.down.tryy(li(j))++;  // TODO uncomment
  const unsigned int li_star =
    arma::accu(model.delta.row(j) == 1u) == 1u ?
    UWORD_NAN :  // spurious column
    arma::find(model.delta.row(j), 2).max();  // row index of second nonzero element in column j
  const arma::uvec li_add_feasible = li_add_exclude(li, j, data.r - S_extend, model);
  const arma::uvec li_add_v = li_add_feasible(arma::find(li_add_feasible < li_star));  // UWORD_NAN is larger than all possible values for li

  if (li_add_v.n_elem > 0) {
    // call compute_delta_move_lead_updown_ratio
    const arma::uvec li_new =
      ([](const arma::uvec& _li,
          const arma::uvec& _li_add_in,
          const unsigned int _j) -> arma::uvec {
       arma::uvec result {_li};
       result(_j) = _li_add_in[0];
       return result;
      })(li, li_add_v, j);
    const double log_ar = compute_delta_move_lead_updown_log_accept_rate(data, model, prior, j, li, li_new);

    // do MH
    if (log_ar >= 0 or std::log(R::unif_rand()) < log_ar) {
      PBFA::mcmcout.mh_lead.down.acc(li(j))++;  // TODO uncomment
      model.delta =
        ([](const arma::umat _delta,
            const unsigned int _li_j,
            const unsigned int _li_new_j,
            const unsigned int _j) -> arma::umat {
         arma::umat result {_delta};
         result(_j, _li_j) = 0u;
         result(_j, _li_new_j) = 1u;
         return result;
        })(model.delta, li(j), li_new(j), j);
      model.li = li_new;
    }
  }
}

void sample_delta_update_row(
    const Data& data,
    Model& model,
    const Prior& prior,
    const unsigned int j,
    const arma::uvec& ichange,
    const double qratio) {
  // not tested explicitly

  const arma::vec odd_ichange =
    compute_delta_odds_likelihood(data, model, prior, j, ichange).odd_lik +
    compute_delta_odds_prior(model.pitau(j));
  const arma::vec log_ar = compute_delta_flip_log_accept_rate(model.delta(arma::uvec{j}, ichange), odd_ichange, qratio);
  // do MH
  const arma::uvec to_update =
    ([](const arma::vec& _log_ar) -> arma::uvec {
      arma::uvec indicator(_log_ar.n_elem);
      std::transform(_log_ar.cbegin(), _log_ar.cend(), indicator.begin(),
          [](const double el) -> arma::uword {
            return el >= 0 or R::unif_rand() <= std::exp(el);
          });
      return arma::find(indicator);
    })(log_ar);
  model.delta(arma::uvec{j}, ichange(to_update)) = 1u - model.delta(arma::uvec{j}, ichange(to_update));
}

void sample_delta_switch(
    const Data& data,
    Model& model,
    const Prior& prior,
    const unsigned int j1,
    const unsigned int j2,
    const arma::uvec& ichange,
    const arma::uvec& li) {
  // not tested explicitly
  PBFA::mcmcout.mh_switch.tryy++;  // TODO uncomment

  const double log_ar = compute_delta_switch_log_accept_rate(data, model, prior, j1, j2, ichange, li);
  // do MH
  if (log_ar >= 0 or R::unif_rand() <= std::exp(log_ar)) {
    PBFA::mcmcout.mh_switch.acc++;  // TODO uncomment
    model.delta(arma::uvec{j1, j2}, ichange) = 1u - model.delta(arma::uvec{j1, j2}, ichange);
    model.li(j1) = arma::as_scalar(arma::find(model.delta.row(j1), 1));
    model.li(j2) = arma::as_scalar(arma::find(model.delta.row(j2), 1));
  }
}

void sample_delta_add_delete_leading_element(
    const Data& data,
    Model& model,
    const Prior& prior,
    const unsigned int j,
    const unsigned int S_extend,
    const arma::uvec& li,
    const double prob_add_no_spurious,
    const bool extended_uniqueness) {
  const arma::uvec li_add_feasible = li_add_exclude(li, j, data.r - S_extend, model);

  // Check if add is possible
  const arma::uvec index_above = li_add_feasible(arma::find(li_add_feasible < li(j)));  // `iup` in Matlab
  const bool possible_add = not index_above.empty();

  // Check if delete is possible
  //   3 elements are enough; `idown(1)` is `li(j)` and we don't care about it
  const arma::uvec index_below = arma::find(model.delta.row(j), 3);  // `idown` in Matlab
  const bool possible_delete =
    index_below.n_elem != 1u and
    li_add_feasible.cend() != std::find(li_add_feasible.cbegin(), li_add_feasible.cend(), index_below(1)) and
    // number of spurious columns == sum(sum(model.delta, 2) == 1), no need for `jall` in Matlab
    (index_below.n_elem > 2u or not extended_uniqueness or arma::accu(arma::sum(model.delta, 1) == 1u) < S_extend);

  double prob_add;
  bool add, del;
  if (possible_add and possible_delete) {
    prob_add = prob_add_no_spurious;
    add = R::unif_rand() <= prob_add;
    del = not add;
  } else if (possible_add) {
    prob_add = 1;
    add = true;
    del = false;
  } else if (possible_delete) {
    prob_add = 0;
    add = false;
    del = true;
  } else {
    add = del = false;
  }

  //add = true;  // TODO delete!
  if (add) {
    PBFA::mcmcout.mh_lead.add.tryy(model.li(j))++;
    const double qratio = std::log(index_above.n_elem * (1. - prob_add_no_spurious) / prob_add),
                 odd_lik = arma::as_scalar(compute_delta_odds_likelihood(data, model, prior, j, {index_above(0)}).odd_lik),
                 log_ar = odd_lik + qratio + compute_delta_odds_prior(model.pitau(j));
    if (log_ar >= 0 or R::unif_rand() <= std::exp(log_ar)) {
      PBFA::mcmcout.mh_lead.add.acc(model.li(j))++;
      model.delta(j, index_above(0)) = 1u;
      model.li(j) = index_above(0);
    }
  } else if (del) {
    PBFA::mcmcout.mh_lead.del.tryy(model.li(j))++;
    const arma::uvec li_add_reverse = li_add_exclude(li, j, data.r - S_extend, model);
    const int n_index_above_reverse = arma::accu(li_add_reverse < index_below(1));
    const double prob_add_reverse = index_below.n_elem == 2 ? 1. : prob_add_no_spurious,
                 qratio = std::log(prob_add_reverse / (1. - prob_add) / n_index_above_reverse),
                 odd_lik = arma::as_scalar(compute_delta_odds_likelihood(data, model, prior, j, {model.li(j)}).odd_lik),
                 log_ar = -(odd_lik + compute_delta_odds_prior(model.pitau(j))) + qratio;
    if (log_ar >= 0 or R::unif_rand() <= std::exp(log_ar)) {
      PBFA::mcmcout.mh_lead.del.acc(model.li(j))++;
      model.delta(j, model.li(j)) = 0u;
      model.li(j) = index_below(1);
    }
  }
}

void sample_delta_move_lead(
    const Data& data,
    Model& model,
    const Prior& prior,
    const unsigned int j,
    const unsigned int S_extend,
    const arma::uvec& li,  // TODO remove these `li` arguments
    const double prob_switch,
    const double prob_move,
    const double prob_add,
    const bool extended_uniqueness) {
  const bool move_lead = R::unif_rand() < prob_move,
             switch_lead = not move_lead and R::unif_rand() < prob_switch / (1 - prob_move),
             add_del = not (move_lead or switch_lead);

  if (move_lead) {
    //const arma::uvec index_below = arma::find(model.delta.row(j), 2);  // this was a misunderstanding
    //if (index_below.n_elem > 1 and index_below(1) > 1u) {
      sample_delta_move_lead_updown(data, model, prior, j, S_extend, model.li);
    //}
  } else if (switch_lead) {
    const arma::uvec j_all_nonzero = arma::find(model.li != UWORD_NAN);
    if (j_all_nonzero.n_elem > 1u) {
      unsigned int j2_index = Rcpp::sample(j_all_nonzero.n_elem - 1, 1, false, ::R_NilValue, false)(0);
      const unsigned int j2 =
        j_all_nonzero(j2_index) >= j ?
        j_all_nonzero(j2_index + 1u) :
        j_all_nonzero(j2_index);
      const unsigned int li_min = std::min(model.li(j), model.li(j2)),
                         li_max = std::max(model.li(j), model.li(j2));
      const arma::uvec ichange = li_min +
        arma::find(model.delta(j, arma::span(li_min, li_max)) != model.delta(j2, arma::span(li_min, li_max)));
      if (not ichange.empty()) {
        sample_delta_switch(data, model, prior, j, j2, ichange, model.li);
      }
    }
  } else if (add_del) {
    sample_delta_add_delete_leading_element(data, model, prior, j, S_extend, model.li, prob_add, extended_uniqueness);
  }
}

Delta_odds_likelihood compute_delta_odds_likelihood(
    const Data& data,
    const Model& model,
    const Prior& prior,
    const unsigned int j,
    const arma::uvec& ichange) {
  const bool joint_sample_delta = ichange.n_elem > 1u;
  arma::vec odd(data.r, arma::fill::zeros);

  if (joint_sample_delta) {
    const arma::urowvec sd = arma::all(model.delta(arma::find(arma::regspace<arma::uvec>(0, model.d - 1u) != j), ichange) == 0, 0);
    {
      const arma::uvec line_z = ichange(arma::find(sd));
      if (line_z.n_elem > 0) {
        const Data dataS {data.sim, data.N, line_z.n_elem, data.y.rows(line_z), data.X.row(j)};
        Model modelS;
        modelS.d = 1;
        modelS.delta.ones(modelS.d, dataS.r);
        const Prior priorS {prior.GLT_RJMCMC, prior.lifix, prior.varsel, prior.triangular,
          prior.a0, prior.b0, prior.bfrac, prior.itype, prior.typetau, prior.type,
          {{prior.par.sigma.c(line_z).t(), prior.par.sigma.C(line_z).t()}, {prior.par.beta.B0inv}}};
        {
          const Facload_one_out post = sample_facload_one(dataS, modelS, priorS);
          odd(line_z) = arma::lgamma(post.sigma.c.t()) - post.sigma.c.t() % arma::log(post.sigma.C.t());
          switch(prior.type) {
            case Prior::Type::FRAC:
              odd(line_z) += 0.5 * dataS.N * prior.bfrac * log(2 * arma::datum::pi) + 0.5 * log(prior.bfrac);
              break;
            case Prior::Type::UNIT:
              odd(line_z) += 0.5 * log(prior.par.beta.B0inv * post.beta.B);
              break;
            default:
              stop("Unknown type of prior");
          }
        }
        modelS.delta.zeros();
        {
          const Facload_one_out post = sample_facload_one(dataS, modelS, priorS);
          odd(line_z) += -arma::lgamma(post.sigma.c.t()) + post.sigma.c.t() % arma::log(post.sigma.C.t());
        }
      }
    }

    {
      const arma::uvec line_nz = ichange(arma::find(1u - sd));
      if (line_nz.n_elem > 0) {
        arma::uvec new_idx(model.d);  // TODO why is this reordering needed?
        for (unsigned int i = 0; i < j; i++) {
          new_idx(i) = i;
        }
        for (unsigned int i = j+1; i < model.d; i++) {
          new_idx(i - 1) = i;
        }
        new_idx.tail(1) = j;
        const Data dataS {data.sim, data.N, line_nz.n_elem, data.y.rows(line_nz), data.X.rows(new_idx)};
        Model modelS;
        modelS.d = model.d;
        modelS.delta = model.delta(new_idx, line_nz);
        modelS.delta.tail_rows(1).fill(1);
        const Prior priorS {prior.GLT_RJMCMC, prior.lifix, prior.varsel, prior.triangular,
          prior.a0, prior.b0, prior.bfrac, prior.itype, prior.typetau, prior.type,
          {{prior.par.sigma.c(line_nz).t(), prior.par.sigma.C(line_nz).t()}, {prior.par.beta.B0inv}}};
        const Facload_joint_out post = posterior_facload_joint(dataS, modelS, priorS);
        odd(line_nz) = post.cT.t() % post.sigmalogcr.t();
        switch(prior.type) {
          case Prior::Type::FRAC:
            odd(line_nz) += 0.5 * log(prior.bfrac);
            break;
          case Prior::Type::UNIT:
            odd(line_nz) += post.betadetl.t() + 0.5 * log(prior.par.beta.B0inv);
            break;
          default:
            stop("Unknown type of prior");
        }
      }
    }
  } else {
    const unsigned int ichang = arma::as_scalar(ichange);
    Model modelnew(model);
    modelnew.delta(j, ichang) = 1;
    const Havard_out havard_out1 = mcfullcondk_havard(data, modelnew, prior, ichang);
    modelnew.delta(j, ichang) = 0;
    const Havard_out havard_out0 = mcfullcondk_havard(data, modelnew, prior, ichang);
    odd(ichang) = havard_out1.marlik - havard_out0.marlik;
  }
  return {odd, odd(ichange)};
}

#undef NOT_YET_IMPLEMENTED

