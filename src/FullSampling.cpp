#include <Rcpp.h>
#include "FFBS.h"
#include "CoalescenceBlocks.h"
#include "FullSampling.h"

Rcpp::List sample_bounded_times_c(Rcpp::NumericVector times,
                                Rcpp::IntegerVector leaves,
                                double ne,
                                double bound,
                                int nsam) {

  // Total the leaves
  int total_leaves = Rcpp::sum(leaves);

  // Storage for samples
  Rcpp::NumericVector sampled_times(nsam * (total_leaves - 1));
  sampled_times.attr("dim") = Rcpp::Dimension(nsam, (total_leaves - 1));

  // Storage for likelihood
  Rcpp::NumericVector likelihood(nsam, 1.0);

  Rcpp::NumericVector forward_probs(total_leaves * (times.size() + 1));
  forward_probs.attr("dim") = Rcpp::Dimension(total_leaves, (times.size() + 1));

  forward_algorithm_c(times, leaves, ne, bound, forward_probs);

  // Counter
  int c, n;

  // Double copy of starting lineages
  double dl;

  // Normalising constant
  double z;

  // Random numbers for sampling
  Rcpp::NumericVector u = Rcpp::runif(nsam * (total_leaves - 1), 0.0, 1.0);
  u.attr("dim") = Rcpp::Dimension(nsam, (total_leaves - 1));

  Rcpp::List back_sampler;
  Rcpp::IntegerVector lineages(times.size() + 1);

  Rcpp::List constraints;
  Rcpp::NumericVector const_lower(total_leaves - 1);
  Rcpp::NumericVector const_upper(total_leaves - 1);
  Rcpp::IntegerVector const_lineages(total_leaves - 1);
  Rcpp::IntegerVector const_events(total_leaves - 1);

  Rcpp::NumericVector single_sample(2);

  for (n = 0; n < nsam; ++n) {

    likelihood(n) *=
      backward_sampler_c(forward_probs, times, leaves, ne, bound, lineages);

    likelihood(n) *=
      constrain_coalescences_c(lineages, times, leaves, ne, bound, const_lower,
                               const_upper, const_lineages, const_events);

    for (c = 0; c < total_leaves - 1; ++c) {

      dl = double(const_lineages(c));

      z = (ne / (dl - 1.0)) * (1.0 - std::exp(((dl - 1.0) / ne) *
        (const_lower(c) - const_upper(c))));

      sampled_times(n, c) = const_upper(c) + (ne / (dl - 1.0)) *
        std::log(1.0 - ((dl - 1.0) / ne) * z * u(n, c));

      likelihood(n) *= (1.0 / z) * std::exp(((dl - 1.0) / ne) *
        (sampled_times(n, c) - const_upper(c)));

    }

  }

  Rcpp::List out = Rcpp::List::create(Rcpp::Named("times") = sampled_times,
                                      Rcpp::Named("likelihood") = likelihood);

  return out;

}


Rcpp::List rejection_bounded_times_c(Rcpp::NumericVector times,
                                   Rcpp::IntegerVector leaves,
                                   double ne,
                                   double bound,
                                   int nsam) {

  // Total the leaves
  int total_leaves = Rcpp::sum(leaves);

  // Storage for samples
  Rcpp::NumericVector sampled_times(nsam * (total_leaves - 1));
  sampled_times.attr("dim") = Rcpp::Dimension(nsam, (total_leaves - 1));

  // Counters
  int c, k, n = 0, total_sims = 0;

  // Exponential inputs
  double scale, dt;

  // Current states
  int current_lineages;
  double current_time;

  // Storage for likelihood
  Rcpp::NumericVector likelihood(nsam, 1.0);

  // Forward filter for normalising constant
  Rcpp::NumericVector forward_probs(total_leaves * (times.size() + 1));
  forward_probs.attr("dim") = Rcpp::Dimension(total_leaves, (times.size() + 1));

  forward_algorithm_c(times, leaves, ne, bound, forward_probs);

  while (n < nsam) {

    likelihood(n) = 1.0;

    c = total_leaves - 2;
    k = times.length() - 1;

    current_lineages = leaves[k];
    current_time = times[k];

    --k;

    while (c >= 0 && current_time > bound) {

      // Check for interrupt
      Rcpp::checkUserInterrupt();

      if (current_lineages > 1) {

        scale = (2.0 * ne) /
          (double(current_lineages) * (double(current_lineages) - 1.0));

        dt = R::rexp(scale);

        if (k < 0 || (current_time - dt) > times[k]) {

          likelihood(n) *= R::dexp(dt, scale, false);

          current_time = current_time - dt;
          --current_lineages;

          sampled_times(n, c) = current_time;
          --c;

        } else {

          likelihood(n) *=
            R::pexp(current_time - times[k], scale, false, false);

          current_time = times[k];
          current_lineages += leaves(k);

          --k;

        }

      } else {

        current_time = times[k];
        current_lineages += leaves(k);

        --k;

      }

    }

    likelihood(n) /= forward_probs(0, 0);

    ++total_sims;

    if (current_time > bound) {

      ++n;

    }

  }

  double ar = double(nsam) / double(total_sims);

  Rcpp::List out = Rcpp::List::create(Rcpp::Named("times") = sampled_times,
                                      Rcpp::Named("likelihood") = likelihood,
                                      Rcpp::Named("acceptance_rate") = ar);

  return out;

}


Rcpp::List sample_topology_c(Rcpp::NumericVector leaf_times,
                             Rcpp::IntegerVector leaves,
                             Rcpp::NumericVector coalescence_times) {

  // Total the leaves
  int total_leaves = Rcpp::sum(leaves);

  // Internal node ancestors
  Rcpp::IntegerVector edge(4 * (total_leaves - 1));
  edge.attr("dim") = Rcpp::Dimension(2 * (total_leaves - 1), 2);
  Rcpp::NumericVector edge_length(2 * (total_leaves - 1));

  // Node ancestors and times
  Rcpp::NumericVector node_times(2 * total_leaves - 1);

  // Counters
  int i, c, k;
  int l_index = total_leaves - 1, n_index = 2 * total_leaves - 2;

  // Current information
  Rcpp::IntegerVector active_nodes(2 * total_leaves - 1);
  int total_active_nodes = 0;

  // Ancestors
  int anc_1, anc_2;

  // Random numbers for sampling
  Rcpp::NumericVector u = Rcpp::runif(2 * (total_leaves - 1), 0.0, 1.0);

  // Cumulative probability
  double sum_prob;

  // Likelihood
  double likelihood = 1.0;

  k = leaf_times.size() - 1;

  for (i = 0; i < leaves(k); ++i) {

    active_nodes(l_index) = 1;
    ++total_active_nodes;

    node_times(l_index) = leaf_times(k);

    --l_index;

  }

  --k;

  c = coalescence_times.size() - 1;

  while (c >= 0) {

    // Check for interrupt
    Rcpp::checkUserInterrupt();

    if (k < 0 || leaf_times(k) < coalescence_times(c)) {

      anc_1 = 0;
      sum_prob = double(active_nodes(anc_1)) / double(total_active_nodes);

      while (sum_prob < u(2 * c + 1)) {

        ++anc_1;
        sum_prob += double(active_nodes(anc_1)) / double(total_active_nodes);

      }

      edge(2 * c + 1, 0) = n_index + 1;
      edge(2 * c + 1, 1) = anc_1 + 1;

      edge_length(2 * c + 1) = node_times(anc_1) - coalescence_times(c);

      likelihood *= 2.0 / double(total_active_nodes);

      active_nodes(anc_1) = 0;
      --total_active_nodes;


      anc_2 = 0;
      sum_prob = double(active_nodes(anc_2)) / double(total_active_nodes);

      while (sum_prob < u(2 * c)) {

        ++anc_2;
        sum_prob += double(active_nodes(anc_2)) / double(total_active_nodes);

      }

      edge(2 * c, 0) = n_index + 1;
      edge(2 * c, 1) = anc_2 + 1;

      edge_length(2 * c) = node_times(anc_2) - coalescence_times(c);

      likelihood *= 1.0 / double(total_active_nodes);

      active_nodes(anc_2) = 0;

      active_nodes(n_index) = 1;
      node_times(n_index) = coalescence_times(c);

      --c;
      --n_index;

    } else {

      for (i = 0; i < leaves(k); ++i) {

        active_nodes(l_index) = 1;
        ++total_active_nodes;

        node_times(l_index) = leaf_times(k);

        --l_index;

      }

      --k;

    }

  }

  Rcpp::List out = Rcpp::List::create(Rcpp::Named("edge") = edge,
                                      Rcpp::Named("edge_length") = edge_length,
                                      Rcpp::Named("likelihood") = likelihood);

  return out;

}

