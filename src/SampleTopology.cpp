#include <Rcpp.h>
#include "SampleTopology.h"

Rcpp::List sample_topology_c(Rcpp::NumericVector leaf_times,
                             Rcpp::IntegerVector leaves,
                             Rcpp::NumericVector coalescence_times) {

  // Total the leaves
  int total_leaves = Rcpp::sum(leaves);

  // Internal node ancestors
  Rcpp::IntegerVector edge(4 * (total_leaves - 1));
  edge.attr("dim") = Rcpp::Dimension(2 * (total_leaves - 1), 2);

  // Node ancestors and times
  Rcpp::IntegerVector node_ancestors(2 * total_leaves - 1);
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

      node_ancestors(anc_1) = n_index + 1;

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

      node_ancestors(anc_2) = n_index + 1;

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

  Rcpp::List out = Rcpp::List::create(Rcpp::Named("times") = node_times,
                                      Rcpp::Named("ancestors") = node_ancestors,
                                      Rcpp::Named("edge") = edge,
                                      Rcpp::Named("likelihood") = likelihood);

  return out;

}
