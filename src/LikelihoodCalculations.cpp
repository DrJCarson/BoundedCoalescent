#include <Rcpp.h>
#include "FFBS.h"

double bounded_times_likelihood_c(Rcpp::NumericVector leaf_times,
                                  Rcpp::IntegerVector leaves,
                                  Rcpp::NumericVector coalescence_times,
                                  double ne,
                                  double bound) {

  // Counters
  int k = leaf_times.size() - 1, c = coalescence_times.size() - 1;

  // Current time
  double current_time = leaf_times(k);

  // Current lineages
  int current_lineages = leaves(k);

  --k;

  // Next event times
  double next_leaf_time;
  if (k >= 0) {

    next_leaf_time = leaf_times(k);

  } else{

    next_leaf_time = bound;

  }

  double next_coalescence_time = coalescence_times(c);

  // Time step
  double dt;

  // Exponential coefficient
  double coef;

  // Likelihood
  double likelihood = 1.0;

  while (current_time > coalescence_times(0)) {

    // Check for interrupt
    Rcpp::checkUserInterrupt();

    if (next_leaf_time > next_coalescence_time) {

      dt = current_time - next_leaf_time;
      coef = (double(current_lineages) * (double(current_lineages) - 1.0)) /
        (2.0 * ne);

      likelihood *= std::exp(-coef * dt);

      current_lineages += leaves(k);
      current_time = next_leaf_time;

      --k;

      if (k >= 0) {

        next_leaf_time = leaf_times(k);

      } else{

        next_leaf_time = bound;

      }

    } else {

      dt = current_time - next_coalescence_time;
      coef = (double(current_lineages) * (double(current_lineages) - 1.0)) /
        (2.0 * ne);

      likelihood *= coef * std::exp(-coef * dt);

      --current_lineages;
      current_time = next_coalescence_time;

      --c;

      if (c >= 0) {

        next_coalescence_time = coalescence_times(c);

      }

    }

  }

  int total_leaves = Rcpp::sum(leaves);

  Rcpp::NumericVector forward_probs(total_leaves * (leaf_times.size() + 1));
  forward_probs.attr("dim") =
    Rcpp::Dimension(total_leaves, (leaf_times.size() + 1));

  forward_algorithm_c(leaf_times, leaves, ne, bound, forward_probs);

  likelihood /= forward_probs(0, 0);


  return likelihood;

}
