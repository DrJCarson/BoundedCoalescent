#include <Rcpp.h>
#include "FFBS.h"
#include "CoalescenceBlocks.h"
#include "RejectionSampling.h"

Rcpp::List rejection_bounded_times(Rcpp::NumericVector times,
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



