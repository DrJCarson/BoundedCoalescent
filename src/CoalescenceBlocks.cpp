#include <Rcpp.h>
#include "FFBS.h"
#include "CoalescenceBlocks.h"

double constrain_coalescences_c(Rcpp::IntegerVector sample,
                                Rcpp::NumericVector times,
                                Rcpp::IntegerVector leaves,
                                double ne,
                                double bound,
                                Rcpp::NumericVector const_lower,
                                Rcpp::NumericVector const_upper,
                                Rcpp::IntegerVector const_lineages,
                                Rcpp::IntegerVector const_events) {

  // Total the leaves
  int total_leaves = Rcpp::sum(leaves);

  // Likelihood increment
  double likelihood = 1.0;

  // Random number for sampling
  double u;

  // Cumulative probability
  double sum_prob;

  // Counters
  int k, m, c = 0;

  // Number of coalescences
  int events, events_lhs, events_rhs;

  // Sub-interval information
  double dt, mid_time, prob_lhs, prob_rhs, prob_norm;

  // Initial constraints, bound interval
  events = sample(1) - sample(0);

  for (m = 0; m < events; ++m) {

    const_lower(c) = bound;
    const_upper(c) = times(0);
    const_lineages(c) = sample(1);
    const_events(c) = events;

    ++c;

  }

  // Initial constraints, leaf intervals
  for (k = 1; k < times.size(); ++k) {

    events = leaves(k - 1) + sample(k + 1) - sample(k);

    for (m = 0; m < events; ++m) {

      const_lower(c) = times(k - 1);
      const_upper(c) = times(k);
      const_lineages(c) = sample(k + 1);
      const_events(c) = events;

      ++c;

    }

  }

  // Separate coalescence events
  for (c = 0; c < (total_leaves - 1); ++c) {

    while (const_events(c) > 1) {

      // Check for interrupt
      Rcpp::checkUserInterrupt();


      events = const_events(c);

      dt = 0.5 * (const_upper(c) - const_lower(c));
      mid_time = const_upper(c) - dt;

      u = R::runif(0.0, 1.0);

      events_lhs = 0;
      events_rhs = events - events_lhs;

      prob_rhs = homochronous_probability(const_lineages(c),
                                          const_lineages(c) - events_rhs,
                                          dt, ne);

      prob_lhs = homochronous_probability(const_lineages(c) - events_rhs,
                                          const_lineages(c) - events,
                                          dt, ne);

      prob_norm = homochronous_probability(const_lineages(c),
                                           const_lineages(c) - events,
                                           2.0 * dt, ne);

      sum_prob = (prob_lhs * prob_rhs) / prob_norm;

      while (u > sum_prob) {

        // Check for interrupt
        Rcpp::checkUserInterrupt();


        ++events_lhs;
        --events_rhs;

        prob_rhs = homochronous_probability(const_lineages(c),
                                            const_lineages(c) - events_rhs,
                                            dt, ne);

        prob_lhs = homochronous_probability(const_lineages(c) - events_rhs,
                                            const_lineages(c) - events,
                                            dt, ne);

        prob_norm = homochronous_probability(const_lineages(c),
                                             const_lineages(c) - events,
                                             2.0 * dt, ne);

        sum_prob += (prob_lhs * prob_rhs) / prob_norm;

      }

      likelihood *= (prob_lhs * prob_rhs) / prob_norm;

      for (m = 0; m < events; ++m) {

        if (m < events_lhs) {

          const_upper(c + m) = mid_time;
          const_lineages(c + m) -= events_rhs;
          const_events(c + m) = events_lhs;

        } else {

          const_lower(c + m) = mid_time;
          const_events(c + m) = events_rhs;

        }

      }

    }

  }

  return likelihood;

}

