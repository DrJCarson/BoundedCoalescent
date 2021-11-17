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

  Rcpp::List back_sampler;
  Rcpp::IntegerVector lineages;

  Rcpp::List constraints;
  Rcpp::NumericVector const_lower, const_upper;
  Rcpp::IntegerVector const_lineages;

  Rcpp::List single_sample;

  for (n = 0; n < nsam; ++n) {

    back_sampler = backward_sampler_c(forward_probs, times, leaves, ne, bound);
    lineages = back_sampler["sample"];
    likelihood(n) *= double(back_sampler["likelihood"]);

    constraints = constrain_coalescences_c(lineages, times, leaves, ne, bound);
    const_lower = constraints["time_lower"];
    const_upper = constraints["time_upper"];
    const_lineages = constraints["lineages"];
    likelihood(n) *= double(constraints["likelihood"]);

    for (c = 0; c < total_leaves - 1; ++c) {

      single_sample = sample_coalescence_time_c(const_lower(c), const_upper(c),
                                                const_lineages(c), ne);



      sampled_times(n, c) = double(single_sample["time"]);
      likelihood(n) *= double(single_sample["likelihood"]);

    }

  }

  Rcpp::List out = Rcpp::List::create(Rcpp::Named("times") = sampled_times,
                                      Rcpp::Named("likelihood") = likelihood);

  return out;

}
