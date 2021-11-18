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

