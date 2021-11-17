#include <Rcpp.h>
#include "FFBS.h"
#include "CoalescenceBlocks.h"

// Determine number of coalescences in each sampling interval
Rcpp::DataFrame block_coalescences_c(Rcpp::IntegerVector sample,
                                Rcpp::NumericVector times,
                                Rcpp::IntegerVector leaves,
                                double bound) {

  // Number of coalescences in each interval
  Rcpp::IntegerVector coalescences(times.size());

  // Interval limits
  Rcpp::NumericVector block_lower(times.size()), block_upper(times.size());

  // Initial number of lineages for each interval
  Rcpp::IntegerVector lineages_upper(times.size());

  // Counter
  int k;

  // Interval at the bound
  coalescences(0) = sample(1) - sample(0);
  block_lower(0) = bound;
  block_upper(0) = times(0);
  lineages_upper(0) = sample(1);

  // Iterate through intervals
  for (k = 1; k < times.size(); ++k) {

    coalescences(k) = leaves(k - 1) + sample(k + 1) - sample(k);
    block_lower(k) = times(k - 1);
    block_upper(k) = times(k);
    lineages_upper(k) = sample(k + 1);

  }

  Rcpp::DataFrame out =
    Rcpp::DataFrame::create(Rcpp::Named("coalescences") = coalescences,
                            Rcpp::Named("time_lower") = block_lower,
                            Rcpp::Named("time_upper") = block_upper,
                            Rcpp::Named("lineages_upper") = lineages_upper);

  return out;

}



Rcpp::List separate_coalescences_c(int coalescences,
                                   double time_lower,
                                   double time_upper,
                                   int lineages_upper,
                                   double ne) {

  // Check for interrupt
  Rcpp::checkUserInterrupt();

  // Counter
  int k;

  // Likelihood increment
  double sub_likelihood = 1.0;

  // Define new partition
  Rcpp::NumericVector sub_times(coalescences);
  sub_times(coalescences - 1) = time_upper;

  // Set leaf additions
  Rcpp::IntegerVector sub_leaves(coalescences);
  sub_leaves(coalescences - 1) = lineages_upper;

  if ((time_upper - time_lower) < (4.0 * ne)) {

    for(k = coalescences - 2; k >= 0; --k){

      sub_times(k) = sub_times(k + 1) -
        (time_upper - time_lower) / double(coalescences);

      sub_leaves(k) = 0;

    }

  } else {

    double lineages_expected = double(lineages_upper);

    for(k = coalescences - 2; k >= 0; --k){

      sub_times(k) = sub_times(k + 1) - ((2.0 * ne) /
        (lineages_expected * (lineages_expected - 1.0)));

      --lineages_expected;

      sub_leaves(k) = 0;

    }

  }


  // Set new bound and size
  double sub_bound = time_lower;
  int sub_bound_size = lineages_upper - coalescences;

  // Forward algorithm for the new partition
  Rcpp::NumericVector sub_probs(lineages_upper * (sub_times.size() + 1));
  sub_probs.attr("dim") =
    Rcpp::Dimension(lineages_upper, (sub_times.size() + 1));

  forward_algorithm_c(sub_times, sub_leaves, ne, sub_bound, sub_probs);

  Rcpp::IntegerVector sub_sample(sub_times.size() + 1);
  sub_likelihood *= backward_sampler_c(sub_probs, sub_times, sub_leaves, ne,
                                       sub_bound, sub_sample, sub_bound_size);

  Rcpp::List sub_const = constrain_coalescences_c(sub_sample, sub_times,
                                                 sub_leaves, ne, sub_bound);

  Rcpp::NumericVector sub_lower = sub_const["time_lower"];
  Rcpp::NumericVector sub_upper = sub_const["time_upper"];
  Rcpp::IntegerVector sub_lineages = sub_const["lineages"];
  sub_likelihood *= double(sub_const["likelihood"]);

  Rcpp::List out =
    Rcpp::List::create(Rcpp::Named("time_lower") = sub_lower,
                       Rcpp::Named("time_upper") = sub_upper,
                       Rcpp::Named("sub_lineages") = sub_lineages,
                       Rcpp::Named("sub_likelihood") = sub_likelihood);

  return out;

}



Rcpp::List constrain_coalescences_c(Rcpp::IntegerVector sample,
                                    Rcpp::NumericVector times,
                                    Rcpp::IntegerVector leaves,
                                    double ne,
                                    double bound) {

  // Total the leaves
  int total_leaves = Rcpp::sum(leaves);

  // Likelihood increment
  double likelihood = 1.0;

  // Interval limits
  Rcpp::NumericVector const_lower(total_leaves - 1);
  Rcpp::NumericVector const_upper(total_leaves - 1);
  Rcpp::IntegerVector const_lineages(total_leaves - 1);

  // Counters
  int k, m, c = 0;

  // Initial blocks
  Rcpp::DataFrame blocks = block_coalescences_c(sample, times, leaves, bound);
  Rcpp::NumericVector coalescences = blocks["coalescences"];
  Rcpp::NumericVector block_lower = blocks["time_lower"];
  Rcpp::NumericVector block_upper = blocks["time_upper"];
  Rcpp::NumericVector lineages_upper = blocks["lineages_upper"];

  // Sub-interval information
  Rcpp::List sub_intervals;
  Rcpp::NumericVector sub_lower, sub_upper;
  Rcpp::IntegerVector sub_lineages;

  // Iterate through blocks
  for (k = 0; k < times.size(); ++k) {

    if (coalescences(k) == 1) {

      const_lower(c) = block_lower(k);
      const_upper(c) = block_upper(k);
      const_lineages(c) = lineages_upper(k);
      ++c;

    } else if(coalescences(k) > 1) {

      sub_intervals = separate_coalescences_c(coalescences(k), block_lower(k),
                                              block_upper(k), lineages_upper(k),
                                              ne);

      sub_lower = sub_intervals["time_lower"];
      sub_upper = sub_intervals["time_upper"];
      sub_lineages = sub_intervals["sub_lineages"];
      likelihood *= double(sub_intervals["sub_likelihood"]);

      for (m = 0; m < coalescences(k); ++m) {

        const_lower(c) = sub_lower(m);
        const_upper(c) = sub_upper(m);
        const_lineages(c) = sub_lineages(m);
        ++c;

      }

    }

  }

  Rcpp::List out =
    Rcpp::List::create(Rcpp::Named("time_lower") = const_lower,
                       Rcpp::Named("time_upper") = const_upper,
                       Rcpp::Named("lineages") = const_lineages,
                       Rcpp::Named("likelihood") = likelihood);

  return out;

}


Rcpp::List sample_coalescence_time_c(double time_lower, double time_upper,
                                     int lineages, double ne) {

  // Random number for sampling
  double u = R::runif(0.0, 1.0);

  //Double copy of starting lineages
  double dl = double(lineages);

  // Normalising constant
  double z = (ne / (dl - 1.0)) *
    (1.0 - std::exp(((dl - 1.0) / ne) * (time_lower - time_upper)));

  // Time of coalescence
  double coalescence_time = time_upper + (ne / (dl - 1.0)) *
    std::log(1.0 - ((dl - 1.0) / ne) * z * u);

  double likelihood = (1.0 / z) *
    std::exp(((dl - 1.0) / ne) * (coalescence_time - time_upper));

  Rcpp::List out = Rcpp::List::create(Rcpp::Named("time") = coalescence_time,
                                      Rcpp::Named("likelihood") = likelihood);

  return out;

}
