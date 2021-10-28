#include <Rcpp.h>
#include "HomochronousProbabilities.h"
#include "ForwardAlgorithm.h"
#include "BackwardSampler.h"
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


