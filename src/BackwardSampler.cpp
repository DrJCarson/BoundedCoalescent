#include <Rcpp.h>
#include "HomochronousProbabilities.h"
#include "ForwardAlgorithm.h"

//' Backward Sampler for the Bounded Coalescent
//'
//' Samples the number of lineages at specified times given forward
//' probabilities
//'
//' @param forward_probs 2D array of probabilities from the forward algorithm.
//' @param times Vector of ordered sampling times for leaves.
//' @param leaves Number of leaves taken at each sampling time.
//' @param ne Effective population size.
//' @param bound Bound time.
//' @param bound_size Number of lineages at the bound (default 1).
//' @export
// [[Rcpp::export]]
Rcpp::List backward_sampler_c(Rcpp::NumericVector forward_probs,
                              Rcpp::NumericVector times,
                              Rcpp::IntegerVector leaves,
                              double ne,
                              double bound,
                              int bound_size = 1) {

  // Sampled lineages including bound
  Rcpp::IntegerVector lineages(times.size() + 1);
  lineages.names() = Rcpp::colnames(forward_probs);
  lineages(0) = bound_size;

  // Likelihood of sample
  double likelihood = 1.0;

  // Total the leaves
  int total_leaves = Rcpp::sum(leaves);

  // Smoothed (sampling) probabilities
  Rcpp::NumericVector smoothed_probs(total_leaves);

  // Counters
  int k, i, j;
  double sum_prob;

  // Time difference and transition probability
  double dt, transition_prob;

  // Random number for sampling
  double u;

  // Initiate backwards sampling
  dt = times(0) - bound;

  // Generate random number for sampling
  u = R::runif(0.0, 1.0);

  // Track cumulative probability for sampling
  sum_prob = 0.0;

  // Evaluate smoothed probabilities, sample, and update likelihood
  for (i = 1; i <= total_leaves; ++i) {

    transition_prob = homochronous_probability(i, lineages(0), dt, ne);

    smoothed_probs(i - 1) = (transition_prob * forward_probs(i - 1, 1)) /
      forward_probs(lineages(0) - 1, 0);

    sum_prob += smoothed_probs(i - 1);

    if (u < sum_prob) {

      lineages(1) = i;
      likelihood *= smoothed_probs(i - 1);
      break;

    }

  }

  // Backwards sampling recursion
  for (k = 1; k < times.size(); ++k) {

    dt = times(k) - times(k - 1);

    // Generate random number for sampling
    u = R::runif(0.0, 1.0);

    // Track cumulative probability for sampling
    sum_prob = 0.0;

    // Evaluate smoothed probabilities, sample, and update likelihood
    for (i = 1; i <= total_leaves; ++i) {

      // Number of lineages before leaves are added
      j = lineages(k) - leaves(k - 1);

      transition_prob = homochronous_probability(i, j, dt, ne);

      smoothed_probs(i - 1) = (transition_prob * forward_probs(i - 1, k + 1)) /
        forward_probs(lineages(k) - 1, k);

      sum_prob += smoothed_probs(i - 1);

      if (u < sum_prob) {

        lineages(k + 1) = i;
        likelihood *= smoothed_probs(i - 1);
        break;

      }

    }

  }

  // Return sample and likelihood
  Rcpp::List out(2);
  out(0) = lineages;
  out(1) = likelihood;
  out.names() = Rcpp::CharacterVector::create("sample", "likelihood");

  return out;

}
