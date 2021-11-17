#include <Rcpp.h>
#include "FFBS.h"

// Coalescent probabilities for the homochronous setting
double homochronous_probability(int i, int j, double dt, double ne) {

  // Counters
  int k, l;

  // Double copies of counters
  double dk, dl;

  // Total and incremental probabilities
  double prob_total, prob_increment;

  // Check that inputs are valid
  if (i <= 0 || j <= 0 || i < j || dt < 0 || ne <= 0) {

    return 0.0;

  }

  if (j == 1) {

    // Initialise total probability
    prob_total = 1.0;

    for (k = 2; k <= i; ++k) {

      dk = double(k);

      // Initialise probability increment
      prob_increment = 1.0;

      for (l = 2; l <= i; ++l) {

        dl = double(l);

        // Calculate coefficients
        if (l != k) {

          prob_increment *= (dl * (dl - 1.0)) /
            (dl * (dl - 1.0) - dk * (dk - 1.0));

        }

      }

      // Calculate exponential
      prob_increment *= std::exp(- ((dk * (dk - 1.0)) / (2.0 * ne)) * dt);

      // Update total probability
      prob_total -= prob_increment;

    }

  } else{

    // Initialise total probability
    prob_total = 0.0;

    for (k = j; k <= i; ++k) {

      dk = double(k);

      // Initialise probability increment
      prob_increment = (dk * (dk - 1.0)) / (double(j) * (double(j) - 1.0));

      for (l = j; l <= i; ++l) {

        dl = double(l);

        // Calculate coefficients
        if (l != k) {

          prob_increment *= (dl * (dl - 1.0)) /
            (dl * (dl - 1.0) - dk * (dk - 1.0));

        }

      }

      // Calculate exponential
      prob_increment *= std::exp(- ((dk * (dk - 1.0)) / (2.0 * ne)) * dt);

      // Update total probability
      prob_total += prob_increment;

    }

  }

  return prob_total;

}


// Forward algorithm for the bounded coalescent
void forward_algorithm_c(Rcpp::NumericVector times,
                         Rcpp::IntegerVector leaves,
                         double ne,
                         double bound,
                         Rcpp::NumericVector forward_probs) {

  // Counters
  int i, j, k = times.size(), sum_leaves = leaves(times.size() - 1);

  // Time difference and transition probability
  double dt, transition_prob;

  // Initiate forward algorithm
  forward_probs(leaves(k - 1) - 1, k) = 1.0;

  // Forward recursion through sampling times
  for (k = times.size() - 1; k > 0; --k) {

    dt = times(k) - times(k - 1);

    for (j = 1; j <= sum_leaves; ++j) {

      for (i = 1; i <= sum_leaves; ++i) {

        transition_prob = homochronous_probability(i, j, dt, ne);

        forward_probs(j + leaves(k - 1) - 1, k) += transition_prob *
          forward_probs(i - 1, k + 1);

      }

    }

    sum_leaves += leaves(k - 1);

  }


  // Bound probabilities
  dt = times(0) - bound;

  for (j = 1; j <= sum_leaves; ++j) {

    for (i = 1; i <= sum_leaves; ++i) {

      transition_prob = homochronous_probability(i, j, dt, ne);

      forward_probs(j - 1, 0) += transition_prob * forward_probs(i - 1, 1);

    }

  }

}


// Backward sampler for the bounded coalescent
Rcpp::List backward_sampler_c(Rcpp::NumericVector forward_probs,
                              Rcpp::NumericVector times,
                              Rcpp::IntegerVector leaves,
                              double ne,
                              double bound,
                              int bound_size) {

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
