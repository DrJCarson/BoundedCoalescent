#include <Rcpp.h>
#include "HomochronousProbabilities.h"

// Construct vector of column names for the forward probabilities
Rcpp::CharacterVector forward_algorithm_colnames(Rcpp::NumericVector times) {

  int c;

  Rcpp::CharacterVector cnames(times.size() + 1, "t");

  for (c = 0; c <= times.size(); ++c) {

    if (c == 0) {

      cnames[c] += "*";

    } else {

      cnames[c] += std::to_string(c);

    }

  }

  return cnames;

}

// Construct vector of row names for the forward probabilities
Rcpp::CharacterVector forward_algorithm_rownames(int total_leaves) {

  int r;

  Rcpp::CharacterVector rnames(total_leaves, "P(");

  for (r = 0; r < total_leaves; ++r) {

    rnames[r] += std::to_string(r + 1) + ")";

  }

  return rnames;

}

//' Forward Algorithm for the Bounded Coalescent
//'
//' Calculate the forward probabilities for the bounded coalescent.
//'
//' @param times Vector of ordered sampling times for leaves.
//' @param leaves Number of leaves taken at each sampling time.
//' @param Ne Effective population size.
//' @param bound Bound time.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector forward_algorithm_c(Rcpp::NumericVector times,
                             Rcpp::IntegerVector leaves,
                             double Ne,
                             double bound) {

  // Total the leaves
  int total_leaves = Rcpp::sum(leaves);

  // Construct array for forward probabilities
  Rcpp::NumericVector forward_probs(total_leaves * (times.size() + 1));
  forward_probs.attr("dim") = Rcpp::Dimension(total_leaves, (times.size() + 1));
  Rcpp::colnames(forward_probs) = forward_algorithm_colnames(times);
  Rcpp::rownames(forward_probs) = forward_algorithm_rownames(total_leaves);

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

      for (i = leaves(k); i <= sum_leaves; ++i) {

        transition_prob = homochronous_probability(i, j, dt, Ne);

        forward_probs(j + leaves(k - 1) - 1, k) += transition_prob *
          forward_probs(i - 1, k + 1);

      }

    }

    sum_leaves += leaves(k - 1);

  }


  // Bound probabilities
  dt = times(0) - bound;

  for (j = 1; j <= sum_leaves; ++j) {

    for (i = leaves(0); i <= sum_leaves; ++i) {

      transition_prob = homochronous_probability(i, j, dt, Ne);

      forward_probs(j - 1, 0) += transition_prob * forward_probs(i - 1, 1);

    }

  }

  return forward_probs;

}
