#include <Rcpp.h>
#include "HomochronousProbabilities.h"

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
