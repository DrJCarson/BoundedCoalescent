#include <Rcpp.h>
#ifndef BACKWARD
#define BACKWARD

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
                              int bound_size = 1);

#endif
