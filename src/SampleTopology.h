#include <Rcpp.h>
#ifndef TOPOLOGY
#define TOPOLOGY

//' Sample topology given coalescent times.
//'
//' @param forward_probs 2D array of probabilities from the forward algorithm.
//' @param times Vector of ordered sampling times for leaves.
//' @param leaves Number of leaves taken at each sampling time.
//' @param ne Effective population size.
//' @param bound Bound time.
//' @param bound_size Number of lineages at the bound (default 1).
//' @export
// [[Rcpp::export]]
Rcpp::List sample_topology_c(Rcpp::NumericVector leaf_times,
                             Rcpp::IntegerVector leaves,
                             Rcpp::NumericVector coalescence_times);


#endif
