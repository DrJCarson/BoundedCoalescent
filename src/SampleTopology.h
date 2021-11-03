#include <Rcpp.h>
#ifndef TOPOLOGY
#define TOPOLOGY

//' Sample topology given coalescent times.
//'
//' @param leaf_times Vector of ordered sampling times for leaves.
//' @param leaves Number of leaves taken at each sampling time.
//' @param coalescence_times Vector of ordered coalescence times.
// [[Rcpp::export]]
Rcpp::List sample_topology_c(Rcpp::NumericVector leaf_times,
                             Rcpp::IntegerVector leaves,
                             Rcpp::NumericVector coalescence_times);


#endif
