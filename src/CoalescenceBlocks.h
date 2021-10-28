#include <Rcpp.h>
#ifndef BLOCKS
#define BLOCKS


//' Determine blocks of coalescent events from backward sample
//'
//' @param sample Sample of number of lineages from backward sampler
//' @param times Vector of ordered sampling times for leaves.
//' @param leaves Number of leaves taken at each sampling time.
//' @param bound Bound time.
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame block_coalescences_c(Rcpp::IntegerVector sample,
                                     Rcpp::NumericVector times,
                                     Rcpp::IntegerVector leaves,
                                     double bound);

#endif
