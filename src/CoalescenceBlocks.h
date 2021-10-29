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


//' Separate coalescences in a block
//'
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame separate_coalescences_c(int coalescences,
                                        double time_lower,
                                        double time_upper,
                                        int lineages_upper,
                                        double ne);


//' Constrain coalescence events to unique time intervals
//'
//' @param sample Sample of number of lineages from backward sampler
//' @param times Vector of ordered sampling times for leaves.
//' @param leaves Number of leaves taken at each sampling time.
//' @param ne Effective population size.
//' @param bound Bound time.
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame constrain_coalescences_c(Rcpp::IntegerVector sample,
                                         Rcpp::NumericVector times,
                                         Rcpp::IntegerVector leaves,
                                         double ne,
                                         double bound);

#endif
