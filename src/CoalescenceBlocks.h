#include <Rcpp.h>
#ifndef BLOCKS
#define BLOCKS


//' Constrain coalescence events to unique time intervals
//'
//' @param sample Sample of number of lineages from backward sampler
//' @param times Vector of ordered sampling times for leaves.
//' @param leaves Number of leaves taken at each sampling time.
//' @param ne Effective population size.
//' @param bound Bound time.
//' @param const_lower Lower bound times for each coalescence event.
//' @param const_upper Upper bound times for each coalescence event.
//' @param const_lineages Number of lineages at the start of the interval.
//' @param const_events Number of coalescence events in the interval.
//' @param norm_tol Threshold to use approximate sampling due to loss of
//' significance.
// [[Rcpp::export]]
double constrain_coalescences_c(Rcpp::IntegerVector sample,
                                Rcpp::NumericVector times,
                                Rcpp::IntegerVector leaves,
                                double ne,
                                double bound,
                                Rcpp::NumericVector const_lower,
                                Rcpp::NumericVector const_upper,
                                Rcpp::IntegerVector const_lineages,
                                Rcpp::IntegerVector const_events,
                                double norm_tol = 1.0e-10);



#endif
