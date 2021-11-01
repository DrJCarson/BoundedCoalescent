#include <Rcpp.h>
#ifndef REJECTIONSAMPLING
#define REJECTIONSAMPLING

//' Rejection sampling coalescence times for the bounded coalescenct
//'
//' @param times Vector of ordered sampling times for leaves.
//' @param leaves Number of leaves taken at each sampling time.
//' @param ne Effective population size.
//' @param bound Bound time.
//' @param nsam Number of samples (deafault 1).
// [[Rcpp::export]]
Rcpp::List rejection_bounded_times(Rcpp::NumericVector times,
                                   Rcpp::IntegerVector leaves,
                                   double ne,
                                   double bound,
                                   int nsam = 1);

#endif
