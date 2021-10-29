#include <Rcpp.h>

#ifndef LIKECALC
#define LIKECALC

//' Calculate the likelihood under the bounded coalescence
//'
//' @param leaf_times Times that leaves are added.
//' @param leaves Number of leaves at each time point.
//' @param coalescence_times Times of coalescences.
//' @param ne Effective population size.
//' @param bound Bound time.
//' @export
// [[Rcpp::export]]
double bounded_times_likelihood_c(Rcpp::NumericVector leaf_times,
                                  Rcpp::IntegerVector leaves,
                                  Rcpp::NumericVector coalescence_times,
                                  double ne,
                                  double bound);



#endif
