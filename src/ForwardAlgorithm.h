#include <Rcpp.h>
#ifndef FORWARD
#define FORWARD

//' Forward Algorithm for the Bounded Coalescent
//'
//' Calculate the forward probabilities for the bounded coalescent.
//'
//' @param times Vector of ordered sampling times for leaves.
//' @param leaves Number of leaves taken at each sampling time.
//' @param ne Effective population size.
//' @param bound Bound time.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector forward_algorithm_c(Rcpp::NumericVector times,
                                        Rcpp::IntegerVector leaves,
                                        double ne,
                                        double bound);

#endif
