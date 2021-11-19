#include <Rcpp.h>
#ifndef FULLSAMPLING
#define FULLSAMPLING

//' Sample coalescence times for the bounded coalescenct
//'
//' @param times Vector of ordered sampling times for leaves.
//' @param leaves Number of leaves taken at each sampling time.
//' @param ne Effective population size.
//' @param bound Bound time.
//' @param nsam Number of samples (deafault 1).
// [[Rcpp::export]]
Rcpp::List sample_bounded_times_c(Rcpp::NumericVector times,
                                Rcpp::IntegerVector leaves,
                                double ne,
                                double bound,
                                int nsam = 1);


//' Rejection sampling coalescence times for the bounded coalescenct
//'
//' @param times Vector of ordered sampling times for leaves.
//' @param leaves Number of leaves taken at each sampling time.
//' @param ne Effective population size.
//' @param bound Bound time.
//' @param nsam Number of samples (deafault 1).
// [[Rcpp::export]]
Rcpp::List rejection_bounded_times_c(Rcpp::NumericVector times,
                                   Rcpp::IntegerVector leaves,
                                   double ne,
                                   double bound,
                                   int nsam = 1);


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
