#include <Rcpp.h>
#ifndef BLOCKS
#define BLOCKS


//' Determine blocks of coalescent events from backward sample
//'
//' @param sample Sample of number of lineages from backward sampler
//' @param times Vector of ordered sampling times for leaves.
//' @param leaves Number of leaves taken at each sampling time.
//' @param bound Bound time.
// [[Rcpp::export]]
Rcpp::DataFrame block_coalescences_c(Rcpp::IntegerVector sample,
                                     Rcpp::NumericVector times,
                                     Rcpp::IntegerVector leaves,
                                     double bound);


//' Separate coalescences in a block
//'
//' @param coalescences Number of coalescences in block.
//' @param time_lower Lower bound time of the block.
//' @param time_upper Upper bound time of the block.
//' @param lineages_upper Number of lineages at the start of the block.
//' @param ne Effective population size.
// [[Rcpp::export]]
Rcpp::List separate_coalescences_c(int coalescences,
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
// [[Rcpp::export]]
Rcpp::List constrain_coalescences_c(Rcpp::IntegerVector sample,
                                         Rcpp::NumericVector times,
                                         Rcpp::IntegerVector leaves,
                                         double ne,
                                         double bound);


//' Constrain coalescence events to unique time intervals
//'
//' @param sample Sample of number of lineages from backward sampler
//' @param times Vector of ordered sampling times for leaves.
//' @param leaves Number of leaves taken at each sampling time.
//' @param ne Effective population size.
//' @param bound Bound time.
// [[Rcpp::export]]
double constrain_coalescences_new_c(Rcpp::IntegerVector sample,
                                    Rcpp::NumericVector times,
                                    Rcpp::IntegerVector leaves,
                                    double ne,
                                    double bound,
                                    Rcpp::NumericVector const_lower,
                                    Rcpp::NumericVector const_upper,
                                    Rcpp::IntegerVector const_lineages,
                                    Rcpp::IntegerVector const_events);

//' Sample a coalescence time within a specified interval
//'
//' @param time_lower Lower bound time of the interval.
//' @param time_upper Upper bound time of the interval.
//' @param lineages Starting number of lineages.
//' @param ne Effective population size.
// [[Rcpp::export]]
Rcpp::List sample_coalescence_time_c(double time_lower, double time_upper,
                                     int lineages, double ne);

#endif
