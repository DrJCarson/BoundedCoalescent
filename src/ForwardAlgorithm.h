#include <Rcpp.h>
#ifndef FORWARD
#define FORWARD

double forward_algorithm_c(Rcpp::NumericVector times,
                           Rcpp::IntegerVector leaves,
                           double Ne,
                           double bound);

#endif
