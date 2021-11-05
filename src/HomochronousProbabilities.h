#include <Rcpp.h>
#ifndef HOMOCHRONOUS
#define HOMOCHRONOUS

//' Calculate coalescent probabilities (homochronous)
//'
//' Calculate the probability of i lineages coalescing down to j <= i lineages
//' in time dt in the homochronous setting.
//'
//' @param i Integer value for the starting number of lineages.
//' @param j Integer value for the final number of lineages.
//' @param dt Time period over which coalescences can occur.
//' @param ne Effective population size.
// [[Rcpp::export]]
double homochronous_probability(int i, int j, double dt, double ne);

#endif
