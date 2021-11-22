#include <Rcpp.h>
#ifndef FFBS
#define FFBS

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


//' Forward Algorithm for the Bounded Coalescent
//'
//' Calculate the forward probabilities for the bounded coalescent.
//'
//' @param times Vector of ordered sampling times for leaves.
//' @param leaves Number of leaves taken at each sampling time.
//' @param ne Effective population size.
//' @param bound Bound time.
//' @param forward_probs Array to store probabilities.
// [[Rcpp::export]]
void forward_algorithm_c(Rcpp::NumericVector times,
                         Rcpp::IntegerVector leaves,
                         double ne,
                         double bound,
                         Rcpp::NumericVector forward_probs);


//' Backward Sampler for the Bounded Coalescent
//'
//' Samples the number of lineages at specified times given forward
//' probabilities
//'
//' @param forward_probs 2D array of probabilities from the forward algorithm.
//' @param times Vector of ordered sampling times for leaves.
//' @param leaves Number of leaves taken at each sampling time.
//' @param ne Effective population size.
//' @param bound Bound time.
//' @param bound_size Number of lineages at the bound (default 1).
// [[Rcpp::export]]
double backward_sampler_c(Rcpp::NumericVector forward_probs,
                          Rcpp::NumericVector times,
                          Rcpp::IntegerVector leaves,
                          double ne,
                          double bound,
                          Rcpp::IntegerVector lineages,
                          int bound_size = 1);

#endif
