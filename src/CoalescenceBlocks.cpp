#include <Rcpp.h>

//' Determine blocks of coalescent events from backward sample
//'
//' @param sample Sample of number of lineages from backward sampler
//' @param times Vector of ordered sampling times for leaves.
//' @param leaves Number of leaves taken at each sampling time.
//' @param ne Effective population size.
//' @param bound Bound time.
//' @export
// [[Rcpp::export]]
Rcpp::List block_coalescences_c(Rcpp::IntegerVector sample,
                                Rcpp::NumericVector times,
                                Rcpp::IntegerVector leaves,
                                double bound) {

  // Number of coalescences in each interval
  Rcpp::IntegerVector coalescences(times.size());

  // Interval limits
  Rcpp::NumericVector time_lower(times.size()), time_upper(times.size());

  // Initial number of lineages for each interval
  Rcpp::IntegerVector lineages_upper(times.size());

  // Counter
  int k;

  // Interval at the bound
  coalescences(0) = sample(1) - sample(0);
  time_lower(0) = bound;
  time_upper(0) = times(0);
  lineages_upper(0) = sample(1);

  // Iterate through intervals
  for (k = 1; k < times.size(); ++k) {

    coalescences(k) = leaves(k - 1) + sample(k + 1) - sample(k);
    time_lower(k) = times(k - 1);
    time_upper(k) = times(k);
    lineages_upper(k) = sample(k + 1);

  }

  // Return block information
  Rcpp::List out(4);
  out(0) = coalescences;
  out(1) = time_lower;
  out(2) = time_upper;
  out(3) = lineages_upper;
  out.names() = Rcpp::CharacterVector::create("coalescences", "time_lower",
            "time_upper", "lineages_upper");

  return out;

}
