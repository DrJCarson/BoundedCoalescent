#' Determine blocks of coalescent events from backward sample
#'
#' @param sample Sample of number of lineages from backward sampler
#' @param times Vector of ordered sampling times for leaves.
#' @param leaves Number of leaves taken at each sampling time.
#' @param ne Effective population size.
#' @param bound Bound time.
#' @export
block_coalescences <- function(sample, times, leaves, bound) {

  # Number of coalescences in the block
  coalescences <- c(0, leaves[1:(length(times) - 1)]) +
    sample[2:(length(times) + 1)] -
    sample[1:length(times)]

  # Lower bound time for the block
  lower_const <- c(bound, times[1:(length(times) - 1)])

  # Upper bound time for the block
  upper_const <- times[1:length(times)]

  # Number of lineages that can coalesce at the start of the block
  initial_lineages <- sample[2:(length(times) + 1)]

  return(data.frame("coalescences" = coalescences,
                    "time_low" = lower_const,
                    "time_upp" = upper_const,
                    "lineages_upp" = initial_lineages))

}
