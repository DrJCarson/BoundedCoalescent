
#' Sample Coalescence Times Under the Bounded Coalescent
#'
#' @param t Vector of leaf sampling times.
#' @param l vector of leaves sampled at each time.
#' @param ne Effective population size.
#' @param b Bound time.
#' @param nsam Number of samples.
#' @export
bounded_times_sample <- function(t, l, ne, b, nsam = 1, method = "direct") {

  if (length(l) != length(t)) {

    stop("t and l are of different lengths.")

  }

  if (!(length(t) > 0)) {

    stop("Need at least one sampling time.")

  }

  if (sum(l) <= 1) {

    stop("Need at least two leaves to coalesce.")

  }

  if (b > min(t)) {

    stop("Bound time is greater than earliest sampling time.")

  }

  if (ne <= 0) {

    stop("Effective population size must be positive.")

  }

  if (nsam <= 0) {

    stop("Sample size must be positive.")

  }

  leaf_order <- order(t)

  ordered_t <- t[leaf_order]
  ordered_l <- l[leaf_order]

  if (method == "direct") {

    full_sample <- sample_bounded_times_c(ordered_t, ordered_l, ne, b, nsam)

  } else if(method == "rejection") {

    full_sample <- rejection_bounded_times(ordered_t, ordered_l, ne, b, nsam)

  }

  return(full_sample)

}


#' Calculate the Likelihood of Coalescence Times Under the Bounded Coalescent
#'
#' @param t Vector of leaf sampling times.
#' @param l vector of leaves sampled at each time.
#' @param c Vector of coalescence times.
#' @param ne Effective population size.
#' @param b Bound time.
#' @export
bounded_times_likelihood <- function(t, l, c, ne, b) {

  if (length(l) != length(t)) {

    stop("t and l are of different lengths.")

  }

  if (!(length(t) > 0)) {

    stop("Need at least one sampling time.")

  }

  if (sum(l) <= 1) {

    stop("Need at least two leaves to coalesce.")

  }

  if (b > min(t)) {

    stop("Bound time is greater than earliest sampling time.")

  }

  if (ne <= 0) {

    stop("Effective population size must be positive.")

  }


  if (length(c) != (sum(l) - 1)) {

    stop("Number of coalescence times inconsistent with number of leaves.")

  }

  leaf_order <- order(t)

  ordered_t <- t[leaf_order]
  ordered_l <- l[leaf_order]

  likelihood <- bounded_times_likelihood_c(ordered_t, ordered_l, sort(c), ne, b)

  return(likelihood)

}
