
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


#' Sample a phylogeny under the bounded coalescent
#'
#' @param t Vector of leaf sampling times.
#' @param l vector of leaves sampled at each time.
#' @param ne Effective population size.
#' @param b Bound time.
#' @param tip_labels Labels for sampled leaves.
#' @param method Sampling method for the coalescence times.
#' @export
bounded_phylo_sample <- function(t, l, ne, b, tip_labels = NA,
                                 method = "direct") {

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

  leaf_order <- order(t)

  ordered_t <- t[leaf_order]
  ordered_l <- l[leaf_order]

  if (method == "direct") {

    times_sample <- sample_bounded_times_c(ordered_t, ordered_l, ne, b)

  } else if(method == "rejection") {

    times_sample <- rejection_bounded_times(ordered_t, ordered_l, ne, b)

  }

  topology_sample <- sample_topology_c(ordered_t,
                                       ordered_l,
                                       times_sample$times)

  edge <- topology_sample$ancestors
  edge_length <- topology_sample$times[topology_sample$ancestors[,2]] -
    topology_sample$times[topology_sample$ancestors[,1]]

  if (is.na(tip_labels)) {

    tip_labels <- character(sum(l))
    tip_index <- 1

    for (i in 1:length(l)) {

      for (j in 1:ordered_l[i]) {

        tip_labels[tip_index] <- paste(i,".",j,sep="")
        tip_index <- tip_index + 1

      }

    }

  }

  phylo_sample <- list(edge = edge,
                       edge.length = edge_length,
                       tip.label = tip_labels,
                       Nnode = as.integer(sum(l) - 1))

  class(phylo_sample) <- "phylo"

  phylo_likelihood <- times_sample$likelihood * topology_sample$likelihood

  return(list(phylo = phylo_sample,
              likelihood = phylo_likelihood))

}
