
#' Sample Coalescence Times Under the Bounded Coalescent
#'
#' @param t Vector of leaf sampling times.
#' @param l vector of leaves sampled at each time.
#' @param ne Effective population size.
#' @param b Bound time.
#' @param nsam Number of samples.
#' @param method Sampling method for the coalescence times.
#' @export
bounded_sample_times <- function(t, l, ne, b, nsam = 1, method = "direct") {

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
#' @param topology Boolean value indicating whether the topology should be
#' included in the likelihood calculation.
#' @export
bounded_likelihood <- function(t, l, c, ne, b, topology = T) {

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
  sorted_c <- sort(c)

  likelihood <-
    bounded_times_likelihood_c(ordered_t, ordered_l, sorted_c, ne, b)

  if (topology) {

    for (i in 1:length(c)) {

      lineages <- sum(l[which(t > sorted_c[i])]) - (length(c) - i)

      likelihood <- likelihood * (2 / (lineages * (lineages - 1)))

    }

  }

  return(likelihood)

}


#' Sample a phylogeny under the bounded coalescent
#'
#' @param t Vector of leaf sampling times.
#' @param l vector of leaves sampled at each time.
#' @param ne Effective population size.
#' @param b Bound time.
#' @param nsam Number of samples.
#' @param tip.label Labels for sampled leaves.
#' @param node.label Labels for nodes.
#' @param method Sampling method for the coalescence times.
#' @export
bounded_sample_phylo <- function(t, l, ne, b, nsam = 1, tip.label, node.label,
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

  if (nsam <= 0) {

    stop("Sample size must be positive.")

  }

  leaf_order <- order(t)

  ordered_t <- t[leaf_order]
  ordered_l <- l[leaf_order]

  if (method == "direct") {

    times_sample <- sample_bounded_times_c(ordered_t, ordered_l, ne, b, nsam)

  } else if(method == "rejection") {

    times_sample <- rejection_bounded_times(ordered_t, ordered_l, ne, b, nsam)

  }

  mphylo <- vector("list", nsam)
  mlikelihood <- numeric(nsam)
  mnodes <- vector("list", nsam)


  if (missing(tip.label)) {

    tip.label <- as.character(1:sum(l))

  }

  if (missing(node.label)) {

    node.label <- as.character(sum(l) + 1:((sum(l) - 1)))

  }

  for (i in 1:nsam) {

    topology_sample <- sample_topology_c(ordered_t,
                                         ordered_l,
                                         times_sample$times[i,])

    nodes <- data.frame(node = 1:(2 * sum(l) - 1),
                             time = topology_sample$times,
                             ancestor = topology_sample$ancestors)

    mnodes[[i]] <- nodes

    edge <- topology_sample$edge
    edge_length <- topology_sample$times[topology_sample$edge[,2]] -
      topology_sample$times[topology_sample$edge[,1]]

    root_time <- min(times_sample$times[i,])

    if (is.finite(b)) {

      root_edge <- root_time - b

    } else {

      root_edge <- NULL

    }

    phylo_sample <- list(edge = edge,
                         edge.length = edge_length,
                         tip.label = tip.label,
                         node.label = node.label,
                         Nnode = as.integer(sum(l) - 1),
                         root.time = root_time,
                         root.edge = root_edge)

    class(phylo_sample) <- 'phylo'

    mphylo[[i]] <- phylo_sample

    phylo_likelihood <- times_sample$likelihood[i] * topology_sample$likelihood

    mlikelihood[i] <- phylo_likelihood

  }

  if (nsam == 1) {

    return(list(phylo = phylo_sample,
                likelihood = phylo_likelihood,
                coalescence_times = c(times_sample$times),
                nodes = nodes
    ))

  } else {

    class(mphylo) <- "multiPhylo"

    return(list(phylo = mphylo,
                likelihood = mlikelihood,
                coalescence_times = times_sample$times,
                nodes = mnodes
    ))


  }


}


#' Calculate the likelihood of a phylogeny under the bounded coalescence
#'
#' @param phy Phylogeny of class 'phylo'.
#' @param ne Effective population size.
#' @param b Bound time.
#' inputs)
#' @param topology Boolean value indicating whether the topology should be
#' included in the likelihood calculation.
#' @export
bounded_likelihood_phylo <- function(phy, ne, b, topology = T) {

  if (class(phy) == "phylo") {

    if (is.null(phy$root.time)) {

      stop("root.time is required.")

    }

    node_times <- phy$root.time+ape::dist.nodes(phy)[ape::Ntip(phy)+1,]

    total_nodes <- length(node_times)

    leaf_indices <- 1:ape::Ntip(phy)
    coalescence_indices <- (ape::Ntip(phy)+1):total_nodes

    t <- node_times[leaf_indices]
    l <- as.integer(rep(1, length(t)))
    c <- node_times[coalescence_indices]

    likelihood <- bounded_likelihood(t, l, c, ne, b, topology)

    return(likelihood)

  } else if (class(phy) == "multiPhylo") {

    likelihood <- numeric(length(phy))

    for (i in 1:length(phy)) {

      if (is.null(phy[[i]]$root.time)) {

        stop("root.time is required.")

      }

      node_times <- phy[[i]]$root.time+
        ape::dist.nodes(phy[[i]])[ape::Ntip(phy[[i]])+1,]

      total_nodes <- length(node_times)

      leaf_indices <- 1:ape::Ntip(phy[[i]])
      coalescence_indices <- (ape::Ntip(phy[[i]])+1):total_nodes

      t <- node_times[leaf_indices]
      l <- as.integer(rep(1, length(t)))
      c <- node_times[coalescence_indices]

      likelihood[i] <- bounded_likelihood(t, l, c, ne, b, topology)

    }

    return(likelihood)

  } else {

    stop("phy is not an allowed class.")

  }

}
