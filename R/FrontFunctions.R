
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
#' @param tip.label Labels for sampled leaves.
#' @param node.label Labels for nodes.
#' @param method Sampling method for the coalescence times.
#' @export
bounded_sample_phylo <- function(t, l, ne, b, tip.label, node.label,
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

  nodes <- data.frame(node = 1:(2 * sum(l) - 1),
                      time = topology_sample$times,
                      ancestor = topology_sample$ancestors)

  edge <- topology_sample$edge
  edge_length <- topology_sample$times[topology_sample$edge[,2]] -
    topology_sample$times[topology_sample$edge[,1]]

  if (missing(tip.label)) {

    tip.label <- as.character(1:sum(l))

  }

  if (missing(node.label)) {

    node.label <- as.character(sum(l) + 1:((sum(l) - 1)))

  }

  phylo_sample <- list(edge = edge,
                       edge.length = edge_length,
                       tip.label = tip.label,
                       node.label = node.label,
                       Nnode = as.integer(sum(l) - 1),
                       root.time = min(times_sample$times))

  class(phylo_sample) <- 'phylo'

  phylo_likelihood <- times_sample$likelihood * topology_sample$likelihood

  return(list(phylo = phylo_sample,
              likelihood = phylo_likelihood,
              coalescence_times = c(times_sample$times),
              nodes = nodes
              ))

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

    node_ancestors <- integer(2 * phy$Nnode + 1)
    node_ancestors[phy$edge[,2]] <- phy$edge[,1]

    node_times <- numeric(2 * phy$Nnode + 1)

    root_nodes <- phy$edge[1, 1]
    node_times[root_nodes] <- phy$root.time

    repeat {

      target_indices <- which(phy$edge[, 1] %in% root_nodes)

      if (length(target_indices) == 0) {
        break
      }

      target_nodes <- phy$edge[target_indices, 2]
      node_times[target_nodes] <- phy$edge.length[target_indices] +
        node_times[node_ancestors[target_nodes]]

      root_nodes <- target_nodes

    }

  } else {

    stop("phy is not an allowed class.")

  }

  total_nodes <- length(node_times)

  leaf_indices <- which(!(1:total_nodes)%in%node_ancestors)
  coalescence_indices <- (1:total_nodes)[-leaf_indices]

  t <- node_times[leaf_indices]
  l <- as.integer(rep(1, length(t)))
  c <- node_times[coalescence_indices]

  likelihood <- bounded_likelihood(t, l, c, ne, b, topology)

  return(likelihood)

}
