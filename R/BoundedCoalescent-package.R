## usethis namespace: start
#' @useDynLib BoundedCoalescent, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("BoundedCoalescent", libpath)
}
