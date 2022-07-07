## usethis namespace: start
#' @importFrom rlang .data
#' @importFrom Rcpp sourceCpp
#' @useDynLib bartcs, .registration = TRUE
## usethis namespace: end
NULL

#' bartcs: Bayesian Nonparametric Adjustment of Confounding
#'
#' Functions in `bartcs` serve one of three purposes.
#' \enumerate{
#'   \item Functions for fitting: `sbart()`, `mbart()`.
#'   \item Functions for summary: `summary()`, `gelman_rubin()`, `plot()`
#'   \item Util function for openMP: `count_omp_thread()`
#' }
#' @docType package
#' @name bartcs-package
NULL
