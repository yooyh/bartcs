## usethis namespace: start
#' @importFrom rlang .data
#' @importFrom Rcpp sourceCpp
#' @useDynLib bartcs, .registration = TRUE
## usethis namespace: end
NULL

#' bartcs: Bayesian Additive Regression Trees for Confounder Selection
#' 
#' Fit Bayesian Regression Additive Trees (BART) models to
#' select relevant confounders among a large set of potential confounders
#' and to estimate average treatment effect. For more information, see
#' Kim et al. (2022).
#'
#' Functions in `bartcs` serve one of three purposes.
#' \enumerate{
#'   \item Functions for fitting: `sbart()` and `mbart()`.
#'   \item Functions for summary: `summary()`, `plot()` and `gelman_rubin()`.
#'   \item Utility function for openMP: `count_omp_thread()`.
#' }
#' 
#' @references
#' Kim, C., Tec, M., & Zigler, C. M. (2022).
#' Bayesian Nonparametric Adjustment of Confounding.
#' *arXiv preprint arXiv:2203.11798*.
#' \doi{10.48550/arXiv.2203.11798}
#' 
#' @docType package
#' @name bartcs-package
NULL
