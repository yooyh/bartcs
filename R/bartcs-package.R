## usethis namespace: start
#' @importFrom rlang .data
#' @importFrom Rcpp sourceCpp
#' @useDynLib bartcs, .registration = TRUE
## usethis namespace: end
NULL

#' bartcs: Bayesian Additive Regression Trees for Confounder Selection
#' 
#' Fit Bayesian Regression Additive Trees (BART) models to
#' select true confounders from a large set of potential confounders and
#' to estimate average treatment effect. For more information, see
#' Kim et al. (2023) \doi{10.1111/biom.13833}.
#'
#' Functions in `bartcs` serve one of three purposes.
#' \enumerate{
#'   \item Functions for fitting: `separate_bart()` and `single_bart()`.
#'   \item Functions for summary: `summary()`, `plot()` and `gelman_rubin()`.
#'   \item Utility function for OpenMP: `count_omp_thread()`.
#' }
#' The code of BART model are based on the 'BART' package by 
#' Sparapani et al. (2021) under the GPL license, with modifications.
#' The modifications from the `BART` package include (but are not limited to):
#' \itemize{
#'   \item Add CHANGE step.
#'   \item Add Single and Separate Model.
#'   \item Add causal effect estimation.
#'   \item Add confounder selection.
#' }
#' 
#' @references
#' Sparapani R, Spanbauer C, McCulloch R (2021).
#' “Nonparametric Machine Learning and Efficient Computation
#' with Bayesian Additive Regression Trees: The BART R Package.”
#' *Journal of Statistical Software*, 97(1), 1–66.
#' \doi{10.18637/jss.v097.i01}
#'
#' Kim, C., Tec, M., & Zigler, C. M. (2023).
#' Bayesian Nonparametric Adjustment of Confounding, *Biometrics*
#' \doi{10.1111/biom.13833}
#' 
#' @docType package
#' @name bartcs-package
NULL
