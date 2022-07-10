#' @name bart
#' @title Fit BART models to select confounders and estimate treatment effect
#'
#' @description
#' Fit Bayesian Regression Additive Trees (BART) models
#'   to select relevant confounders among a large set of potential confounders
#'   and to estimate average treatment effect \eqn{(Y(1) - Y(0))}.
#'
#' @usage
#' sbart(
#'   Y, trt, X,
#'   trt_treated     = 1,
#'   trt_control     = 0,
#'   num_tree        = 50,
#'   num_chain       = 4,
#'   num_burn_in     = 100,
#'   num_thin        = 0,
#'   num_post_sample = 100,
#'   step_prob       = c(0.28, 0.28, 0.44),
#'   alpha           = 0.95,
#'   beta            = 2,
#'   nu              = 3,
#'   q               = 0.95,
#'   dir_alpha       = 5,
#'   boot_size       = NULL,
#'   parallel        = NULL,
#'   verbose         = TRUE
#' )
#'
#' mbart(
#'   Y, trt, X,
#'   trt_treated     = 1,
#'   trt_control     = 0,
#'   num_tree        = 50,
#'   num_chain       = 4,
#'   num_burn_in     = 100,
#'   num_thin        = 0,
#'   num_post_sample = 100,
#'   step_prob       = c(0.28, 0.28, 0.44),
#'   alpha           = 0.95,
#'   beta            = 2,
#'   nu              = 3,
#'   q               = 0.95,
#'   dir_alpha       = 5,
#'   boot_size       = NULL,
#'   parallel        = NULL,
#'   verbose         = TRUE
#' )
#'
#' @param Y Outcome variable.
#' @param trt Treatment variable.
#' @param X Potential confounders.
#' @param trt_treated Value of `trt` for treated group.
#' @param trt_control Value of `trt` for control group.
#' @param num_tree Number of trees in BART model.
#' @param num_chain Number of MCMC chains.
#'   Need to set `num_chain > 1` for Gelman-Rubin diagnostic.
#' @param num_burn_in Number of MCMC samples to be discarded per chain
#'   as initial burn-in periods.
#' @param num_thin Number of thinning per chain.
#'   One in every `num_thin` samples are selected.
#' @param num_post_sample Final number of posterior samples per chain.
#'   Number of MCMC iterations per chain is
#'   `burn_in + num_thin * num_post_sample`.
#' @param step_prob A vector of tree alteration probabilities (GROW, PRUNE, CHANGE).
#'   Each alteration is proposed to change the tree structure.
#'   Default setting is `(0.28, 0.28, 0.44)`.
#' @param alpha,beta Hyperparameters for tree regularization prior.
#'   A terminal node of depth `d` will split with
#'   probability of `alpha * (1 + d)^(-beta)`.
#' @param nu,q Values to calibrate hyperparameter of sigma prior.
#'   Default setting is `(nu, q) = (3, 0.95)` from Chipman et al. (2010).
#' @param dir_alpha Hyperparameter of Dirichlet prior for selection probabilities.
#' @param boot_size Number of bootstrap sample size.
#'   Bootstrap samples will be used to compute potential outcomes \eqn{Y(1)} and \eqn{Y(0)}.
#' @param parallel If `TRUE`, model fitting will be parallelized
#'   with respect to `n = nrow(X)`.
#'   Parallelization is recommended for very high `n` only.
#' @param verbose If `TRUE`, message will be printed during training.
#'   If `FALSE`, message will be suppressed.
#'
#' @details
#' `sbart()` and `mbart()` fit an exposure model and outcome model(s)
#' for estimating treatment effect with adjustment of confounders
#' in the presence of a large set of potential confounders (Kim et al. 2022).
#'
#' The exposure model \eqn{E[A|X]} and the outcome model(s) \eqn{E[Y|A,X]} are
#' linked together with a common Dirichlet prior that accrues
#' posterior selection probability to confounders (\eqn{X}) on the basis of
#' association with both the exposure (\eqn{A}) and the outcome (\eqn{Y}).
#'
#' There is a distinction between fitting each outcome model for the treated and control groups and
#' fitting a single outcome model for the entire sample.
#' \itemize{
#'   \item `sbart()` specifies two **"separate"** outcome models
#'     for two binary treatment levels.
#'     Thus, it fits three models:
#'       one exposure model and two separate outcome models for \eqn{A = 0, 1}.
#'
#'   \item `mbart()` specifies a single **"marginal"** outcome models.
#'     Thus, it fits two models:
#'       one exposure model and one outcome model for the entire sample.
#' }
#'
#' All inferences are made with outcome model(s).
#'
#' @references
#' Chipman, H. A., George, E. I., & McCulloch, R. E. (2010).
#' BART: Bayesian additive regression trees.
#' *The Annals of Applied Statistics, 4*(1), 266-298.
#' \doi{10.1214/09-AOAS285}
#'
#' Kim, C., Tec, M., & Zigler, C. M. (2022).
#' Bayesian Nonparametric Adjustment of Confounding.
#' *arXiv preprint arXiv:2203.11798*.
#' \doi{10.48550/arXiv.2203.11798}
#'
#' @return
#' \code{bartcs} object is a list with following components.
#'
#' \item{ATE}{Aggregated posterior samples of average treatment effect \eqn{(Y(1) - Y(0))}.}
#' \item{Y1}{Aggregated posterior samples of potential outcome \eqn{Y(1)}.}
#' \item{Y0}{Aggregated posterior samples of potential outcome \eqn{Y(0)}.}
#' \item{var_prob}{Aggregated posterior inclusion probability of each variable.}
#' \item{chains}{A list of results from each MCMC chain.
#'   Each list element consists of followings.}
#'   \itemize{
#'     \item \code{ATE}        Posterior sample of average treatment effect \eqn{(Y(1) - Y(0))}.
#'     \item \code{Y1}         Posterior sample of potential outcome \eqn{Y(1)}.
#'     \item \code{Y0}         Posterior sample of potential outcome \eqn{Y(0)}.
#'     \item \code{var_prob}   Posterior inclusion probability of each variable.
#'     \item \code{var_count}  Number of selection of each variable in each MCMC iteration.
#'       Its dimension is `num_post_sample * ncol(X)`.
#'     \item \code{sigma2_out} Posterior sample of `sigma2` in the outcome model.
#'     \item \code{dir_alpha}  Posterior sample of `dir_alpha.`
#'   }
#' \item{model}{`sbart` or `mbart`.}
#' \item{label}{Column names of `X`.}
#' \item{params}{Parameters used in the model.}
#'
NULL
