#' @name bart
#' @title Fit BART models to select confounders and estimate treatment effect
#'
#' @description
#' Fit Bayesian Regression Additive Trees (BART) models
#'   to select relevant confounders among a large set of potential confounders
#'   and to estimate average treatment effect \eqn{E[Y(1) - Y(0)]}.
#'
#' @usage
#' separate_bart(
#'   Y, trt, X,
#'   trt_treated     = 1,
#'   trt_control     = 0,
#'   num_tree        = 50,
#'   num_chain       = 4,
#'   num_burn_in     = 100,
#'   num_thin        = 1,
#'   num_post_sample = 100,
#'   step_prob       = c(0.28, 0.28, 0.44),
#'   alpha           = 0.95,
#'   beta            = 2,
#'   nu              = 3,
#'   q               = 0.95,
#'   dir_alpha       = 5,
#'   parallel        = FALSE,
#'   verbose         = TRUE
#' )
#'
#' single_bart(
#'   Y, trt, X,
#'   trt_treated     = 1,
#'   trt_control     = 0,
#'   num_tree        = 50,
#'   num_chain       = 4,
#'   num_burn_in     = 100,
#'   num_thin        = 1,
#'   num_post_sample = 100,
#'   step_prob       = c(0.28, 0.28, 0.44),
#'   alpha           = 0.95,
#'   beta            = 2,
#'   nu              = 3,
#'   q               = 0.95,
#'   dir_alpha       = 5,
#'   parallel        = FALSE,
#'   verbose         = TRUE
#' )
#'
#' @param Y A vector of outcome values.
#' @param trt A vector of treatment values. Binary treatment works for both model
#'   and continuous treatment works for single_bart(). For binary treatment,
#'   use 1 to indicate the treated group and 0 for the control group.
#' @param X A matrix of potential confounders.
#' @param trt_treated Value of `trt` for the treated group.
#'   The default value is set to 1.
#' @param trt_control Value of `trt` for the control group.
#'   The default value is set to 0.
#' @param num_tree Number of trees in BART model. The default value is set to 100.
#' @param num_chain Number of MCMC chains.
#'   Need to set `num_chain > 1` for the Gelman-Rubin diagnostic.
#'   The default value is set to 4.
#' @param num_burn_in Number of MCMC samples to be discarded per chain
#'   as initial burn-in periods.
#'   The default value is set to 100.
#' @param num_thin Number of thinning per chain.
#'   One in every `num_thin` samples are selected.
#'   The default value is set to 1.
#' @param num_post_sample Final number of posterior samples per chain.
#'   Number of MCMC iterations per chain is
#'   `burn_in + num_thin * num_post_sample`.
#'   The default value is set to 100.
#' @param step_prob A vector of tree alteration probabilities (GROW, PRUNE, CHANGE).
#'   Each alteration is proposed to change the tree structure. \cr
#'   The default setting is `(0.28, 0.28, 0.44)`.
#' @param alpha,beta Hyperparameters for tree regularization prior.
#'   A terminal node of depth `d` will split with
#'   probability of `alpha * (1 + d)^(-beta)`. \cr
#'   The default setting is
#'   `(alpha, beta) = (0.95, 2)` from Chipman et al. (2010).
#' @param nu,q Values to calibrate hyperparameter of sigma prior. \cr
#'   The default setting is `(nu, q) = (3, 0.95)` from Chipman et al. (2010).
#' @param dir_alpha Hyperparameter of Dirichlet prior for selection probabilities.
#'   The default value is 5.
#' @param parallel If `TRUE`, model fitting will be parallelized
#'   with respect to `N = nrow(X)`.
#'   Parallelization is recommended for very high `n` only.
#'   The default setting is FALSE.
#' @param verbose If `TRUE`, message will be printed during training.
#'   If `FALSE`, message will be suppressed.
#'
#' @details
#' `separate_bart()` and `single_bart()` fit an exposure model and outcome model(s)
#' for estimating treatment effect with adjustment of confounders
#' in the presence of a large set of potential confounders (Kim et al. 2023).
#'
#' The exposure model \eqn{E[A|X]} and the outcome model(s) \eqn{E[Y|A,X]} are
#' linked together with a common Dirichlet prior that accrues
#' posterior selection probabilities to the corresponding confounders (\eqn{X}) 
#' on the basis of association with both the exposure (\eqn{A}) and the 
#' outcome (\eqn{Y}).
#'
#' There is a distinction between fitting separate outcome models for the treated 
#' and control groups and fitting a single outcome model for both groups.
#' \itemize{
#'   \item `separate_bart()` specifies two **"separate"** outcome models
#'     for two binary treatment levels.
#'     Thus, it fits three models:
#'       one exposure model and two separate outcome models for \eqn{A = 0, 1}.
#'
#'   \item `single_bart()` specifies one **"single"** outcome model.
#'     Thus, it fits two models:
#'       one exposure model and one outcome model for the entire sample.
#' }
#'
#' All inferences are made with outcome model(s).
#'
#' @references
#' Chipman, H. A., George, E. I., & McCulloch, R. E. (2010).
#' BART: Bayesian additive regression trees.
#' *The Annals of Applied Statistics*, 4(1), 266-298.
#' \doi{10.1214/09-AOAS285}
#'
#' Kim, C., Tec, M., & Zigler, C. M. (2023).
#' Bayesian Nonparametric Adjustment of Confounding, *Biometrics*
#' \doi{10.1111/biom.13833}
#'
#' @return
#' A `bartcs` object. A list object contains the following components.
#' 
#' \item{mcmc_outcome}{A `mcmc.list` object from \pkg{coda} package.
#'   `mcmc_outcome` contains the following items.}
#'   \itemize{
#'     \item `ATE`        Posterior sample of average treatment effect \eqn{E[Y(1) - Y(0)]}.
#'     \item `Y1`         Posterior sample of potential outcome \eqn{E[Y(1)]}.
#'     \item `Y0`         Posterior sample of potential outcome \eqn{E[Y(0)]}.
#'   }
#' \item{mcmc_param}{A `mcmc.list` object from \pkg{coda} package.
#'   `mcmc_param` contains the following items.}
#'   \itemize{
#'     \item `dir_alpha`  Posterior sample of `dir_alpha.`
#'     \item `sigma2_out` Posterior sample of `sigma2` in the outcome model.
#'   }
#' \item{var_prob}{Aggregated posterior inclusion probability of each variable.}
#' \item{var_count}{Number of selection of each variable in each MCMC iteration.
#'   Its dimension is `num_post_sample * ncol(X)`.}
#' \item{chains}{A list of results from each MCMC chain.}
#' \item{model}{`separate` or `single`.}
#' \item{label}{Column names of `X`.}
#' \item{params}{Parameters used in the model.}
#'
#' @examples
#' data(ihdp, package = "bartcs")
#' single_bart(
#'   Y               = ihdp$y_factual,
#'   trt             = ihdp$treatment,
#'   X               = ihdp[, 6:30],
#'   num_tree        = 10,
#'   num_chain       = 2,
#'   num_post_sample = 20,
#'   num_burn_in     = 10,
#'   verbose         = FALSE
#' )
#' separate_bart(
#'   Y               = ihdp$y_factual,
#'   trt             = ihdp$treatment,
#'   X               = ihdp[, 6:30],
#'   num_tree        = 10,
#'   num_chain       = 2,
#'   num_post_sample = 20,
#'   num_burn_in     = 10,
#'   verbose         = FALSE
#' )
NULL
