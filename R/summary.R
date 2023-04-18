#' Summary for `bartcs` object
#'
#' @description
#' Provide summary for `bartcs` object.
#'
#' @param object A `bartcs` object.
#' @param ... Additional arguments. Not yet supported.
#'
#' @details
#' `summary()` provides 95% posterior credible interval for both
#' aggregated outcome and individual outcomes from each MCMC chain.
#'
#' @return
#' Provide list with the following components
#'
#' \item{model}{`separate_bart` or `single_bart`.}
#' \item{trt_value}{Treatment values for each treatment group:
#'   `trt_treated` for the treatment group and `trt_control` for
#' the control group.}
#' \item{tree_params}{Parameters for the tree structure.}
#' \item{chain_params}{Parameters for MCMC chains.}
#' \item{outcome}{Summary of outcomes from the model. This includes
#'   both aggregated outcome and individual outcomes from each MCMC chain.}
#'
#' @examples
#' data(ihdp, package = "bartcs")
#' x <- single_bart(
#'   Y               = ihdp$y_factual,
#'   trt             = ihdp$treatment,
#'   X               = ihdp[, 6:30],
#'   num_tree        = 10,
#'   num_chain       = 2,
#'   num_post_sample = 20,
#'   num_burn_in     = 10,
#'   verbose         = FALSE
#' )
#' summary(x)
#'
#' @exportS3Method
summary.bartcs <- function(object, ...) {
  estimand  <- c("ATE", "Y1", "Y0")
  num_chain <- object$params$num_chain

  res <- list()
  res$model     <- object$model
  res$trt_value <- c(object$params$trt_treated, object$params$trt_control)

  res$tree_params  <- list(
    num_tree  = object$params$num_tree,
    step_prob = object$params$step_prob,
    alpha     = object$params$alpha,
    beta      = object$params$beta,
    nu        = object$params$nu,
    q         = object$params$q
  )

  res$chain_params <- list(
    num_chain       = object$params$num_chain,
    num_chain_iter  = object$params$num_chain_iter,
    num_post_sample = object$params$num_post_sample,
    num_burn_in     = object$params$num_burn_in,
    num_thin        = object$params$num_thin
  )

  outcome <- data.frame(
    estimand    = character(),
    chain       = factor(levels = c(seq_len(num_chain), "agg")),
    `2.5%`      = numeric(),
    `1Q`        = numeric(),
    mean        = numeric(),
    median      = numeric(),
    `3Q`        = numeric(),
    `97.5%`     = numeric(),
    check.names = FALSE
  )
  outcome[1:(3 * (num_chain + 1)), 1] <- rep(estimand, each = num_chain + 1)
  outcome[1:(3 * (num_chain + 1)), 2] <- rep(c(seq_len(num_chain), "agg"), 3)
  for (i in seq_along(estimand)) {
    idx <- (i - 1) * (num_chain + 1)
    est <- estimand[i]
    outcome[(idx + 1):(idx + num_chain), 3:ncol(outcome)] <- t(vapply(
      seq_len(num_chain),
      function(chain_idx) {t(c(
        stats::quantile(object$mcmc_outcome[[chain_idx]][, est],
                        probs = c(0.025, 0.25)),
        mean(object$mcmc_outcome[[chain_idx]][, i]),
        stats::quantile(object$mcmc_outcome[[chain_idx]][, est],
                        probs = c(0.5, 0.75, 0.975))
      ))},
      numeric(6)
    ))
    mat <- do.call("rbind", object$mcmc_outcome)
    outcome[idx + num_chain + 1, 3:ncol(outcome)] <- c(
      stats::quantile(mat[, est], probs = c(0.025, 0.25)),
      mean(mat[, est]),
      stats::quantile(mat[, est], probs = c(0.5, 0.75, 0.975))
    )
  }
  res$outcome <- outcome

  class(res) <- "bartcs_summary"
  res
}

#' @exportS3Method
print.bartcs_summary <- function(x, ...) {
  width = 6
  cat(
    "`bartcs` fit by `", x$model, "_bart()`", "\n",
    "\n", sep = ""
  )
  cat(
    "Treatment Value\n",
    "  Treated group    : ", format(x$trt_value[1], width = width), "\n",
    "  Control group    : ", format(x$trt_value[2], width = width), "\n",
    "\n", sep = ""
  )

  # Tree Summary
  cat(
    "Tree Parameters\n",
    "  Number of Tree   : ",   format(x$tree_params$num_tree,     width = width), "\t",
    "\tValue  of alpha    : ", format(x$tree_params$alpha,        width = width), "\n",
    "  Prob.  of Grow   : ",   format(x$tree_params$step_prob[1], width = width), "\t",
    "\tValue  of beta     : ", format(x$tree_params$beta,         width = width), "\n",
    "  Prob.  of Prune  : ",   format(x$tree_params$step_prob[2], width = width), "\t",
    "\tValue  of nu       : ", format(x$tree_params$nu,           width = width), "\n",
    "  Prob.  of Change : ",   format(x$tree_params$step_prob[3], width = width), "\t",
    "\tValue  of q        : ", format(x$tree_params$q,            width = width), "\n",
    "\n", sep = ""
  )

  # Chain Summary
  cat(
    "Chain Parameters\n",
    "  Number of Chains : ",   format(x$chain_params$num_chain,       width = width), "\t",
    "\tNumber of burn-in  : ", format(x$chain_params$num_burn_in,     width = width), "\n",
    "  Number of Iter   : ",   format(x$chain_params$num_chain_iter,  width = width),  "\t",
    "\tNumber of thinning : ", format(x$chain_params$num_thin,        width = width), "\n",
    "  Number of Sample : ",   format(x$chain_params$num_post_sample, width = width), "\n",
    "\n", sep = ""
  )

  cat(
    "Outcome \n"
  )
  print(x$outcome, row.names = FALSE)
}
