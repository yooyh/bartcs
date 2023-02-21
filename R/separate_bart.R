#' @rdname bart
#' @usage NULL
#' @export
separate_bart <- function(
  Y, trt, X,
  trt_treated     = 1,
  trt_control     = 0,
  num_tree        = 50,
  num_chain       = 4,
  num_burn_in     = 100,
  num_thin        = 1,
  num_post_sample = 100,
  step_prob       = c(0.28, 0.28, 0.44),
  alpha           = 0.95,
  beta            = 2,
  nu              = 3,
  q               = 0.95,
  dir_alpha       = 5,
  parallel        = FALSE,
  verbose         = TRUE
) {

  # ---- check input ----
  check_input(
    Y, trt, X, trt_treated, trt_control,
    num_tree, num_chain,
    num_burn_in, num_thin, num_post_sample,
    step_prob, alpha, beta, nu, q,
    dir_alpha, verbose
  )

  if (length(unique(trt)) != 2)
    stop(
      "`trt` must be binary vector.\n",
      "  For non-binary `trt`, try `single_bart()` instead."
    )
  binary_trt <- isTRUE(
    all.equal(sort(unique(trt)), sort(c(trt_treated, trt_control)))
  )
  if (!binary_trt)
    stop(
      "`trt_treated` or `trt_control` must be value of `trt`.\n",
      "  For other values, try `single_bart()` instead."
    )


  # ---- data preprocessing ----
  N <- nrow(X)
  P <- ncol(X)

  # check for factor variable then change it to dummy variables
  if (sum(vapply(X, is.factor, TRUE))) {
    X <- fct_to_dummy(X)
    P <- ncol(X)
  }

  # convert to numeric vector and matrix
  if (!is.numeric(Y))
    Y <- as.numeric(Y)
  if (!is.numeric(trt))
    trt <- as.numeric(trt)
  if (!is.matrix(X))
    X <- as.matrix(X)

  # shift and rescale to [-0.5, 0.5]
  Y_max <- max(Y)
  Y_min <- min(Y)
  Y     <- (Y - Y_min) / (Y_max - Y_min) - 0.5

  # assign variable names if there are no name
  if (is.null(colnames(X)))
    colnames(X) <- paste0("X", seq_len(P))


  # ---- specific preprocessing step for separate model ----
  # separate data with respect to treatment
  Y_treated <- Y[trt == trt_treated]
  Y_control <- Y[trt == trt_control]
  X_treated <- X[trt == trt_treated, ]
  X_control <- X[trt == trt_control, ]


  # ---- calculate lambda before MCMC iterations ----
  sigma2_out1 <- stats::var(Y_treated)
  sigma2_out0 <- stats::var(Y_control)
  f <- function(lambda) {
    invgamma::qinvgamma(
      q, nu / 2, rate = lambda * nu / 2, lower.tail = TRUE, log.p = FALSE
    ) - sqrt(sigma2_out1)
  }
  lambda_out1 <- rootSolve::uniroot.all(f, c(0.1^5, 10))

  f <- function(lambda) {
    invgamma::qinvgamma(
      q, nu / 2, rate = lambda * nu / 2, lower.tail = TRUE, log.p = FALSE
    ) - sqrt(sigma2_out0)
  }
  lambda_out0 <- rootSolve::uniroot.all(f, c(0.1^5, 10))


  # ---- run MCMC and save result of each chain ----
  chains         <- list()
  num_chain_iter <- num_burn_in + num_thin * num_post_sample
  if (verbose) {
    cat(
      "Fitting ", num_chain, " chains with ", num_chain_iter, " iters each...",
      "\n\n",
      sep = ""
    )
  }

  # Call Rcpp
  chains <- cseparate_bart(
    X, Y_treated, X_treated, Y_control, X_control, trt, Y_min, Y_max, 
    step_prob, num_chain, num_chain_iter, num_burn_in, num_thin, num_post_sample,
    num_tree, alpha, beta, nu, lambda_out1, lambda_out0, 
    dir_alpha, sigma2_out1, sigma2_out0, parallel, verbose
  )


  # ---- post processing ----
  names(chains) <- paste0("chain", seq_len(num_chain))

  # merge result
  ATE <- Y1 <- Y0 <- NULL
  var_prob <- vector(mode = "numeric", length = P)
  for (chain_idx in seq_len(num_chain)) {
    ATE <- c(ATE, chains[[chain_idx]]$ATE)
    Y1  <- c(Y1,  chains[[chain_idx]]$Y1)
    Y0  <- c(Y0,  chains[[chain_idx]]$Y0)
    var_prob <- var_prob + chains[[chain_idx]]$var_prob
  }
  var_prob        <- var_prob / num_chain
  names(var_prob) <- colnames(X)

  # return as bartcs object
  structure(
    list(
      ATE = ATE, Y1 = Y1, Y0 = Y0,
      var_prob = var_prob,
      chains   = chains,
      model    = "separate",
      label    = colnames(X),
      params   = list(
        trt_treated     = trt_treated,
        trt_control     = trt_control,
        num_tree        = num_tree,
        num_chain_iter  = num_chain_iter,
        num_chain       = num_chain,
        num_burn_in     = num_burn_in,
        num_thin        = num_thin,
        num_post_sample = num_post_sample,
        step_prob       = step_prob,
        alpha           = alpha,
        beta            = beta,
        nu              = nu,
        q               = q
      )
    ),
    class = "bartcs"
  )
}
