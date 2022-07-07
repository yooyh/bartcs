#' @rdname bart
#' @usage NULL
#' @export
mbart <- function(
  Y, trt, X,
  trt_treated     = 1,
  trt_control     = 0,
  num_tree        = 50,
  num_chain       = 4,
  num_burn_in     = 100,
  num_thin        = 0,
  num_post_sample = 100,
  step_prob       = c(0.28, 0.28, 0.44),
  alpha           = 0.95,
  beta            = 2,
  nu              = 3,
  q               = 0.95,
  dir_alpha       = 5,
  boot_size       = NULL,
  parallel        = NULL,
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


  # ---- data preprocessing ----
  n <- nrow(X)
  p <- ncol(X)

  # check for factor variable then change it to dummy variables
  if (sum(vapply(X, is.factor, TRUE))) {
    X <- fct_to_dummy(X)
    p <- ncol(X)
  }

  # convert to numeric vector and matrix
  if (!is.numeric(Y))
    Y <- as.numeric(Y)
  if (!is.numeric(trt))
    trt <- as.numeric(trt)
  if (!is.matrix(X))
    X <- as.matrix(X)

  # shift and rescale to [-0.5, 0.5]
  Y_mean <- mean(Y)
  Y      <- Y - Y_mean
  Y_max  <- max(Y)
  Y_min  <- min(Y)
  Y      <- (Y - Y_min) / (Y_max - Y_min) - 0.5

  # initialize bootstrap sample size and parallel
  if (is.null(boot_size))
    boot_size <- 2 * n
  if (is.null(parallel))
    parallel <- ifelse(n < 15e4, TRUE, FALSE)

  # assign variable names if there are no name
  if (is.null(colnames(X)))
    colnames(X) <- paste0("X", seq_len(p))


  # ---- mbart specific preprocessing step ----
  # check whether it is binary treatment
  is_binary_trt <- isTRUE(
    all.equal(sort(unique(trt)), sort(c(trt_treated, trt_control)))
  )


  # ---- calculate lambda before MCMC iterations ----
  sigma2_exp <- ifelse(is_binary_trt, 1, stats::var(Y))
  sigma2_out <- stats::var(Y)
  if (is_binary_trt) {
    lambda_exp <- 0 # arbitrary value
  } else {
    f <- function(lambda) {
      invgamma::qinvgamma(
        q, nu / 2, rate = lambda * nu / 2, lower.tail = TRUE, log.p = FALSE
      ) - sqrt(sigma2_exp)
    }
    lambda_exp <- rootSolve::uniroot.all(f, c(0.1^5, 10))
  }

  f <- function(lambda) {
    invgamma::qinvgamma(
      q, nu / 2, rate = lambda * nu / 2, lower.tail = TRUE, log.p = FALSE
    ) - sqrt(sigma2_out)
  }
  lambda_out <- rootSolve::uniroot.all(f, c(0.1^5, 10))


  # ---- run MCMC and save result of each chain ----
  chains         <- list()
  num_chain_iter <- num_burn_in + (num_thin + 1) * num_post_sample
  if (verbose) {
    cat(
      "\n",
      "Fitting ", num_chain, " chains with ", num_chain_iter, " iters each...",
      "\n\n",
      sep = ""
    )
  }

  for (chain_idx in seq_len(num_chain)) {
    # placeholder for MCMC samples
    # Y1 and Y0 are potential outcomes. ATE = Y1 - Y0
    Y1  <- vector(mode = "numeric", length = num_post_sample)
    Y0  <- vector(mode = "numeric", length = num_post_sample)

    # placeholder for inclusion probabilities and initialize var_prob
    var_count <- matrix(0, nrow = num_post_sample, ncol = p + 1)
    var_prob  <- MCMCpack::rdirichlet(1, rep(dir_alpha, p + 1))

    # placeholder for parameters
    sigma2_out_hist    <- vector(mode = "numeric", length = num_chain_iter + 1)
    dir_alpha_hist     <- vector(mode = "numeric", length = num_chain_iter + 1)
    sigma2_out_hist[1] <- stats::var(Y)
    dir_alpha_hist[1]  <- dir_alpha

    # call Rcpp implementation
    fit_mbart(
      Y1, Y0, var_count, var_prob,
      sigma2_exp, sigma2_out_hist, dir_alpha_hist,
      Y, trt, X, trt_treated, trt_control,
      chain_idx, num_chain,
      num_chain_iter, num_burn_in, num_thin, num_post_sample,
      num_tree, step_prob, alpha, beta,
      nu, lambda_exp, lambda_out,
      boot_size, is_binary_trt, parallel, verbose
    )

    # post-processing for each MCMC chain
    var_prob <- colMeans(ifelse(var_count > 1, 1, 0))

    # rescale result
    Y1  <- (Y1  + 0.5) * (Y_max - Y_min) + Y_min + Y_mean
    Y0  <- (Y0  + 0.5) * (Y_max - Y_min) + Y_min + Y_mean
    ATE <-  Y1 - Y0

    chains[[chain_idx]] <- list(
      ATE = ATE, Y1 = Y1, Y0 = Y0,
      var_count  = var_count,
      var_prob   = var_prob,
      sigma2_out = sigma2_out_hist,
      dir_alpha  = dir_alpha_hist
    )
  }


  # ---- post processing ----
  names(chains) <- paste0("chain", seq_len(num_chain))

  # merge result
  ATE <- Y1 <- Y0 <- NULL
  var_prob <- vector(mode = "numeric", length = p + 1)
  for (chain_idx in seq_len(num_chain)) {
    ATE <- c(ATE, chains[[chain_idx]]$ATE)
    Y1  <- c(Y1,  chains[[chain_idx]]$Y1)
    Y0  <- c(Y0,  chains[[chain_idx]]$Y0)
    var_prob <- var_prob + chains[[chain_idx]]$var_prob
  }
  var_prob        <- var_prob / num_chain
  names(var_prob) <- c(colnames(X), "trt")

  cat("\n")

  # return as bartcs object
  structure(
    list(
      ATE = ATE, Y1 = Y1, Y0 = Y0,
      var_prob = var_prob,
      chains   = chains,
      model    = "mbart",
      label    = c(colnames(X), "trt"),
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
