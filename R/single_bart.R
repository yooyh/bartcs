#' @rdname bart
#' @usage NULL
#' @export
single_bart <- function(
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


  # ---- specific preprocessing step for single model ----
  # check whether it is binary treatment
  binary_trt <- isTRUE(
    all.equal(sort(unique(trt)), sort(c(trt_treated, trt_control)))
  )


  # ---- calculate lambda before MCMC iterations ----
  sigma2_exp <- ifelse(binary_trt, 1, stats::var(Y))
  sigma2_out <- stats::var(Y)
  if (binary_trt) {
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
  num_chain_iter <- num_burn_in + num_thin * num_post_sample
  if (verbose) {
    cat(
      "Fitting ", num_chain, " chains with ", num_chain_iter, " iters each...",
      "\n\n",
      sep = ""
    )
  }

  # Call Rcpp
  chains <- csingle_bart(
    Y, X, trt, trt_treated, trt_control, Y_min, Y_max,
    step_prob, num_chain, num_chain_iter, num_burn_in, num_thin, num_post_sample,
    num_tree, alpha, beta, nu, lambda_exp, lambda_out,
    dir_alpha, sigma2_exp, sigma2_out, binary_trt, parallel, verbose
  )


  # ---- post processing ----
  names(chains) <- paste0("chain", seq_len(num_chain))

  # merge result
  ATE <- Y1 <- Y0 <- NULL
  var_prob <- vector(mode = "numeric", length = P + 1)
  for (chain_idx in seq_len(num_chain)) {
    ATE <- c(ATE, chains[[chain_idx]]$ATE)
    Y1  <- c(Y1,  chains[[chain_idx]]$Y1)
    Y0  <- c(Y0,  chains[[chain_idx]]$Y0)
    var_prob <- var_prob + chains[[chain_idx]]$var_prob
  }
  var_prob        <- var_prob / num_chain
  names(var_prob) <- c("trt", colnames(X))

  # return as bartcs object
  structure(
    list(
      ATE = ATE, Y1 = Y1, Y0 = Y0,
      var_prob = var_prob,
      chains   = chains,
      model    = "single",
      label    = c("trt", colnames(X)),
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
