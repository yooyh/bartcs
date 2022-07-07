# utils.R

# change factor variables to dummy variables
fct_to_dummy <- function(X) {
  fct_idx <- (seq_len(ncol(X)))[vapply(X, is.factor, TRUE)]

  for (i in fct_idx) {
    fct_levels <- levels(X[[i]])[-1]
    for (level in fct_levels) {
      var_name <- paste0(colnames(X)[i], "_", level, sep = "")
      X[[var_name]] <- ifelse(X[[i]] == level, 1, 0)
    }
  }
  X[, -fct_idx]
}

# test inputs for sbart and mbart
check_input <- function(
  Y, trt, X,
  trt_treated     = 1,
  trt_control     = 0,
  num_tree        = 50,
  num_chain       = 4,
  num_burn_in     = 250,
  num_thin        = 0,
  num_post_sample = 1000,
  step_prob       = c(0.28, 0.28, 0.44),
  alpha           = 0.95,
  beta            = 2,
  nu              = 3,
  q               = 0.95,
  dir_alpha       = 5,
  verbose         = TRUE
) {

  # type check ----
  if (!is.logical(verbose))
    stop("`verbose` must be boolean.")

  if (!is.logical(Y)   && !is.numeric(Y))
    stop("`Y` must be a numeric vector.")

  if (!is.logical(trt) && !is.numeric(trt))
    stop("`trt` must be a numeric vector.")

  if (!is.matrix(X)) {
    if (is.data.frame(X)) {
      ncol_numeric <- sum(vapply(X, is.numeric, TRUE))
      ncol_logical <- sum(vapply(X, is.logical, TRUE))
      ncol_factor  <- sum(vapply(X, is.factor, TRUE))
      if (ncol_numeric + ncol_logical + ncol_factor != ncol(X))
        stop("Columns of `X` must be numeric or logical.")
    } else {
      stop("`X` must be matrix or data.frame.")
    }
  }

  if (!is.numeric(trt_treated) && is.logical(trt_treated)) 
    stop("`trt_treated` must be numeric.")
  if (!is.numeric(trt_control) && is.logical(trt_control))
    stop("`trt_control` must be numeric.")

  if (!is.numeric(num_tree))
    stop("`num_tree` must be numeric.")
  if (!is.numeric(num_chain))
    stop("`num_chain` must be numeric.")
  if (!is.numeric(num_burn_in))
    stop("`num_burn_in` must be numeric.")
  if (!is.numeric(num_thin))
    stop("`num_thin` must be numeric.")
  if (!is.numeric(num_post_sample))
    stop("`num_post_sample` must be numeric.")

  if (!is.numeric(step_prob))
    stop("`step_prob` must be numeric.")
  if (!is.numeric(alpha))
    stop("`alpha` must be numeric.")
  if (!is.numeric(beta))
    stop("`beta` must be numeric.")
  if (!is.numeric(nu))
    stop("`nu` must be numeric.")
  if (!is.numeric(q))
    stop("`q` must be numeric.")
  if (!is.numeric(dir_alpha))
    stop("`dir_alpha` must be numeric.")


  # value check
  if (num_tree <= 0)
    stop("`num_tree` must be greater than 0.")

  if (num_chain <= 0)
    stop("`num_tree` must be greater than 0.")

  if (num_burn_in < 0)
    stop("`num_burn_in` must be non-negative.")

  if (num_thin < 0)
    stop("`num_thin` must be non-negative.")

  if (num_post_sample < 0)
    stop("`num_post_sample` must be greater than 0.")

  if (alpha < 0 || alpha > 1)
    stop("`alpha` must be between 0 and 1.")

  if (beta < 0)
    stop("`beta` must be greater than 0.")

  if (q < 0 || q > 1)
    stop("`q` must be between 0 and 1.")

  if (dir_alpha <= 0)
    stop("`dir_alpha` must be greater than 0.")


  # check dimension
  n <- nrow(X)
  p <- ncol(X)
  if (length(Y) != n || length(trt) != n)
    stop("Length of `Y`, number of rows of `X` and length of `trt` must match.")
  
  if (length(unique(trt)) < 2)
    stop("`trt` must have at least 2 unique values.")

  if (length(step_prob) != 3)
    stop("Length of `step_prob` must be 3.")
}
