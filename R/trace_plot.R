trace_plot <- function(x, parameter) {
  if (length(parameter) != 1)
    stop("`trace_plot` can take single parameter.\n")

  if (x$model == "single" && parameter %in% c("sigma2_out1", "sigma2_out0")) {
    message(
      "Current parameter is", parameter, ". Did you mean `sigma2_out`?, \n",
      "  Proceeding with `sigma2_out`..."
    )
    parameter <- "sigma2_out"
  }
  if (x$model == "separate" && parameter == "sigma2_out") {
    message(
      "There are `sigma2_out1` and `sigma2_out0`.\n",
      "  Proceeding with `sigma2_out1`..."
    )
    parameter <- "sigma2_out1"
  }
  if (parameter == "alpha")
    parameter <- "dir_alpha"

  outcome <- c("ATE", "Y1", "Y0")
  if (x$model == "separate")
    params <- c("sigma2_out1", "sigma2_out0", "dir_alpha")
  else if (x$model == "single")
    params <- c("sigma2_out", "dir_alpha")


  if (parameter %in% outcome) {
    num_post_sample <- x$params$num_post_sample
    num_chain       <- x$params$num_chain

    # gather data
    df <- data.frame(
      iter   = rep(seq_len(num_post_sample), num_chain),
      chains = rep(paste0("chain", seq_len(num_chain)), each = num_post_sample)
    )
    df[[parameter]] <- x[[parameter]]

    # draw plot
    if (parameter == "ATE") {
      title <- "Traceplot of ATE"
      ylab  <- "ATE"
    } else {
      title <- paste0("Traceplot of Potential Outcome ", parameter, sep = "")
      ylab  <- parameter
    }


    res <- ggplot2::ggplot(data = df) +
      ggplot2::xlim(0, num_post_sample) +
      ggplot2::labs(y = ylab, x = "Iteration", title = title) +
      ggplot2::geom_path(
        mapping = ggplot2::aes(x     = .data$iter,
                               y     = .data[[parameter]],
                               color = .data$chains)
      )

    return(res)

  } else if (parameter %in% params) {
    num_chain_iter <- x$params$num_chain_iter + 1
    num_chain      <- x$params$num_chain
    num_burn_in    <- x$params$num_burn_in

    # gather data
    df <- data.frame(
      iter   = rep(seq_len(num_chain_iter), num_chain),
      chains = rep(paste0("chain", seq_len(num_chain)), each = num_chain_iter)
    )
    df[[parameter]] <- as.vector(vapply(
      seq_len(num_chain),
      function(chain_idx) x$chains[[chain_idx]][[parameter]],
      numeric(num_chain_iter)
    ))

    ylab  <- ifelse(parameter == "dir_alpha", "Alpha", "Sigma 2")
    title <- NULL
    if (parameter == "sigma2_out")
      title <- "Traceplot of Sigma2 of Outcome Model"
    else if (parameter == "sigma2_out1")
      title <- "Traceplot of Sigma2 of Outcome Model 1"
    else if (parameter == "sigma2_out0")
      title <- "Traceplot of Sigma2 of Outcome Model 0"
    else if (parameter == "dir_alpha")
      title <- "Traceplot of Alpha"


    res <- ggplot2::ggplot(data = df) +
      ggplot2::xlim(0, num_chain_iter) +
      ggplot2::labs(y = ylab, x = "Iteration", title = title) +
      ggplot2::geom_vline(xintercept = num_burn_in, linetype = "dashed") +
      ggplot2::geom_path(
        mapping = ggplot2::aes(x     = .data$iter,
                               y     = .data[[parameter]],
                               color = .data$chains)
      )

    return(res)

  } else {
    stop(
      "There is no parameter named ", parameter, ". \n",
      "  Please try one of following parameters\n",
      "  * `ATE`, `Y1`, `Y0`, `sigma2_out`, `sigma2_out1`, `sigma2_out0` or `dir_alpha`."
    )
  }
}
