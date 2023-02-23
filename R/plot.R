#' Draw plot for `bartcs` object
#'
#' @description
#' Two options are available:
#'   posterior inclusion probability (PIP) plot and trace plot.
#'
#' @param x A `bartcs` object.
#' @param method "`pip`" for posterior inclusion probability plot
#'   or "`trace`" for trace plot.
#' @param parameter Parameter for traceplot.
#' @param ... Additional arguments for PIP plot.
#'   Check `?ggcharts::bar_chart` for possible arguments.
#'
#' @details
#' ## PIP plot
#' When a posterior sample is sampled during training,
#' `separate_bart()` or `single_bart()` also counts
#' which variables are included in the model and
#' compute PIP for each variable.
#' For `bartcs` object `x`,
#' this is stored in `x$var_count` and `x$var_prob` respectively.
#' `plot(method = "pip")` uses this information and
#' draws plot using `ggcharts::bar_chart()`.
#'
#' ## Traceplot
#' Parameters are recorded for each MCMC iterations.
#' Parameters include "`ATE`", "`Y1`", "`Y0`", "`dir_alpha`",
#' and either "`sigma2_out`" from `single_bart()`
#' or "`sigma2_out1`" and "`sigma2_out0`" from \cr 
#' `separate_bart()`.
#' Vertical line indicates burn-in.
#'
#' @return
#' A `ggplot` object of either PIP plot or trace plot.
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
#'
#' # PIP plot
#' plot(x, method = "pip")
#' plot(x, method = "pip", top_n = 10)
#' plot(x, method = "pip", threshold = 0.5)
#' # Check `?ggcharts::bar_chart` for other possible arguments.
#'
#' # trace plot
#' plot(x, method = "trace")
#' plot(x, method = "trace", "Y1")
#' plot(x, method = "trace", "dir_alpha")
#'
#' @exportS3Method
plot.bartcs <- function(x, method = NULL, parameter = NULL, ...) {
  if (is.null(method))
    stop(
      "You must choose method.\n",
      "  * For PIP plot, set method = \"pip\".\n",
      "  * For trace plot, set method = \"trace\".\n",
      "\nTry `?plot.bartcs` for more detail."
    )

  if (method == "pip") {
    if (!is.null(parameter))
      message("Argument `parameter` not used with `method = `pip`.")

    pip_plot(x, ...)

  } else if (method == "trace") {
    if (is.null(parameter))
      parameter <- "ATE"

    trace_plot(x, parameter)
  }
}