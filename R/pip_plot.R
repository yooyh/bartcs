pip_plot <- function(x, ...) {
  label    <- x$label
  var_prob <- x$var_prob

  df <- data.frame(label = label, var_prob = var_prob)

  if (x$model == "separate")
    title <- "PIP (Separate Model)"
  else if (x$model == "single")
    title <- "PIP (Single Model)"

  ggcharts::bar_chart(df, label, var_prob, ...) +
    ggplot2::labs(
      y     = "Probability",
      x     = "Variables",
      title = title
    )
}
