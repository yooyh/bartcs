#' @exportS3Method
print.bartcs <- function(x, ...) {
  cat(
    "`bartcs` fit by `", x$model, "()`",
    "\n\n", sep = ""
  )

  estimand <- c("ATE", "Y1", "Y0")
  df <- data.frame(
    mean = vapply(estimand, function(est) mean(x[[est]]), numeric(1)),
    ci   = t(vapply(
      estimand,
      function(est) stats::quantile(x[[est]], probs = c(0.025, 0.975)),
      numeric(2)
    ))
  )
  colnames(df) <- c("mean", "2.5%", "97.5%")
  print(df)
  cat("\n")
}
