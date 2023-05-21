#' @exportS3Method
print.bartcs <- function(x, ...) {
  cat(
    "`bartcs` fit by `", x$model, "_bart()`",
    "\n\n", sep = ""
  )

  mat <- do.call("rbind", x$mcmc_list[, 1:3])
  df <- data.frame(
    mu = apply(mat, 2, mean),
    ci = t(apply(mat, 2, stats::quantile, c(0.025, 0.975)))
  )
  colnames(df) <- c("mean", "2.5%", "97.5%")
  print(df)
}
