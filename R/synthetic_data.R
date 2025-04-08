#' Synthetic dataset for simulation
#'
#' @description
#' Create synthetic dataset for simulation.
#'
#' @param N Number of observations for dataset. The default value is set to 300.
#' @param P Number of potential confounders for dataset. 
#' Need to set X > 7 for data generation. The default value is set to 100.
#' @param seed Seed value for simulation. The default value is set to 42.
#' @param binary_trt Whether the treatment is binary. The default value is set to TRUE.
#'
#' @details
#' `synthetic_data()` generates synthetic dataset for Scenario 1 from 
#' Kim et al. (2023). Among possible confounders, X1 - X5 are true confounders.
#'
#' @references
#' Kim, C., Tec, M., & Zigler, C. M. (2023).
#' Bayesian Nonparametric Adjustment of Confounding, *Biometrics*
#' \doi{10.1111/biom.13833}
#'
#' @return
#' Provide list with the following components
#'
#' \item{Y}{A vector of outcome values.}
#' \item{Trt}{A vector of binary treatment values.}
#' \item{X}{A matrix of potential confounders.}
#'
#' @examples
#' synthetic_data()
#'
#' @export
synthetic_data <- function(N=300, P=100, seed=42, binary_trt=TRUE) {
  set.seed(seed)
  cov <- list()
  for (i in 1:P) {
    cov[[i]] <- stats::rnorm(N, 0, 1)
  }
  X <- do.call(cbind, cov)
  h1 <- ifelse(X[, 1] < 0, 1, -1)
  h2 <- ifelse(X[, 2] < 0, -1, 1)

  if (binary_trt) {
    prob <- stats::pnorm(0.5 + h1 + h2 - 0.5 * abs(X[, 3] - 1) + 1.5 * X[, 4] * X[, 5])
    Trt <- stats::rbinom(N, 1, prob)
    mu1 <- 1 * h1 + 1.5 * h2 - 1 + 2 * abs(X[, 3] + 1) + 2 * X[, 4] + exp(0.5 * X[, 5]) -
      0.5 * 1 * abs(X[, 6]) - 1 * 1 * abs(X[, 7] + 1)
    mu0 <- 1 * h1 + 1.5 * h2 - 0 + 2 * abs(X[, 3] + 1) + 2 * X[, 4] + exp(0.5 * X[, 5]) -
      0.5 * 0 * abs(X[, 6]) - 1 * 0 * abs(X[, 7] + 1)
    Y1 <- stats::rnorm(N, mu1, 0.3)
    Y0 <- stats::rnorm(N, mu0, 0.3)
    Y <- Trt * Y1 + (1 - Trt) * Y0
  } else {
    mu_trt <- 0.5 + h1 + h2 - 0.5 * abs(X[, 3] - 1) + 0.5 * X[, 4] * X[, 5]
    Trt <- stats::rnorm(N, mu_trt, 0.3)
    mu_y <- 1 * h1 + 1 * h2 - Trt + 1 * abs(X[, 3] + 1) + 1 * X[, 4] + exp(0.5 * X[, 5]) -
      0.5 * Trt * abs(X[, 6]) - 0.5 * Trt * abs(X[, 7] + 1)
    Y <- stats::rnorm(N, mu_y, 0.3)
  }
  list(Y = Y, Trt = Trt, X = X)
}
