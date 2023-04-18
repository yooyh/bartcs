# install package from source
# install.packages("devtools")
devtools::install()

# load package and data
library(bartcs)
# install.packages("Matching")
data("lalonde", package = "Matching")

?plot.bartcs
count_omp_thread()

Y   <- lalonde$re78
trt <- lalonde$treat
X   <- subset(lalonde, select = -c(re78, treat))

# fit model ----
# right call
# call with data.frame
single_bart(Y = Y, trt = trt, X = X,
      num_burn_in = 0, num_post_sample = 10)

# call with tibble
single_bart(Y = Y, trt = trt, X = tibble::as_tibble(X),
      num_burn_in = 0, num_post_sample = 10)

# call with matrix
single_bart(Y = Y, trt = trt, X = as.matrix(X),
      num_burn_in = 0, num_post_sample = 10)

# errors
# error types: wrong types
single_bart(Y = "Y", trt = trt, X = X)
single_bart(Y = Y, trt = "trt", X = X)
single_bart(Y = Y, trt = trt, X = "X")
single_bart(Y = Y, trt = trt, X = data.frame(foo = c("a", "b"), bar = 1:2))
single_bart(Y = Y, trt = trt, X = X, verbose = -1)
single_bart(Y = Y, trt = trt, X = X, num_tree = "many")

# error types: wrong value
single_bart(Y = Y, trt = trt, X = X, num_tree = -100)             # wrong value
single_bart(Y = Y, trt = rep(1, 100), X = X)                      # wrong dimension
single_bart(Y = Y, trt = rep(1, length(Y)), X = X)                # wrong trt
single_bart(Y = Y, trt = trt, X = X, step_prob = c(1, 1, 1, 1))   # wrong step_prob


# print and summary ----
# normal case
fit <- single_bart(Y = Y, trt = trt, X = X)
fit
summary(fit)

# when number of MCMC chain = 1 -> message and NA
fit2 <- single_bart(Y = Y, trt = trt, X = X, num_chain = 1)
summary(fit2)

# data visualization ----
plot(fit) # error -> choose method

# pip plot
plot(fit, method = "pip")                                   # works fine
plot(fit, method = "pip", "bar")                            # works with message
plot(fit, method = "pip", top_n = 10)                       # works fine
plot(fit, method = "pip", "bar", top_n = 10)                # works with message

# trace plot
plot(fit, method = "trace")
plot(fit, method = "trace", "Y1")
plot(fit, method = "trace", "dir_alpha")
plot(fit, method = "trace", "foo")
plot(fit, method = "trace", "sigma2_out1")                  # auto-switch to sigma2_out

# plot when number of MCMC chain = 1 -> works fine
plot(fit2, "pip")
plot(fit2, "trace")
plot(fit2, "trace", "Y1")
plot(fit2, "trace", "dir_alpha")
