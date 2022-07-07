# install package from source
# install.packages("devtools")
devtools::install()

# load package and data
library(bartcs)
data("lalonde", package = "Matching")

?sbart
?plot.bartcs
count_omp_thread()

Y   <- lalonde$re78
trt <- lalonde$treat
X   <- subset(lalonde, select = -c(re78, treat))

# fit model ----
# right call
# call with data.frame
mbart(Y = Y, trt = trt, X = X,
      num_burn_in = 0, num_post_sample = 10)

# call with tibble
mbart(Y = Y, trt = trt, X = tibble::as_tibble(X),
      num_burn_in = 0, num_post_sample = 10)

# call with matrix
mbart(Y = Y, trt = trt, X = as.matrix(X),
      num_burn_in = 0, num_post_sample = 10)

# errors
# error types: wrong types
mbart(Y = "Y", trt = trt, X = X)
mbart(Y = Y, trt = "trt", X = X)
mbart(Y = Y, trt = trt, X = "X")
mbart(Y = Y, trt = trt, X = data.frame(foo = c("a", "b"), bar = 1:2))
mbart(Y = Y, trt = trt, X = X, verbose = -1)
mbart(Y = Y, trt = trt, X = X, num_tree = "many")

# error types: wrong value
mbart(Y = Y, trt = trt, X = X, num_tree = -100)             # wrong value
mbart(Y = Y, trt = rep(1, 100), X = X)                      # wrong dimension
mbart(Y = Y, trt = rep(1, length(Y)), X = X)                # wrong trt
mbart(Y = Y, trt = trt, X = X, step_prob = c(1, 1, 1, 1))   # wrong step_prob


# print and summary ----
# normal case
fit <- mbart(Y = Y, trt = trt, X = X)
fit
summary(fit)
gelman_rubin(fit)

# when number of MCMC chain = 1 -> message and NA
fit2 <- mbart(Y = Y, trt = trt, X = X, num_chain = 1)
gelman_rubin(fit2)
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
