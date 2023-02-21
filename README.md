
# bartcs

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/bartcs)](https://CRAN.R-project.org/package=bartcs)
[![R-CMD-check](https://github.com/yooyh/bartcs/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yooyh/bartcs/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Fit Bayesian Additive Regression Trees (BART) models to select confounders and estimate treatment effect.

## Installation

You can install from CRAN with:

```r
install.packages("bartcs")
```

Or you can install the development version from
[GitHub](https://github.com/) with:
``` r
remotes::install_github("yooyh/bartcs")
```

### For **Mac** Users Only

bartcs supports multi-threading by OpenMP. If you

- Have OpenMP installed
- And want to use OpenMP multi-threading

then install package from source with:

```r
# build from source
install.packages("bartcs", type = "source")

# count OpemMP thread
bartcs::count_omp_thread() # this should be greater than 1
```