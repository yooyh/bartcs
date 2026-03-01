
# bartcs

<!-- badges: start -->
[![CRAN](https://www.r-pkg.org/badges/version/bartcs)](https://cran.r-project.org/package=bartcs)
[![R-CMD-check](https://github.com/yooyh/bartcs/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yooyh/bartcs/actions/workflows/R-CMD-check.yaml)
[![R Journal](https://img.shields.io/badge/The%20R%20Journal-2025-blue)](https://journal.r-project.org/articles/RJ-2025-019/)
<!-- badges: end -->

Fit Bayesian Additive Regression Trees (BART) models to select confounders and estimate treatment effect.

## Publication

This package is described in The R Journal:

Yoo, Y. and Kim, C. (2025).
*bartcs: An R Package for Bayesian Nonparametric Adjustment for Confounding*.
The R Journal, 17(2).
https://journal.r-project.org/articles/RJ-2025-019/

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
