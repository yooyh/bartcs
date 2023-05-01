# bartcs 1.1.0

* Change output format to `mcmc.list` from `coda` package.

# bartcs 1.0.0

## Changes

* Changed function names.
    * `mbart()` -> `single_bart()`
    * `sbart()` -> `separate_bart()`
* Removed bootstrapping in prediction.
* `trace_plot()` can take both `dir_alpha` and `alpha` as input.
* Set default `num_thin` to 1.
* Fixed post-processing of `var_prob`.
* Fixed scaling of `Y`.

## Improvements

* Improved tree computation by referencing [BART](https://CRAN.R-project.org/package=BART) package. 
* Fix computation process in MCMC algorithm to reflect the changes in the original paper. 

# bartcs 0.1.2

* All OpenMP pragmas now have #ifdef _OPENMP conditions as recommended.

# bartcs 0.1.1

* Deleted commented examples and changed them to executable examples.

# bartcs 0.1.0

* Add package docs
