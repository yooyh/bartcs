#include <Rcpp.h>
#include "BartTree.h"

using namespace Rcpp;
using namespace std;

void BartTree::updateLatentVariable(NumericVector& latent_variable, const bool is_binary_trt)
{
    if (is_binary_trt)
    {
        NumericVector fitted_values = rowSums(leaf_values_);
        double Ystar;
        for (int i = 0; i < fitted_values.length(); i++)
        {
            Ystar = R::rnorm(fitted_values(i), 1);
            latent_variable(i) = trt_(i) * max(Ystar, 0.0) + (1 - trt_(i)) * min(Ystar, 0.0);
        }
    }
    else
    {
        NumericVector fitted_values = rowSums(leaf_values_);
        for (int i = 0; i < fitted_values.length(); i++)
        {
            latent_variable(i) = R::rnorm(fitted_values(i), sigma2_);
        }
    }
}

void BartTree::updateResidual(const NumericVector& latent_variable, const int t)
{
    NumericVector fitted_values = rowSums(leaf_values_);
    for (int i = 0; i < fitted_values.length(); i++)
    {
        residual_(i) = latent_variable(i) - (fitted_values(i) - leaf_values_(i, t));
    }
}

void BartTree::updateSigma2(
    const Function&      rinvgamma,
    const NumericVector& Y,
    const double         nu,
    const double         lambda
) {
    const int     NUM_OBS       = X_.nrow();
    NumericVector fitted_values = rowSums(leaf_values_);

    const double  SHAPE         = nu / 2        + NUM_OBS / 2;
    const double  SCALE         = nu*lambda / 2 + sum(pow(Y - fitted_values, 2)) / 2;

    NumericVector sigma2_temp   = rinvgamma(1, SHAPE, SCALE);
    sigma2_ = sigma2_temp(0);
}
