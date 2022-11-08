#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

void updateDirAlpha(double& dir_alpha, const NumericVector& var_prob)
{
    const int NUM_VAR        = var_prob.length(); // is X.ncol() + 1 in mar model

    double    prop_dir_alpha = max(R::rnorm(dir_alpha, 0.1), pow(0.1, 10));
    NumericVector log_var_prob(NUM_VAR);
    {
        const double LB     = pow(0.1, 300);
        const double LOG_LB = -300*log(10);
        for (int i = 0; i < NUM_VAR; i++)
        {
            if(var_prob(i) > LB) 
            {
                log_var_prob(i) = log(var_prob(i));
            }
            else
            {
                log_var_prob(i) = LOG_LB;
            }
        }
    }

    double prop_likelihood =
        sum(log_var_prob * (rep(prop_dir_alpha/NUM_VAR, NUM_VAR) - 1))
        + lgamma(sum(rep(prop_dir_alpha/NUM_VAR, NUM_VAR)))
        - sum(lgamma(rep(prop_dir_alpha/NUM_VAR, NUM_VAR)))
        ;

    double likelihood =
        sum(log_var_prob * (rep(dir_alpha/NUM_VAR, NUM_VAR) - 1))
        + lgamma(sum(rep(dir_alpha/NUM_VAR, NUM_VAR)))
        - sum(lgamma(rep(dir_alpha/NUM_VAR, NUM_VAR)))
        ;

    double ratio =
        prop_likelihood
        - 0.5 * log(prop_dir_alpha)
        - 1.5 * log(NUM_VAR + prop_dir_alpha)
        - likelihood
        + 0.5 * log(dir_alpha)
        + 1.5 * log(NUM_VAR + dir_alpha)
        ;

    if (ratio > log(R::runif(0,1)))
    {
        dir_alpha = prop_dir_alpha;
    }
}

// update var_prob for marginal model
void updateVarProb(
    NumericVector&       var_prob,
    NumericVector&       var_count_exp,
    const NumericVector& var_count_out,
    const Function&      rdirichlet,
    const NumericVector& post_dir_alpha
) {
    // this p is not X.ncol() in marginal model !!
    const int NUM_VAR = var_prob.length();
    const int TRT_IDX = NUM_VAR - 1;

    // Assign var_count_out(TRT_IDX) to var_count_exp(TRT_IDX)
    // 
    // During probability update, we assume TRT counts of exposure model equal
    // TRT counts of outcome model since exposure model does not include TRT 
    // in the model.
    // i.e. var_count_exp(TRT_IDX) = var_count_out(TRT_IDX)
    var_count_exp.push_back(var_count_out(TRT_IDX));

    NumericVector alpha         = post_dir_alpha + var_count_exp + var_count_out;
    NumericVector prop_var_prob = rdirichlet(1, alpha);

    // change probability into log form
    NumericVector log_var_prob(NUM_VAR), log_prop_var_prob(NUM_VAR);
    {
        const double LB     = pow(0.1, 300);
        const double LOG_LB = -300 * log(10);
        for (int i = 0; i < NUM_VAR; i++)
        {
            // logarize var_prob
            if (var_prob(i) > LB)
            {
                log_var_prob(i) = log(var_prob(i));
            }
            else
            {
                log_var_prob(i) = LOG_LB;
            }
            // logarize prop_var_prob
            if (prop_var_prob(i) > LB)
            {
                log_prop_var_prob(i) = log(prop_var_prob(i));
            }
            else
            {
                log_prop_var_prob(i) = LOG_LB;
            }
        }
    }

    double prop_dir_lik = sum(var_count_exp) * (- log(1 - prop_var_prob(TRT_IDX)));
    double dir_lik      = sum(var_count_exp) * (- log(1 - var_prob(TRT_IDX)));
    for (int i = 0; i < NUM_VAR; i++)
    {
        double temp   = var_count_exp(i) + var_count_out(i) + post_dir_alpha(i);
        prop_dir_lik += temp * log_prop_var_prob(i);
        dir_lik      += temp * log_var_prob(i);
    }

    double ratio =
        prop_dir_lik
        + sum((post_dir_alpha + var_count_exp + var_count_out - 1.0) * log_var_prob)
        - dir_lik
        - sum((post_dir_alpha + var_count_exp + var_count_out - 1.0) * log_prop_var_prob)
        ;

    if (ratio > log(R::runif(0,1)))
    {
        var_prob = clone(prop_var_prob);
    }
}
