#include <Rcpp.h>
#include "BartTree.h"
#include "ProgressBar.h"

using namespace Rcpp;
using namespace R;
using namespace std;

// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
void fit_sbart(
    NumericVector&       Y1,
    NumericVector&       Y0,
    NumericMatrix&       var_count,
    NumericVector&       var_prob,
    NumericVector&       sigma2_out1_hist,
    NumericVector&       sigma2_out0_hist,
    NumericVector&       dir_alpha_hist,
    const NumericVector& Y_treated,
    const NumericVector& Y_control,
    const NumericVector& trt,
    const NumericMatrix& X,
    const NumericMatrix& X_treated,
    const NumericMatrix& X_control,
    const int            chain_idx,
    const int            num_chain,
    const int            total_iter,
    const int            num_burn_in,
    const int            num_thin,
    const int            num_post_sample,
    const int            num_tree,
    const NumericVector  step_prob,
    const double         alpha,
    const double         beta,
    const double         nu,
    const double         lambda_out1,
    const double         lambda_out0,
    const int            boot_size,
    const bool           is_binary_trt,
    const bool           parallel,
    const bool           verbose
) {

    const int NUM_OBS = X.nrow();
    const int NUM_VAR = X.ncol();

    // unique value of potential confounders
    vector<NumericVector> Xcut(NUM_VAR), Xcut1(NUM_VAR), Xcut0(NUM_VAR);
    for (int j = 0; j < NUM_VAR; j++)
    {
        NumericVector temp;
        temp = unique(X(_, j));
        temp.sort();
        Xcut[j] = clone(temp);

        temp = unique(X_treated(_, j));
        temp.sort();
        Xcut1[j] = clone(temp);

        temp = unique(X_control(_, j));
        temp.sort();
        Xcut0[j] = clone(temp);
    }

    NumericVector latent_variable;
    {
        const double MEAN = R::qnorm(mean(trt), 0, 1, true, false);
        const double SD   = 1.0;
        latent_variable = Rcpp::rnorm(NUM_OBS, MEAN, SD);
    }

    // initialize parameters
    double sigma2_exp  = 1; // probit model of binary treatment
    double sigma2_out1 = sigma2_out1_hist(0);
    double sigma2_out0 = sigma2_out0_hist(0);
    double dir_alpha   = dir_alpha_hist(0);

    // sigma_mu based on min/max of Y, Y (trt=1) and Y (trt=0)
    const double sigma_mu_exp  = max(
        pow(min(latent_variable) / (-2*sqrt(num_tree)), 2),
        pow(max(latent_variable) / ( 2*sqrt(num_tree)), 2)
    );
    const double sigma_mu_out1 = max(
        pow(min(Y_treated)       / (-2*sqrt(num_tree)), 2),
        pow(max(Y_treated)       / ( 2*sqrt(num_tree)), 2)
    );
    const double sigma_mu_out0 = max(
        pow(min(Y_control)       / (-2*sqrt(num_tree)), 2),
        pow(max(Y_control)       / ( 2*sqrt(num_tree)), 2)
    );

    // Initial values of R
    NumericVector residual  = clone(latent_variable);
    NumericVector residual1 = clone(Y_treated);
    NumericVector residual0 = clone(Y_control);

    // Initial values for the selection probabilities
    NumericVector post_dir_alpha = rep(1.0, NUM_VAR);

    // Obtaining namespace of MCMCpack package
    // Then pick up rinvgamma() and rdirichlet() function from MCMCpack package
    Environment MCMCpack   = Environment::namespace_env("MCMCpack");
    Function    rinvgamma  = MCMCpack["rinvgamma"];
    Function    rdirichlet = MCMCpack["rdirichlet"];

    // initialize model
    BartTree exposure = BartTree(
        residual, var_prob, sigma2_exp,                // mutable variables
        0, trt, X, Xcut, step_prob, num_tree,          // const variables 1
        alpha, beta, sigma_mu_exp,  parallel           // const variables 2
    );
    BartTree outcome1 = BartTree(
        residual1, var_prob, sigma2_out1,              // mutable variables
        0, trt, X_treated, Xcut1, step_prob, num_tree, // const variables 1
        alpha, beta, sigma_mu_out1, parallel           // const variables 2
    );
    BartTree outcome0 = BartTree(
        residual0, var_prob, sigma2_out0,              // mutable variables
        0, trt, X_control, Xcut0, step_prob, num_tree, // const variables 1
        alpha, beta, sigma_mu_out0, parallel           // const variables 2
    );

    IntegerVector boot_idx = sample(NUM_OBS, boot_size, true) - 1;

    int thin_count = 0, res_idx = 0;
    
    // init ProgressBar
    ProgressBar progress_bar = ProgressBar(chain_idx, num_chain, total_iter, 40);

    // run MCMC
    for (int iter = 1; iter < total_iter + 1; iter++)
    {
        // check for break from R
        Rcpp::checkUserInterrupt();
        if (verbose)
            progress_bar.print(iter);

        // update latent_variable
        exposure.updateLatentVariable(latent_variable, is_binary_trt);

        // update tree
        exposure.step(latent_variable, is_binary_trt);
        outcome1.step(Y_treated,       is_binary_trt);
        outcome0.step(Y_control,       is_binary_trt);

        // update sigma
        outcome1.updateSigma2(rinvgamma, Y_treated, nu, lambda_out1);
        outcome0.updateSigma2(rinvgamma, Y_control, nu, lambda_out0);
        sigma2_out1_hist(iter) = outcome1.getSigma2();
        sigma2_out0_hist(iter) = outcome0.getSigma2();

        // count included
        NumericVector var_count_exp  = exposure.countSelectedVariables();
        NumericVector var_count_out1 = outcome1.countSelectedVariables();
        NumericVector var_count_out0 = outcome0.countSelectedVariables();

        // MH for dir_alpha
        if (iter < total_iter / 10)
        {
            // warm up
            post_dir_alpha = rep(1.0, NUM_VAR) + var_count_exp + var_count_out1 + var_count_out0;
        }
        else
        {
            // use MH algorithm to update dir_alpha
            exposure.updateDirAlpha(dir_alpha);

            // then update post_dir_alpha
            post_dir_alpha = rep(dir_alpha / NUM_VAR, NUM_VAR) + var_count_exp + var_count_out1 + var_count_out0;
        }
        dir_alpha_hist(iter) = dir_alpha;

        // update confounder_prob
        var_prob = rdirichlet(1, post_dir_alpha);

        // sample E[Y(1) - Y(0)]
        if (iter > num_burn_in)
        {
            if (thin_count == num_thin)
            {
                // record selected variables
                var_count(res_idx, _) = var_count_out1 + var_count_out0;

                // predict effect and potential outcomes
                Y1(res_idx)  = outcome1.predict(boot_idx, X);
                Y0(res_idx)  = outcome0.predict(boot_idx, X);
                if (res_idx == num_post_sample)
                    break;
                res_idx++;
                thin_count = 0;
            }
            else
            {
                thin_count++;
            }
        }
    // end of an MCMC iteration
    } 
}
