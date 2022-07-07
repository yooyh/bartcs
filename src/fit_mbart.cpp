#include <Rcpp.h>
#include "BartTree.h"
#include "ProgressBar.h"

using namespace Rcpp;
using namespace R;
using namespace std;

// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
void fit_mbart(
        NumericVector&       Y1,
        NumericVector&       Y0,
        NumericMatrix&       var_count,
        NumericVector&       var_prob,
        double&              sigma2_exp,
        NumericVector&       sigma2_out_hist,
        NumericVector&       dir_alpha_hist,
        const NumericVector& Y,
        const NumericVector& trt,
        const NumericMatrix& X,
        const double         trt_treated,
        const double         trt_control,
        const int            chain_idx,
        const int            num_chain,
        const int            total_iter,
        const int            num_burn_in,
        const int            num_thin,
        const int            num_post_sample,
        const int            num_tree,
        const NumericVector& step_prob,
        const double         alpha,
        const double         beta,
        const double         nu,
        const double         lambda_exp,
        const double         lambda_out,
        const int            boot_size,
        const bool           is_binary_trt,
        const bool           parallel,
        const bool           verbose
) {
    
    const int NUM_OBS = X.nrow();
    const int NUM_VAR = X.ncol();
    const int TRT_IDX = X.ncol();
    
    // unique value of potential confounders
    vector<NumericVector> Xcut (NUM_VAR + 1);
    for (int j = 0; j < NUM_VAR + 1; j++)
    {
        NumericVector temp;
        if (j == TRT_IDX)
        {
            // keep unique values of trt in Xcut
            temp = unique(trt);
            temp.sort();
            Xcut[j] = clone(temp);
        } 
        else 
        {
            temp = unique(X(_, j));
            temp.sort();
            Xcut[j] = clone(temp);
        }
    }
    
    NumericVector latent_variable;
    {
        const double MEAN = R::qnorm(mean(trt), 0, 1, true, false);
        const double SD   = 1.0;
        latent_variable   = Rcpp::rnorm(NUM_OBS, MEAN, SD);
    }
    
    double sigma2_out = sigma2_out_hist(0);
    double dir_alpha  = dir_alpha_hist(0);
    
    // sigma_mu based on min/max of Y, Y (Tr=1) and Y (Tr=0)    
    const double sigma_mu_exp = max(
        pow(min(latent_variable) / (-2*sqrt(num_tree)), 2),
        pow(max(latent_variable) / ( 2*sqrt(num_tree)), 2)
    );
    const double sigma_mu_out = max(
        pow(min(Y) / (-2*sqrt(num_tree)), 2),
        pow(max(Y) / ( 2*sqrt(num_tree)), 2)
    );
    
    // Initial values of R
    NumericVector residual_exp = clone(latent_variable);
    NumericVector residual_out = clone(Y);
    
    // Initial values for the selection probabilities
    NumericVector post_dir_alpha = rep(1.0, NUM_VAR + 1);
    
    // Obtaining namespace of MCMCpack package
    // Then pick up rinvgamma() and rdirichlet() function from MCMCpack package
    Environment MCMCpack   = Environment::namespace_env("MCMCpack");
    Function    rinvgamma  = MCMCpack["rinvgamma"];
    Function    rdirichlet = MCMCpack["rdirichlet"];
    
    // initialize model
    BartTree exposure = BartTree(
        residual_exp, var_prob, sigma2_exp,   // mutable variables
        1, trt, X, Xcut, step_prob, num_tree, // const variables 1
        alpha, beta, sigma_mu_exp, parallel   // const variables 2
    );
    BartTree outcome  = BartTree(
        residual_out, var_prob, sigma2_out,   // mutable variables
        2, trt, X, Xcut, step_prob, num_tree, // const variables 1
        alpha, beta, sigma_mu_out, parallel   // const variables 2
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
         outcome.step(Y,               is_binary_trt);
        
        // update sigma
        if (!is_binary_trt)
            exposure.updateSigma2(rinvgamma, Y, nu, lambda_exp);
        outcome.updateSigma2(rinvgamma,  Y, nu, lambda_out);
        sigma2_out_hist(iter) = outcome.getSigma2();

        // count included
        NumericVector var_count_exp = exposure.countSelectedVariables();
        NumericVector var_count_out =  outcome.countSelectedVariables();
        
        // MH to update dir_alpha
        exposure.updateDirAlpha(dir_alpha);
        dir_alpha_hist(iter) = dir_alpha;
        
        // then update post_dir_alpha
        post_dir_alpha = rep(dir_alpha / (NUM_VAR + 1), NUM_VAR + 1);
        
        // MH algorithm to update inclusion probabilities
        exposure.updateVarProb(
            rdirichlet, post_dir_alpha, var_count_exp, var_count_out
        );
        
        // sample E[Y(1) - Y(0)]
        if (iter > num_burn_in) 
        {
            if (thin_count == num_thin)
            {
                // record selected variables
                var_count(res_idx, _) = var_count_out;
                
                // predict effect and potential outcomes
                Y1(res_idx) = outcome.predict(boot_idx, trt_treated);
                Y0(res_idx) = outcome.predict(boot_idx, trt_control);
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
