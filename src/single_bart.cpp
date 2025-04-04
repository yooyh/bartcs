#include "Models.h"
#include "ProgressBar.h"

using namespace std;

// [[Rcpp::export]]
Rcpp::List csingle_bart(
    const Rcpp::NumericVector& Y_src,
    const Rcpp::NumericMatrix& X_src,
    const Rcpp::NumericVector& TRT_src,
    const double trt_treated, 
    const double trt_control,
    const double Y_min, 
    const double Y_max,
    const Rcpp::NumericVector& step_prob,
    const int num_chain,
    const int num_chain_iter,
    const int num_burn_in,
    const int num_thin,
    const int num_post_sample,
    const int num_tree,
    const double alpha,
    const double beta,
    const double nu,
    const double lambda_exp,
    const double lambda_out,
    const double initial_dir_alpha,
    const double initial_sigma2_exp,
    const double initial_sigma2_out,
    const bool binary_trt,
    const bool parallel,
    const bool verbose
) {
    // preprocessing
    const int N = X_src.nrow(), P = X_src.ncol();

    vector<double> Z (N), Y (N); // latent variable
    copy(Y_src.begin(), Y_src.end(), Y.begin());
    vector<vector<double>> X (N, vector<double> (P + 1, 0.0)), Xcut;
    get_data(X_src, TRT_src, X, Xcut);

    // Obtaining namespace of MCMCpack package
    // Then pick up rinvgamma() and rdirichlet() function from MCMCpack package
    Environment MCMCpack   = Environment::namespace_env("MCMCpack");
    Function    rinvgamma  = MCMCpack["rinvgamma"];
    Function    rdirichlet = MCMCpack["rdirichlet"];


    // initialize model
    Rcpp::NumericVector var_prob_exp (P,     0.0);
    Rcpp::NumericVector var_prob_out (P + 1, 0.0);
    auto exposure = SingleModel(X, Xcut, N, P,     num_tree, step_prob, alpha, beta, var_prob_exp, parallel);
    auto outcome  = SingleModel(X, Xcut, N, P + 1, num_tree, step_prob, alpha, beta, var_prob_out, parallel);

    Rcpp::List chains (num_chain);
    for (int chain = 0; chain < num_chain; chain++)
    {
        // init ProgressBar
        ProgressBar pb = ProgressBar(chain, num_chain, num_chain_iter, 40);

        // init
        double dir_alpha = initial_dir_alpha;
        var_prob_out = rdirichlet(1, rep(dir_alpha, P + 1));
        normalize(var_prob_exp, var_prob_out);

        int thin_count = 0, res_id = 0;
        init_Z(Z, TRT_src, binary_trt);
        exposure.init(Z, initial_sigma2_exp);
        outcome.init(Y,  initial_sigma2_out);

        // create placeholder
        Rcpp::NumericVector POY1 (num_post_sample, 0.0), POY0 (num_post_sample, 0.0);
        Rcpp::NumericMatrix y1_sample (N, num_post_sample), y0_sample (N, num_post_sample);
        Rcpp::NumericVector dir_alpha_hist (num_post_sample, 0.0);
        Rcpp::NumericVector sigma2_hist    (num_post_sample, 0.0);
        Rcpp::IntegerMatrix var_count      (num_post_sample, P + 1);
        Rcpp::NumericVector post_dir_alpha (P + 1, dir_alpha / (P + 1));

        for (int iter = 1; iter <= num_chain_iter; iter++)
        {
            // check for break from R
            Rcpp::checkUserInterrupt();
            if (verbose) pb.print(iter);

            // update Z
            exposure.update_Z(Z, TRT_src, binary_trt);

            // update tree
            exposure.draw(Z);
            outcome.draw(Y);

            // update sigma
            if (!binary_trt) exposure.draw_sigma2(rinvgamma, Y, nu, lambda_exp);
            outcome.draw_sigma2(rinvgamma, Y, nu, lambda_out);

            // update dir_alpha and var_prob
            mh_dir_alpha(var_prob_out, dir_alpha, post_dir_alpha);
            mh_var_prob(rdirichlet, post_dir_alpha, exposure.var_count(), outcome.var_count(), var_prob_out);
            normalize(var_prob_exp, var_prob_out);

            if (iter <= num_burn_in) continue;
            thin_count++;
            if (thin_count == num_thin)
            {
                // record sigma and dir_alpha
                dir_alpha_hist[res_id] = dir_alpha;
                sigma2_hist[res_id]    = outcome.sigma2();

                // record selected variable
                auto cnt = outcome.var_count();
                for (int j = 0; j < P + 1; j++)
                {
                    // var_count(res_id, j) = cnt[j];
                    int c = (j == P) ? 0 : j + 1;
                    var_count[res_id + num_post_sample * c] = cnt[j];
                }
                
                // predict effect and potential outcomes
                outcome.predict(POY1, y1_sample, res_id, trt_treated);
                outcome.predict(POY0, y0_sample, res_id, trt_control);
                if (res_id == num_post_sample)
                    break;
                res_id++;
                thin_count = 0;
            }
        }
        // post processing
        POY1 = (POY1 + 0.5) * (Y_max - Y_min) + Y_min;
        POY0 = (POY0 + 0.5) * (Y_max - Y_min) + Y_min;
        Rcpp::NumericVector ATE = POY1 - POY0;
        Rcpp::NumericVector PIP (P + 1, 0.0);
        for (int j = 0; j < P + 1; j++)
        {
            for (int i = 0; i < num_post_sample; i++)
                PIP[j] += var_count[i + j * num_post_sample] > 0 ? 1.0 : 0.0;
            PIP[j] /= num_post_sample;
        }

        // gather result
        chains[chain] = Rcpp::List::create(
            Named("ATE")        = ATE,
            Named("Y1")         = POY1,
            Named("Y0")         = POY0,
            Named("Y1_sample")  = y1_sample,
            Named("Y0_sample")  = y0_sample,
            Named("var_count")  = var_count,
            Named("var_prob")   = PIP,
            Named("sigma2_out") = sigma2_hist,
            Named("dir_alpha")  = dir_alpha_hist
        );
    }
    return chains;
}