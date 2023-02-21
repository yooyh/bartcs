#include "Models.h"
#include "ProgressBar.h"

using namespace std;

// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
Rcpp::List cseparate_bart(
    const Rcpp::NumericMatrix& X_src,
    const Rcpp::NumericVector& Y1_src,
    const Rcpp::NumericMatrix& X1_src,
    const Rcpp::NumericVector& Y0_src,
    const Rcpp::NumericMatrix& X0_src,
    const Rcpp::NumericVector& TRT_src,
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
    const double lambda_out1,
    const double lambda_out0,
    const double initial_dir_alpha,
    const double initial_sigma2_out1,
    const double initial_sigma2_out0,
    const bool parallel,
    const bool verbose
) {
    // preprocessing
    const int N  = X_src.nrow(),  P  = X_src.ncol();
    const int N1 = X1_src.nrow(), N0 = X0_src.nrow();

    vector<double> Z (N), Y1 (N1), Y0 (N0); // Z : latent variable
    copy(Y1_src.begin(), Y1_src.end(), Y1.begin());
    copy(Y0_src.begin(), Y0_src.end(), Y0.begin());
    vector<vector<double>> X  (N,  vector<double> (P, 0.0));
    vector<vector<double>> X1 (N1, vector<double> (P, 0.0));
    vector<vector<double>> X0 (N0, vector<double> (P, 0.0));
    vector<vector<double>> Xcut, Xcut1, Xcut0;
    get_data(X_src,  X,  Xcut);
    get_data(X1_src, X1, Xcut1);
    get_data(X0_src, X0, Xcut0);

    // Obtaining namespace of MCMCpack package
    // Then pick up rinvgamma() and rdirichlet() function from MCMCpack package
    Environment MCMCpack   = Environment::namespace_env("MCMCpack");
    Function    rinvgamma  = MCMCpack["rinvgamma"];
    Function    rdirichlet = MCMCpack["rdirichlet"];


    // initialize model
    Rcpp::NumericVector var_prob (P);
    auto exposure = SeparateModel(X,  Xcut,  num_tree, step_prob, alpha, beta, var_prob, parallel);
    auto outcome1 = SeparateModel(X1, Xcut1, num_tree, step_prob, alpha, beta, var_prob, parallel);
    auto outcome0 = SeparateModel(X0, Xcut0, num_tree, step_prob, alpha, beta, var_prob, parallel);

    Rcpp::List chains (num_chain);
    for (int chain = 0; chain < num_chain; chain++)
    {
        // init ProgressBar
        ProgressBar pb = ProgressBar(chain, num_chain, num_chain_iter, 40);

        // init
        double dir_alpha = initial_dir_alpha;
        var_prob = rdirichlet(1, rep(dir_alpha, P));
        int thin_count = 0, res_id = 0;
        init_Z(Z, TRT_src, true);
        exposure.init(Z,  1); 
        outcome1.init(Y1, initial_sigma2_out1); 
        outcome0.init(Y0, initial_sigma2_out0);

        // create placeholder
        Rcpp::NumericVector POY1 (num_post_sample, 0.0), POY0 (num_post_sample, 0.0);
        Rcpp::NumericVector sigma2_out1_hist (num_chain_iter + 1, 0.0);
        Rcpp::NumericVector sigma2_out0_hist (num_chain_iter + 1, 0.0);
        Rcpp::NumericVector dir_alpha_hist   (num_chain_iter + 1, 0.0);
        Rcpp::IntegerMatrix var_count        (num_post_sample, P);
        Rcpp::NumericVector post_dir_alpha   (P, 1.0);
        sigma2_out1_hist[0] = initial_sigma2_out1;
        sigma2_out0_hist[0] = initial_sigma2_out0;
        dir_alpha_hist[0]   = initial_dir_alpha;

        for (int iter = 1; iter <= num_chain_iter; iter++)
        {
            // check for break from R
            Rcpp::checkUserInterrupt();
            if (verbose) pb.print(iter);

            // update Z
            exposure.update_Z(Z, TRT_src, true);

            // update tree
            exposure.draw(Z);
            outcome1.draw(Y1);
            outcome0.draw(Y0);

            // update sigma
            outcome1.draw_sigma2(rinvgamma, Y1, nu, lambda_out1);
            outcome0.draw_sigma2(rinvgamma, Y0, nu, lambda_out0);
            sigma2_out1_hist[iter] = outcome1.sigma2();
            sigma2_out0_hist[iter] = outcome0.sigma2();

            // update dir_alpha and var_prob
            mh_dir_alpha(iter, num_chain_iter, var_prob,
                         exposure.var_count(), outcome1.var_count(), outcome0.var_count(), 
                         dir_alpha, post_dir_alpha);
            dir_alpha_hist[iter] = dir_alpha;
            var_prob = rdirichlet(1, post_dir_alpha);

            if (iter <= num_burn_in) continue;
            thin_count++;
            if (thin_count == num_thin)
            {
                // record selected variable
                auto cnt1 = outcome1.var_count(), cnt0 = outcome0.var_count();
                for (int j = 0; j < P; j++)
                    var_count[res_id + num_post_sample * j] = cnt1[j] + cnt0[j];
                
                // predict effect and potential outcomes
                outcome1.predict(POY1, res_id, X);
                outcome0.predict(POY0, res_id, X);
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
        Rcpp::NumericVector PIP (P, 0.0);
        for (int j = 0; j < P; j++)
        {
            for (int i = 0; i < num_post_sample; i++)
                PIP[j] += (double) var_count[i + j * num_post_sample] > 0 ? 1 : 0;
            PIP[j] /= num_post_sample;
        }

        // gather result
        chains[chain] = Rcpp::List::create(
            Named("ATE")         = ATE,
            Named("Y1")          = POY1,
            Named("Y0")          = POY0,
            Named("var_count")   = var_count,
            Named("var_prob")    = PIP,
            Named("sigma2_out1") = sigma2_out1_hist,
            Named("sigma2_out0") = sigma2_out0_hist,
            Named("alpha")       = dir_alpha_hist
        );
    }
    return chains;
}