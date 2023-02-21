#include "util.h"

// util functions for both model
void get_Xcut(const vector<vector<double>>& X, vector<vector<double>>& Xcut)
{
    const int N = X.size(), P = X[0].size();
    Xcut.resize(P);
    set<double> unique_value;
    for (int j = 0; j < P; j++)
    {
        unique_value.clear();
        for (int i = 0; i < N; i++)
            unique_value.insert(X[i][j]);
        Xcut[j].resize(unique_value.size());
        copy(unique_value.begin(), unique_value.end(), Xcut[j].begin());
    }
}

void get_data(const Rcpp::NumericMatrix& X_src, vector<vector<double>>& X, vector<vector<double>>& Xcut)
{
    const int N = X_src.nrow(), P = X_src.ncol();
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < P; j++)
            X[i][j] = X_src[i + j * N];
    }
    get_Xcut(X, Xcut);
}
void get_data(
    const Rcpp::NumericMatrix& X_src, const Rcpp::NumericVector& TRT_src,
    vector<vector<double>>& X, vector<vector<double>>& Xcut
) {
    const int N = X_src.nrow(), P = X_src.ncol();
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < P; j++)
            X[i][j] = X_src[i + j * N];
        X[i][P] = TRT_src[i];
    }
    get_Xcut(X, Xcut);
}
void init_Z(vector<double>& Z, const Rcpp::NumericVector& TRT_src, const bool binary_trt)
{
    int N = TRT_src.size();
    const double MEAN = binary_trt ? R::qnorm(Rcpp::mean(TRT_src), 0, 1, true, false) : 0.0;
    const double SD   = 1.0;
    auto         tmp  = Rcpp::rnorm(N, MEAN, SD);
    copy(tmp.begin(), tmp.end(), Z.begin());
}
void draw_dir_alpha(const Rcpp::NumericVector& var_prob, double& dir_alpha) 
{
    int P = var_prob.size();
    double prop_dir_alpha = max(R::rnorm(dir_alpha, 0.1), pow(0.1, 10));
    vector<double> log_var_prob (P);
    {
        const double LB     = pow(0.1, 300);
        const double LOG_LB = -300 * log(10);
        for (int j = 0; j < P; j++)
        {
            if(var_prob[j] > LB)
                log_var_prob[j] = log(var_prob[j]);
            else
                log_var_prob[j] = LOG_LB;
        }
    }

    Rcpp::NumericVector p_dir_alpha_vec (P, prop_dir_alpha / P);
    Rcpp::NumericVector dir_alpha_vec   (P, dir_alpha / P);

    double prop_likelihood = lgamma(prop_dir_alpha) - Rcpp::sum(Rcpp::lgamma(p_dir_alpha_vec));
    for (int j = 0; j < P; j++)
        prop_likelihood += log_var_prob[j] * prop_dir_alpha / P - 1.0;
    
    double likelihood = lgamma(dir_alpha) - Rcpp::sum(Rcpp::lgamma(dir_alpha_vec));
    for (int j = 0; j < P; j++)
        likelihood += log_var_prob[j] * dir_alpha / P - 1.0;

    double ratio = prop_likelihood - 0.5 * log(prop_dir_alpha) - 1.5 * log(P + prop_dir_alpha)
                - likelihood    + 0.5 * log(dir_alpha)      + 1.5 * log(P + dir_alpha);

    if (ratio > log(R::runif(0,1)))
        dir_alpha = prop_dir_alpha;
}

// util functions for separate model
void mh_dir_alpha(
    const int iter, const int total_iter,
    const Rcpp::NumericVector& var_prob,
    const vector<int>& var_count_exp, 
    const vector<int>& var_count_out1, 
    const vector<int>& var_count_out0,
    double& dir_alpha, Rcpp::NumericVector& post_dir_alpha
) {
    int P = var_count_exp.size();
    if (iter < total_iter / 10)
    {
        // warm up
        for (int j = 0; j < P; j++)
        {
            double total_count = var_count_exp[j] + var_count_out1[j] + var_count_out0[j];
            post_dir_alpha[j] = 1.0 + total_count;
        }
    }
    else
    {
        // draw alpha
        draw_dir_alpha(var_prob, dir_alpha);

        // upate dir_alpha
        for (int j = 0; j < P; j++)
        {
            double total_count = var_count_exp[j] + var_count_out1[j] + var_count_out0[j];
            post_dir_alpha[j] = dir_alpha / P + total_count;
        }
    }
}

// util functions for single model
void normalize(Rcpp::NumericVector& var_prob_exp, const Rcpp::NumericVector& var_prob)
{
    int P = var_prob_exp.length();
    double sum_prob = sum(var_prob) - var_prob[P];
    for (int j = 0; j < P; j++)
        var_prob_exp[j] = var_prob[j] / sum_prob;
}

void mh_dir_alpha(
    const Rcpp::NumericVector& var_prob,
    double& dir_alpha, Rcpp::NumericVector& post_dir_alpha
) {
    int P = var_prob.size();
    draw_dir_alpha(var_prob, dir_alpha);
    
    // update post dir alpha
    for (auto& n : post_dir_alpha)
        n = dir_alpha / P;
}

void mh_var_prob(
    const Rcpp::Function& rdirichlet,
    const Rcpp::NumericVector& post_dir_alpha,
    const vector<int>& var_count_exp, 
    const vector<int>& var_count_out, 
    Rcpp::NumericVector& var_prob
) {
    int P = var_count_exp.size(); // var_count_out.size() = P + 1
    int TRT_ID = P;
    Rcpp::NumericVector alpha (P + 1);
    for (int j = 0; j < P; j++)
        alpha[j] = post_dir_alpha[j] + (double) (var_count_exp[j] + var_count_out[j]);
    alpha[TRT_ID] = post_dir_alpha[TRT_ID] + (double) (2 * var_count_out[TRT_ID]);

    // propose new var prob
    Rcpp::NumericVector p_var_prob = rdirichlet(1, alpha);

    // calculate mh ratio
    vector<double> log_var_prob(P + 1), log_prop_var_prob(P + 1);
    {
        const double LB     = pow(0.1, 300);
        const double LOG_LB = -300 * log(10);
        for (int i = 0; i < P + 1; i++)
        {
            // logarize var_prob
            if (var_prob[i] > LB)
                log_var_prob[i] = log(var_prob[i]);
            else
                log_var_prob[i] = LOG_LB;
            
            // logarize prop_var_prob
            if (p_var_prob[i] > LB)
                log_prop_var_prob[i] = log(p_var_prob[i]);
            else
                log_prop_var_prob[i] = LOG_LB;
        }
    }

    //double sum_exp = sum(var_count_exp);
    //var_count_exp.push_back(var_count_out(TRT_IDX));
    double sum_exp = 0.0;
    for (auto n : var_count_exp)
        sum_exp += n;
    double prop_dir_lik = sum_exp * (- log(1 - p_var_prob[TRT_ID]));
    double dir_lik      = sum_exp * (- log(1 -   var_prob[TRT_ID]));
    for (int j = 0; j < P; j++)
    {
        double temp   = var_count_exp[j] + var_count_out[j] + post_dir_alpha[j] - 1.0;
        prop_dir_lik += temp * log_prop_var_prob[j];
        dir_lik      += temp * log_var_prob[j];
    }
    prop_dir_lik += (var_count_out[TRT_ID] + post_dir_alpha[TRT_ID] - 1.0) * log_prop_var_prob[TRT_ID];
    dir_lik      += (var_count_out[TRT_ID] + post_dir_alpha[TRT_ID] - 1.0) * log_var_prob[TRT_ID];

    double ratio = prop_dir_lik - dir_lik;
    for (int j = 0; j < P; j++)
        ratio += (post_dir_alpha[j] + var_count_exp[j] + var_count_out[j] - 1.0) 
        * (log_var_prob[j] - log_prop_var_prob[j]);
    ratio += (post_dir_alpha[TRT_ID] + 2 * var_count_out[TRT_ID] - 1.0)
        * (log_var_prob[TRT_ID] - log_prop_var_prob[TRT_ID]);
        
    if (ratio > log(R::runif(0,1)))
        var_prob = clone(p_var_prob);
}
