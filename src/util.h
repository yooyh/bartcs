#pragma once
#include "Node.h"

// util functions for both model
void get_Xcut(const vector<vector<double>>& X, vector<vector<double>>& Xcut);
void get_data(const Rcpp::NumericMatrix& X_src, vector<vector<double>>& X, vector<vector<double>>& Xcut);
void get_data(const Rcpp::NumericMatrix& X_src, const Rcpp::NumericVector& TRT_src, 
              vector<vector<double>>& X, vector<vector<double>>& Xcut);
void init_Z(vector<double>& Z, const Rcpp::NumericVector& TRT_src, const bool binary_trt);
void draw_dir_alpha(const Rcpp::NumericVector& var_prob, double& dir_alpha);

// util functions for separate model
void mh_dir_alpha(
    const int iter, const int total_iter,
    const Rcpp::NumericVector& var_prob,
    const vector<int>& var_count_exp, 
    const vector<int>& var_count_out1, 
    const vector<int>& var_count_out0,
    double& dir_alpha, Rcpp::NumericVector& post_dir_alpha
);

// util functions for single model
void normalize(Rcpp::NumericVector& var_prob_exp, const Rcpp::NumericVector& var_prob);
void mh_dir_alpha(const Rcpp::NumericVector& var_prob, double& dir_alpha, Rcpp::NumericVector& post_dir_alpha);
void mh_var_prob(
    const Rcpp::Function& rdirichlet,
    const Rcpp::NumericVector& post_dir_alpha,
    const vector<int>& var_count_exp, 
    const vector<int>& var_count_out, 
    Rcpp::NumericVector& var_prob
);
