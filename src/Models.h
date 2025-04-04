#pragma once
#include "BART.h"
#include "util.h"

using namespace std;

class SeparateModel: public BART
{
public:
    SeparateModel(
        const vector<vector<double>>& X, const vector<vector<double>>& Xcut, const int n_tree, 
        const Rcpp::NumericVector& step_prob, const double alpha, const double beta, Rcpp::NumericVector& var_prob, 
        const bool parallel
    );
    void predict(Rcpp::NumericVector& outcome, Rcpp::NumericMatrix& outcome_sample, const int id, const vector<vector<double>>& full_X) const;
};

class SingleModel: public BART
{
private:
public:
    SingleModel(
        const vector<vector<double>>& X, const vector<vector<double>>& Xcut, const int N, const int P,
        const int n_tree, const Rcpp::NumericVector& step_prob,
        const double alpha, const double beta, Rcpp::NumericVector& var_prob, 
        const bool parallel
    );
    void set_prob(const Rcpp::NumericVector& prob);
    void predict(Rcpp::NumericVector& outcome, Rcpp::NumericMatrix& outcome_sample, const int id, const double trt_value) const;
};
