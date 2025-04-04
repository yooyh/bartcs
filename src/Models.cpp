#include "Models.h"

// method for separate model
SeparateModel::SeparateModel(
    const vector<vector<double>>& X, const vector<vector<double>>& Xcut, const int n_tree, 
    const Rcpp::NumericVector& step_prob, const double alpha, const double beta, Rcpp::NumericVector& var_prob, 
    const bool parallel
) : BART(X, Xcut, X.size(), Xcut.size(), n_tree, step_prob, alpha, beta, var_prob, parallel) {}

void SeparateModel::predict(Rcpp::NumericVector& outcome, Rcpp::NumericMatrix& outcome_sample, const int id, const vector<vector<double>>& full_X) const
{
    int N = full_X.size();
    double res = 0.0;

    #ifdef _OPENMP
        #pragma omp parallel for reduction(+ : res) if (parallel_)
    #endif
    for (const auto& x : full_X)
    for (int i = 0; i < N; i++)
    {
        double sum_mu = 0.0;
        for (auto& tree : tree_)
            sum_mu += tree.assigned_node(Xcut_, x)->mu();
        outcome_sample(i, id) = sum_mu;
        res += sum_mu;
    }
    outcome(id) = res / N;
}

// method for single model
SingleModel::SingleModel(
    const vector<vector<double>>& X, const vector<vector<double>>& Xcut, const int N, const int P,
    const int n_tree, const Rcpp::NumericVector& step_prob, 
    const double alpha, const double beta, Rcpp::NumericVector& var_prob, 
    const bool parallel
) : BART(X, Xcut, N, P, n_tree, step_prob, alpha, beta, var_prob, parallel) {}

void SingleModel::set_prob(const Rcpp::NumericVector& prob)
{
    int P = prob_.size();
    double sum_prob = sum(prob) - prob(P);

    #ifdef _OPENMP
        #pragma omp parallel for if (parallel_)
    #endif
    for (int i = 0; i < P; i++)
        prob_(i) = prob(i) / sum_prob;
}

void SingleModel::predict(Rcpp::NumericVector& outcome, Rcpp::NumericMatrix& outcome_sample, const int id, const double trt_value) const
{
    int N = X_.size();
    int TRT_id = P - 1;
    double res = 0.0;
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+ : res) if (parallel_)
    #endif
    for (int i = 0; i < N; i++)
    {
        double sum_mu = 0.0;
        for (const auto& tree : tree_)
        {
            auto node = &tree;
            while (!node->is_terminal())
            {
                int var = node->var();
                int cut = node->cut();
                double x = (var == TRT_id) ? trt_value : X_[i][var];
                if (x < Xcut_[var][cut])
                    node = node->left();
                else
                    node = node->right();
            }
            sum_mu += node->mu();
        }
        outcome_sample(i, id) = sum_mu;
        res += sum_mu;
    }
    outcome(id) = res / N;
}