#pragma once
#include "Node.h"

using namespace std;

class BART
{
protected:
    // data
    const vector<vector<double>>& X_;
    const vector<vector<double>>& Xcut_;
    const int N, P;

    // parameter
    const int n_tree_;
    vector<Node> tree_;
    const Rcpp::NumericVector& step_prob_;
    const double alpha_, beta_;
    double sigma_mu_, sigma2_;
    Rcpp::NumericVector& prob_;

    // other
    vector<double> fitted_;
    vector<double> fit_tmp_;
    vector<double> residual_;
    vector<int> var_count_;
    const bool parallel_;

public:
    BART(
        const vector<vector<double>>& X, const vector<vector<double>>& Xcut, const int N, const int P,
        const int n_tree, const Rcpp::NumericVector& step_prob,
        const double alpha, const double beta, Rcpp::NumericVector& var_prob, const bool parallel
    );
    void init(const vector<double>& Y, double sigma2);

    // get, set
    const double sigma2() const {return sigma2_;}
    const vector<int>& var_count() const {return var_count_;}
    double get_sigma_mu(const vector<double>& Y, const int n_tree) const;

    // draw and step
    void draw(const vector<double>& Y);
    void step(Node& tree);
    void grow(Node& tree);
    void prune(Node& tree);
    void change(Node& tree);
    void draw_mu(Node& tree);
    void draw_sigma2(const Rcpp::Function& rinvgamma, const vector<double>& Y, const double nu, const double lambda);
    void update_Z(vector<double>& Z, const Rcpp::NumericVector& TRT, const bool binary_trt);

    // util function
    void get_SS_grow(
        const Node& tree, const Node* prop_node, const int var, const int cut,
        int& nl, int& nr, double& rl, double& rr, int& n_unique
    ) const;
    void get_SS_prune(
        const Node& tree, const Node* prop_node, const int var, const int cut,
        int& nl, int& nr, double& rl, double& rr, int& n_unique
    ) const;
    void get_SS_change(
        const Node& tree, const Node* prop_node,
        const int cvar, const int ccut, int& cnl, int& cnr, double& crl, double& crr,
        const int pvar, const int pcut, int& pnl, int& pnr, double& prl, double& prr
    ) const;
    void get_ratio(
        const int& n_unique, const int& n_terminal, const int& n_singly, 
        const int depth, const double& log_prob,
        const int& nl, const int& nr, const double& rl, const double& rr,
        double& ratio
    ) const;
    void fit(const Node& tree, vector<double>& res) const;
    void get_vars(const Node* node, vector<int>& vars) const;
};
