// BartTree.h
#pragma once

#include <Rcpp.h>
#include "BartNode.h"

using namespace Rcpp;
using namespace std;

// BART model
class BartTree
{
    vector<BartNode*>            root_nodes_;
    vector<vector<BartNode*>>    assigned_nodes_;
    NumericVector&               residual_;
    NumericMatrix                leaf_values_;
    NumericVector&               var_prob_;
    double                       sigma2_;

    // const variables
    const int                    model_;     // 0 : sep,  1 : mar_exp, 2 : mar_out
    const NumericVector&         trt_;
    const NumericMatrix&         X_;
    const vector<NumericVector>& Xcut_;
    const NumericVector&         step_prob_; // 0 : GROW, 1 : PRUNE,   2 : CHANGE
    const int                    num_tree_;
    const double                 alpha_;
    const double                 beta_;
    const double                 sigma_mu_;
    const bool                   parallel_;

public:
    ~BartTree()
    {
        for (int t = 0; t < num_tree_; t++)
            delete root_nodes_[t];
    };
    BartTree(
        NumericVector&               residual,
        NumericVector&               var_prob,
        double&                      sigma2,
        const int                    model,
        const NumericVector&         trt,
        const NumericMatrix&         X,
        const vector<NumericVector>& Xcut,
        const NumericVector&         step_prob,
        const int                    num_tree,
        const double                 alpha,
        const double                 beta,
        const double                 sigma_mu,
        const bool                   parallel
    ) :
        residual_(residual),
        var_prob_(var_prob),
        sigma2_(sigma2),
        model_(model),
        trt_(trt),
        X_(X),
        Xcut_(Xcut),
        step_prob_(step_prob),
        num_tree_(num_tree),
        alpha_(alpha),
        beta_(beta),
        sigma_mu_(sigma_mu),
        parallel_(parallel)
    {
        // initialize leaf_values_, root_nodes_ and assigned_nodes_
        leaf_values_ = NumericMatrix(X.nrow(), num_tree);

        root_nodes_.resize(num_tree);
        for (int t = 0; t < num_tree; t++)
            root_nodes_[t] = new BartNode();

        assigned_nodes_.resize(num_tree);
        for (int t = 0; t < num_tree; t++)
            assigned_nodes_[t] = vector<BartNode*> (X.nrow(), root_nodes_[t]);
    };

    inline bool isSeparateModel() const { return (model_ == 0); };
    inline bool isMarginalModel() const { return (model_); };

    // getters
    inline double        getSigma2()              const { return sigma2_; };
    inline NumericVector getFittedValues()        const { return rowSums(leaf_values_); };
    inline BartNode*     getRootNode(const int t) const { return root_nodes_[t]; };
    inline double        getXValue(const int obs_idx, const int var_idx) const
    {
        return (var_idx == X_.ncol()) ? trt_(obs_idx) : X_(obs_idx, var_idx);
    };
    vector<BartNode*>    getTerminalNodes(const int t) const
    {
        vector<BartNode*> res;
        root_nodes_[t]->getTerminalNodes(res);
        return res;
    };
    vector<BartNode*>    getSinglyNodes(const int t)   const
    {
        vector<BartNode*> res;
        root_nodes_[t]->getSinglyNodes(res);
        return res;
    };

    // step related methods
    void step(const NumericVector& latent_variable, const bool is_binary_trt);
    void grow(const int t);
    void prune(const int t);
    void change(const int t);
    void drawLeafValue(const int t);

    // predict methods
    double predict(const NumericMatrix& X);
    double predict(const double trt_value);

    // update methods
    void updateFlag(
        LogicalVector&        obs_flag,
        LogicalVector&        var_flag,
        const int             baseline_idx,
        const BartNode* const prop_node,
        const int             t,
        const bool            is_terminal
    );
    void updateLatentVariable(NumericVector& latent_variable, const bool is_binary_trt);
    void updateResidual(const NumericVector& latent_variable, const int t);
    void updateSigma2(
        const Function&      rinvgamma,
        const NumericVector& Y,
        const double         nu,
        const double         lambda
    );
    void updateDirAlpha(double& dir_alpha);
    void updateVarProb(
        NumericVector&       var_count_exp,
        const NumericVector& var_count_out,
        const Function&      rdirichlet,
        const NumericVector& post_dir_alpha
    );

    // util methods
    double        findMinValue(const LogicalVector& obs_flag, const int prop_var_idx);
    double        findMaxValue(const LogicalVector& obs_flag, const int prop_var_idx);
    double        proposeRule( const LogicalVector& obs_flag, const int prop_var_idx, const int num_uniques);
    int           findCutIdx(const int prop_var_idx, const int num_uniques, const double rule);
    int           countUniqueValues(const BartNode* prop_node, const int prop_var_idx, const int t, const bool is_grow);
    int           countSinglyFromProposedTree(const BartNode* prop_node, const int t);
    NumericVector countSelectedVariables();
};

