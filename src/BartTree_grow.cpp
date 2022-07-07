#include <Rcpp.h>
#include "BartTree.h"

using namespace Rcpp;
using namespace std;

void BartTree::grow(const int t)
{
    const int NUM_OBS      = X.nrow();
    const int NUM_VAR      = X.ncol();
    const int TRT_IDX      = X.ncol();

    BartNode* prop_node    = nullptr;
    int       prop_var_idx = -1;
    int       prop_cut_idx = -1;

    int num_terminal_nodes, num_uniques, num_new_singly;
    double rule, log_prop_prob;

    if (root_nodes_[t]->isTerminal())
    {
        // there is only one node, so grow from root node
        prop_node          = root_nodes_[t];
        num_terminal_nodes = 1;

        switch (model)
        {
            // 0 : sep model
            // 1 : marginal exposure model
            case 0: // sep model
                prop_var_idx = sample(NUM_VAR, 1, false, var_prob_)(0) - 1;
                break;

            case 1: 
                // mar exposure model
                // in marginal model, last index of var_prob indicate TRT.
                // In marginal exposure model, we do not use TRT.
                do
                {
                    prop_var_idx = sample(NUM_VAR + 1, 1, false, var_prob_)(0) - 1;
                } 
                while (prop_var_idx == TRT_IDX);
                break;

            case 2: 
                // marginal outcome model
                // in marginal model, last index of var_prob indicate TRT.
                // In marginal outcome model, we can select TRT.
                prop_var_idx = sample(NUM_VAR + 1, 1, false, var_prob_)(0) - 1;
                break;
        }

        num_uniques    = Xcut[prop_var_idx].length();
        num_new_singly = 1; // number of singly nodes in new tree
        log_prop_prob  = log(var_prob_(prop_var_idx));

        double min_cutpoint = min(Xcut[prop_var_idx]);
        switch (num_uniques)
        {
            case 1: 
                // cannot grow -> all values are same
                return;
            
            case 2:
                // binary variable
                rule         = max(Xcut[prop_var_idx]);
                prop_cut_idx = 1;
                break;

            default:
                // other types of variables
                do
                {
                    // select rule but it should not be min(Xcut)
                    prop_cut_idx = sample(Xcut[prop_var_idx].length(), 1)(0) - 1;
                    rule         = Xcut[prop_var_idx](prop_cut_idx);
                }
                while (rule == min_cutpoint);
        }
    }
    else
    {
        // propose a terminal node for grow
        vector<BartNode*> terminal_nodes = getTerminalNodes(t);
        num_terminal_nodes = terminal_nodes.size();
        prop_node = terminal_nodes[sample(num_terminal_nodes, 1)(0) - 1];

        // propose variable and observation for new rule
        // we will find variable with enough unique observations (at least 2)
        LogicalVector obs_flag = rep(false, NUM_OBS);
        LogicalVector var_flag = rep(false, var_prob_.length());

        // find observation with proposed node and set it as baseline observation
        // if other observation has different value with baseline observation
        // then it has at least two unique observation
        vector<BartNode*>::iterator begin = assigned_nodes_[t].begin();
        vector<BartNode*>::iterator end   = assigned_nodes_[t].end();
        int baseline_idx = find(begin, end, prop_node) - begin;

        updateFlag(obs_flag, var_flag, baseline_idx, prop_node, t, true);

        if ((sum(obs_flag) == 1) || (sum(var_flag) == 0))
        {
            // there is no predictor with unique values
            return;
        }

        // sample predictor with unique values
        NumericVector flagged_var_prob = clone(var_prob_);
        for (int i = 0; i < flagged_var_prob.length(); i++)
        {
            if (!var_flag(i)) flagged_var_prob(i) = 0;
        }
        prop_var_idx = sample(var_prob_.length(), 1, false, flagged_var_prob)(0) - 1;

        // count unique values of chosen predictor
        num_uniques  = countUniqueValues(prop_var_idx);
        if (num_uniques == 1)
            return;

        // sample cut point with sampled predictor
        rule           = proposeRule(obs_flag, prop_var_idx, num_uniques);
        prop_cut_idx   = findCutIdx(prop_var_idx, num_uniques, rule);
        
        // count number of singly node with new node
        num_new_singly = countSinglyFromProposedTree(prop_node, t);  
        log_prop_prob  = log(var_prob_(prop_var_idx)) - log(sum(flagged_var_prob));
    } // end of proposing new node


    // calculate likelihood of current tree
    int    num_left      = 0,   num_right      = 0;
    double residual_left = 0.0, residual_right = 0.0;
    #pragma omp parallel for reduction(+ : num_left, num_right, residual_left, residual_right) if (parallel)
    for (int i = 0; i < NUM_OBS; i++)
    {
        const BartNode* assigned_node = assigned_nodes_[t][i];
        if (assigned_node == prop_node)
        {
            double value = getXValue(i, prop_var_idx);
            if (value < rule)
            {
                num_left       += 1;
                residual_left  += residual_(i);
            }
            else
            {
                num_right      += 1;
                residual_right += residual_(i);
            }
        }
    }

    double transition = 
        log(step_prob(1))      + log(num_terminal_nodes) - log_prop_prob
        + log(num_uniques - 1) - log(step_prob(0))       - log(num_new_singly)
    ;

    double likelihood = 
        0.5 *   log(sigma2_) 
        + 0.5 * log(sigma2_ + sigma_mu * (num_left + num_right))
        - 0.5 * log(sigma2_ + sigma_mu * num_left)
        - 0.5 * log(sigma2_ + sigma_mu * num_right)
        + (sigma_mu / (2 * sigma2_)
        *(pow(residual_left,                  2) / (sigma2_ + sigma_mu *  num_left)
        + pow(residual_right,                 2) / (sigma2_ + sigma_mu *  num_right)
        - pow(residual_left + residual_right, 2) / (sigma2_ + sigma_mu * (num_left + num_right))))
    ;

    int depth = prop_node->getDepth();
    double structure = (
        log(alpha) + 2*log(1 - alpha / pow(2 + depth, beta))
        - log(pow(1 + depth, beta) - alpha)
        + log_prop_prob
        - log(num_uniques - 1)
    );

    double ratio = transition + likelihood + structure;

    // accept or reject tree
    if (ratio > log(R::runif(0, 1)))
    {
        prop_node->grow(prop_var_idx, prop_cut_idx);

        // update assigned nodes
        #pragma omp parallel for if (parallel)
        for (int i = 0; i < NUM_OBS; i++)
        {
            const BartNode* assigned_node = assigned_nodes_[t][i];
            if (assigned_node == prop_node)
            {
                double value = getXValue(i, prop_var_idx);
                if (value < rule)
                    assigned_nodes_[t][i] = prop_node->getChildLeft();
                else
                    assigned_nodes_[t][i] = prop_node->getChildRight();
            }
        }
    }
}
