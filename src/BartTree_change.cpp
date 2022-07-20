#include <Rcpp.h>
#include "BartTree.h"

using namespace Rcpp;
using namespace std;

void BartTree::change(const int t)
{
    const int NUM_OBS = X.nrow();
    // const int NUM_VAR = X.ncol();
    // const int TRT_IDX = X.ncol();

    // propose a singly node to change
    vector<BartNode*> singly_nodes = getSinglyNodes(t);
    BartNode*         prop_node    = singly_nodes[sample(singly_nodes.size(), 1)(0) - 1];

    // propose predictor and observation for new rule
    // we will use obs_flag and var_flag
    // to find predictor with enough unique observations (at least 2)
    LogicalVector obs_flag = rep(false, NUM_OBS);
    LogicalVector var_flag = rep(false, var_prob_.length());

    // find observation with proposed node and set it as baseline observation
    // if other observation has different value with baseline observation
    // then it has at least two unique observation
    vector<BartNode*>::iterator begin = assigned_nodes_[t].begin();
    vector<BartNode*>::iterator end   = assigned_nodes_[t].end();
    int baseline_idx = find_if(begin, end,
        [&](BartNode* node)->bool {return node->getParent() == prop_node;}
    ) - begin;

    updateFlag(obs_flag, var_flag, baseline_idx, prop_node, t, false);

    if ((sum(obs_flag) == 1) || (sum(var_flag) == 0))
    {
        // there is no predictor with unique values
        return;
    }

    NumericVector flagged_var_prob = clone(var_prob_);
    for (int i = 0; i < flagged_var_prob.length(); i++)
    {
        if (!var_flag(i)) flagged_var_prob(i) = 0;
    }
    int prop_var_idx = sample(var_prob_.length(), 1, false, flagged_var_prob)(0) - 1;

    // count unique values of chosen predictor
    int num_uniques  = countUniqueValues(prop_var_idx);
    if (num_uniques == 1)
        return;

    // sample cut point with sampled predictor
    double rule          = proposeRule(obs_flag, prop_var_idx, num_uniques);
    int    prop_cut_idx  = findCutIdx(prop_var_idx, num_uniques, rule);

    int    num_left      = 0,   num_right      = 0,   prop_num_left      = 0,   prop_num_right      = 0;
    double residual_left = 0.0, residual_right = 0.0, prop_residual_left = 0.0, prop_residual_right = 0.0;
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+ : num_left, num_right, residual_left, residual_right, prop_num_left, prop_num_right, prop_residual_left, prop_residual_right) if (parallel)
    #endif
    for (int i = 0; i < NUM_OBS; i++)
    {
        const BartNode* assigned_node = assigned_nodes_[t][i];

        if (assigned_node->getParent() == prop_node)
        {
            // compute current likelihood
            if (assigned_node == prop_node->getChildLeft())
            {
                num_left       += 1;
                residual_left  += residual_(i);
            }
            else
            {
                num_right      += 1;
                residual_right += residual_(i);
            }
            // compute new likelihood
            double value = getXValue(i, prop_var_idx);
            if (value < rule) {
                prop_num_left       += 1;
                prop_residual_left  += residual_(i);
            }
            else
            {
                prop_num_right      += 1;
                prop_residual_right += residual_(i);
            }
        }
    }

    double likelihood = 
        0.5 *   log(sigma2_ / sigma_mu + num_left)
        + 0.5 * log(sigma2_ / sigma_mu + num_right)
        - 0.5 * log(sigma2_ / sigma_mu + prop_num_left)
        - 0.5 * log(sigma2_ / sigma_mu + prop_num_right)
        + (0.5 / sigma2_
        *(pow(prop_residual_left,  2) / (sigma2_ / sigma_mu + prop_num_left)
        + pow(prop_residual_right, 2) / (sigma2_ / sigma_mu + prop_num_right)
        - pow(residual_left,       2) / (sigma2_ / sigma_mu + num_left)
        - pow(residual_right,      2) / (sigma2_ / sigma_mu + num_right)))
    ;

    if (likelihood > log(R::runif(0, 1)))
    {
        prop_node->change(prop_var_idx, prop_cut_idx);

        // update assigned nodes
        #ifdef _OPENMP
            #pragma omp parallel for if (parallel)
        #endif
        for (int i = 0; i < NUM_OBS; i++) 
        {
            const BartNode* assigned_node = assigned_nodes_[t][i];
            if (assigned_node->getParent() == prop_node)
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


