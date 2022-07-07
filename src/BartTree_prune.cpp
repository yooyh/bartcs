#include <Rcpp.h>
#include "BartTree.h"

using namespace Rcpp;
using namespace std;

void BartTree::prune(const int t)
{
    const int NUM_OBS = X.nrow();
    // const int NUM_VAR = X.ncol();
    // const int TRT_IDX = X.ncol();

    // propose a singly node to prune
    vector<BartNode*> singly_nodes = getSinglyNodes(t);
    BartNode*         prop_node    = singly_nodes[sample(singly_nodes.size(), 1)(0) - 1];
    int               var_idx      = prop_node->getVarIdx();

    // find variables with unique values
    // and calculate likelihood of current tree
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

    int    num_left      = 0,   num_right      = 0;
    double residual_left = 0.0, residual_right = 0.0;
    #pragma omp parallel for reduction(+ : num_left, num_right, residual_left, residual_right) if (parallel)
    for (int i = 0; i < NUM_OBS; i++)
    {
        const BartNode* assigned_node = assigned_nodes_[t][i];
        if (assigned_node->getParent() == prop_node)
        {
            if (assigned_node == prop_node->getChildLeft())
            {
                num_left      += 1;
                residual_left += residual_(i);
            }
            else
            {
                num_right      += 1;
                residual_right += residual_(i);
            }
        }
    }

    if ((sum(obs_flag) == 1) || (sum(var_flag) == 0))
    {
        // there is no predictor with unique values
        return;
    }

    // find leaf node for computation of transition
    vector<BartNode*> terminal_nodes = getTerminalNodes(t);

    NumericVector flagged_var_prob   = var_prob_[var_flag];
    double        log_prop_prob      = log(var_prob_(var_idx)) - log(sum(flagged_var_prob));
    int           num_uniques        = countUniqueValues(var_idx);

    double transition = 
        log(step_prob(0))  - log(terminal_nodes.size() - 1) + log_prop_prob
        - log(num_uniques) - log(step_prob(1))              + log(singly_nodes.size())
    ;

    double likelihood = 
        - 0.5 * log(sigma2_)
        - 0.5 * log(sigma2_ + sigma_mu * (num_left + num_right))
        + 0.5 * log(sigma2_ + sigma_mu * num_left)
        + 0.5 * log(sigma2_ + sigma_mu * num_right)
        + (sigma_mu / (2*sigma2_)
        * (- pow(residual_left,                   2) / (sigma2_ + sigma_mu *  num_left)
           - pow(residual_right,                  2) / (sigma2_ + sigma_mu *  num_right)
           + pow(residual_left + residual_right,  2) / (sigma2_ + sigma_mu * (num_left + num_right))))
    ;

    int depth = prop_node->getDepth();
    double structure = 
        - log(alpha) - 2*log(1 - alpha / pow(2 + depth, beta))
        + log(pow(1 + depth, beta) - alpha)
        - log_prop_prob
        + log(num_uniques)
    ;

    double ratio = transition + likelihood + structure;

    if (ratio > log(R::runif(0, 1)))
    {
        // update assigned nodes
        #pragma omp parallel for if (parallel)
        for (int i = 0; i < NUM_OBS; i++)
        {
            const BartNode* assigned_node = assigned_nodes_[t][i];
            if (assigned_node->getParent() == prop_node)
                assigned_nodes_[t][i] = prop_node;
        }

        prop_node->prune();
    }
}
