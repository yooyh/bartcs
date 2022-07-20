#include <Rcpp.h>
#include "BartTree.h"

using namespace Rcpp;
using namespace std;

double BartTree::findMinValue(const LogicalVector& obs_flag, const int prop_var_idx)
{
    const int NUM_OBS    = X.nrow();
    const int TRT_IDX    = X.ncol();
    double    min_value  = __DBL_MAX__;
    if ((model == 2) && (prop_var_idx == TRT_IDX))
    {
        #ifdef _OPENMP
            #pragma omp parallel for if (parallel)
        #endif
        for (int i = 0; i < NUM_OBS; i++)
        {
            if (obs_flag(i))
            {
                if (trt(i) < min_value)
                {
                    #ifdef _OPENMP
                        #pragma omp critical
                    #endif
                    min_value = trt(i);
                }
            }
        }
    }
    else
    {
        #ifdef _OPENMP
            #pragma omp parallel for if (parallel)
        #endif
        for (int i = 0; i < NUM_OBS; i++)
        {
            if (obs_flag(i))
            {
                if (X(i, prop_var_idx) < min_value)
                {
                    #ifdef _OPENMP
                        #pragma omp critical
                    #endif
                    min_value = X(i, prop_var_idx);
                }
            }
        }
    }
    return min_value;
}

double BartTree::findMaxValue(const LogicalVector& obs_flag, const int prop_var_idx)
{
    const int NUM_OBS    = X.nrow();
    const int TRT_IDX    = X.ncol();
    double    max_value  = __DBL_MIN__;
    if ((model == 2) && (prop_var_idx == TRT_IDX))
    {
        #ifdef _OPENMP
            #pragma omp parallel for if (parallel)
        #endif
        for (int i = 0; i < NUM_OBS; i++)
        {
            if (obs_flag(i))
            {
                if (trt(i) > max_value)
                {
                    #ifdef _OPENMP
                        #pragma omp critical
                    #endif
                    max_value = trt(i);
                }
            }
        }
    }
    else
    {
        #ifdef _OPENMP
            #pragma omp parallel for if (parallel)
        #endif
        for (int i = 0; i < NUM_OBS; i++)
        {
            if (obs_flag(i))
            {
                if (X(i, prop_var_idx) > max_value)
                {
                    #ifdef _OPENMP
                        #pragma omp critical
                    #endif
                    max_value = X(i, prop_var_idx);
                }
            }
        }
    }
    return max_value;
}

double BartTree::proposeRule(const LogicalVector& obs_flag, const int prop_var_idx, const int num_uniques)
{
    const int    NUM_OBS   = X.nrow();
    const int    TRT_IDX   = X.ncol();

    const double MIN_VALUE = findMinValue(obs_flag, prop_var_idx);

    double rule;
    switch (num_uniques)
    {   
        case 2:
            // binary variables
            rule = findMaxValue(obs_flag, prop_var_idx);
            break;
        
        default:
            // non-binary variables
            int prop_obs_idx;
            if ((model == 2) && (prop_var_idx == TRT_IDX))
            {
                // marginal outcome model and trt is chosen
                do
                {
                    prop_obs_idx = sample(NUM_OBS, 1)(0) - 1;
                }
                while (!obs_flag(prop_obs_idx) || (trt(prop_obs_idx) == MIN_VALUE));
                rule = trt(prop_obs_idx);
            }
            else
            {
                do
                {
                    prop_obs_idx = sample(NUM_OBS, 1)(0) - 1;
                }
                while (!obs_flag(prop_obs_idx) || (X(prop_obs_idx, prop_var_idx) == MIN_VALUE));
                rule = X(prop_obs_idx, prop_var_idx);
            }
    }
    return rule;
}

int BartTree::findCutIdx(const int prop_var_idx, const int num_uniques, const double rule)
{
    if (num_uniques == 2)
        return 1;

    const int NUM_OBS = Xcut[prop_var_idx].length();

    int cut_idx = -1;
    volatile bool found = false;
    #ifdef _OPENMP
        #pragma omp parallel for shared(found) if (parallel)
    #endif
    for (int i = 0; i < NUM_OBS; i++)
    {
        if (found) continue;
        if (Xcut[prop_var_idx](i) == rule)
        {
            cut_idx = i;
            found   = true;
        }
    }
    return cut_idx;
}

int BartTree::countUniqueValues(const int prop_var_idx)
{
    const int NUM_OBS = X.nrow();
    map<double, bool> map_uniques;
    if (prop_var_idx == X.ncol())
    {
        // is treatment
        for (int i = 0; i < NUM_OBS; i++)
        {
            map_uniques[trt(i)] = true;
        }
    }
    else
    {
        for (int i = 0; i < NUM_OBS; i++)
        {
            map_uniques[X(i, prop_var_idx)] = true;
        }
    }
    
    return map_uniques.size();
}

int BartTree::countSinglyFromProposedTree(const BartNode* prop_node, const int t)
{
    vector<BartNode*> singly_nodes = getSinglyNodes(t);
    int res = singly_nodes.size();
    if (!prop_node->getParent()->isSingly())
        res++;
    return res;
}

NumericVector BartTree::countSelectedVariables()
{
    NumericVector res(var_prob_.length());
    for (int t = 0; t < num_tree; t++)
    {
        root_nodes_[t]->countSelectedVariables(res);
    }
    return res;
}

void BartTree::updateFlag(
    LogicalVector&        obs_flag, 
    LogicalVector&        var_flag, 
    const int             baseline_idx, 
    const BartNode* const prop_node, 
    const int             t,
    const bool            is_terminal
) {
    const int NUM_OBS = X.nrow();
    const int NUM_VAR = X.ncol();
    const int TRT_IDX = X.ncol();

    volatile bool found = false;
    #ifdef _OPENMP
        #pragma omp parallel for shared(found) if (parallel)
    #endif
    for (int i = baseline_idx; i < NUM_OBS; i++)
    {
        const BartNode* assigned_node = assigned_nodes_[t][i];
        bool cond_terminal = (assigned_node == prop_node);
        bool cond_singly   = (assigned_node->getParent() == prop_node);
        if (is_terminal ? cond_terminal : cond_singly)
        {
            obs_flag(i) = true; // flag obs with prop_node
            if (found)
            {
                // skip remaining part because all predictors are with unique values 
                continue;
            }

            switch (model)
            {
                case 2:
                    // marginal outcome model
                    if (var_flag(TRT_IDX)) continue;
                    if (trt(i) != trt(baseline_idx))
                    {
                        #ifdef _OPENMP
                            #pragma omp atomic
                        #endif
                        var_flag(TRT_IDX) |= true;
                    }

                default:
                    for (int j = 0; j < NUM_VAR; j++)
                    {
                        if (var_flag(j)) continue;
                        if (X(i, j) != X(baseline_idx, j))
                        {
                            #ifdef _OPENMP
                                #pragma omp atomic
                            #endif
                            var_flag(j) |= true;
                        }
                    }
            } // end switch

            if (sum(var_flag) == var_flag.length())
            {
                // all predictors has unique values
                #ifdef _OPENMP
                    #pragma omp atomic
                #endif
                found |= true;
            }
        }
    } // end for i
}
