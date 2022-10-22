#include <Rcpp.h>
#include "BartTree.h"

using namespace Rcpp;
using namespace std;

double BartTree::predict(const NumericMatrix& X) 
{
    // predict function for sep model
    const int NUM_BOOT = X.nrow();
    double    res      = 0.0;

    #ifdef _OPENMP
        #pragma omp parallel for reduction(+ : res) if (parallel_)
    #endif
    for (int i = 0; i < NUM_BOOT; i++)
    {
        double temp = 0.0;
        for (int t = 0; t < num_tree_; t++)
        {
            const BartNode* current_node = root_nodes_[t];
            while (!current_node->isTerminal())
            {
                double value = X(i, current_node->getVarIdx());
                if (value < Xcut_[current_node->getVarIdx()][current_node->getCutIdx()])
                    current_node = current_node->getChildLeft();
                else
                    current_node = current_node->getChildRight();
            }
            temp += current_node->getLeafValue();
        }
        res += temp;
    }
    res /= NUM_BOOT;
    return res;
}

double BartTree::predict(const double trt_value) 
{
    // predict function for marginal model
    const int NUM_OBS = X_.nrow();
    const int TRT_IDX = X_.ncol();

    double    res = 0.0;

    #ifdef _OPENMP
        #pragma omp parallel for reduction(+ : res) if (parallel_)
    #endif
    for (int i = 0; i < NUM_OBS; i++)
    {
        double temp = 0.0;
        for (int t = 0; t < num_tree_; t++)
        {
            const BartNode* current_node = root_nodes_[t];
            while (!current_node->isTerminal())
            {
                double value;
                if (current_node->getVarIdx() == TRT_IDX)
                    value = trt_value;
                else
                    value = X_(i, current_node->getVarIdx());

                if (value < Xcut_[current_node->getVarIdx()][current_node->getCutIdx()])
                    current_node = current_node->getChildLeft();
                else
                    current_node = current_node->getChildRight();
            }
            temp += current_node->getLeafValue();
        }
        res += temp;
    }
    res /= NUM_OBS;
    return res;
}
