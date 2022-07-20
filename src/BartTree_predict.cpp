#include <Rcpp.h>
#include "BartTree.h"

using namespace Rcpp;
using namespace std;

double BartTree::predict(
    const IntegerVector& boot_idx,
    const NumericMatrix& boot_X
) {
    // predict function for sep model
    const int NUM_BOOT = boot_idx.length();
    double    res      = 0.0;

    #ifdef _OPENMP
        #pragma omp parallel for reduction(+ : res) if (parallel)
    #endif
    for (int i = 0; i < NUM_BOOT; i++)
    {
        double temp = 0.0;
        for (int t = 0; t < num_tree; t++)
        {
            const BartNode* current_node = root_nodes_[t];
            while (!current_node->isTerminal())
            {
                double value = boot_X(boot_idx(i), current_node->getVarIdx());
                if (value < Xcut[current_node->getVarIdx()][current_node->getCutIdx()])
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

double BartTree::predict(
        const IntegerVector& boot_idx,
        const double         trt_value
) {
    // predict function for marginal model
    const int NUM_BOOT = boot_idx.length();
    const int TRT_IDX  = X.ncol();

    double    res = 0.0;

    #ifdef _OPENMP
        #pragma omp parallel for reduction(+ : res) if (parallel)
    #endif
    for (int i = 0; i < NUM_BOOT; i++)
    {
        double temp = 0.0;
        for (int t = 0; t < num_tree; t++)
        {
            const BartNode* current_node = root_nodes_[t];
            while (!current_node->isTerminal())
            {
                double value;
                if (current_node->getVarIdx() == TRT_IDX)
                    value = trt_value;
                else
                    value = X(boot_idx(i), current_node->getVarIdx());

                if (value < Xcut[current_node->getVarIdx()][current_node->getCutIdx()])
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

