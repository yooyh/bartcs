#include <Rcpp.h>
#include "BartTree.h"

using namespace Rcpp;
using namespace std;

void BartTree::step(
    const NumericVector& latent_variable, 
    const bool           is_binary_trt
) {
    for (int t = 0; t < num_tree; t++)
    {
        updateResidual(latent_variable, t);
        // select node for the step
        if (root_nodes_[t]->isTerminal())
        {
            // root node has no child node -> can only grow
            grow(t);
        }
        else
        {
            int step = sample(3, 1, false, step_prob)(0);
            switch (step)
            {
                case 1:
                    grow(t);
                    break;

                case 2:
                    prune(t);
                    break;

                case 3:
                    change(t);
                    break;

                default: {};
            }
        }
        drawLeafValue(t);
    }
}
