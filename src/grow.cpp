#include <Rcpp.h>
#include "BART.h"

using namespace std;

void BART::grow(Node& tree)
{
    vector<Node*> tnodes;
    tree.get_terminal_nodes(tnodes);

    // draw node to grow
    int node_id = Rcpp::sample(tnodes.size(), 1)(0) - 1;
    Node* prop_node = tnodes[node_id];

    // draw variable
    vector<int> vars; // variable id can split
    get_vars(prop_node, vars);
    if (vars.size() == 0) return;

    Rcpp::NumericVector flagged_prob (vars.size());
    for (auto i = 0u; i < vars.size(); i++)
        flagged_prob(i) = prob_(vars[i]);
    
    int var = Rcpp::sample(vars.size(), 1, false, flagged_prob)(0) - 1;

    int L = 0, U = Xcut_[var].size() - 1;     // for cutpoint draw
    tree.find_region(var, &L, &U);
    int cut = L + Rcpp::sample(U - L, 1)(0);  // avoid min value

    // compute metropolis ratio
    double log_prob   = log(flagged_prob(var)) - log(sum(flagged_prob));  // prob of var
    int    depth      = prop_node->depth();                               // depth of proposed node
    int    n_terminal = tnodes.size();
    int    n_singly   = int();                                            // number of singly node in new tree
    {
        auto parent = prop_node->parent();
        if (parent && parent->is_singly())
            n_singly = tree.singly_size();
        else 
            n_singly = tree.singly_size() + 1;
    }

    int n_unique; // number of possible cut to select
    int nl, nr; double rl, rr;
    get_SS_grow(tree, prop_node, var, cut, nl, nr, rl, rr, n_unique);
    if (nl == 0 || nr == 0 || n_unique < 2) return;

    double ratio = double();
    get_ratio(n_unique - 1, n_terminal, n_singly, depth, log_prob, nl, nr, rl, rr, ratio);

    if (ratio > log(R::runif(0, 1)))
    {
        prop_node->grow(var, cut);
        var_count_[var]++;
    }
}