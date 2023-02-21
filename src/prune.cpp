#include <Rcpp.h>
#include "BART.h"

using namespace std;

void BART::prune(Node& tree)
{
    // collect singly nodes
    vector<Node*> snodes;
    tree.get_singly_nodes(snodes);
    int node_id = Rcpp::sample(snodes.size(), 1)(0) - 1;
    auto prop_node = snodes[node_id];
    int var = prop_node->var();
    int cut = prop_node->cut();

    // draw variable
    vector<int> vars; // variable id can split
    get_vars(prop_node, vars);

    // compute metropolis ratio
    int depth = prop_node->depth(); // depth of proposed node
    double log_prob = double();
    {
        double sum_prob = 0.0;
        for (auto i : vars)
            sum_prob += prob_(i);
        log_prob = log(prob_(var)) - log(sum_prob);
    }

    int L = 0, U = Xcut_[var].size() - 1; // for cutpoint draw
    tree.find_region(var, &L, &U);
    int n_singly   = snodes.size();
    int n_terminal = tree.terminal_size();

    int n_unique; // number of possible cut to select
    int nl, nr; double rl, rr;
    get_SS_prune(tree, prop_node, var, cut, nl, nr, rl, rr, n_unique);

    double ratio = double();
    get_ratio(n_unique, n_terminal - 1, n_singly, depth, log_prob, nl, nr, rl, rr, ratio);
    ratio *= -1;

    if (ratio > log(R::runif(0, 1)))
    {
        prop_node->prune();
        var_count_[var]--;
    }
}