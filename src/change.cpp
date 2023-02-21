#include <Rcpp.h>
#include "BART.h"

using namespace std;

void BART::change(Node& tree)
{
    // collect singly nodes
    vector<Node*> snodes;
    tree.get_singly_nodes(snodes);
    int node_id = Rcpp::sample(snodes.size(), 1)(0) - 1;
    auto prop_node = snodes[node_id];

    // draw variable
    vector<int> vars; // variable id can split
    get_vars(prop_node, vars);

    Rcpp::NumericVector flagged_prob (vars.size());
    for (int i = 0; i < vars.size(); i++)
        flagged_prob(i) = prob_(vars[i]);
    
    int var = Rcpp::sample(vars.size(), 1, false, flagged_prob)(0) - 1;

    int L = 0, U = Xcut_[var].size() - 1; // for cutpoint draw
    tree.find_region(var, &L, &U);
    int cut = L + Rcpp::sample(U - L, 1)(0); // avoid min value

    // compute metropolis ratio
    int cnl, cnr; double crl, crr; // current likelihood
    int pnl, pnr; double prl, prr; // proposed likelihood
    get_SS_change(tree, prop_node, 
                  prop_node->var(), prop_node->cut(), cnl, cnr, crl, crr, 
                  var, cut, pnl, pnr, prl, prr);
    if (pnl == 0 || pnr == 0) return;

    double ratio = 
          0.5 * log(sigma2_ / sigma_mu_ + cnl)
        + 0.5 * log(sigma2_ / sigma_mu_ + cnr)
        - 0.5 * log(sigma2_ / sigma_mu_ + pnl)
        - 0.5 * log(sigma2_ / sigma_mu_ + pnr)
        + (0.5 / sigma2_
        *(pow(prl, 2) / (sigma2_ / sigma_mu_ + pnl)
        + pow(prr, 2) / (sigma2_ / sigma_mu_ + pnr)
        - pow(crl, 2) / (sigma2_ / sigma_mu_ + cnl)
        - pow(crr, 2) / (sigma2_ / sigma_mu_ + cnr)));
    
    if (ratio > log(R::runif(0, 1)))
    {
        var_count_[prop_node->var()]--;
        prop_node->change(var, cut);
        var_count_[var]++;
    }
}