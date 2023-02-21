#include "BART.h"

BART::BART(
    const vector<vector<double>>& X, const vector<vector<double>>& Xcut, const int N, const int P,
    const int n_tree, const Rcpp::NumericVector& step_prob,
    const double alpha, const double beta, Rcpp::NumericVector& var_prob, const bool parallel
) : X_(X), Xcut_(Xcut), N(N), P(P), n_tree_(n_tree), step_prob_(step_prob), alpha_(alpha), beta_(beta), prob_(var_prob), parallel_(parallel) {
    tree_.resize(n_tree_);
    for (int t = 0; t < n_tree_; t++)
        tree_[t] = Node();

    fitted_.resize(N);
    fit_tmp_.resize(N);
    var_count_.resize(P);
    residual_.resize(N);
}

double BART::get_sigma_mu(const vector<double>& Y, const int n_tree) const
{
    auto min_max = minmax_element(Y.begin(), Y.end());
    return max(pow(*min_max.first  / (-2*sqrt(n_tree_)), 2), pow(*min_max.second / ( 2*sqrt(n_tree_)), 2));
}

void BART::init(const vector<double>& Y, double sigma2)
{
    sigma2_   = sigma2;
    sigma_mu_ = get_sigma_mu(Y, n_tree_);

    // reset tree
    for (auto& tree : tree_)
        tree.reset();

    // reset value to 0
    #ifdef _OPENMP
        #pragma omp parallel for if (parallel_)
    #endif
    for (int i = 0; i < N; i++)
    {
        fitted_[i] = 0.0;
        fit_tmp_[i] = 0.0;
        residual_[i] = Y[i];
    }
    for (auto& n : var_count_) n = 0;
}

// draw and step function
void BART::draw(const vector<double>& Y)
{
    for (auto& tree : tree_)
    {
        fit(tree, fit_tmp_);
        #ifdef _OPENMP
            #pragma omp parallel for if (parallel_)
        #endif
        for (int i = 0; i < N; i++)
        {
            fitted_[i]  -= fit_tmp_[i];
            residual_[i] = Y[i] - fitted_[i];
        }
        step(tree);
        draw_mu(tree);
        fit(tree, fit_tmp_);
        #ifdef _OPENMP
            #pragma omp parallel for if (parallel_)
        #endif
        for (int i = 0; i < N; i++)
            fitted_[i] += fit_tmp_[i];
    }
}

void BART::step(Node& tree)
{
    if (tree.is_terminal())
    {
        grow(tree);
        return;
    }
    int step = Rcpp::sample(3, 1, false, step_prob_)(0);
    switch (step)
    {
        case 1: // grow
            grow(tree);
            return;
        
        case 2: // prune:
            prune(tree);
            return;
        
        case 3: // change
            change(tree);
            return;
    }       
}

// util functions
void BART::get_SS(
    const Node& tree, const Node* prop_node, const int var, const int cut,
    int& nl, int& nr, double& rl, double& rr
) const {
    nl = 0; rl = 0.0;
    nr = 0; rr = 0.0;
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+ : nl, nr, rl, rr) if (parallel_)
    #endif
    for (int i = 0; i < N; i++)
    {
        auto node = tree.assigned_node(Xcut_, X_[i]);
        if (node != prop_node) continue;
        if (X_[i][var] < Xcut_[var][cut])
        {
            nl++;
            rl += residual_[i];
        }
        else
        {
            nr++;
            rr += residual_[i];
        }
    }
}
void BART::get_SS_grow(
    const Node& tree, const Node* prop_node, const int var, const int cut,
    int& nl, int& nr, double& rl, double& rr, int& n_unique
) const {
    nl = 0; rl = 0.0;
    nr = 0; rr = 0.0;
    unordered_set<double> set;
    set.reserve(N);
    if (parallel_)
    {
        #ifdef _OPENMP
            #pragma omp parallel if (parallel_)
        #endif
        {
            unordered_set<double> local_set;
            #ifdef _OPENMP
                #pragma omp for reduction(+ : nl, nr, rl, rr)
            #endif
            for (int i = 0; i < N; i++)
            {
                auto node = tree.assigned_node(Xcut_, X_[i]);
                if (node != prop_node) continue;
                local_set.insert(X_[i][var]);
                if (X_[i][var] < Xcut_[var][cut])
                {
                    nl++;
                    rl += residual_[i];
                }
                else
                {
                    nr++;
                    rr += residual_[i];
                }
            }
            #ifdef _OPENMP
                #pragma omp critical
            #endif
            {
                set.insert(local_set.begin(), local_set.end());
            }
        }
    }
    else
    {
        for (int i = 0; i < N; i++)
        {
            auto node = tree.assigned_node(Xcut_, X_[i]);
            if (node != prop_node) continue;
            set.insert(X_[i][var]);
            if (X_[i][var] < Xcut_[var][cut])
            {
                nl++;
                rl += residual_[i];
            }
            else
            {
                nr++;
                rr += residual_[i];
            }
        }
    }
    n_unique = set.size();
}
void BART::get_SS_prune(
    const Node& tree, const Node* prop_node, const int var, const int cut,
    int& nl, int& nr, double& rl, double& rr, int& n_unique
) const {
    nl = 0; rl = 0.0;
    nr = 0; rr = 0.0;
    unordered_set<double> set;
    set.reserve(N);
    if (parallel_)
    {
        #ifdef _OPENMP
            #pragma omp parallel if (parallel_)
        #endif
        {
            unordered_set<double> local_set;
            #ifdef _OPENMP
                #pragma omp for reduction(+ : nl, nr, rl, rr)
            #endif
            for (int i = 0; i < N; i++)
            {
                auto node = tree.assigned_node(Xcut_, X_[i]);
                if (node->parent() != prop_node) continue;
                local_set.insert(X_[i][var]);
                if (X_[i][var] < Xcut_[var][cut])
                {
                    nl++;
                    rl += residual_[i];
                }
                else
                {
                    nr++;
                    rr += residual_[i];
                }
            }
            #ifdef _OPENMP
                #pragma omp critical
            #endif
            {
                set.insert(local_set.begin(), local_set.end());
            }
        }
    }
    else
    {
        for (int i = 0; i < N; i++)
        {
            auto node = tree.assigned_node(Xcut_, X_[i]);
            if (node->parent() != prop_node) continue;
            set.insert(X_[i][var]);
            if (X_[i][var] < Xcut_[var][cut])
            {
                nl++;
                rl += residual_[i];
            }
            else
            {
                nr++;
                rr += residual_[i];
            }
        }
    }
    n_unique = set.size();
}
void BART::get_SS_change(
    const Node& tree, const Node* prop_node,
    const int cvar, const int ccut, int& cnl, int& cnr, double& crl, double& crr,
    const int pvar, const int pcut, int& pnl, int& pnr, double& prl, double& prr
) const {
    cnl = 0; cnr = 0; crl = 0.0; crr = 0.0;
    pnl = 0; pnr = 0; prl = 0.0; prr = 0.0;
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+ : cnl, cnr, crl, crr, pnl, pnr, prl, prr) if (parallel_)
    #endif
    for (int i = 0; i < N; i++)
    {
        auto node = tree.assigned_node(Xcut_, X_[i]);
        if (node->parent() != prop_node) continue;
        if (X_[i][cvar] < Xcut_[cvar][ccut])
        {
            cnl++;
            crl += residual_[i];
        }
        else
        {
            cnr++;
            crr += residual_[i];
        }
        if (X_[i][pvar] < Xcut_[pvar][pcut])
        {
            pnl++;
            prl += residual_[i];
        }
        else
        {
            pnr++;
            prr += residual_[i];
        }
    }
}

// get ratio for grow
void BART::get_ratio(
    const int& n_unique, const int& n_terminal, const int& n_singly, 
    const int depth, const double& log_prob,
    const int& nl, const int& nr, const double& rl, const double& rr,
    double& ratio
) const {
    double transition = log(step_prob_(1)) + log(n_terminal) - log_prob
        + log(n_unique) - log(step_prob_(0)) - log(n_singly);
    
    double likelihood = 
          0.5 * log(sigma2_) 
        + 0.5 * log(sigma2_ + sigma_mu_ * (nl + nr))
        - 0.5 * log(sigma2_ + sigma_mu_ * nl)
        - 0.5 * log(sigma2_ + sigma_mu_ * nr)
        + (sigma_mu_ / (2 * sigma2_)
        *(pow(rl,      2) / (sigma2_ + sigma_mu_ *  nl)
        + pow(rr,      2) / (sigma2_ + sigma_mu_ *  nr)
        - pow(rl + rr, 2) / (sigma2_ + sigma_mu_ * (nl + nr))));

    double structure = 
          log(alpha_) + 2*log(1 - alpha_ / pow(2 + depth, beta_))
        - log(pow(1 + depth, beta_) - alpha_)
        + log_prob - log(n_unique);
    
    ratio = transition + likelihood + structure;
}

// draw mu from each terminal nodes
void BART::draw_mu(Node& tree)
{
    vector<Node*> tnodes;
    tree.get_terminal_nodes(tnodes);

    vector<int>    NN (tnodes.size(), 0);
    vector<double> RR (tnodes.size(), 0.0);
    unordered_map<const Node*, int> node2id;
    for (int i = 0; i < tnodes.size(); i++)
        node2id[tnodes[i]] = i;
    
    for (int i = 0; i < N; i++)
    {
        auto node = tree.assigned_node(Xcut_, X_[i]);
        int id = node2id[node];
        NN[id]++;
        RR[id] += residual_[i];
    }

    for (int i = 0; i < tnodes.size(); i++)
        tnodes[i]->draw_mu(NN[i], RR[i], sigma2_, sigma_mu_);
}

void BART::draw_sigma2(
    const Rcpp::Function& rinvgamma,
    const vector<double>& Y,
    const double          nu,
    const double          lambda
) {
    double sum = 0.0;
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+ : sum) if (parallel_)
    #endif
    for (int i = 0; i < N; i++)
        sum += pow(Y[i] - fitted_[i], 2);
    
    double shape = nu / 2 + N / 2;
    double scale = nu * lambda / 2 + sum / 2;
    Rcpp::NumericVector tmp = rinvgamma(1, shape, scale);
    sigma2_ = tmp(0);
}

//fit tree and save at res
void BART::fit(
    const Node& tree,
    vector<double>& res
) const {
    #ifdef _OPENMP
        #pragma omp parallel for if (parallel_)
    #endif
    for (int i = 0; i < N; i++)
    {
        auto node = tree.assigned_node(Xcut_, X_[i]);
        res[i] = node->mu();
    }
}

//find variables n can split on, put their indices in vars
void BART::get_vars(const Node* node, vector<int>& vars) const
{
    vars.clear();
    int L, U;
    for (int v = 0; v < P; v++) {//try each variable
        L = 0; 
        U = Xcut_[v].size() - 1;
        node->find_region(v, &L, &U);
        if (U >= L) 
            vars.push_back(v);
    }
}

void BART::update_Z(vector<double>& Z, const Rcpp::NumericVector& TRT, const bool binary_trt)
{
    if (binary_trt)
    {
        #ifdef _OPENMP
            #pragma omp parallel for if (parallel_)
        #endif
        for (int i = 0; i < N; i++)
        {
            double Ystar = R::rnorm(fitted_[i], 1);
            Z[i] = TRT[i] * max(Ystar, 0.0) + (1 - TRT[i]) * min(Ystar, 0.0);
        }
    }
    else
    {
        #ifdef _OPENMP
            #pragma omp parallel for if (parallel_)
        #endif
        for (int i = 0; i < N; i++)
            Z[i] = R::rnorm(fitted_[i], sigma2_);
    }
}