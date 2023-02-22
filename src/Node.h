#pragma once

#include <Rcpp.h>

using namespace std;

class Node 
{
private:
    int    var_;
    int    cut_;
    double mu_;
    Node*  parent_;
    Node*  left_;
    Node*  right_;

public:
    Node() : var_(-1), cut_(-1), mu_(0.0), parent_(nullptr), left_(nullptr), right_(nullptr) {}
    ~Node();
    void reset();
    
    // setter
    void draw_mu(const int& n, const double& r, const double& sigma2, const double& sigma_mu);

    // getter
    const int    var() const {return var_;}
    const int    cut() const {return cut_;}
    const double mu()     const {return mu_;}
    const Node*  parent() const {return parent_;}
    const Node*  left()   const {return left_;}
    const Node*  right()  const {return right_;}

    // node functions
    int  depth() const; // depth of node
    bool is_terminal() const;
    bool is_singly() const;
    void grow(int var, int cut);
    void prune();
    void change(int var, int cut);

    // tree functions
    int  terminal_size() const;  // number of terminal nodes
    int  singly_size() const;    // number of singly nodes
    void get_terminal_nodes(vector<Node*>& tnodes);
    void get_singly_nodes(vector<Node*>& snodes);
    void find_region(int var, int* L, int* U) const;
    
    Node* assigned_node(const vector<vector<double>>& Xcut, const vector<double>& x);
    const Node* assigned_node(const vector<vector<double>>& Xcut, const vector<double>& x) const;
};
