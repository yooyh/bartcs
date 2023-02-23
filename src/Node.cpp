#include "Node.h"

using namespace std;

Node::~Node()
{
    if (left_)  delete left_;
    if (right_) delete right_;
}
void Node::reset()
{
    if (left_)  delete left_;
    if (right_) delete right_;
    var_    = -1;
    cut_    = -1;
    mu_     = 0.0;
    parent_ = nullptr;
    left_   = nullptr;
    right_  = nullptr;
}

void Node::draw_mu(const int& n, const double& r, const double& sigma2, const double& sigma_mu)
{
    double VAR = 1 / (1/sigma_mu + n / sigma2);
    double MEAN = VAR * r / sigma2;
    mu_ = R::rnorm(MEAN, sqrt(VAR));
}

// node functions
int Node::depth() const
{
    if (!parent_) 
        return 0;
    return (1 + parent_->depth());
}
bool Node::is_terminal() const {return !left_;}
bool Node::is_singly() const
{
    if (left_ == nullptr || right_ == nullptr) 
        return false;
    else if (left_->is_terminal() && right_->is_terminal())
        return true;
    else
        return false;
}
void Node::grow(int var, int cut)
{
    this  ->left_   = new Node();
    this  ->right_  = new Node();
    this  ->var_    = var;
    this  ->cut_    = cut;
    left_ ->parent_ = this;
    right_->parent_ = this;
}
void Node::prune()
{
    delete this->left_;
    delete this->right_;
    this->left_  = nullptr;
    this->right_ = nullptr;
    this->var_   = -1;
    this->cut_   = -1;
}
void Node::change(int var, int cut)
{
    this->var_ = var;
    this->cut_ = cut;
}

// tree functions
int Node::terminal_size() const
{
    if (is_terminal()) return 1;
    return left_->terminal_size() + right_->terminal_size();
}
int Node::singly_size() const
{
    if (is_terminal()) return 0;
    if (left_->left_ || right_->left_)
        return left_->singly_size() + right_->singly_size();
    return 1;
}
void Node::get_terminal_nodes(vector<Node*>& tnodes)
{
    if (is_terminal())
    {
        tnodes.push_back(this);
        return;
    }
    left_ ->get_terminal_nodes(tnodes);
    right_->get_terminal_nodes(tnodes);
}
void Node::get_singly_nodes(vector<Node*>& snodes)
{
    if (is_singly()) 
    {
        snodes.push_back(this);
        return;
    }
    if (!is_terminal())
    {
        left_ ->get_singly_nodes(snodes);
        right_->get_singly_nodes(snodes);
    }   
}
void Node::find_region(int var, int* L, int* U) const
{
    // modified tree::rg() from BART package
    if (!parent_) return;
    if (parent_->var_ == var)        // does my parent use var?
    {
        if (this == parent_->left_)  // am I left child?
        {
            if (parent_->cut_ <= *U)
                *U = (parent_->cut_) - 1;
            parent_->find_region(var, L, U);
        }
        else
        {
            if (parent_->cut_ >= *L)
                *L = (parent_->cut_) + 1;
            parent_->find_region(var, L, U);
        }
    }
    else
    {
        parent_->find_region(var, L, U);
    }
}

Node* Node::assigned_node(const vector<vector<double>>& Xcut, const vector<double>& x)
{
    // modified tree::bn() from BART package
    if (is_terminal()) return this;
    if (x[var_] < Xcut[var_][cut_])
        return left_ ->assigned_node(Xcut, x);
    else
        return right_->assigned_node(Xcut, x);
}
const Node* Node::assigned_node(const vector<vector<double>>& Xcut, const vector<double>& x) const
{
    // modified tree::bn() from BART package
    if (is_terminal()) return this;
    if (x[var_] < Xcut[var_][cut_])
        return left_->assigned_node(Xcut, x);
    else
        return right_->assigned_node(Xcut, x);
}
