// BartNode.h
#pragma once

#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// Node object for BART
class BartNode
{
    double    leaf_value_  = -1;
    int       var_idx_     = -1;
    int       cut_idx_     = -1;
    BartNode* parent_      = nullptr;
    BartNode* child_left_  = nullptr;
    BartNode* child_right_ = nullptr;

public:
    BartNode() {};
    ~BartNode()
    {
        if (!isTerminal())
        {
            if (child_left_)  delete child_left_;
            if (child_right_) delete child_right_;
        }
    }

    inline double    getLeafValue()  const { return leaf_value_;   }
    inline int       getVarIdx()     const { return var_idx_; }
    inline int       getCutIdx()     const { return cut_idx_; }
    inline BartNode* getParent()     const { return parent_; }
    inline BartNode* getChildLeft()  const { return child_left_; }
    inline BartNode* getChildRight() const { return child_right_; }

    inline void setLeafValue(double leaf_value) { leaf_value_ = leaf_value; }

    inline bool isRoot()     const { return parent_ == nullptr; }
    inline bool isTerminal() const 
    {
        return (child_left_ == nullptr) && (child_right_ == nullptr);
    }
    inline bool isSingly()   const
    {
        if (child_left_ == nullptr || child_right_ == nullptr)
        {
            return false;
        }
        else
        {
            if (child_left_->isTerminal() && child_right_->isTerminal())
                return true;
            else
                return false;
        }
    }

    int getDepth() const
    {
        int depth = 0;
        const BartNode* currentNode = this;
        while (!currentNode->isRoot())
        {
            depth++;
            currentNode = currentNode->parent_;
        }
        return depth;
    }

    void grow(int prop_var_idx, int prop_cut_idx)
    {
        // give two child node
        var_idx_              = prop_var_idx;
        cut_idx_              = prop_cut_idx;
        child_left_           = new BartNode();
        child_right_          = new BartNode();
        child_left_->parent_  = this;
        child_right_->parent_ = this;
    }

    void change(int prop_var_idx, int prop_cut_idx)
    {
        // change split criterion
        var_idx_ = prop_var_idx;
        cut_idx_ = prop_cut_idx;
    }

    void prune()
    {
        // remove child node and make terminal node
        delete this->child_left_;
        delete this->child_right_;
        this->child_left_  = nullptr;
        this->child_right_ = nullptr;
        this->var_idx_     = -1;
        this->cut_idx_     = -1;
    }

    void countSelectedVariables(NumericVector& res) const
    {
        if (isTerminal())
        {
            return;
        }
        else
        {
            res(var_idx_) += 1;
            if (child_left_)  child_left_->countSelectedVariables(res);
            if (child_right_) child_right_->countSelectedVariables(res);
        }
    }

    void getTerminalNodes(vector<BartNode*> &res)
    {
        // get vector of leaf nodes
        if (isTerminal())
        {
            res.push_back(this);
        }
        else
        {
            if (child_left_)  child_left_->getTerminalNodes(res);
            if (child_right_) child_right_->getTerminalNodes(res);
        }
    }

    void getSinglyNodes(vector<BartNode*> &res)
    {
        // get vector of singly nodes
        if (isSingly())
        {
            res.push_back(this);
        }
        else
        {
            if (child_left_)  child_left_->getSinglyNodes(res);
            if (child_right_) child_right_->getSinglyNodes(res);
        }
    }
};

