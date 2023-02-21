#pragma once

#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

class ProgressBar
{
    const int chain_idx_; 
    const int chain_pad_;
    const int iter_pad_;
    const int total_iter_;
    const int bar_width_;

public:
    ProgressBar(
        const int chain_idx, 
        const int num_chain, 
        const int total_iter,
        const int bar_width
    ) :
        chain_idx_  (chain_idx),
        chain_pad_  (1 + (int) log10(num_chain)),
        iter_pad_   (1 + log10(total_iter)),
        total_iter_ (total_iter),
        bar_width_  (bar_width) {};

    inline void print(const int current_iter) const
    {
        if (current_iter == total_iter_)
        {
            Rcout << "\rChain " << setw(chain_pad_) << chain_idx_ << ": 100% [";
            for (int i = 0; i < bar_width_; i++)
                Rcout << "=";
            Rcout << "] " << current_iter << "/" << total_iter_ << "\n";
        }
        else
        {
            int progress    = 100        * current_iter / total_iter_;
            int current_pos = bar_width_ * current_iter / total_iter_;

            Rcout << "\r";
            Rcout << "Chain " << setw(chain_pad_) << chain_idx_ << ": ";
            Rcout << setw(3) << progress << "% [";
            for (int i = 0; i < bar_width_; i++)
            {
                if (i < current_pos) 
                    Rcout << "=";
                else if (i == current_pos) 
                    Rcout << ">";
                else 
                    Rcout << " ";
            }
            Rcout << "] "  << setw(iter_pad_);
            Rcout << current_iter << "/" << total_iter_ ;
        }
    }
};
