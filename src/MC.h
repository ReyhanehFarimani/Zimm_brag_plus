#pragma once
#include <random>
#include "input.h"
#include "chain.h"

// Owns the RNG and drives one MC sweep = N attempts of each move type
// (state flip, position). Orientation moves come later with sterics.
class MC {
public:
    MC(const Input& in, Chain& chain);
    void sweep();
    void reset_acceptance();

    double acc_state() const { return try_state_ ? double(acc_state_) / try_state_ : 0.0; }
    double acc_pos()   const { return try_pos_   ? double(acc_pos_)   / try_pos_   : 0.0; }
    double acc_hinge() const { return try_hinge_ ? double(acc_hinge_) / try_hinge_ : 0.0; }
    double acc_pivot() const { return try_pivot_ ? double(acc_pivot_) / try_pivot_ : 0.0; }

private:
    const Input& in_;
    Chain&       chain_;
    std::mt19937_64 rng_;

    long try_state_ = 0, acc_state_ = 0;
    long try_pos_   = 0, acc_pos_   = 0;
    long try_hinge_ = 0, acc_hinge_ = 0;
    long try_pivot_ = 0, acc_pivot_ = 0;
};
