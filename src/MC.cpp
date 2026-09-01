#include "MC.h"

#include "update_position.h"
#include "update_state.h"
#include "update_hinge.h"
#include "update_pivot.h"

MC::MC(const Input& in, Chain& chain)
    : in_(in), chain_(chain), rng_(in.seed) {}

// One sweep = N attempted state moves + N attempted position moves, each on a random residue.
void MC::sweep() {
    const int N = chain_.N();
    std::uniform_int_distribution<int> pick(0, N - 1);
    for (int k = 0; k < N; ++k) {
        ++try_state_;
        if (try_state_move(chain_, pick(rng_), in_, rng_)) ++acc_state_;

        ++try_pos_;
        if (try_position_move(chain_, pick(rng_), in_, rng_)) ++acc_pos_;
    }
    for (int k = 0; k < in_.n_hinge; ++k) {
        ++try_hinge_;
        if (try_hinge_move(chain_, pick(rng_), in_, rng_)) ++acc_hinge_;
    }
    if (N > 2) {
        std::uniform_int_distribution<int> pick_inner(1, N - 2);
        for (int k = 0; k < in_.n_pivot; ++k) {
            ++try_pivot_;
            if (try_pivot_move(chain_, pick_inner(rng_), in_, rng_)) ++acc_pivot_;
        }
    }
}

void MC::reset_acceptance() {
    try_state_ = acc_state_ = 0;
    try_pos_   = acc_pos_   = 0;
    try_hinge_ = acc_hinge_ = 0;
    try_pivot_ = acc_pivot_ = 0;
}
