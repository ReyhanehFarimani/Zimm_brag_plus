#include "MC.h"

#include "update_position.h"

MC::MC(const Input& in, Chain& chain)
    : in_(in), chain_(chain), rng_(in.seed) {}

// One sweep = N attempted state moves + N attempted position moves, each on a random residue.
void MC::sweep() {
    const int N = chain_.N();
    std::uniform_int_distribution<int> pick(0, N - 1);
    for (int k = 0; k < N; ++k) {
        ++try_state_;   // TODO: try_state_move()

        ++try_pos_;
        if (try_position_move(chain_, pick(rng_), in_, rng_)) ++acc_pos_;
    }
}

void MC::reset_acceptance() {
    try_state_ = acc_state_ = 0;
    try_pos_   = acc_pos_   = 0;
}
