#pragma once
#include <random>
#include "chain.h"

// Metropolis single-bead displacement.
// Bead i is moved by a uniform random vector in [-max_disp, max_disp]^3; the energy
// change comes from the (at most two) bonds attached to i. Returns true if accepted.
bool try_position_move(Chain& chain, int i, const Input& in, std::mt19937_64& rng);
