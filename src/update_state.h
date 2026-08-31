#pragma once
#include <random>
#include "chain.h"

// Metropolis spin move on residue i: the spin (L=-1, C=0, R=+1) is raised or lowered by one
// with equal probability. C -> R or L; R -> C; L -> C. A step outside [-1, 1] is rejected
// (keeps the proposal symmetric). With n_states = 2 the spin is flipped R <-> L directly. The energy change includes the state energies, the bond
// energies of the two bonds touching i and the bends at i-1, i, i+1, since bond and bend
// parameters depend on the pair / triple classes.
bool try_state_move(Chain& chain, int i, const Input& in, std::mt19937_64& rng);
