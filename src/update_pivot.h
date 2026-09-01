#pragma once
#include <random>
#include "chain.h"

// Pivot move: rotate every bead beyond bead i rigidly about a random axis through bead i by a random angle
// in [-max_rot, max_rot]. Changes only the bend at i and the non-bonded pairs across the cut, which are
// evaluated in full (valid with any interaction). Returns true if accepted.
bool try_pivot_move(Chain& chain, int i, const Input& in, std::mt19937_64& rng);
