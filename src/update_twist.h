#pragma once
#include <random>
#include "chain.h"
#include "input.h"

// Twist move: the whole same-handed run containing helical residue i is rotated about its own
// axis by one angle delta in [-twist_step, twist_step]; Metropolis on the non-bonded pairs touching
// the run.  Only meaningful when the pair potential depends on the registries (nb_hh = db).
bool try_twist_move(Chain& chain, int i, const Input& in, std::mt19937_64& rng);
