#pragma once
#include <random>
#include "chain.h"
#include "input.h"

// Twist move: the registry of helical residue i is rotated about its own axis by delta in
// [-twist_step, twist_step] degrees; Metropolis on the non-bonded pairs of i and the twist term
// of its two junctions.  Only meaningful when the pair potential depends on the registries
// (nb_hh = db).
bool try_twist_move(Chain& chain, int i, const Input& in, std::mt19937_64& rng);
