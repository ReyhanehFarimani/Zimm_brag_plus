#pragma once
#include <random>
#include "chain.h"

// Combined spin-flip + hinge-resampling move (for chains without non-bonded interactions).
// The spin of residue i is changed as in the state move; every hinge in {i-1, i, i+1} whose bend parameters
// change gets a new bond direction drawn from its NEW hinge Boltzmann distribution p(theta) ~ sin(theta) e^{-U},
// and the chain downstream of that bond is rotated rigidly. Because the angular proposal is the exact hinge
// distribution, the acceptance reduces to  exp(-dE_spin - dE_bond) * prod (Z_hinge,new / Z_hinge,old),
// i.e. the bending energy barrier disappears (only the free-energy cost remains).
bool try_hinge_move(Chain& chain, int i, const Input& in, std::mt19937_64& rng);
