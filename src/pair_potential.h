#pragma once
#include "chain.h"

// Harmonic bond between residues i and i+1:
//   U = 1/2 k (r - r0)^2, with (r0, k) taken from the pair class of the bond.
double bond_energy(const Chain& chain, int i);            // i in [0, N-2]
double bond_energy_of(const Chain& chain, int i, double r); // same bond, given length r

// Sum over all N-1 bonds.
double total_bond_energy(const Chain& chain);
