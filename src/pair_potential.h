#pragma once
#include "chain.h"

// Harmonic bond between residues i and i+1:
//   U = 1/2 k (r - r0)^2, with (r0, k) taken from the pair class of the bond.
double bond_energy(const Chain& chain, int i);            // i in [0, N-2]
double bond_energy_of(const Chain& chain, int i, double r); // same bond, given length r

// Sum over all N-1 bonds.
double total_bond_energy(const Chain& chain);

// Harmonic bending at residue i. theta is the valence angle at the bead (straight chain: theta = pi):
//   U = 1/2 kappa (theta - theta0)^2, with (kappa, theta0) from the triple class of (i-1, i, i+1).
double bend_energy_angle(const BendParams& p, double theta);
double bend_energy(const Chain& chain, int i);                               // i in [1, N-2]
double bend_energy_of(const Chain& chain, int i, const Vec3& a, const Vec3& b); // same bend, given bond vectors a = bond(i-1), b = bond(i)

// Sum over all N-2 bends.
double total_bend_energy(const Chain& chain);

// ln of the angular partition function of one hinge with parameters p:  ln[ 2pi int sin(theta) e^{-U} dtheta ]
// (4 pi for a free hinge). Used by the hinge move.
double hinge_log_weight(const BendParams& p);

// Nearest-neighbour state energy of the pair (i, i+1): -J0 (HH), +J1 (CH), 0 (CC), +J2 (RL).
double state_energy(const Chain& chain, int i);             // i in [0, N-2]
double total_state_energy(const Chain& chain);
