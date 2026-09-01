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

// ---- non-bonded (steric) interaction, all pairs with |i-j| > 1 ----
// isotropic pair energy (gauss / wca) as a function of the squared distance
double nb_iso_energy(const Input& in, double r2);
// purely repulsive Gay-Berne between rods with axes u1, u2 (unit) separated by rvec (from 1 to 2)
double gb_energy(const Input& in, const Vec3& u1, const Vec3& u2, const Vec3& rvec);
// full pair energy of beads a and b of the chain (states, positions and, for helix-helix, tangents)
double nb_pair_energy(const Chain& chain, int a, int b);
// true if the pair energy depends on the tangents (helix-helix Gay-Berne active)
bool nb_anisotropic(const Input& in);
// energy of all pairs involving bead i (|i-j| > 1)
double nb_bead_energy(const Chain& chain, int i);
// energy of all pairs that can change when bead i moves: pairs involving beads i-1, i, i+1 when the
// interaction depends on the tangents, pairs involving bead i otherwise
double nb_local_energy(const Chain& chain, int i);
// energy of all pairs (a, b) with a <= i <= b and b - a > 1: the pairs that change when the chain beyond
// bead i is rotated rigidly about it (tangents of the tail rotate rigidly, the tangent at i does not)
double nb_pivot_energy(const Chain& chain, int i);
double total_nb_energy(const Chain& chain);
