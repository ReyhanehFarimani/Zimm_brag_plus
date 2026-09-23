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
// On-site energy of residue i: E_helix if it is R or L, 0 if coil.
double site_energy(const Chain& chain, int i);
// all pair terms plus all site terms
double total_state_energy(const Chain& chain);

// ---- non-bonded (steric) interaction, all pairs with |i-j| > 1 ----
// isotropic pair energy (gauss / wca) as a function of the squared distance
double nb_iso_energy(const Input& in, double r2);
// purely repulsive Gay-Berne between rods with axes u1, u2 (unit) separated by rvec (from 1 to 2)
double gb_energy(const Input& in, const Vec3& u1, const Vec3& u2, const Vec3& rvec);
// measured helix-helix radial pair potential (nb_hh = fit): erfc radial sectors x quadratic
// angular maps, handedness-resolved; hand1/hand2 = spin (+1 R, -1 L) of the two residues
double helixfit_energy(const Input& in, const Vec3& u1, const Vec3& u2, const Vec3& rvec,
                       int hand1, int hand2);
// tabulated helix-helix potential (nb_hh = db): nearest node of the six-argument table in
// (r, betaA, betaB, psi, alphaA, alphaB); m1, m2 = the registry directions of the two residues
double db_energy(const Input& in, const Vec3& u1, const Vec3& u2, const Vec3& rvec,
                 int hand1, int hand2, const Vec3& m1, const Vec3& m2);
// full pair energy of beads a and b of the chain (states, positions and, for helix-helix, tangents)
double nb_pair_energy(const Chain& chain, int a, int b);
// energy of every pair (k, j) with k in [lo, hi], |k-j| > 1; pairs inside the range counted once
// (the pairs that change when the states or registries of [lo, hi] change; geometry unchanged)
double nb_range_energy(const Chain& chain, int lo, int hi, bool use_list = true);
// true if the pair energy depends on the tangents (helix-helix Gay-Berne or fit active)
bool nb_anisotropic(const Input& in);
// true if the pair energy depends on the helix handedness (nb_hh = fit): state moves that
// change R <-> L, including domain flips, must then recompute the non-bonded energy
bool nb_chiral(const Input& in);
// largest distance at which any pair energy can be non-zero (iso cutoff, Gay-Berne reach, fit cutoff)
double nb_cutoff_max(const Input& in);

// ---- Verlet neighbour list (Chain::nl): changes the cost only, never the energies ----
// true if the list is switched on; rebuilds it first when it is dirty
bool nl_ready(const Chain& chain);
// true if the list stays valid with bead i at position p (within skin/2 of its reference position)
bool nl_covers(const Chain& chain, int i, const Vec3& p);
inline void nl_invalidate(const Chain& chain) { chain.nl.dirty = true; }
// consistency check: every pair closer than the largest cutoff is in the list (always true if off)
bool nl_verify(const Chain& chain);

// energy of all pairs involving bead i (|i-j| > 1)
double nb_bead_energy(const Chain& chain, int i);
// energy of all pairs that can change when bead i moves: pairs involving beads i-1, i, i+1 when the
// interaction depends on the tangents, pairs involving bead i otherwise.
// use_list = false forces the all-pairs loop (for a trial position the list does not cover)
double nb_local_energy(const Chain& chain, int i, bool use_list = true);
// energy of all pairs (a, b) with a <= i <= b and b - a > 1: the pairs that change when the chain beyond
// bead i is rotated rigidly about it (tangents of the tail rotate rigidly, the tangent at i does not)
double nb_pivot_energy(const Chain& chain, int i);
double total_nb_energy(const Chain& chain);
// star core (n_arms > 0): harmonic wall at core_radius + 1/2 for every residue and the graft tether of the
// first residue of each arm; zero for a linear chain
double core_energy_at(const Chain& chain, int i, const Vec3& p);
double core_energy(const Chain& chain, int i);
double core_range_energy(const Chain& chain, int lo, int hi);
double total_core_energy(const Chain& chain);
// debug: the energy change the last move computed for its acceptance test (set by every move)
extern double g_last_dE;
