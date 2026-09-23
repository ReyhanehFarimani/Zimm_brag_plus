#pragma once
#include <random>
#include <vector>
#include "chain.h"

// Registry of a residue: a unit direction m_i perpendicular to its axis (the local tangent), fixed
// to the rod as a material direction.  EVERY residue carries one (a coil residue's is a dummy that
// enters no energy); a helical residue's registry enters the tabulated pair potential (nb_hh = db)
// and the TWIST term between sequential same-handed residues,
//
//     E_twist(i, i+1) = 1/2 kappa_t (dalpha_i - s dalpha_0)^2,   s = +1 (R.R), -1 (L.L),
//
// dalpha_i = signed rotation of the registry from residue i to i+1 about their common axis
// n = t_i + t_{i+1} (right-handed about the chain direction), wrapped to (-pi, pi].  The measured
// values (grant repo, prelim/cg/twist_potentials.py, K* = 20 graded runs; 2026-09-23):
//   theta0 = 100: dalpha_0 = -135.3 deg (ideal screw, 7 x turn per bead), sigma = 42.5 deg
//   theta0 =  45: dalpha_0 =  +25.3 deg,                                   sigma = 22.3 deg
// kappa_t = kT / sigma^2 reproduces the measured variance; the fitted well positions are -134.6 and
// +21.3 deg.  Opposite-hand neighbours (R.L) and helix-coil neighbours have no twist term.
//
// Registries move with the geometry: a rigid rotation of a segment (pivot, hinge) rotates its
// registries with it; a residue whose tangent changes has its registry carried by the minimal
// rotation old tangent -> new tangent (a bijection of the circle that is its own inverse under the
// reverse move, so every geometric move keeps the plain Metropolis rule).  A state change at i
// redraws m_i uniformly in both directions (symmetric proposal); a whole-run handedness flip
// negates the run's junction twists (an involution) so the twist energy is invariant and the flip
// is judged by the pair potential alone.

// minimal rotation taking u_old to u_new, applied to m (keeps m perpendicular to the axis)
Vec3 registry_transport(const Vec3& m, const Vec3& u_old, const Vec3& u_new);
// a unit vector perpendicular to u at azimuth phi from an arbitrary but fixed reference
Vec3 registry_from_angle(const Vec3& u, double phi);
// signed angle (right-handed about u) from a to b, both perpendicular to u
double registry_angle(const Vec3& u, const Vec3& a, const Vec3& b);

// registry twist between residues i and i+1 (rad, in (-pi, pi]); geometry only, no state test
double registry_twist(const Chain& chain, int i);
// twist energy of the junction (i, i+1): 0 unless both are helical and same-handed
double twist_energy(const Chain& chain, int i);
// twist energy of every junction touching residues [lo, hi], i.e. junctions lo-1 .. hi
double twist_range_energy(const Chain& chain, int lo, int hi);
double total_twist_energy(const Chain& chain);

// initial registries: uniform, every residue
void registry_init(Chain& chain, std::mt19937_64& rng);
// redraw m_i uniformly (state change at i); returns the old value for the undo
Vec3 registry_resample(Chain& chain, int i, std::mt19937_64& rng);
// negate the junction twists of residues [a, b] (whole-run handedness flip), keeping m_a;
// an involution: applying it twice restores every registry
void registry_negate_twists(Chain& chain, int a, int b);
// consistency check (debug): every registry is a unit vector perpendicular to its tangent;
// returns the number of violations
int registry_check(const Chain& chain, double tol = 1e-6);
