#pragma once
#include <random>
#include <vector>
#include "chain.h"

// Registry of a helical residue: a unit direction m_i perpendicular to its axis (the local
// tangent), fixed to the rod.  A contiguous run of same-handed helical residues shares one
// registry: it is created by copying (transported) from the neighbour, twisted as a whole by the
// twist move, and split/merged with the rules below, whose proposal factors cancel under a uniform
// measure for the registries so that every move keeps the plain Metropolis rule.
//   on a split (a residue of a run turns coil or flips hand): the RIGHT part draws a fresh uniform twist
//   on a merge (a residue joins two same-handed runs): the RIGHT run is twisted to match the junction
//   a residue joining one run copies its neighbour's registry (left first); joining none: uniform.

// minimal rotation taking u_old to u_new, applied to m (keeps m perpendicular to the axis)
Vec3 registry_transport(const Vec3& m, const Vec3& u_old, const Vec3& u_new);
// a unit vector perpendicular to u at azimuth phi from an arbitrary but fixed reference
Vec3 registry_from_angle(const Vec3& u, double phi);
// signed angle (right-handed about u) from a to b, both perpendicular to u
double registry_angle(const Vec3& u, const Vec3& a, const Vec3& b);

struct RegistryUndo {
    int lo = 0, hi = -1;
    std::vector<Vec3> m;
};
// index range [i, hi] whose registries change when residue i goes from s_old to s_new (states as
// they are NOW, i.e. before the change); hi = i when only residue i is affected
int registry_change_hi(const Chain& chain, int i, State s_old, State s_new);
// apply the registry rules for the state change at i (call AFTER chain.state[i] = s_new);
// what is overwritten goes into `undo`
void registry_on_state_change(Chain& chain, int i, State s_old, State s_new,
                              std::mt19937_64& rng, RegistryUndo& undo);
void registry_undo(Chain& chain, const RegistryUndo& undo);
// initial registries: one uniform twist per run, transported along it
void registry_init(Chain& chain, std::mt19937_64& rng);
// One registry per run: the LEFTMOST residue's m is the anchor, the others are derived by
// parallel transport along the run's current geometry.  After a geometry change at residue k
// (its tangent moved), rederive [k, hi of its run]: if k is the anchor its m must already have
// been transported by the caller.  Returns hi (the last residue whose registry may have changed).
int registry_rederive(Chain& chain, int k);
// run bounds of a helical residue
void registry_run(const Chain& chain, int i, int& lo, int& hi);
// consistency check (debug): every helical m is perpendicular to its tangent and every run's
// registries agree under transport; returns the number of violations
int registry_check(const Chain& chain, double tol = 1e-6);
