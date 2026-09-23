#pragma once
#include <cmath>
#include <vector>
#include "input.h"

struct Vec3 {
    double x, y, z;
    Vec3() : x(0.0), y(0.0), z(0.0) {}
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
};
inline Vec3   operator+(Vec3 a, Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3   operator-(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3   operator*(double s, Vec3 a) { return {s * a.x, s * a.y, s * a.z}; }
inline double dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline double norm2(Vec3 a) { return dot(a, a); }
inline double norm(Vec3 a) { return std::sqrt(norm2(a)); }
inline Vec3   cross(Vec3 a, Vec3 b) { return Vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x); }
// rotate v about the unit axis k by angle a (Rodrigues)
inline Vec3 rotate(Vec3 v, Vec3 k, double a) {
    const double c = std::cos(a), s = std::sin(a);
    return c * v + s * cross(k, v) + (dot(k, v) * (1.0 - c)) * k;
}

// Per-residue state (C / R / L) and bead position. No orientations for now.
class Chain {
public:
    explicit Chain(const Input& in);
    void init();                       // states per Input::init_state, positions per Input::init

    int N() const { return N_; }
    const Input& input() const { return in_; }

    // Which parameter set applies.
    //   bond i (between i and i+1): pair class of (state[i], state[i+1])
    //   bend at i (i-1, i, i+1):     triple class of (state[i-1], state[i], state[i+1])
    Pair   bond_class(int i) const { return pair_class(state[i], state[i + 1]); }
    Triple bend_class(int i) const {
        if (pair_keyed_bends_) {
            switch (pair_class(state[i - 1], state[i])) {
            case Pair::HH: return Triple::HHH;
            case Pair::RL: return Triple::RLR;
            case Pair::CH: return Triple::CHC;
            default:       return Triple::CCC;
            }
        }
        return triple_class(state[i - 1], state[i], state[i + 1]);
    }
    const BondParams& bond_par(int i) const { return in_.p(bond_class(i)); }
    const BendParams& bend_par(int i) const { return in_.p(bend_class(i)); }

    // geometry helpers
    Vec3   bond(int i) const { return pos[i + 1] - pos[i]; }   // i in [0, N-2]
    double bend_angle(int i) const;    // valence angle at bead i between (pos[i-1]-pos[i]) and (pos[i+1]-pos[i]);
                                       // straight chain = pi. i in [1, N-2]
    // Local tangent of the chain at bead i (unit vector): the bisector of the two adjacent bond directions,
    // t_i ~ b(i-1)/|b(i-1)| + b(i)/|b(i)|  (one bond at the chain ends). Used as the axis of a helical
    // residue when it is represented as a rod; it depends only on positions, so it needs no extra degrees
    // of freedom or moves.
    Vec3 tangent(int i) const;
    // Orientation of the rod at bead i as a unit quaternion (x, y, z, w): the minimal rotation taking the
    // z axis onto tangent(i). OVITO convention (Orientation property), cylinder/ellipsoid axis along z.
    void quaternion(int i, double q[4]) const;
    double end_to_end() const;
    double rg2() const;                // squared radius of gyration

    // state statistics
    int    count(State s) const;
    int    n_helix() const { return count(State::R) + count(State::L); }
    double helicity() const { return double(n_helix()) / N_; }

    std::vector<State> state;
    std::vector<Vec3>  pos;
    // registry direction of each residue (unit, perpendicular to tangent(i)); meaningful for helical
    // residues only, kept by registry.cpp (creation/split/merge rules) and moved with the rod by the
    // geometry moves.  Used by the tabulated pair potential (nb_hh = db).
    std::vector<Vec3>  reg;

    // Verlet neighbour list for the non-bonded energy (built and used in pair_potential.cpp).
    // Invariant while !dirty: every bead is within skin/2 of its position at the last build (ref), so
    // every pair closer than the largest cutoff is listed. A move that would break the invariant takes
    // the all-pairs path and, if accepted, marks the list dirty (it is rebuilt on its next use).
    struct NeighbourList {
        bool   on = false;
        bool   dirty = true;
        double rc_max2 = 0.0;                      // (largest pair cutoff)^2
        double r_list2 = 0.0;                      // (largest pair cutoff + skin)^2
        double half_skin2 = 0.0;                   // (skin / 2)^2
        long   n_build = 0;
        std::vector<Vec3> ref;
        std::vector<std::vector<int> > nbrs;       // per bead, ascending, only |i-j| > 1
    };
    mutable NeighbourList nl;

private:
    int N_;
    const Input& in_;
    bool pair_keyed_bends_;
};
