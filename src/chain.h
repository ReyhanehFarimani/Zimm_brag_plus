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

// Per-residue state (C / R / L) and bead position. No orientations for now.
class Chain {
public:
    explicit Chain(const Input& in);
    void init();                       // all coil; positions per Input::init ("rod" or "walk")

    int N() const { return N_; }

    // Which parameter set applies.
    //   bond i (between i and i+1): pair class of (state[i], state[i+1])
    //   bend at i (i-1, i, i+1):     triple class of (state[i-1], state[i], state[i+1])
    Pair   bond_class(int i) const { return pair_class(state[i], state[i + 1]); }
    Triple bend_class(int i) const { return triple_class(state[i - 1], state[i], state[i + 1]); }
    const BondParams& bond_par(int i) const { return in_.p(bond_class(i)); }
    const BendParams& bend_par(int i) const { return in_.p(bend_class(i)); }

    // geometry helpers
    Vec3   bond(int i) const { return pos[i + 1] - pos[i]; }   // i in [0, N-2]
    double bend_angle(int i) const;    // angle between bond(i-1) and bond(i), i in [1, N-2]
    double end_to_end() const;
    double rg2() const;                // squared radius of gyration

    // state statistics
    int    count(State s) const;
    int    n_helix() const { return count(State::R) + count(State::L); }
    double helicity() const { return double(n_helix()) / N_; }

    std::vector<State> state;
    std::vector<Vec3>  pos;

private:
    int N_;
    const Input& in_;
};
