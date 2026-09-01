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
    bool pair_keyed_bends_;
};
