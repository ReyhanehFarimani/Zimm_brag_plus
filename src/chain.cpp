#include "chain.h"

#include <algorithm>
#include <cmath>
#include <random>

Chain::Chain(const Input& in)
    : state(in.N, State::Coil), pos(in.N), N_(in.N), in_(in) {}

void Chain::init() {
    // every residue coil; bonds have the coil-coil length
    const double b = in_.p(Pair::CC).bond_len;
    std::fill(state.begin(), state.end(), State::Coil);

    if (in_.init == "walk") {
        // freely-jointed random walk: uniformly random bond directions.
        // Own generator (offset from the MC seed) so different seeds give different starts.
        std::mt19937_64 rng(in_.seed ^ 0x9E3779B97F4A7C15ULL);
        std::uniform_real_distribution<double> unif(0.0, 1.0);
        pos[0] = Vec3(0.0, 0.0, 0.0);
        for (int i = 1; i < N_; ++i) {
            const double cz  = 2.0 * unif(rng) - 1.0;           // uniform on the sphere
            const double sz  = std::sqrt(std::max(0.0, 1.0 - cz * cz));
            const double phi = 2.0 * M_PI * unif(rng);
            pos[i] = pos[i - 1] + b * Vec3(sz * std::cos(phi), sz * std::sin(phi), cz);
        }
    } else {
        // straight rod along z
        for (int i = 0; i < N_; ++i) pos[i] = Vec3(0.0, 0.0, i * b);
    }
}

double Chain::bend_angle(int i) const {
    const Vec3 a = bond(i - 1), b = bond(i);
    double c = dot(a, b) / (norm(a) * norm(b));
    c = std::max(-1.0, std::min(1.0, c));   // guard against round-off
    return std::acos(c);
}

double Chain::end_to_end() const {
    return norm(pos[N_ - 1] - pos[0]);
}

double Chain::rg2() const {
    Vec3 com;
    for (const auto& p : pos) com = com + p;
    com = (1.0 / N_) * com;
    double r2 = 0.0;
    for (const auto& p : pos) r2 += norm2(p - com);
    return r2 / N_;
}

int Chain::count(State s) const {
    return static_cast<int>(std::count(state.begin(), state.end(), s));
}
