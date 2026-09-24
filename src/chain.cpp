#include "chain.h"
#include "registry.h"
#include "pair_potential.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

Chain::Chain(const Input& in)
    : state(in.n_arms > 0 ? in.n_arms * in.N : in.N, State::Coil), pos(in.n_arms > 0 ? in.n_arms * in.N : in.N),
      N_(in.n_arms > 0 ? in.n_arms * in.N : in.N), n_arms_(in.n_arms), arm_len_(in.N), in_(in),
      pair_keyed_bends_(in.bend_key == "pair") {
    box_ = in.box;
    const double rc = nb_cutoff_max(in);
    nl.rc_max2    = rc * rc;
    nl.on         = in.nb_type != "none" && in.nl_skin > 0.0;
    nl.r_list2    = (rc + in.nl_skin) * (rc + in.nl_skin);
    nl.half_skin2 = 0.25 * in.nl_skin * in.nl_skin;
}

void Chain::init() {
    // Own generator (offset from the MC seed) so different seeds give different starts.
    std::mt19937_64 rng(in_.seed ^ 0x9E3779B97F4A7C15ULL);
    std::uniform_real_distribution<double> unif(0.0, 1.0);

    if (!in_.state_pattern.empty()) {
        for (int i = 0; i < N_; ++i) {
            const char ch = in_.state_pattern[i];
            state[i] = ch == 'R' ? State::R : (ch == 'L' ? State::L : State::Coil);
        }
    } else if (in_.init_state == "random") {
        for (int i = 0; i < N_; ++i) state[i] = unif(rng) < 0.5 ? State::R : State::L;
    } else {
        std::fill(state.begin(), state.end(), State::Coil);
    }

    // bonds start at the coil-coil length
    const double b = in_.p(Pair::CC).bond_len;
    if (n_arms_ > 0 && box_ > 0.0) {
        // periodic box: every chain a straight rod of random orientation from a random point of the box
        // (overlaps are relaxed by the MC: the coil core is a soft Gaussian and the states start as coil)
        for (int k = 0; k < n_arms_; ++k) {
            const Vec3 p0(box_ * unif(rng), box_ * unif(rng), box_ * unif(rng));
            const double cz = 2.0 * unif(rng) - 1.0, sz = std::sqrt(std::max(0.0, 1.0 - cz * cz)), phi = 2.0 * M_PI * unif(rng);
            const Vec3 u(sz * std::cos(phi), sz * std::sin(phi), cz);
            for (int m = 0; m < arm_len_; ++m) pos[k * arm_len_ + m] = p0 + (m * b) * u;
        }
    } else if (n_arms_ > 0) {
        // star: graft sites = Fibonacci points on the sphere of radius core_radius + 1/2; each arm starts at
        // its site and runs radially outward as a straight chain (no overlaps, no core penetration)
        graft.resize(n_arms_);
        const double R0 = in_.core_radius + 0.5, ga = M_PI * (3.0 - std::sqrt(5.0));
        for (int k = 0; k < n_arms_; ++k) {
            const double z = 1.0 - 2.0 * (k + 0.5) / n_arms_, r = std::sqrt(std::max(0.0, 1.0 - z * z)), ph = ga * k;
            const Vec3 u(r * std::cos(ph), r * std::sin(ph), z);
            graft[k] = R0 * u;
            for (int m = 0; m < arm_len_; ++m) pos[k * arm_len_ + m] = (R0 + m * b) * u;
        }
    } else if (in_.init == "walk") {
        // freely-jointed random walk: uniformly random bond directions
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
    nl.dirty = true;
    registry_init(*this, rng);          // one uniform registry per residue
    if (!in_.restart_file.empty()) load_restart();
}

// Restart from a frame of a trajectory written by Logger::dump_config (extended xyz: "Properties=species:S:1:pos:R:3:
// ...:registry:R:3:..." in the comment line, the star core as an extra first particle X). States, positions and
// registries are taken over; the neighbour list is rebuilt on first use.
void Chain::load_restart() {
    std::ifstream f(in_.restart_file);
    if (!f) { std::fprintf(stderr, "restart_file %s: cannot open\n", in_.restart_file.c_str()); std::exit(EXIT_FAILURE); }
    std::string line; std::vector<std::string> frame; std::string header; int want = in_.restart_frame, idx = 0, got = -1; bool have = false;
    while (std::getline(f, line)) {
        const int n = std::atoi(line.c_str()); std::string hdr; std::getline(f, hdr);
        std::vector<std::string> rows; rows.reserve(n);
        for (int k = 0; k < n; ++k) { std::getline(f, line); rows.push_back(line); }
        if (want < 0 || idx == want) { frame.swap(rows); header = hdr; have = true; got = idx; }
        if (idx == want) break;
        ++idx;
    }
    if (!have) { std::fprintf(stderr, "restart_file %s: frame %d not found (%d frames)\n", in_.restart_file.c_str(), want, idx); std::exit(EXIT_FAILURE); }
    // column layout from the Properties string
    int col = 0, c_pos = -1, c_reg = -1, n_reg = 0;
    const size_t p0 = header.find("Properties="); std::string props = p0 == std::string::npos ? "" : header.substr(p0 + 11);
    props = props.substr(0, props.find(' '));
    std::istringstream ps(props); std::string tok; std::vector<std::string> t;
    while (std::getline(ps, tok, ':')) t.push_back(tok);
    for (size_t k = 0; k + 2 < t.size(); k += 3) {
        const int cnt = std::atoi(t[k + 2].c_str());
        if (t[k] == "pos") c_pos = col;
        if (t[k] == "registry") { c_reg = col; n_reg = cnt; }
        col += cnt;
    }
    if (c_pos < 0) { std::fprintf(stderr, "restart_file: no pos column\n"); std::exit(EXIT_FAILURE); }
    int i = 0;
    for (const std::string& r : frame) {
        std::istringstream ss(r); std::vector<std::string> w; while (ss >> tok) w.push_back(tok);
        if (w.empty() || w[0] == "X") continue;                       // the star core
        if (i >= N_) { std::fprintf(stderr, "restart_file: more residues than N = %d\n", N_); std::exit(EXIT_FAILURE); }
        state[i] = w[0] == "R" ? State::R : (w[0] == "L" ? State::L : State::Coil);
        pos[i] = Vec3(std::atof(w[c_pos].c_str()), std::atof(w[c_pos + 1].c_str()), std::atof(w[c_pos + 2].c_str()));
        if (c_reg >= 0 && n_reg == 3 && i < (int)reg.size()) reg[i] = Vec3(std::atof(w[c_reg].c_str()), std::atof(w[c_reg + 1].c_str()), std::atof(w[c_reg + 2].c_str()));
        ++i;
    }
    if (i != N_) { std::fprintf(stderr, "restart_file: %d residues, N = %d\n", i, N_); std::exit(EXIT_FAILURE); }
    const size_t sp = header.find("sweep="); const long sw = sp == std::string::npos ? -1 : std::atol(header.c_str() + sp + 6);
    std::printf("# restart: %d residues from %s frame %d of %d (sweep %ld)%s\n", i, in_.restart_file.c_str(), got, idx, sw, c_reg >= 0 ? ", registries taken over" : "");
    nl.dirty = true;
}

double Chain::bend_angle(int i) const {
    const Vec3 a = bond(i - 1), b = bond(i);
    double c = -dot(a, b) / (norm(a) * norm(b));   // minus: angle between (pos[i-1]-pos[i]) and (pos[i+1]-pos[i])
    c = std::max(-1.0, std::min(1.0, c));         // guard against round-off
    return std::acos(c);
}

double Chain::end_to_end() const {
    if (box_ > 0.0 && n_arms_ > 0) {                       // box: mean end-to-end distance of the chains
        double s = 0.0;
        for (int k = 0; k < n_arms_; ++k) s += norm(pos[k * arm_len_ + arm_len_ - 1] - pos[k * arm_len_]);
        return s / n_arms_;
    }
    return norm(pos[N_ - 1] - pos[0]);
}

double Chain::rg2() const {
    if (box_ > 0.0 && n_arms_ > 0) {                       // box: mean squared radius of gyration of the chains
        double s = 0.0;
        for (int k = 0; k < n_arms_; ++k) {
            Vec3 c; for (int m = 0; m < arm_len_; ++m) c = c + pos[k * arm_len_ + m];
            c = (1.0 / arm_len_) * c;
            for (int m = 0; m < arm_len_; ++m) s += norm2(pos[k * arm_len_ + m] - c);
        }
        return s / N_;
    }
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

Vec3 Chain::tangent(int i) const {
    Vec3 t;
    if (has_bond(i - 1)) { const Vec3 b = bond(i - 1); t = t + (1.0 / norm(b)) * b; }
    if (has_bond(i))     { const Vec3 b = bond(i);     t = t + (1.0 / norm(b)) * b; }
    const double n = norm(t);
    // the two bonds are exactly antiparallel (folded back): tangent undefined, fall back to the outgoing bond
    if (n < 1e-12) { const Vec3 b = has_bond(i) ? bond(i) : bond(i - 1); return (1.0 / norm(b)) * b; }
    return (1.0 / n) * t;
}

void Chain::quaternion(int i, double q[4]) const {
    const Vec3 t = tangent(i);
    const double c = std::max(-1.0, std::min(1.0, t.z));       // cos(angle between z and t)
    if (c < -1.0 + 1e-12) { q[0] = 1.0; q[1] = q[2] = q[3] = 0.0; return; }   // t = -z: rotate by pi about x
    // axis = z x t = (-t.y, t.x, 0);  q = (axis_unit sin(a/2), cos(a/2))  with  |z x t| = sin a
    const double s = std::sqrt(2.0 * (1.0 + c));               // = 2 cos(a/2)
    q[0] = -t.y / s; q[1] = t.x / s; q[2] = 0.0; q[3] = 0.5 * s;
}
