// helix_pair_db.hpp -- nearest-node reader for the pair-potential table written by
// prelim/cg/helix_pair_db_full.py of the grant repo ("HELIX_PAIR_DB v2"): a text header
// (one line per item, ends with "END"), then ONE float32 array U in C order with the axis
// order  theta0 x type x r x betaA x betaB x psi x alphaA x alphaB.  Every axis is uniform
// (value = start + step*i); psi and the alphas are periodic (360), r and the betas clamped.
// Angles in degrees, r in the campaign unit a, U in kT.  U >= ucap means overlap.
#pragma once
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct HelixPairDB {
    enum Type { RR = 0, LL = 1, RL = 2, LR = 3 };
    struct Axis { int n; double start, step; bool periodic; };
    Axis ax[6];                 // r, betaA, betaB, psi, alphaA, alphaB
    int shape[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<int> theta0;    // the theta0 values in the file, in index order
    double ucap = 100.0;
    std::string level;
    std::vector<float> u;

    void load(const std::string& path) {
        std::ifstream f(path, std::ios::binary);
        if (!f) throw std::runtime_error("cannot open " + path);
        std::string line; int nax = 0;
        while (std::getline(f, line)) {
            if (line == "END") break;
            std::istringstream ss(line); std::string key; ss >> key;
            if (key == "shape") { for (int i = 0; i < 8; ++i) ss >> shape[i]; }
            else if (key == "theta0") { int t; while (ss >> t) theta0.push_back(t); }
            else if (key == "ucap") { ss >> ucap; }
            else if (key == "level") { std::getline(ss, level); }
            else if (key == "axis") {
                std::string name, kw, per; Axis a;
                ss >> name >> kw >> a.n >> kw >> a.start >> kw >> a.step >> per;
                a.periodic = (per == "periodic");
                if (nax < 6) ax[nax++] = a;
            }
        }
        if (nax != 6 || shape[0] <= 0) throw std::runtime_error("bad header in " + path);
        size_t total = 1; for (int i = 0; i < 8; ++i) total *= (size_t)shape[i];
        u.resize(total);
        f.read(reinterpret_cast<char*>(u.data()), total * sizeof(float));
        if (!f) throw std::runtime_error("short read in " + path);
    }
    int theta_index(int th) const {
        for (size_t i = 0; i < theta0.size(); ++i) if (theta0[i] == th) return (int)i;
        return -1;
    }
    double r_max() const { return ax[0].start + ax[0].step * (ax[0].n - 1); }

    int idx(int k, double v) const {
        const Axis& a = ax[k];
        long i = (long)std::floor((v - a.start) / a.step + 0.5);   // nearest node, consistent ties
        if (a.periodic) { const long n = a.n; i %= n; if (i < 0) i += n; }
        else { if (i < 0) i = 0; if (i >= a.n) i = a.n - 1; }
        return (int)i;
    }
    double U(int ti, int type, double r, double bA, double bB, double psi, double aA, double aB) const {
        size_t o = (size_t)ti;
        o = o * shape[1] + type;
        o = o * shape[2] + idx(0, r);
        o = o * shape[3] + idx(1, bA);
        o = o * shape[4] + idx(2, bB);
        o = o * shape[5] + idx(3, psi);
        o = o * shape[6] + idx(4, aA);
        o = o * shape[7] + idx(5, aB);
        return u[o];
    }
};
