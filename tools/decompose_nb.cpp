// Decompose the non-bonded energy of dumped configurations with the SAME potential code as the MC:
//   E_nb = E_iso (coil-coil and coil-helix gaussian core) + E_hh,achiral (A Ia) + E_hh,chiral (eps_s C Ic)
// The chiral part is obtained exactly as E_nb(eps_s) - E_nb(eps_s = 0) on the same configuration.
// Usage: decompose_nb <input.dat> <conf.xyz>      (not part of the zimm build: compile separately, see below)
//   g++ -std=c++11 -O2 -Isrc tools/decompose_nb.cpp src/pair_potential.cpp src/chain.cpp src/input.cpp -o tools/decompose_nb
// Output, one line per frame:
//   sweep  n_helix  n_hh_pairs(inside the fit cutoff)  E_iso  E_hh_achiral  E_chiral  sum|chiral per pair|  E_chiral_homo  E_chiral_hetero
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include "pair_potential.h"

int main(int argc, char** argv) {
    if (argc < 3) { std::fprintf(stderr, "usage: decompose_nb input.dat conf.xyz\n"); return 1; }
    Input in; if (!read_input(argv[1], in)) return 1;
    const double eps = in.hf_eps_s;
    Chain chain(in);
    std::ifstream f(argv[2]); std::string line;
    const double rc_fit = 3.0 * in.hf_len;
    while (std::getline(f, line)) {
        const int n = std::atoi(line.c_str()); if (n != chain.N()) break;
        std::getline(f, line); long sweep = 0; { const auto p = line.find("sweep="); if (p != std::string::npos) sweep = std::atol(line.c_str() + p + 6); }
        bool ok = true;
        for (int i = 0; i < n; ++i) {
            if (!std::getline(f, line)) { ok = false; break; }
            std::istringstream ss(line); std::string sp; double x, y, z; ss >> sp >> x >> y >> z;
            chain.state[i] = sp == "R" ? State::R : (sp == "L" ? State::L : State::Coil); chain.pos[i] = Vec3(x, y, z);
        }
        if (!ok) break;                                                  // half-written last frame of a running job
        double e_iso = 0, e_a = 0, e_c = 0, e_abs = 0, e_homo = 0, e_het = 0; int npair = 0;
        for (int a = 0; a < n; ++a) for (int b = a + 2; b < n; ++b) {
            const bool hh = is_helix(chain.state[a]) && is_helix(chain.state[b]);
            in.hf_eps_s = 0.0; const double e0 = nb_pair_energy(chain, a, b);
            if (!hh) { e_iso += e0; continue; }
            in.hf_eps_s = eps; const double c = nb_pair_energy(chain, a, b) - e0;
            e_a += e0; e_c += c; e_abs += std::fabs(c);
            (chain.state[a] == chain.state[b] ? e_homo : e_het) += c;
            if (norm(chain.pos[b] - chain.pos[a]) < rc_fit) ++npair;
        }
        std::printf("%ld %d %d %.6f %.6f %.6f %.6f %.6f %.6f\n", sweep, chain.n_helix(), npair, e_iso, e_a, e_c, e_abs, e_homo, e_het);
    }
    return 0;
}
