#include "update_pivot.h"

#include <cmath>
#include <algorithm>
#include <vector>
#include "pair_potential.h"
#include "registry.h"

bool try_pivot_move(Chain& chain, int i, const Input& in, std::mt19937_64& rng) {
    const int N = chain.N();
    if (i < 1 || i > N - 2) return false;
    std::uniform_real_distribution<double> unif(0.0, 1.0);

    // random unit axis and angle
    const double cz = 2.0 * unif(rng) - 1.0, sz = std::sqrt(std::max(0.0, 1.0 - cz * cz)), ph = 2.0 * M_PI * unif(rng);
    const Vec3 axis(sz * std::cos(ph), sz * std::sin(ph), cz);
    const double ang = in.max_rot * (2.0 * unif(rng) - 1.0);

    // registries (nb_hh = db): the tail rotates rigidly (registries with it, transport-consistent);
    // residue i's tangent changes, so its run is rederived from i on -- those residues' pairs with the
    // rest of the tail change too and are added to the pivot energy (all-pairs, counted once)
    const bool reg = in.nb_hh == "db";
    int hi_run = i;
    if (reg && is_helix(chain.state[i])) { int lo; registry_run(chain, i, lo, hi_run); }
    auto e_all = [&]() {
        double e = bend_energy(chain, i) + nb_pivot_energy(chain, i);
        if (reg && hi_run > i) {
            for (int a = i + 1; a <= hi_run; ++a)
                for (int b = i + 1; b < N; ++b) {
                    if (b <= a + 1 && b >= a - 1) continue;
                    if (b <= hi_run && b < a) continue;                 // in-run pairs once
                    e += nb_pair_energy(chain, a, b);
                }
        }
        return e;
    };
    const double e_old = e_all();

    std::vector<Vec3> old_tail(chain.pos.begin() + i + 1, chain.pos.end());
    const Vec3 pivot = chain.pos[i];
    std::vector<Vec3> old_reg;
    const Vec3 u_i_old = chain.tangent(i);
    if (reg) old_reg.assign(chain.reg.begin() + i, chain.reg.end());
    for (int k = i + 1; k < N; ++k) chain.pos[k] = pivot + rotate(chain.pos[k] - pivot, axis, ang);
    if (reg) {
        for (int k = i + 1; k < N; ++k) chain.reg[k] = rotate(chain.reg[k], axis, ang);
        if (is_helix(chain.state[i])) {
            int lo, hi; registry_run(chain, i, lo, hi);
            if (i == lo) chain.reg[i] = registry_transport(old_reg[0], u_i_old, chain.tangent(i));
            registry_rederive(chain, std::max(i, lo + 1));
        }
    }

    const double dE = e_all() - e_old;
    g_last_dE = dE;
    // a pivot moves the whole tail by arbitrary distances: its energy is an all-pairs sum (nb_pivot_energy
    // does not use the neighbour list) and an accepted pivot invalidates the list
    if (dE <= 0.0 || unif(rng) < std::exp(-dE / in.kT)) { nl_invalidate(chain); return true; }
    std::copy(old_tail.begin(), old_tail.end(), chain.pos.begin() + i + 1);   // rejected: restore
    if (reg) std::copy(old_reg.begin(), old_reg.end(), chain.reg.begin() + i);
    return false;
}
