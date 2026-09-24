#include "update_pivot.h"

#include <cmath>
#include <algorithm>
#include <vector>
#include "pair_potential.h"
#include "registry.h"

bool try_pivot_move(Chain& chain, int i, const Input& in, std::mt19937_64& rng, bool torsion) {
    const int N = chain.N();
    // the tail i+1 .. last residue of i's arm rotates rigidly about bead i (a star arm ends at arm_last)
    const int t1 = chain.arm_last(i);
    if (i + 1 > t1) return false;
    std::uniform_real_distribution<double> unif(0.0, 1.0);

    // random unit axis and angle
    Vec3 axis;
    double ang;
    if (torsion) {                                                  // about the bond i -> i+1
        axis = chain.pos[i + 1] - chain.pos[i];
        const double L = norm(axis);
        if (L < 1e-9) return false;
        axis = (1.0 / L) * axis;
        ang = in.torsion_rot * (2.0 * unif(rng) - 1.0);
    } else {
        const double cz = 2.0 * unif(rng) - 1.0, sz = std::sqrt(std::max(0.0, 1.0 - cz * cz)), ph = 2.0 * M_PI * unif(rng);
        axis = Vec3(sz * std::cos(ph), sz * std::sin(ph), cz);
        ang = in.max_rot * (2.0 * unif(rng) - 1.0);
    }

    // registries (nb_hh = db): the tail rotates rigidly (registries with it); residue i's tangent
    // changes, its registry is carried by the minimal rotation; the twist term of the junctions
    // (i-1, i) and (i, i+1) changes, the tail's junctions are rigidly rotated (unchanged)
    const bool reg = in.nb_hh == "db";
    auto e_all = [&]() {
        double e = bend_energy(chain, i) + nb_pivot_energy(chain, i) + core_range_energy(chain, i + 1, t1);
        if (reg) e += twist_range_energy(chain, i, i);
        return e;
    };
    const double e_old = e_all();

    std::vector<Vec3> old_tail(chain.pos.begin() + i + 1, chain.pos.begin() + t1 + 1);
    const Vec3 pivot = chain.pos[i];
    std::vector<Vec3> old_reg;
    const Vec3 u_i_old = chain.tangent(i);
    if (reg) old_reg.assign(chain.reg.begin() + i, chain.reg.begin() + t1 + 1);
    for (int k = i + 1; k <= t1; ++k) chain.pos[k] = pivot + rotate(chain.pos[k] - pivot, axis, ang);
    if (reg) {
        for (int k = i + 1; k <= t1; ++k) chain.reg[k] = rotate(chain.reg[k], axis, ang);
        chain.reg[i] = registry_transport(old_reg[0], u_i_old, chain.tangent(i));
    }

    const double dE = e_all() - e_old;
    g_last_dE = dE;
    // a pivot moves the whole tail by arbitrary distances: its energy is an all-pairs sum (nb_pivot_energy
    // does not use the neighbour list) and an accepted pivot invalidates the list
    if (dE <= 0.0 || unif(rng) < std::exp(-dE / in.kT)) {
        if (chain.nl.on && !chain.nl.dirty && t1 - i <= chain.N() / 8) { for (int k = i + 1; k <= t1; ++k) nl_update_bead(chain, k); }
        else nl_invalidate(chain);                                     // long tails: the full rebuild is cheaper
        return true;
    }
    std::copy(old_tail.begin(), old_tail.end(), chain.pos.begin() + i + 1);   // rejected: restore
    if (reg) std::copy(old_reg.begin(), old_reg.end(), chain.reg.begin() + i);
    return false;
}
