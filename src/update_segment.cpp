#include "update_segment.h"

#include <algorithm>
#include <cmath>
#include <vector>
#include "pair_potential.h"
#include "registry.h"

bool try_segment_move(Chain& chain, int i, const Input& in, std::mt19937_64& rng) {
    std::uniform_int_distribution<int> klen(1, std::max(1, in.seg_len_max));
    const int k = klen(rng), j = i + k + 1;
    if (i < 0 || j >= chain.N() || j > chain.arm_last(i)) return false;       // both anchors in i's arm
    Vec3 axis = chain.pos[j] - chain.pos[i];
    const double L = norm(axis);
    if (L < 1e-9) return false;
    axis = (1.0 / L) * axis;
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    const double ang = in.seg_rot * (2.0 * unif(rng) - 1.0);
    const bool reg = in.nb_hh == "db";

    std::vector<Vec3> old_pos(chain.pos.begin() + i + 1, chain.pos.begin() + j), new_pos(k);
    const Vec3 p0 = chain.pos[i];
    for (int s = 0; s < k; ++s) new_pos[s] = p0 + rotate(old_pos[s] - p0, axis, ang);

    // neighbour list: re-reference the moving beads at their current positions, then the list covers the
    // trial positions if every displacement is below skin/2; otherwise fall back to the all-pairs loops
    bool use_list = chain.nl.on;
    if (use_list) {
        for (int s = 0; s < k; ++s) {
            const int b = i + 1 + s;
            if (!nl_covers(chain, b, new_pos[s])) nl_update_bead(chain, b);
            if (!nl_covers(chain, b, new_pos[s])) use_list = false;
        }
    }
    auto e_all = [&]() {
        double e = 0.0;
        if (chain.has_bend(i)) e += bend_energy(chain, i);
        if (chain.has_bend(j)) e += bend_energy(chain, j);
        e += nb_range_energy(chain, i, j, use_list) + core_range_energy(chain, i + 1, j - 1);
        if (reg) e += twist_range_energy(chain, i, j);
        return e;
    };
    const double e_old = e_all();

    std::vector<Vec3> old_reg; Vec3 u_i_old, u_j_old;
    if (reg) { old_reg.assign(chain.reg.begin() + i, chain.reg.begin() + j + 1); u_i_old = chain.tangent(i); u_j_old = chain.tangent(j); }
    for (int s = 0; s < k; ++s) chain.pos[i + 1 + s] = new_pos[s];
    if (reg) {
        for (int s = 1; s <= k; ++s) chain.reg[i + s] = rotate(chain.reg[i + s], axis, ang);   // rigid with the segment
        chain.reg[i] = registry_transport(old_reg[0], u_i_old, chain.tangent(i));            // anchors: minimal rotation
        chain.reg[j] = registry_transport(old_reg[k + 1], u_j_old, chain.tangent(j));
    }
    const double dE = e_all() - e_old;
    g_last_dE = dE;
    if (dE <= 0.0 || unif(rng) < std::exp(-dE / in.kT)) {
        if (!use_list && chain.nl.on && !chain.nl.dirty)
            for (int s = 0; s < k; ++s) nl_update_bead(chain, i + 1 + s);                 // beads left their skin
        return true;
    }
    std::copy(old_pos.begin(), old_pos.end(), chain.pos.begin() + i + 1);                  // rejected: restore
    if (reg) std::copy(old_reg.begin(), old_reg.end(), chain.reg.begin() + i);
    return false;
}
