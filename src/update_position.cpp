#include "update_position.h"

#include <algorithm>
#include <cmath>
#include <vector>
#include "pair_potential.h"
#include "registry.h"

namespace {

// bond vector j evaluated with bead i at position p (all other beads as stored)
Vec3 bond_with(const Chain& chain, int j, int i, const Vec3& p) {
    const Vec3& a = (j == i)     ? p : chain.pos[j];
    const Vec3& b = (j + 1 == i) ? p : chain.pos[j + 1];
    return b - a;
}

// energy of everything that depends on bead i, evaluated with bead i at position p:
// the bonds (i-1,i), (i,i+1) and the bends at i-1, i, i+1
double local_energy(const Chain& chain, int i, const Vec3& p) {
    double e = core_energy_at(chain, i, p);
    if (chain.has_bond(i - 1)) e += bond_energy_of(chain, i - 1, norm(bond_with(chain, i - 1, i, p)));
    if (chain.has_bond(i))     e += bond_energy_of(chain, i,     norm(bond_with(chain, i,     i, p)));
    for (int k = i - 1; k <= i + 1; ++k) {
        if (!chain.has_bend(k)) continue;
        e += bend_energy_of(chain, k, bond_with(chain, k - 1, i, p), bond_with(chain, k, i, p));
    }
    return e;
}

} // namespace

bool try_position_move(Chain& chain, int i, const Input& in, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> disp(-in.max_disp, in.max_disp);
    std::uniform_real_distribution<double> unif(0.0, 1.0);

    const Vec3 old_p = chain.pos[i];
    const Vec3 new_p = old_p + Vec3(disp(rng), disp(rng), disp(rng));

    double dE = local_energy(chain, i, new_p) - local_energy(chain, i, old_p);
    // non-bonded part: apply the move tentatively (the tangents of i-1, i, i+1 depend on pos[i])
    if (in.nb_type[0] != 'n') {
        // the neighbour list is exact for both configurations only if the trial position keeps bead i
        // within skin/2 of its reference; otherwise this move uses the all-pairs loop and, if accepted,
        // the list is rebuilt
        const bool use_list = nl_covers(chain, i, new_p);
        const double e_nb_old = nb_local_energy(chain, i, use_list);
        // registries (nb_hh = db): the tangents of i-1, i, i+1 change with pos[i]; their registries
        // are carried by the minimal rotation old -> new tangent, and the twist term of the
        // junctions i-2 .. i+1 changes with the geometry
        const bool reg = in.nb_hh == "db";
        const int N = chain.N();
        const int lo_e = std::max(0, i - 1), hi_e = std::min(N - 1, i + 1);
        Vec3 m_old[3], u_old[3];
        if (reg) {
            for (int k = 0; k < 3; ++k) { const int j = i - 1 + k;
                if (j >= 0 && j < N) { u_old[k] = chain.tangent(j); m_old[k] = chain.reg[j]; } }
        }
        const double e_nb_old_r = reg ? nb_range_energy(chain, lo_e, hi_e, use_list) + twist_range_energy(chain, lo_e, hi_e) : e_nb_old;
        chain.pos[i] = new_p;
        if (reg) {
            for (int k = 0; k < 3; ++k) {
                const int j = i - 1 + k;
                if (j < 0 || j >= N) continue;
                chain.reg[j] = registry_transport(m_old[k], u_old[k], chain.tangent(j));
            }
        }
        dE += (reg ? nb_range_energy(chain, lo_e, hi_e, use_list) + twist_range_energy(chain, lo_e, hi_e)
                   : nb_local_energy(chain, i, use_list)) - e_nb_old_r;
        g_last_dE = dE;
        if (dE <= 0.0 || unif(rng) < std::exp(-dE / in.kT)) {
            if (!use_list) nl_invalidate(chain);
            return true;
        }
        chain.pos[i] = old_p;
        if (reg) for (int k = 0; k < 3; ++k) { const int j = i - 1 + k; if (j >= 0 && j < N) chain.reg[j] = m_old[k]; }
        return false;
    }

    if (dE <= 0.0 || unif(rng) < std::exp(-dE / in.kT)) {
        chain.pos[i] = new_p;
        return true;
    }
    return false;
}
