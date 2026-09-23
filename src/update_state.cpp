#include "update_state.h"

#include <cmath>
#include "pair_potential.h"
#include "registry.h"

namespace {

// everything whose parameters depend on state[i], at the current geometry:
// state + bond energy of the two bonds touching i, and the bends at i-1, i, i+1
double local_energy(const Chain& chain, int i) {
    double e = site_energy(chain, i);
    if (chain.has_bond(i - 1)) e += state_energy(chain, i - 1) + bond_energy(chain, i - 1);
    if (chain.has_bond(i))     e += state_energy(chain, i)     + bond_energy(chain, i);
    for (int k = i - 1; k <= i + 1; ++k)
        if (chain.has_bend(k)) e += bend_energy(chain, k);
    return e;
}

} // namespace

bool try_state_move(Chain& chain, int i, const Input& in, std::mt19937_64& rng) {
    std::uniform_int_distribution<int>     coin(0, 1);
    std::uniform_real_distribution<double> unif(0.0, 1.0);

    int q_new;
    if (in.n_states == 2) {
        q_new = -spin(chain.state[i]);                   // two-state model: flip R <-> L directly
    } else {
        q_new = spin(chain.state[i]) + (coin(rng) ? 1 : -1);
        if (q_new < -1 || q_new > 1) return false;       // R -> "+2" or L -> "-2": rejected
    }

    const State old_s = chain.state[i], new_s = from_spin(q_new);
    // the non-bonded pair type of (i, j) changes with the state of i (helix-helix vs isotropic);
    // with registries (nb_hh = db) m_i is redrawn uniformly (symmetric proposal) and the twist
    // term of the junctions (i-1, i), (i, i+1) changes with the states
    const bool reg = in.nb_hh == "db";
    const double e_old = local_energy(chain, i)
        + (nb_anisotropic(in) ? nb_bead_energy(chain, i) : 0.0)
        + (reg ? twist_range_energy(chain, i, i) : 0.0);
    chain.state[i] = new_s;
    Vec3 m_old;
    if (reg) m_old = registry_resample(chain, i, rng);
    const double dE = local_energy(chain, i)
        + (nb_anisotropic(in) ? nb_bead_energy(chain, i) : 0.0)
        + (reg ? twist_range_energy(chain, i, i) : 0.0) - e_old;

    g_last_dE = dE;
    if (dE <= 0.0 || unif(rng) < std::exp(-dE / in.kT)) return true;
    chain.state[i] = old_s;                              // rejected: restore
    if (reg) chain.reg[i] = m_old;
    return false;
}
