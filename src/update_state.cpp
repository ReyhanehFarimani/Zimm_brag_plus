#include "update_state.h"

#include <cmath>
#include "pair_potential.h"

namespace {

// everything whose parameters depend on state[i], at the current geometry:
// state + bond energy of the two bonds touching i, and the bends at i-1, i, i+1
double local_energy(const Chain& chain, int i) {
    const int N = chain.N();
    double e = 0.0;
    if (i > 0)     e += state_energy(chain, i - 1) + bond_energy(chain, i - 1);
    if (i < N - 1) e += state_energy(chain, i)     + bond_energy(chain, i);
    for (int k = i - 1; k <= i + 1; ++k)
        if (k >= 1 && k <= N - 2) e += bend_energy(chain, k);
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

    const State old_s = chain.state[i];
    const double e_old = local_energy(chain, i);
    chain.state[i] = from_spin(q_new);
    const double dE = local_energy(chain, i) - e_old;

    if (dE <= 0.0 || unif(rng) < std::exp(-dE / in.kT)) return true;
    chain.state[i] = old_s;                              // rejected: restore
    return false;
}
