#include "update_twist.h"

#include <cmath>
#include "pair_potential.h"
#include "registry.h"

bool try_twist_move(Chain& chain, int i, const Input& in, std::mt19937_64& rng) {
    if (!is_helix(chain.state[i])) return false;
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    const double delta = in.twist_step * M_PI / 180.0 * (2.0 * unif(rng) - 1.0);
    const double e_old = nb_range_energy(chain, i, i) + twist_range_energy(chain, i, i);
    const Vec3 old_m = chain.reg[i];
    chain.reg[i] = rotate(chain.reg[i], chain.tangent(i), delta);
    const double dE = nb_range_energy(chain, i, i) + twist_range_energy(chain, i, i) - e_old;
    g_last_dE = dE;
    if (dE <= 0.0 || unif(rng) < std::exp(-dE / in.kT)) return true;
    chain.reg[i] = old_m;                                             // rejected: restore
    return false;
}
