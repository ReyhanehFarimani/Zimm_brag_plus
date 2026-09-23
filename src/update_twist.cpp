#include "update_twist.h"

#include <cmath>
#include <vector>
#include "pair_potential.h"

bool try_twist_move(Chain& chain, int i, const Input& in, std::mt19937_64& rng) {
    if (!is_helix(chain.state[i])) return false;
    int lo = i, hi = i;
    const State s = chain.state[i];
    while (lo > 0 && chain.state[lo - 1] == s) --lo;
    while (hi + 1 < chain.N() && chain.state[hi + 1] == s) ++hi;
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    const double delta = in.twist_step * (2.0 * unif(rng) - 1.0);

    const double e_old = nb_range_energy(chain, lo, hi);
    std::vector<Vec3> old_m(chain.reg.begin() + lo, chain.reg.begin() + hi + 1);
    for (int k = lo; k <= hi; ++k) chain.reg[k] = rotate(chain.reg[k], chain.tangent(k), delta);
    const double dE = nb_range_energy(chain, lo, hi) - e_old;
    g_last_dE = dE;
    if (dE <= 0.0 || unif(rng) < std::exp(-dE / in.kT)) return true;
    for (int k = lo; k <= hi; ++k) chain.reg[k] = old_m[k - lo];      // rejected: restore
    return false;
}
