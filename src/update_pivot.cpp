#include "update_pivot.h"

#include <cmath>
#include <vector>
#include "pair_potential.h"

bool try_pivot_move(Chain& chain, int i, const Input& in, std::mt19937_64& rng) {
    const int N = chain.N();
    if (i < 1 || i > N - 2) return false;
    std::uniform_real_distribution<double> unif(0.0, 1.0);

    // random unit axis and angle
    const double cz = 2.0 * unif(rng) - 1.0, sz = std::sqrt(std::max(0.0, 1.0 - cz * cz)), ph = 2.0 * M_PI * unif(rng);
    const Vec3 axis(sz * std::cos(ph), sz * std::sin(ph), cz);
    const double ang = in.max_rot * (2.0 * unif(rng) - 1.0);

    const double e_old = bend_energy(chain, i) + nb_pivot_energy(chain, i);

    std::vector<Vec3> old_tail(chain.pos.begin() + i + 1, chain.pos.end());
    const Vec3 pivot = chain.pos[i];
    for (int k = i + 1; k < N; ++k) chain.pos[k] = pivot + rotate(chain.pos[k] - pivot, axis, ang);

    const double dE = bend_energy(chain, i) + nb_pivot_energy(chain, i) - e_old;
    if (dE <= 0.0 || unif(rng) < std::exp(-dE / in.kT)) return true;
    std::copy(old_tail.begin(), old_tail.end(), chain.pos.begin() + i + 1);   // rejected: restore
    return false;
}
