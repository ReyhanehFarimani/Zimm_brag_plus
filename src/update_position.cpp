#include "update_position.h"

#include <cmath>
#include "pair_potential.h"

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
    const int N = chain.N();
    double e = 0.0;
    if (i > 0)     e += bond_energy_of(chain, i - 1, norm(bond_with(chain, i - 1, i, p)));
    if (i < N - 1) e += bond_energy_of(chain, i,     norm(bond_with(chain, i,     i, p)));
    for (int k = i - 1; k <= i + 1; ++k) {
        if (k < 1 || k > N - 2) continue;
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

    const double dE = local_energy(chain, i, new_p) - local_energy(chain, i, old_p);

    if (dE <= 0.0 || unif(rng) < std::exp(-dE / in.kT)) {
        chain.pos[i] = new_p;
        return true;
    }
    return false;
}
