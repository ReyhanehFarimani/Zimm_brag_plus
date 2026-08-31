#include "update_position.h"

#include <cmath>
#include "pair_potential.h"

namespace {

// bond energy of the bonds touching bead i, evaluated with bead i at position p
double local_bond_energy(const Chain& chain, int i, const Vec3& p) {
    double e = 0.0;
    if (i > 0)             e += bond_energy_of(chain, i - 1, norm(p - chain.pos[i - 1]));
    if (i < chain.N() - 1) e += bond_energy_of(chain, i,     norm(chain.pos[i + 1] - p));
    return e;
}

} // namespace

bool try_position_move(Chain& chain, int i, const Input& in, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> disp(-in.max_disp, in.max_disp);
    std::uniform_real_distribution<double> unif(0.0, 1.0);

    const Vec3 old_p = chain.pos[i];
    const Vec3 new_p = old_p + Vec3(disp(rng), disp(rng), disp(rng));

    const double dE = local_bond_energy(chain, i, new_p) - local_bond_energy(chain, i, old_p);

    if (dE <= 0.0 || unif(rng) < std::exp(-dE / in.kT)) {
        chain.pos[i] = new_p;
        return true;
    }
    return false;
}
