#include "pair_potential.h"

double bond_energy_of(const Chain& chain, int i, double r) {
    const BondParams& p = chain.bond_par(i);
    const double dr = r - p.bond_len;
    return 0.5 * p.k_bond * dr * dr;
}

double bond_energy(const Chain& chain, int i) {
    return bond_energy_of(chain, i, norm(chain.bond(i)));
}

double total_bond_energy(const Chain& chain) {
    double e = 0.0;
    for (int i = 0; i < chain.N() - 1; ++i) e += bond_energy(chain, i);
    return e;
}
