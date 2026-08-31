#include "pair_potential.h"

#include <algorithm>
#include <cmath>

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

double bend_energy_angle(const BendParams& p, double theta) {
    const double dth = theta - p.theta0;
    return 0.5 * p.kappa * dth * dth;
}

double bend_energy_of(const Chain& chain, int i, const Vec3& a, const Vec3& b) {
    // valence angle at bead i: a = bond(i-1) points into the bead, b = bond(i) out of it -> cos = -a.b
    double c = -dot(a, b) / std::sqrt(norm2(a) * norm2(b));
    c = std::max(-1.0, std::min(1.0, c));
    return bend_energy_angle(chain.bend_par(i), std::acos(c));
}

double hinge_log_weight(const BendParams& p) {
    // ln of  2 pi * int_0^pi sin(theta) exp(-U(theta)) dtheta  (Simpson), U from bend_energy_angle
    const int n = 4000;
    const double h = M_PI / n;
    double sum = 0.0;
    for (int k = 0; k <= n; ++k) {
        const double th = k * h;
        const double f  = std::sin(th) * std::exp(-bend_energy_angle(p, th));
        sum += f * ((k == 0 || k == n) ? 1.0 : (k % 2 ? 4.0 : 2.0));
    }
    return std::log(2.0 * M_PI * sum * h / 3.0);
}

double bend_energy(const Chain& chain, int i) {
    return bend_energy_of(chain, i, chain.bond(i - 1), chain.bond(i));
}

double total_bend_energy(const Chain& chain) {
    double e = 0.0;
    for (int i = 1; i < chain.N() - 1; ++i) e += bend_energy(chain, i);
    return e;
}

double state_energy(const Chain& chain, int i) {
    return chain.input().pair_energy(chain.bond_class(i));
}

double total_state_energy(const Chain& chain) {
    double e = 0.0;
    for (int i = 0; i < chain.N() - 1; ++i) e += state_energy(chain, i);
    return e;
}
