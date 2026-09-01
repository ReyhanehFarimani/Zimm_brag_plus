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

double nb_iso_energy(const Input& in, double r2) {
    if (in.nb_type[0] == 'n') return 0.0;
    if (r2 >= in.nb_rcut * in.nb_rcut) return 0.0;
    if (in.nb_type[0] == 'g') {                                   // gauss
        return in.nb_A * std::exp(-r2 / (in.nb_sigma * in.nb_sigma));
    }
    const double s2 = in.nb_sigma * in.nb_sigma / r2;             // wca
    const double s6 = s2 * s2 * s2;
    return 4.0 * in.nb_A * (s6 * s6 - s6) + in.nb_A;
}

double gb_energy(const Input& in, const Vec3& u1, const Vec3& u2, const Vec3& rvec) {
    const double D = 2.0 * in.rod_r, L = in.rod_L, s0 = D;
    const double k = L / D, chi = (k * k - 1.0) / (k * k + 1.0);
    const double r = norm(rvec);
    if (r <= 0.0) return 1e30;
    // quick reject: the largest contact distance is L (end to end), the truncation adds (2^(1/6)-1) s0
    if (r > L + 0.13 * s0) return 0.0;
    const Vec3 rh = (1.0 / r) * rvec;
    const double a = dot(rh, u1), b = dot(rh, u2), c = dot(u1, u2);
    const double ap = a + b, am = a - b;
    const double sig = s0 / std::sqrt(1.0 - 0.5 * chi * (ap * ap / (1.0 + chi * c) + am * am / (1.0 - chi * c)));
    const double rho = r - sig + s0;
    const double rc  = std::pow(2.0, 1.0 / 6.0) * s0;
    if (rho >= rc) return 0.0;
    double eps = in.gb_eps0;
    if (in.gb_aniso_eps) {
        const double chip = (std::pow(in.gb_kappa_p, 1.0 / in.gb_mu) - 1.0) / (std::pow(in.gb_kappa_p, 1.0 / in.gb_mu) + 1.0);
        const double e1 = std::pow(1.0 - chi * chi * c * c, -0.5 * in.gb_nu);
        const double e2 = std::pow(1.0 - 0.5 * chip * (ap * ap / (1.0 + chip * c) + am * am / (1.0 - chip * c)), in.gb_mu);
        eps *= e1 * e2;
    }
    if (rho <= 0.0) return 1e30;                                  // deep overlap
    const double q2 = (s0 / rho) * (s0 / rho), q6 = q2 * q2 * q2;
    return 4.0 * eps * (q6 * q6 - q6) + eps;
}

bool nb_anisotropic(const Input& in) { return in.nb_type[0] != 'n' && in.nb_hh == "gb"; }

double nb_pair_energy(const Chain& chain, int a, int b) {
    const Input& in = chain.input();
    if (in.nb_hh == "gb" && is_helix(chain.state[a]) && is_helix(chain.state[b]))
        return gb_energy(in, chain.tangent(a), chain.tangent(b), chain.pos[b] - chain.pos[a]);
    return nb_iso_energy(in, norm2(chain.pos[a] - chain.pos[b]));
}

double nb_bead_energy(const Chain& chain, int i) {
    if (chain.input().nb_type[0] == 'n') return 0.0;
    double e = 0.0;
    const int N = chain.N();
    for (int j = 0; j < N; ++j)
        if (j < i - 1 || j > i + 1) e += nb_pair_energy(chain, i, j);
    return e;
}

double nb_local_energy(const Chain& chain, int i) {
    const Input& in = chain.input();
    if (in.nb_type[0] == 'n') return 0.0;
    if (!nb_anisotropic(in)) return nb_bead_energy(chain, i);
    const int N = chain.N();
    const int lo = std::max(0, i - 1), hi = std::min(N - 1, i + 1);
    double e = 0.0;
    for (int a = lo; a <= hi; ++a) {
        for (int b = 0; b < N; ++b) {
            if (b >= lo && b <= hi) { if (b > a + 1) e += nb_pair_energy(chain, a, b); continue; }   // inside the window: once
            if (b < a - 1 || b > a + 1) e += nb_pair_energy(chain, a, b);
        }
    }
    return e;
}

double nb_pivot_energy(const Chain& chain, int i) {
    if (chain.input().nb_type[0] == 'n') return 0.0;
    double e = 0.0;
    const int N = chain.N();
    for (int a = 0; a <= i; ++a)
        for (int b = std::max(i, a + 2); b < N; ++b)
            e += nb_pair_energy(chain, a, b);
    return e;
}

double total_nb_energy(const Chain& chain) {
    if (chain.input().nb_type[0] == 'n') return 0.0;
    double e = 0.0;
    const int N = chain.N();
    for (int a = 0; a < N; ++a)
        for (int b = a + 2; b < N; ++b)
            e += nb_pair_energy(chain, a, b);
    return e;
}
