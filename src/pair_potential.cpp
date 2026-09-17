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

double site_energy(const Chain& chain, int i) {
    return chain.input().site_energy(chain.state[i]);
}

double total_state_energy(const Chain& chain) {
    double e = 0.0;
    for (int i = 0; i < chain.N() - 1; ++i) e += state_energy(chain, i);
    for (int i = 0; i < chain.N(); ++i) e += site_energy(chain, i);
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

// ---------------------------------------------------------------------------
// nb_hh = fit: the measured helix-helix radial pair potential (umbrella
// campaign 2026-09).  U(r, Omega) = A(Omega) Ia(r) + eps_s C(Omega) Ic(r):
//   A = a0 + a1 cos psi + a2 cos 2psi           (achiral sector)
//   C = b1 sin psi + b2 sin 2psi                (chiral sector, mirror-odd)
//   coefficient(e1, e2) = quadratic in S = e1+e2, D = e1-e2
//     with e1 = rh.u1, e2 = rh.u2, psi = atan2((u1 x u2).rh, u1.u2 - e1 e2),
//   I(r) = sqrt(pi/2) sig erfc(r / sqrt2 sig) exp(r0^2 / 2 sig^2), r0 = 1.7 a
// so the repulsion-positive radial force is exactly -dU/dr.  Units inside:
// kT and campaign length a; code lengths are divided by in.hf_len.
// Fitted domain: e in [-0.5, 0.5], r in [1.5, 2.5] a; the theta0 = 100
// near-planar locking pocket is NOT in this smooth form.

struct HelixFitCoef { double a0[6], a1[6], a2[6], b1[6], b2[6]; };

// homochiral (joint R.R/L.L mirror fit) and R.L, theta0 = 100 deg
static const HelixFitCoef HF_HOMO_100 = {
    {+2.347, +1.035, -0.009, +1.277, +2.733, -0.631},
    {+0.256, -1.228, -0.562, +0.646, +0.708, +1.327},
    {+0.985, +1.151, -0.328, -1.476, +0.242, -0.437},
    {-0.025, +0.136, -0.038, -0.313, +0.175, +0.237},
    {+0.039, -0.315, +0.090, +0.572, +0.221, -0.228}};
static const HelixFitCoef HF_RL_100 = {
    {+2.509, 0.0, -0.172, +2.529, +2.702, 0.0},
    {+0.049, 0.0, -0.415, -0.764, +1.023, 0.0},
    {+1.129, 0.0, -0.455, +0.039, +0.371, 0.0},
    {0.0, +0.007, 0.0, 0.0, 0.0, -0.052},
    {0.0, +0.097, 0.0, 0.0, 0.0, +0.319}};
static const HelixFitCoef HF_HOMO_45 = {
    {+7.205, -0.246, -0.613, +1.816, +1.625, +0.672},
    {+0.198, -1.875, +0.553, +1.555, +0.959, -1.467},
    {+0.047, -0.774, +0.778, +2.336, +1.985, -1.450},
    {-0.159, +0.861, +0.141, -0.898, +0.129, -0.434},
    {-0.477, -0.066, +0.329, -0.606, -0.458, -0.478}};
static const HelixFitCoef HF_RL_45 = {
    {+7.130, 0.0, -0.521, +1.649, +1.740, 0.0},
    {-0.123, 0.0, +0.242, -0.260, +1.101, 0.0},
    {-0.065, 0.0, +0.564, +1.216, +1.765, 0.0},
    {0.0, -0.062, 0.0, 0.0, 0.0, -0.265},
    {0.0, +0.631, 0.0, 0.0, 0.0, +1.075}};

namespace {

const double HF_R0   = 1.7;   // reference shell [a]
const double HF_RCUT = 3.0;   // cutoff [a]; tail there ~1e-2 kT

// family MEDIANS of the per-window radial widths (robust)
double hf_sig_a(int th) { return th == 45 ? 0.87 : 0.93; }
double hf_sig_c(int th) { return th == 45 ? 0.51 : 0.70; }

double hf_quad(const double c[6], double S, double D) {
    return c[0] + c[1] * S + c[2] * D + c[3] * S * S + c[4] * D * D + c[5] * S * D;
}

double hf_tail(double r, double sig) {
    // integral_r^inf exp(-(x^2 - r0^2) / 2 sig^2) dx
    return std::sqrt(M_PI / 2.0) * sig * std::exp(HF_R0 * HF_R0 / (2.0 * sig * sig))
         * std::erfc(r / (std::sqrt(2.0) * sig));
}

} // namespace

// hand1, hand2 = +1 (R) / -1 (L): the spins of the two helical residues
double helixfit_energy(const Input& in, const Vec3& u1, const Vec3& u2,
                       const Vec3& rvec, int hand1, int hand2) {
    const double r = norm(rvec) / in.hf_len;     // -> campaign units [a]
    if (r <= 0.0) return 1e30;
    if (r >= HF_RCUT) return 0.0;
    const Vec3 rh = (1.0 / norm(rvec)) * rvec;
    double e1 = dot(rh, u1), e2 = dot(rh, u2);
    const double x = dot(u1, u2);
    const double P = dot(cross(u1, u2), rh);     // pseudoscalar
    double psi = std::atan2(P, x - e1 * e2);     // COLVARS psi
    int h1 = hand1, h2 = hand2;
    if (e1 + e2 < 0.0) {                         // exchange 1 <-> 2
        const double t = e1; e1 = -e2; e2 = -t;
        psi = -psi; std::swap(h1, h2);
    }
    double chir = 1.0;
    const bool homo = (h1 == h2);
    if (homo) chir = (h1 > 0 ? 1.0 : -1.0);      // L.L mirrors R.R
    else if (h1 < 0) psi = -psi;                 // L.R mirrors R.L
    const HelixFitCoef& c = (in.hf_theta0 == 45) ? (homo ? HF_HOMO_45 : HF_RL_45)
                                                 : (homo ? HF_HOMO_100 : HF_RL_100);
    const double S = e1 + e2, D = e1 - e2;
    const double A = hf_quad(c.a0, S, D)
                   + hf_quad(c.a1, S, D) * std::cos(psi)
                   + hf_quad(c.a2, S, D) * std::cos(2.0 * psi);
    const double C = hf_quad(c.b1, S, D) * std::sin(psi)
                   + hf_quad(c.b2, S, D) * std::sin(2.0 * psi);
    return in.hf_scale * (A * hf_tail(r, hf_sig_a(in.hf_theta0))
                          + in.hf_eps_s * chir * C * hf_tail(r, hf_sig_c(in.hf_theta0)));
}

bool nb_anisotropic(const Input& in) {
    return in.nb_type[0] != 'n' && (in.nb_hh == "gb" || in.nb_hh == "fit");
}

// true if the pair energy depends on the helix handedness (R vs L)
bool nb_chiral(const Input& in) { return in.nb_type[0] != 'n' && in.nb_hh == "fit"; }

// helix-helix: Gay-Berne rods (nb_hh = gb) or the measured potential
// (nb_hh = fit); coil-coil and helix-coil: the isotropic core
double nb_pair_energy(const Chain& chain, int a, int b) {
    const Input& in = chain.input();
    if (is_helix(chain.state[a]) && is_helix(chain.state[b])) {
        if (in.nb_hh == "gb")
            return gb_energy(in, chain.tangent(a), chain.tangent(b), chain.pos[b] - chain.pos[a]);
        if (in.nb_hh == "fit")
            return helixfit_energy(in, chain.tangent(a), chain.tangent(b), chain.pos[b] - chain.pos[a],
                                   spin(chain.state[a]), spin(chain.state[b]));
    }
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
