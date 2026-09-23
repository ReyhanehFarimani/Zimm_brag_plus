#include "registry.h"

#include <algorithm>
#include <cmath>

Vec3 registry_transport(const Vec3& m, const Vec3& u_old, const Vec3& u_new) {
    Vec3 axis = cross(u_old, u_new);
    const double s = norm(axis), c = std::max(-1.0, std::min(1.0, dot(u_old, u_new)));
    Vec3 r;
    if (s < 1e-12) {
        r = (c > 0.0) ? m : (-1.0) * m;               // parallel: nothing; antiparallel: any pi turn
    } else {
        r = rotate(m, (1.0 / s) * axis, std::atan2(s, c));
    }
    r = r - dot(r, u_new) * u_new;                    // re-orthogonalise (numerical drift)
    const double n = norm(r);
    return n > 1e-12 ? (1.0 / n) * r : registry_from_angle(u_new, 0.0);
}

Vec3 registry_from_angle(const Vec3& u, double phi) {
    Vec3 t = std::fabs(u.x) < 0.9 ? Vec3(1, 0, 0) : Vec3(0, 1, 0);
    Vec3 e1 = t - dot(t, u) * u; e1 = (1.0 / norm(e1)) * e1;
    const Vec3 e2 = cross(u, e1);
    return std::cos(phi) * e1 + std::sin(phi) * e2;
}

double registry_angle(const Vec3& u, const Vec3& a, const Vec3& b) {
    return std::atan2(dot(cross(a, b), u), dot(a, b));
}

namespace {

// common axis of the junction (i, i+1) and the two registries projected onto its plane
inline void junction(const Chain& chain, int i, Vec3& n, Vec3& p1, Vec3& p2) {
    n = chain.tangent(i) + chain.tangent(i + 1);
    const double nn = norm(n);
    n = nn > 1e-12 ? (1.0 / nn) * n : chain.tangent(i);     // antiparallel tangents: degenerate, use t_i
    p1 = chain.reg[i]     - dot(chain.reg[i], n) * n;
    p2 = chain.reg[i + 1] - dot(chain.reg[i + 1], n) * n;
}

inline double wrap_pi(double a) {
    while (a > M_PI)  a -= 2.0 * M_PI;
    while (a <= -M_PI) a += 2.0 * M_PI;
    return a;
}

} // namespace

double registry_twist(const Chain& chain, int i) {
    Vec3 n, p1, p2;
    junction(chain, i, n, p1, p2);
    return registry_angle(n, p1, p2);
}

double twist_energy(const Chain& chain, int i) {
    const State a = chain.state[i], b = chain.state[i + 1];
    if (!is_helix(a) || a != b) return 0.0;
    const Input& in = chain.input();
    if (!in.twist_on) return 0.0;
    const double d = wrap_pi(registry_twist(chain, i) - spin(a) * in.twist_alpha0_rad);
    return 0.5 * in.twist_kappa * d * d;
}

double twist_range_energy(const Chain& chain, int lo, int hi) {
    double e = 0.0;
    for (int j = std::max(0, lo - 1); j <= std::min(chain.N() - 2, hi); ++j) e += twist_energy(chain, j);
    return e;
}

double total_twist_energy(const Chain& chain) {
    return twist_range_energy(chain, 1, chain.N() - 2);
}

void registry_init(Chain& chain, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    const int N = chain.N();
    chain.reg.assign(N, Vec3(1, 0, 0));
    for (int i = 0; i < N; ++i) chain.reg[i] = registry_from_angle(chain.tangent(i), 2.0 * M_PI * unif(rng));
}

Vec3 registry_resample(Chain& chain, int i, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    const Vec3 old = chain.reg[i];
    chain.reg[i] = registry_from_angle(chain.tangent(i), 2.0 * M_PI * unif(rng));
    return old;
}

void registry_negate_twists(Chain& chain, int a, int b) {
    // m_a kept; m_{j+1} is rebuilt from m_j so that the junction twist is minus the old one.  The old
    // twists are read first (they are a function of the old registries only).
    std::vector<double> tw;
    for (int j = a; j < b; ++j) tw.push_back(registry_twist(chain, j));
    for (int j = a; j < b; ++j) {
        Vec3 n, p1, p2;
        junction(chain, j, n, p1, p2);
        // rotate the projected m_j about n by -twist, then carry it onto the plane of t_{j+1}
        const double np1 = norm(p1);
        Vec3 q = np1 > 1e-12 ? rotate((1.0 / np1) * p1, n, -tw[j - a]) : chain.reg[j];
        chain.reg[j + 1] = registry_transport(q, n, chain.tangent(j + 1));
    }
}

int registry_check(const Chain& chain, double tol) {
    int bad = 0;
    for (int i = 0; i < chain.N(); ++i) {
        const Vec3 u = chain.tangent(i);
        if (std::fabs(dot(u, chain.reg[i])) > tol || std::fabs(norm(chain.reg[i]) - 1.0) > tol) ++bad;
    }
    return bad;
}
