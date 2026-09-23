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

// maximal run of residues with state s starting at i and going right
int run_right(const Chain& chain, int i, State s) {
    int hi = i;
    while (hi + 1 < chain.N() && chain.state[hi + 1] == s) ++hi;
    return hi;
}

void twist_run(Chain& chain, int lo, int hi, double delta) {
    for (int k = lo; k <= hi; ++k) chain.reg[k] = rotate(chain.reg[k], chain.tangent(k), delta);
}

} // namespace

int registry_change_hi(const Chain& chain, int i, State s_old, State s_new) {
    const int N = chain.N();
    int hi = i;
    if (is_helix(s_old) && i + 1 < N && chain.state[i + 1] == s_old)       // split off the right part
        hi = run_right(chain, i + 1, s_old);
    if (is_helix(s_new) && i > 0 && chain.state[i - 1] == s_new &&
        i + 1 < N && chain.state[i + 1] == s_new)                           // merge: right run re-twisted
        hi = std::max(hi, run_right(chain, i + 1, s_new));
    return hi;
}

void registry_on_state_change(Chain& chain, int i, State s_old, State s_new,
                              std::mt19937_64& rng, RegistryUndo& undo) {
    const int N = chain.N();
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    undo.lo = i;
    undo.hi = registry_change_hi(chain, i, s_old, s_new);
    undo.m.assign(chain.reg.begin() + undo.lo, chain.reg.begin() + undo.hi + 1);
    // 1. the run that continued to the right of the OLD state is split off: fresh twist
    if (is_helix(s_old) && i + 1 < N && chain.state[i + 1] == s_old)
        twist_run(chain, i + 1, run_right(chain, i + 1, s_old), 2.0 * M_PI * unif(rng));
    // 2. the NEW state joins / creates a run
    if (!is_helix(s_new)) return;
    const Vec3 u = chain.tangent(i);
    const bool left  = i > 0 && chain.state[i - 1] == s_new;
    const bool right = i + 1 < N && chain.state[i + 1] == s_new;
    if (left) {
        chain.reg[i] = registry_transport(chain.reg[i - 1], chain.tangent(i - 1), u);
        if (right) {                                                        // merge: twist the right run to match
            const Vec3 u1 = chain.tangent(i + 1);
            const Vec3 want = registry_transport(chain.reg[i], u, u1);
            twist_run(chain, i + 1, run_right(chain, i + 1, s_new),
                      registry_angle(u1, chain.reg[i + 1], want));
        }
    } else if (right) {
        chain.reg[i] = registry_transport(chain.reg[i + 1], chain.tangent(i + 1), u);
    } else {
        chain.reg[i] = registry_from_angle(u, 2.0 * M_PI * unif(rng));
    }
}

void registry_undo(Chain& chain, const RegistryUndo& undo) {
    for (int k = undo.lo; k <= undo.hi; ++k) chain.reg[k] = undo.m[k - undo.lo];
}

void registry_init(Chain& chain, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    const int N = chain.N();
    chain.reg.assign(N, Vec3(1, 0, 0));
    int i = 0;
    while (i < N) {
        if (!is_helix(chain.state[i])) { ++i; continue; }
        const int hi = run_right(chain, i, chain.state[i]);
        chain.reg[i] = registry_from_angle(chain.tangent(i), 2.0 * M_PI * unif(rng));
        for (int k = i + 1; k <= hi; ++k)
            chain.reg[k] = registry_transport(chain.reg[k - 1], chain.tangent(k - 1), chain.tangent(k));
        i = hi + 1;
    }
}

void registry_run(const Chain& chain, int i, int& lo, int& hi) {
    const State s = chain.state[i];
    lo = hi = i;
    while (lo > 0 && chain.state[lo - 1] == s) --lo;
    while (hi + 1 < chain.N() && chain.state[hi + 1] == s) ++hi;
}

int registry_rederive(Chain& chain, int k) {
    if (!is_helix(chain.state[k])) return k;
    int lo, hi;
    registry_run(chain, k, lo, hi);
    for (int j = std::max(k, lo + 1); j <= hi; ++j)
        chain.reg[j] = registry_transport(chain.reg[j - 1], chain.tangent(j - 1), chain.tangent(j));
    return hi;
}

int registry_check(const Chain& chain, double tol) {
    int bad = 0;
    const int N = chain.N();
    for (int i = 0; i < N; ++i) {
        if (!is_helix(chain.state[i])) continue;
        const Vec3 u = chain.tangent(i);
        if (std::fabs(dot(u, chain.reg[i])) > tol || std::fabs(norm(chain.reg[i]) - 1.0) > tol) ++bad;
        if (i + 1 < N && chain.state[i + 1] == chain.state[i]) {
            const Vec3 t = registry_transport(chain.reg[i], u, chain.tangent(i + 1));
            if (norm(t - chain.reg[i + 1]) > tol) ++bad;
        }
    }
    return bad;
}
