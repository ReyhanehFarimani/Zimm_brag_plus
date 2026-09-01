#include "update_hinge.h"

#include <cmath>
#include <map>
#include "pair_potential.h"

namespace {

// ln Z_hinge per bend-parameter set, computed once (parameters are fixed for a run)
double log_weight(const BendParams& p) {
    static std::map<std::pair<double, double>, double> cache;
    const auto key = std::make_pair(p.kappa, p.theta0);
    auto it = cache.find(key);
    if (it == cache.end()) it = cache.emplace(key, hinge_log_weight(p)).first;
    return it->second;
}

// draw a valence angle from p(theta) ~ sin(theta) exp(-U(theta)) on [0, pi]
double draw_angle(const BendParams& p, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    if (p.kappa <= 0.0) return std::acos(1.0 - 2.0 * unif(rng));           // free hinge: uniform on the sphere
    std::normal_distribution<double> gauss(p.theta0, 1.0 / std::sqrt(p.kappa));
    for (int it = 0; it < 100000; ++it) {                                   // Gaussian proposal, accept with sin(theta)
        const double th = gauss(rng);
        if (th < 0.0 || th > M_PI) continue;
        if (unif(rng) < std::sin(th)) return th;
    }
    return p.theta0;                                                        // (unreachable in practice)
}

// Redraw bond j (j >= 1) so that the valence angle at bead j is theta, azimuth random; rotate beads j+1.. rigidly.
void resample_bond(Chain& chain, int j, double theta, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    const Vec3 u = (1.0 / norm(chain.bond(j - 1))) * chain.bond(j - 1);     // incoming direction
    // orthonormal frame (e1, e2) perpendicular to u
    Vec3 t = std::fabs(u.x) < 0.9 ? Vec3(1, 0, 0) : Vec3(0, 1, 0);
    Vec3 e1 = t - dot(t, u) * u; e1 = (1.0 / norm(e1)) * e1;
    const Vec3 e2 = cross(u, e1);
    const double alpha = M_PI - theta;                                      // angle between consecutive bond vectors
    const double phi   = 2.0 * M_PI * unif(rng);
    const Vec3 v_new = std::cos(alpha) * u + std::sin(alpha) * (std::cos(phi) * e1 + std::sin(phi) * e2);

    const Vec3 b_old = chain.bond(j);
    const double len = norm(b_old);
    const Vec3 v_old = (1.0 / len) * b_old;
    // rotation taking v_old -> v_new (minimal rotation), applied to every bead downstream of bead j
    Vec3 axis = cross(v_old, v_new);
    const double sa = norm(axis), ca = std::max(-1.0, std::min(1.0, dot(v_old, v_new)));
    const Vec3 pivot = chain.pos[j];
    if (sa < 1e-12) {
        if (ca > 0) return;                                                 // already aligned
        axis = e1;                                                          // antiparallel: rotate by pi about e1
    } else {
        axis = (1.0 / sa) * axis;
    }
    const double ang = std::atan2(sa, ca);
    for (int k = j + 1; k < chain.N(); ++k)
        chain.pos[k] = pivot + rotate(chain.pos[k] - pivot, axis, ang);
}

} // namespace

bool try_hinge_move(Chain& chain, int i, const Input& in, std::mt19937_64& rng) {
    std::uniform_int_distribution<int>     coin(0, 1);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    const int N = chain.N();

    // proposed new state (same proposal rules as the plain state move)
    int q_new;
    if (in.n_states == 2) {
        q_new = -spin(chain.state[i]);
    } else {
        q_new = spin(chain.state[i]) + (coin(rng) ? 1 : -1);
        if (q_new < -1 || q_new > 1) return false;
    }
    const State old_s = chain.state[i];

    // spin + bond energy before (bond parameters may depend on the pair class)
    double e_old = 0.0;
    if (i > 0)     e_old += state_energy(chain, i - 1) + bond_energy(chain, i - 1);
    if (i < N - 1) e_old += state_energy(chain, i)     + bond_energy(chain, i);
    // hinges in the window and their old parameters
    BendParams old_par[3]; bool has[3] = {false, false, false};
    for (int m = 0; m < 3; ++m) {
        const int j = i - 1 + m;
        if (j >= 1 && j <= N - 2) { has[m] = true; old_par[m] = chain.bend_par(j); }
    }

    chain.state[i] = from_spin(q_new);

    double e_new = 0.0;
    if (i > 0)     e_new += state_energy(chain, i - 1) + bond_energy(chain, i - 1);
    if (i < N - 1) e_new += state_energy(chain, i)     + bond_energy(chain, i);
    double log_ratio = -(e_new - e_old) / in.kT;
    for (int m = 0; m < 3; ++m) {
        if (!has[m]) continue;
        const BendParams& np = chain.bend_par(i - 1 + m);
        if (np.kappa == old_par[m].kappa && np.theta0 == old_par[m].theta0) continue;   // unchanged hinge
        log_ratio += log_weight(np) - log_weight(old_par[m]);
    }

    if (log_ratio < 0.0 && unif(rng) >= std::exp(log_ratio)) {
        chain.state[i] = old_s;                                             // rejected
        return false;
    }
    // accepted: redraw every changed hinge from its new distribution (upstream first)
    for (int m = 0; m < 3; ++m) {
        if (!has[m]) continue;
        const int j = i - 1 + m;
        const BendParams& np = chain.bend_par(j);
        if (np.kappa == old_par[m].kappa && np.theta0 == old_par[m].theta0) continue;
        resample_bond(chain, j, draw_angle(np, rng), rng);
    }
    return true;
}
