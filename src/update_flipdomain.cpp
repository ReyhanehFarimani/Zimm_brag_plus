#include "update_flipdomain.h"

#include <vector>

#include <algorithm>
#include <cmath>
#include "pair_potential.h"
#include "registry.h"

namespace {
// flip the handedness of every residue in [a, b] (all helical): the shared core of the single-site, whole-domain
// and segment flips. i_single >= 0 marks the single-site move (its registry is redrawn), otherwise the junction
// twists inside [a, b] are negated (an involution).
bool flip_range(Chain& chain, int a, int b, int i_single, const Input& in, std::mt19937_64& rng) {
    const int N = chain.N();
    const int i = i_single >= 0 ? i_single : a;
    // registries (nb_hh = db): a single-site flip redraws m_i uniformly (symmetric proposal); a
    // whole-domain flip negates the domain's junction twists (an involution, so the twist term is
    // invariant and the flip is judged by the pair potential); the twist term of the edge
    // junctions (a-1, a), (b, b+1) and, for a single site, of (i-1, i), (i, i+1) is included
    const bool reg = in.nb_hh == "db";
    const bool mixed = i_single < 0;                              // a segment may hold both hands: flip each

    // energies that can change: the edge bonds (a-1, a) and (b, b+1) and the bends at a-1, a, b, b+1;
    // with a handedness-dependent non-bonded potential (nb_hh = fit) also every non-bonded pair
    // involving a domain bead (in-domain pairs counted once)
    auto local = [&]() {
        double e = 0.0;
        if (a > 0)     e += state_energy(chain, a - 1) + bond_energy(chain, a - 1);
        if (b < N - 1) e += state_energy(chain, b)     + bond_energy(chain, b);
        int js[4] = {a - 1, a, b, b + 1};
        for (int k = 0; k < 4; ++k) {
            const int j = js[k];
            if (j < 1 || j > N - 2) continue;
            if (k > 0 && j == js[k - 1]) continue;                 // avoid double counting when the domain is short
            e += bend_energy(chain, j);
        }
        if (reg) e += twist_range_energy(chain, a, b);
        if (nb_chiral(in)) {
            if (reg) { e += nb_range_energy(chain, a, b); return e; }
            const bool list = nl_ready(chain);                     // the geometry does not change in this move
            for (int k = a; k <= b; ++k) {
                if (list) {
                    for (int j : chain.nl.nbrs[k]) {
                        if (j >= a && j <= b) { if (j > k + 1) e += nb_pair_energy(chain, k, j); continue; }
                        e += nb_pair_energy(chain, k, j);
                    }
                    continue;
                }
                for (int j = 0; j < N; ++j) {
                    if (j >= a && j <= b) { if (j > k + 1) e += nb_pair_energy(chain, k, j); continue; }
                    if (j < k - 1 || j > k + 1) e += nb_pair_energy(chain, k, j);
                }
            }
        }
        return e;
    };
    const double e_old = local();
    std::vector<Vec3> reg_old;
    if (reg) reg_old.assign(chain.reg.begin() + a, chain.reg.begin() + b + 1);
    std::vector<State> old_state(chain.state.begin() + a, chain.state.begin() + b + 1);
    for (int k = a; k <= b; ++k) chain.state[k] = chain.state[k] == State::R ? State::L : State::R;
    if (reg) {
        if (i_single >= 0) registry_resample(chain, i, rng);
        else               registry_negate_twists(chain, a, b);
    }
    (void)mixed;
    const double dE = local() - e_old;

    std::uniform_real_distribution<double> unif(0.0, 1.0);
    g_last_dE = dE;
    if (dE <= 0.0 || unif(rng) < std::exp(-dE / in.kT)) return true;
    std::copy(old_state.begin(), old_state.end(), chain.state.begin() + a);   // rejected: restore
    if (reg) std::copy(reg_old.begin(), reg_old.end(), chain.reg.begin() + a);
    return false;
}
}  // namespace

bool try_domain_flip(Chain& chain, int i, const Input& in, std::mt19937_64& rng) {
    const State s = chain.state[i];
    if (!is_helix(s)) return false;
    std::uniform_int_distribution<int> coin(0, 1);
    int a = i, b = i;
    if (coin(rng)) {                                              // (b) whole maximal domain
        while (chain.has_bond(a - 1) && chain.state[a - 1] == s) --a;
        while (chain.has_bond(b) && chain.state[b + 1] == s) ++b;
        // a flip that merges with a helical neighbour cannot be reversed by the same move: reject
        if (chain.has_bond(a - 1) && is_helix(chain.state[a - 1])) return false;
        if (chain.has_bond(b) && is_helix(chain.state[b + 1])) return false;
        return flip_range(chain, a, b, -1, in, rng);
    }
    return flip_range(chain, i, i, i, in, rng);                   // (a) single site
}

// Segment flip (2026-09-24, user: "flip can be not for full chain but segments as well"): the residues
// i .. i + len - 1 of i's arm, len uniform in 1..segflip_len_max, all helical (any mix of hands), flip their
// handedness; the segment is chosen independently of the domain structure, so the reverse move is proposed with
// the same probability and no merge rule is needed.
bool try_segment_flip(Chain& chain, int i, const Input& in, std::mt19937_64& rng) {
    std::uniform_int_distribution<int> lenu(1, std::max(1, in.segflip_len_max));
    const int len = lenu(rng), b = std::min(i + len - 1, chain.arm_last(i));
    for (int k = i; k <= b; ++k) if (!is_helix(chain.state[k])) return false;
    return flip_range(chain, i, b, -1, in, rng);
}
