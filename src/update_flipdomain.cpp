#include "update_flipdomain.h"

#include <algorithm>
#include <cmath>
#include "pair_potential.h"

bool try_domain_flip(Chain& chain, int i, const Input& in, std::mt19937_64& rng) {
    const State s = chain.state[i];
    if (!is_helix(s)) return false;
    const int N = chain.N();
    std::uniform_int_distribution<int> coin(0, 1);
    int a = i, b = i;
    if (coin(rng)) {                                              // (b) whole maximal domain
        while (a > 0 && chain.state[a - 1] == s) --a;
        while (b < N - 1 && chain.state[b + 1] == s) ++b;
        // a flip that merges with a helical neighbour cannot be reversed by the same move: reject
        if (a > 0 && is_helix(chain.state[a - 1])) return false;
        if (b < N - 1 && is_helix(chain.state[b + 1])) return false;
    }                                                             // else (a) single site: a = b = i

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
        if (nb_chiral(in)) {
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
    const State t = (s == State::R) ? State::L : State::R;
    for (int k = a; k <= b; ++k) chain.state[k] = t;
    const double dE = local() - e_old;

    std::uniform_real_distribution<double> unif(0.0, 1.0);
    if (dE <= 0.0 || unif(rng) < std::exp(-dE / in.kT)) return true;
    for (int k = a; k <= b; ++k) chain.state[k] = s;              // rejected: restore
    return false;
}
