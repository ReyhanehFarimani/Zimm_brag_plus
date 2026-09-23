// Consistency test of the move energies with registries (nb_hh = db): for every accepted move the
// change of the TOTAL energy (all terms recomputed from scratch) must equal the dE the move used.
//   g++ -std=c++11 -O2 -Isrc tools/test_moves.cpp src/{chain,input,pair_potential,registry,update_*,logging,MC}.cpp -o tools/test_moves
//   tools/test_moves runs/db_test/mech.dat
#include <cmath>
#include <cstdio>
#include <random>
#include <cstdlib>
#include "input.h"
#include "chain.h"
#include "pair_potential.h"
#include "registry.h"
#include "update_state.h"
#include "update_position.h"
#include "update_pivot.h"
#include "update_flipdomain.h"
#include "update_twist.h"

static double total(const Chain& c) {
    return total_state_energy(c) + total_bond_energy(c) + total_bend_energy(c) + total_nb_energy(c);
}

int main(int argc, char** argv) {
    Input in;
    if (!read_input(argc > 1 ? argv[1] : "input.dat", in)) return 1;
    in.nl_skin = (argc > 2) ? atof(argv[2]) : 1.0;
    Chain chain(in); chain.init();
    std::mt19937_64 rng(7);
    std::uniform_int_distribution<int> pick(0, in.N - 1), pick_inner(1, in.N - 2);
    const char* names[5] = {"state", "position", "pivot", "flip", "twist"};
    long n_acc[5] = {0, 0, 0, 0, 0}, n_bad[5] = {0, 0, 0, 0, 0};
    double worst[5] = {0, 0, 0, 0, 0};
    for (int it = 0; it < 6000; ++it) {
        const int kind = it % 5;
        const double e0 = total(chain);
        bool acc = false;
        switch (kind) {
        case 0: acc = try_state_move(chain, pick(rng), in, rng); break;
        case 1: acc = try_position_move(chain, pick(rng), in, rng); break;
        case 2: acc = try_pivot_move(chain, pick_inner(rng), in, rng); break;
        case 3: acc = try_domain_flip(chain, pick(rng), in, rng); break;
        case 4: acc = try_twist_move(chain, pick(rng), in, rng); break;
        }
        if (!acc) continue;
        ++n_acc[kind];
        const double d = total(chain) - e0 - g_last_dE;
        if (std::fabs(d) > 1e-8 * (1.0 + std::fabs(g_last_dE))) { ++n_bad[kind]; if (std::fabs(d) > worst[kind]) worst[kind] = std::fabs(d); }
        if (registry_check(chain)) { std::printf("registry invariant violated after %s\n", names[kind]); return 2; }
    }
    for (int k = 0; k < 5; ++k)
        std::printf("%-9s accepted %5ld  dE mismatches %ld  (worst %.2e)\n", names[k], n_acc[k], n_bad[k], worst[k]);
    return 0;
}
