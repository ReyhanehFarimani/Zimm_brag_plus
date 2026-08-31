// Zimm-Bragg helix-coil chain with sterics, bonding and state-dependent bending.
// Monte Carlo driver: read input, build chain, run sweeps, log observables.

#include <cstdio>
#include <cstdlib>
#include <string>

#include "input.h"
#include "chain.h"
#include "MC.h"
#include "logging.h"

int main(int argc, char** argv) {
    // Usage: ./zimm [input_file]
    const std::string input_file = (argc > 1) ? argv[1] : "input.dat";

    Input in;
    if (!read_input(input_file, in)) {
        std::fprintf(stderr, "error: could not read or parse input file '%s'\n", input_file.c_str());
        return EXIT_FAILURE;
    }
    print_input(in);

    // Chain holds per-residue state (helix/coil), positions and orientations.
    Chain chain(in);
    chain.init();

    // MC owns the RNG and the three move types (state flip, position, orientation).
    MC mc(in, chain);

    Logger log(in);
    log.header();

    // Equilibration: no output.
    for (long sweep = 0; sweep < in.n_equil; ++sweep) {
        mc.sweep();
    }
    mc.reset_acceptance();

    // Production.
    for (long sweep = 0; sweep < in.n_sweeps; ++sweep) {
        mc.sweep();
        if (sweep % in.log_every == 0) {
            log.sample(sweep, chain, mc);
        }
        if (in.dump_every > 0 && sweep % in.dump_every == 0) {
            log.dump_config(sweep, chain);
        }
    }

    log.summary(chain, mc);
    return EXIT_SUCCESS;
}
