#include "logging.h"
#include "pair_potential.h"

#include <cstdlib>
#include <string>

Logger::Logger(const Input& in) : in_(in) {
    const std::string obs_name = in_.out_prefix + "_obs.dat";
    obs_ = std::fopen(obs_name.c_str(), "w");
    if (!obs_) {
        std::fprintf(stderr, "error: cannot open '%s' for writing\n", obs_name.c_str());
        std::exit(EXIT_FAILURE);
    }
    if (in_.dump_every > 0) {
        const std::string conf_name = in_.out_prefix + "_conf.xyz";
        conf_ = std::fopen(conf_name.c_str(), "w");
        if (!conf_) {
            std::fprintf(stderr, "error: cannot open '%s' for writing\n", conf_name.c_str());
            std::exit(EXIT_FAILURE);
        }
    }
}

Logger::~Logger() {
    if (obs_)  std::fclose(obs_);
    if (conf_) std::fclose(conf_);
}

void Logger::header() {
    std::fprintf(obs_, "# %10s %6s %6s %10s %12s %12s %12s %12s %12s %10s %10s %10s\n",
                 "sweep", "n_R", "n_L", "helicity", "E_state", "E_bond", "E_bend", "Ree", "Rg2", "acc_state", "acc_pos", "acc_hinge");
    std::fflush(obs_);
}

void Logger::sample(long sweep, const Chain& chain, const MC& mc) {
    const int    nR  = chain.count(State::R);
    const int    nL  = chain.count(State::L);
    const double th  = chain.helicity();
    const double es  = total_state_energy(chain);
    const double eb  = total_bond_energy(chain);
    const double ek  = total_bend_energy(chain);
    const double ree = chain.end_to_end();
    const double rg2 = chain.rg2();

    std::fprintf(obs_, "  %10ld %6d %6d %10.6f %12.6f %12.6f %12.6f %12.6f %12.6f %10.4f %10.4f %10.4f\n",
                 sweep, nR, nL, th, es, eb, ek, ree, rg2, mc.acc_state(), mc.acc_pos(), mc.acc_hinge());

    std::fflush(obs_);   // keep the file usable if a long run is killed

    ++n_samples_;
    sum_helix_ += th;
    sum_estate_ += es;
    sum_ebond_ += eb;
    sum_ebend_ += ek;
    sum_ree_   += ree;
    sum_rg2_   += rg2;
}

void Logger::dump_config(long sweep, const Chain& chain) {
    if (!conf_) return;
    // extended-xyz style frame: element = C (coil), R or L (helix sense)
    std::fprintf(conf_, "%d\nsweep=%ld\n", chain.N(), sweep);
    for (int i = 0; i < chain.N(); ++i) {
        std::fprintf(conf_, "%s %14.8f %14.8f %14.8f\n",
                     state_name(chain.state[i]), chain.pos[i].x, chain.pos[i].y, chain.pos[i].z);
    }
    std::fflush(conf_);
}

void Logger::summary(const Chain& chain, const MC& mc) {
    std::fflush(obs_);
    std::printf("# ---- summary ----\n");
    std::printf("# samples        = %ld\n", n_samples_);
    if (n_samples_ > 0) {
        std::printf("# <helicity>     = %.6f\n", sum_helix_ / n_samples_);
        std::printf("# <E_state>      = %.6f\n", sum_estate_ / n_samples_);
        std::printf("# <E_bond>       = %.6f\n", sum_ebond_ / n_samples_);
        std::printf("# <E_bend>       = %.6f\n", sum_ebend_ / n_samples_);
        std::printf("# <Ree>          = %.6f\n", sum_ree_   / n_samples_);
        std::printf("# <Rg2>          = %.6f\n", sum_rg2_   / n_samples_);
    }
    std::printf("# final n_R/n_L  = %d / %d  (N = %d)\n",
                chain.count(State::R), chain.count(State::L), chain.N());
    std::printf("# acc_state      = %.4f\n", mc.acc_state());
    std::printf("# acc_pos        = %.4f\n", mc.acc_pos());
    std::printf("# acc_hinge      = %.4f\n", mc.acc_hinge());
    std::printf("# -----------------\n");
}
