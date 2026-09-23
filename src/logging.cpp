#include "logging.h"
#include "pair_potential.h"
#include "registry.h"

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
    std::fprintf(obs_, "# %10s %6s %6s %10s %12s %12s %12s %12s %12s %12s %12s %10s %10s %10s %10s %10s\n",
                 "sweep", "n_R", "n_L", "helicity", "E_state", "E_bond", "E_bend", "E_nb", "E_twist", "Ree", "Rg2", "acc_state", "acc_pos", "acc_hinge", "acc_pivot", "acc_flip");
    if (in_.n_arms > 0) std::fprintf(obs_, "# (star: E_bond includes the core wall and graft tethers; Ree = first-to-last residue of arm 0)\n");
    std::fflush(obs_);
}

void Logger::sample(long sweep, const Chain& chain, const MC& mc) {
    // guard: the neighbour list must contain every interacting pair (E_nb below is an all-pairs sum)
    if (!nl_verify(chain)) {
        std::fprintf(stderr, "error: neighbour list inconsistent at sweep %ld\n", sweep);
        std::exit(EXIT_FAILURE);
    }
    const int    nR  = chain.count(State::R);
    const int    nL  = chain.count(State::L);
    const double th  = chain.helicity();
    const double es  = total_state_energy(chain);
    const double eb  = total_bond_energy(chain) + total_core_energy(chain);
    const double ek  = total_bend_energy(chain);
    const double en  = total_nb_energy(chain);
    const double et  = total_twist_energy(chain);
    const double ree = chain.end_to_end();
    const double rg2 = chain.rg2();

    std::fprintf(obs_, "  %10ld %6d %6d %10.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
                 sweep, nR, nL, th, es, eb, ek, en, et, ree, rg2, mc.acc_state(), mc.acc_pos(), mc.acc_hinge(), mc.acc_pivot(), mc.acc_flip());

    std::fflush(obs_);   // keep the file usable if a long run is killed

    ++n_samples_;
    sum_helix_ += th;
    sum_estate_ += es;
    sum_ebond_ += eb;
    sum_ebend_ += ek;
    sum_enb_   += en;
    sum_etwist_ += et;
    sum_ree_   += ree;
    sum_rg2_   += rg2;
}

void Logger::dump_config(long sweep, const Chain& chain) {
    if (!conf_) return;
    // extended-xyz frame (OVITO-readable): species = C (coil), R or L (helix sense); position;
    // orientation = unit quaternion (x, y, z, w) rotating the z axis onto the local tangent -- for EVERY
    // residue, coils included (user 2026-09-23; before, coils carried the identity); tangent = the
    // tangent itself; registry = the registry direction m_i (unit, perpendicular to the tangent; a
    // coil's is the dummy that enters no energy); rod_L, rod_D = rod length and diameter (helix), 0 and
    // sphere diameter (coil); spin (-1/0/+1)
    const int extra = chain.n_arms() > 0 ? 1 : 0;                    // the core as one more particle (species X)
    std::fprintf(conf_, "%d\nProperties=species:S:1:pos:R:3:orientation:R:4:tangent:R:3:registry:R:3:rod_L:R:1:rod_D:R:1:spin:I:1 sweep=%ld\n",
                 chain.N() + extra, sweep);
    if (extra) std::fprintf(conf_, "X %14.8f %14.8f %14.8f %10.6f %10.6f %10.6f %10.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %6.3f %6.3f %2d\n",
                            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 2.0 * in_.core_radius, 0);
    for (int i = 0; i < chain.N(); ++i) {
        const bool helix = is_helix(chain.state[i]);
        double q[4] = {0.0, 0.0, 0.0, 1.0};
        chain.quaternion(i, q);
        const Vec3 t = chain.tangent(i);
        const Vec3 m = (i < (int)chain.reg.size()) ? chain.reg[i] : Vec3(0, 0, 0);
        std::fprintf(conf_, "%s %14.8f %14.8f %14.8f %10.6f %10.6f %10.6f %10.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %6.3f %6.3f %2d\n",
                     state_name(chain.state[i]), chain.pos[i].x, chain.pos[i].y, chain.pos[i].z,
                     q[0], q[1], q[2], q[3], t.x, t.y, t.z, m.x, m.y, m.z,
                     helix ? in_.rod_L : 0.0, 2.0 * in_.rod_r, spin(chain.state[i]));
    }
    std::fflush(conf_);
}

void Logger::summary(const Chain& chain, const MC& mc) {
    std::fflush(obs_);
    std::printf("# ---- summary ----\n");
    std::printf("# samples        = %ld\n", n_samples_);
    if (in_.n_twist > 0) std::printf("# acc_twist      = %.4f\n", mc.acc_twist());
    if (n_samples_ > 0) {
        std::printf("# <helicity>     = %.6f\n", sum_helix_ / n_samples_);
        std::printf("# <E_state>      = %.6f\n", sum_estate_ / n_samples_);
        std::printf("# <E_bond>       = %.6f\n", sum_ebond_ / n_samples_);
        std::printf("# <E_bend>       = %.6f\n", sum_ebend_ / n_samples_);
        std::printf("# <E_nb>         = %.6f\n", sum_enb_   / n_samples_);
        if (in_.nb_hh == "db") std::printf("# <E_twist>      = %.6f\n", sum_etwist_ / n_samples_);
        std::printf("# <Ree>          = %.6f\n", sum_ree_   / n_samples_);
        std::printf("# <Rg2>          = %.6f\n", sum_rg2_   / n_samples_);
    }
    std::printf("# final n_R/n_L  = %d / %d  (N = %d)\n",
                chain.count(State::R), chain.count(State::L), chain.N());
    std::printf("# acc_state      = %.4f\n", mc.acc_state());
    std::printf("# acc_pos        = %.4f\n", mc.acc_pos());
    std::printf("# acc_hinge      = %.4f\n", mc.acc_hinge());
    std::printf("# acc_pivot      = %.4f\n", mc.acc_pivot());
    std::printf("# acc_flip       = %.4f\n", mc.acc_flip());
    std::printf("# nl_builds      = %ld\n", chain.nl.n_build);
    std::printf("# -----------------\n");
}
