#include "input.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

namespace {

// strip leading/trailing whitespace
std::string trim(const std::string& s) {
    const char* ws = " \t\r\n";
    const auto b = s.find_first_not_of(ws);
    if (b == std::string::npos) return "";
    const auto e = s.find_last_not_of(ws);
    return s.substr(b, e - b + 1);
}

int pair_index(const std::string& suf) {
    for (int i = 0; i < N_PAIRS; ++i)
        if (suf == pair_name(static_cast<Pair>(i))) return i;
    return -1;
}
int triple_index(const std::string& suf) {
    for (int i = 0; i < N_TRIPLES; ++i)
        if (suf == triple_name(static_cast<Triple>(i))) return i;
    return -1;
}

// parse val into field; returns false on failure
template <typename T>
bool parse(const std::string& val, T& field) {
    std::istringstream ss(val);
    return static_cast<bool>(ss >> field);
}

// assign one key/value pair; returns false on unknown key or bad value
bool assign(Input& in, const std::string& key, const std::string& val) {
    // keyed parameters: <base>_<suffix>
    const auto us = key.rfind('_');
    if (us != std::string::npos) {
        const std::string base = key.substr(0, us);
        const std::string suf  = key.substr(us + 1);
        if (base == "bond_len" || base == "k_bond") {
            const int i = pair_index(suf);
            if (i < 0) return false;
            return base == "bond_len" ? parse(val, in.bond[i].bond_len) : parse(val, in.bond[i].k_bond);
        }
        if (base == "kappa" || base == "theta0") {
            const int i = triple_index(suf);
            if (i < 0) return false;
            if (base == "kappa") return parse(val, in.bend[i].kappa);
            // theta0 accepts radians, or degrees with a "deg" suffix:  theta0_HHH = 135 deg
            if (!parse(val, in.bend[i].theta0)) return false;
            if (val.find("deg") != std::string::npos) in.bend[i].theta0 *= M_PI / 180.0;
            return true;
        }
    }

    if (key == "N")          return parse(val, in.N);
    if (key == "init")       { in.init = val; return true; }
    if (key == "init_state") { in.init_state = val; return true; }
    if (key == "n_states")   return parse(val, in.n_states);
    if (key == "bend_key")   { in.bend_key = val; return true; }
    if (key == "J0")         return parse(val, in.J0);
    if (key == "J1")         return parse(val, in.J1);
    if (key == "J2")         return parse(val, in.J2);
    if (key == "kT")         return parse(val, in.kT);
    if (key == "n_equil")    return parse(val, in.n_equil);
    if (key == "n_sweeps")   return parse(val, in.n_sweeps);
    if (key == "max_disp")   return parse(val, in.max_disp);
    if (key == "n_hinge")    return parse(val, in.n_hinge);
    if (key == "n_pivot")    return parse(val, in.n_pivot);
    if (key == "max_rot")    return parse(val, in.max_rot);
    if (key == "nb_type")    { in.nb_type = val; return true; }
    if (key == "nb_A")       return parse(val, in.nb_A);
    if (key == "nb_sigma")   return parse(val, in.nb_sigma);
    if (key == "nb_rcut")    return parse(val, in.nb_rcut);
    if (key == "nb_hh")      { in.nb_hh = val; return true; }
    if (key == "gb_eps0")    return parse(val, in.gb_eps0);
    if (key == "gb_aniso_eps") return parse(val, in.gb_aniso_eps);
    if (key == "gb_mu")      return parse(val, in.gb_mu);
    if (key == "gb_nu")      return parse(val, in.gb_nu);
    if (key == "gb_kappa_p") return parse(val, in.gb_kappa_p);
    if (key == "rod_L")      return parse(val, in.rod_L);
    if (key == "rod_r")      return parse(val, in.rod_r);
    if (key == "seed")       return parse(val, in.seed);
    if (key == "log_every")  return parse(val, in.log_every);
    if (key == "dump_every") return parse(val, in.dump_every);
    if (key == "out_prefix") { in.out_prefix = val; return true; }
    return false;
}

} // namespace

bool read_input(const std::string& filename, Input& in) {
    std::ifstream f(filename);
    if (!f) return false;

    std::string line;
    int lineno = 0;
    bool ok = true;
    while (std::getline(f, line)) {
        ++lineno;
        // drop comments
        const auto hash = line.find('#');
        if (hash != std::string::npos) line.erase(hash);
        line = trim(line);
        if (line.empty()) continue;

        const auto eq = line.find('=');
        if (eq == std::string::npos) {
            std::fprintf(stderr, "input: line %d: expected 'key = value'\n", lineno);
            ok = false;
            continue;
        }
        const std::string key = trim(line.substr(0, eq));
        const std::string val = trim(line.substr(eq + 1));
        if (!assign(in, key, val)) {
            std::fprintf(stderr, "input: line %d: unknown key or bad value '%s = %s'\n",
                         lineno, key.c_str(), val.c_str());
            ok = false;
        }
    }

    // basic sanity checks
    if (in.N < 2)          { std::fprintf(stderr, "input: N must be >= 2\n");        ok = false; }
    if (in.kT <= 0)        { std::fprintf(stderr, "input: kT must be > 0\n");        ok = false; }
    if (in.init != "rod" && in.init != "walk") {
        std::fprintf(stderr, "input: init must be 'rod' or 'walk'\n"); ok = false;
    }
    if (in.init_state != "coil" && in.init_state != "random") {
        std::fprintf(stderr, "input: init_state must be 'coil' or 'random'\n"); ok = false;
    }
    if (in.n_states != 2 && in.n_states != 3) {
        std::fprintf(stderr, "input: n_states must be 2 or 3\n"); ok = false;
    }
    if (in.n_states == 2 && in.init_state != "random") {
        std::fprintf(stderr, "input: n_states = 2 requires init_state = random\n"); ok = false;
    }
    if (in.nb_type != "none" && in.nb_type != "gauss" && in.nb_type != "wca") {
        std::fprintf(stderr, "input: nb_type must be none, gauss or wca\n"); ok = false;
    }
    if (in.nb_hh != "same" && in.nb_hh != "gb") {
        std::fprintf(stderr, "input: nb_hh must be same or gb\n"); ok = false;
    }
    if (in.nb_hh == "gb" && in.nb_type == "none") {
        std::fprintf(stderr, "input: nb_hh = gb needs nb_type = gauss or wca for the coil pairs\n"); ok = false;
    }
    if (in.nb_type != "none" && in.n_hinge > 0) {
        std::fprintf(stderr, "input: the hinge move (n_hinge > 0) is only valid without non-bonded interactions\n"); ok = false;
    }
    if (in.nb_rcut <= 0.0) in.nb_rcut = (in.nb_type == "wca") ? std::pow(2.0, 1.0 / 6.0) * in.nb_sigma : 3.5 * in.nb_sigma;
    if (in.bend_key != "triple" && in.bend_key != "pair") {
        std::fprintf(stderr, "input: bend_key must be 'triple' or 'pair'\n"); ok = false;
    }
    if (in.log_every <= 0) { std::fprintf(stderr, "input: log_every must be > 0\n"); ok = false; }
    for (int c = 0; c < N_PAIRS; ++c) {
        if (in.bond[c].bond_len <= 0) {
            std::fprintf(stderr, "input: bond_len_%s must be > 0\n", pair_name(static_cast<Pair>(c)));
            ok = false;
        }
    }
    return ok;
}

void print_input(const Input& in) {
    std::printf("# ---- input ----\n");
    std::printf("# N           = %d\n", in.N);
    std::printf("# init        = %s\n", in.init.c_str());
    std::printf("# init_state  = %s\n", in.init_state.c_str());
    std::printf("# n_states    = %d\n", in.n_states);
    std::printf("# bend_key    = %s\n", in.bend_key.c_str());
    std::printf("# %-6s %10s %10s\n", "bond", "bond_len", "k_bond");
    for (int c = 0; c < N_PAIRS; ++c) {
        const BondParams& b = in.bond[c];
        std::printf("# %-6s %10g %10g\n", pair_name(static_cast<Pair>(c)), b.bond_len, b.k_bond);
    }
    std::printf("# %-6s %10s %10s\n", "bend", "kappa", "theta0");
    for (int t = 0; t < N_TRIPLES; ++t) {
        const BendParams& b = in.bend[t];
        std::printf("# %-6s %10g %10g\n", triple_name(static_cast<Triple>(t)), b.kappa, b.theta0);
    }
    std::printf("# J0 (E_HH=-J0) = %g\n", in.J0);
    std::printf("# J1 (E_CH=+J1) = %g\n", in.J1);
    std::printf("# J2 (E_RL=+J2) = %g\n", in.J2);
    std::printf("# kT          = %g\n",   in.kT);
    std::printf("# n_equil     = %ld\n",  in.n_equil);
    std::printf("# n_sweeps    = %ld\n",  in.n_sweeps);
    std::printf("# max_disp    = %g\n",   in.max_disp);
    std::printf("# n_hinge     = %d\n",   in.n_hinge);
    std::printf("# n_pivot     = %d  (max_rot = %g)\n", in.n_pivot, in.max_rot);
    std::printf("# nb_type     = %s  (A = %g, sigma = %g, rcut = %g)\n", in.nb_type.c_str(), in.nb_A, in.nb_sigma, in.nb_rcut);
    std::printf("# nb_hh       = %s  (eps0 = %g, aniso_eps = %d, mu = %g, nu = %g, kappa' = %g)\n",
                in.nb_hh.c_str(), in.gb_eps0, in.gb_aniso_eps, in.gb_mu, in.gb_nu, in.gb_kappa_p);
    std::printf("# rod         = L %g, r %g  (D = %g, L/D = %g)\n", in.rod_L, in.rod_r, 2 * in.rod_r, in.rod_L / (2 * in.rod_r));
    std::printf("# seed        = %lu\n",  in.seed);
    std::printf("# log_every   = %ld\n",  in.log_every);
    std::printf("# dump_every  = %ld\n",  in.dump_every);
    std::printf("# out_prefix  = %s\n",   in.out_prefix.c_str());
    std::printf("# ---------------\n");
}
