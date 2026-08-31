#include "input.h"

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
            return base == "kappa" ? parse(val, in.bend[i].kappa) : parse(val, in.bend[i].theta0);
        }
    }

    if (key == "N")          return parse(val, in.N);
    if (key == "init")       { in.init = val; return true; }
    if (key == "s_R")        return parse(val, in.s_R);
    if (key == "s_L")        return parse(val, in.s_L);
    if (key == "sigma_zb")   return parse(val, in.sigma_zb);
    if (key == "kT")         return parse(val, in.kT);
    if (key == "n_equil")    return parse(val, in.n_equil);
    if (key == "n_sweeps")   return parse(val, in.n_sweeps);
    if (key == "max_disp")   return parse(val, in.max_disp);
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
    std::printf("# s_R         = %g\n",   in.s_R);
    std::printf("# s_L         = %g\n",   in.s_L);
    std::printf("# sigma_zb    = %g\n",   in.sigma_zb);
    std::printf("# kT          = %g\n",   in.kT);
    std::printf("# n_equil     = %ld\n",  in.n_equil);
    std::printf("# n_sweeps    = %ld\n",  in.n_sweeps);
    std::printf("# max_disp    = %g\n",   in.max_disp);
    std::printf("# seed        = %lu\n",  in.seed);
    std::printf("# log_every   = %ld\n",  in.log_every);
    std::printf("# dump_every  = %ld\n",  in.dump_every);
    std::printf("# out_prefix  = %s\n",   in.out_prefix.c_str());
    std::printf("# ---------------\n");
}
