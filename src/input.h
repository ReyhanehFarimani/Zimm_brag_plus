#pragma once
#include <string>

// Residue state: coil, right-handed helix, left-handed helix.
enum class State : int { Coil = 0, R = 1, L = 2 };
inline const char* state_name(State s) { return s == State::Coil ? "C" : (s == State::R ? "R" : "L"); }
inline bool is_helix(State s) { return s != State::Coil; }

// ---------------------------------------------------------------------------
// Pair class of two bonded residues; bond parameters are keyed by this.
//   CC : coil - coil
//   CH : coil - helix (either sense, either order)
//   HH : helix - helix, same sense (RR or LL)
//   RL : helix - helix, opposite sense (RL or LR)
enum class Pair : int { CC = 0, CH = 1, HH = 2, RL = 3 };
constexpr int N_PAIRS = 4;
inline const char* pair_name(Pair p) {
    static const char* names[N_PAIRS] = {"CC", "CH", "HH", "RL"};
    return names[static_cast<int>(p)];
}
inline Pair pair_class(State a, State b) {
    if (!is_helix(a) && !is_helix(b)) return Pair::CC;
    if (!is_helix(a) || !is_helix(b)) return Pair::CH;
    return a == b ? Pair::HH : Pair::RL;
}

// ---------------------------------------------------------------------------
// Triple class of three consecutive residues; bend parameters are keyed by this.
// Classes are symmetric under reversal (abc == cba) and mirror (R <-> L).
// 'H' is used when all helix residues in the triple share one sense; explicit
// R/L is used only when two senses are present.
//   CCC
//   CCH : CCR CCL RCC LCC
//   CHC : CRC CLC
//   CHH : CRR CLL RRC LLC
//   CRL : CRL CLR LRC RLC
//   HCH : RCR LCL
//   RCL : RCL LCR
//   HHH : RRR LLL
//   RRL : RRL LLR LRR RLL
//   RLR : RLR LRL
enum class Triple : int { CCC = 0, CCH, CHC, CHH, CRL, HCH, RCL, HHH, RRL, RLR };
constexpr int N_TRIPLES = 10;
inline const char* triple_name(Triple t) {
    static const char* names[N_TRIPLES] =
        {"CCC", "CCH", "CHC", "CHH", "CRL", "HCH", "RCL", "HHH", "RRL", "RLR"};
    return names[static_cast<int>(t)];
}
inline Triple triple_class(State a, State b, State c) {
    const int nh = is_helix(a) + is_helix(b) + is_helix(c);
    switch (nh) {
    case 0:  return Triple::CCC;
    case 1:  return is_helix(b) ? Triple::CHC : Triple::CCH;
    case 2:
        if (!is_helix(b)) return a == c ? Triple::HCH : Triple::RCL;      // H C H
        return (is_helix(a) ? a : c) == b ? Triple::CHH : Triple::CRL;    // C H H / H H C
    default:
        if (a == b && b == c) return Triple::HHH;
        return a == c ? Triple::RLR : Triple::RRL;
    }
}

// ---------------------------------------------------------------------------
struct BondParams {
    double bond_len = 1.0;   // equilibrium bond length
    double k_bond   = 100.0; // bond spring constant
};
struct BendParams {
    double kappa  = 1.0;     // bending stiffness
    double theta0 = 0.0;     // preferred bend angle (0 = straight)
};

// All run parameters. Filled by read_input(); main only threads this through.
// Input file format: one "key = value" per line, '#' starts a comment.
// Bond keys carry a pair suffix   (bond_len_CC, k_bond_HH, ...).
// Bend keys carry a triple suffix (kappa_CCC, theta0_HHH, ...).
struct Input {
    int N = 50;                       // number of residues
    std::string init = "rod";         // initial conformation: "rod" (straight) or "walk" (random walk)
    BondParams bond[N_PAIRS];         // indexed by Pair
    BendParams bend[N_TRIPLES];       // indexed by Triple

    // Zimm-Bragg: propagation for each helix sense, common nucleation penalty
    double s_R      = 1.0;
    double s_L      = 1.0;
    double sigma_zb = 1e-3;

    // thermodynamics / MC
    double kT          = 1.0;
    long   n_equil     = 1000;
    long   n_sweeps    = 10000;
    double max_disp    = 0.1;         // trial displacement amplitude
    unsigned long seed = 12345;

    // output
    long        log_every  = 100;
    long        dump_every = 0;       // 0 = never dump configurations
    std::string out_prefix = "run";

    const BondParams& p(Pair c)   const { return bond[static_cast<int>(c)]; }
    const BendParams& p(Triple t) const { return bend[static_cast<int>(t)]; }
};

bool read_input(const std::string& filename, Input& in);
void print_input(const Input& in);
