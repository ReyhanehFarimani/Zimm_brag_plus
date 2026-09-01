#pragma once
#include <string>

// Residue state: coil, right-handed helix, left-handed helix.
enum class State : int { Coil = 0, R = 1, L = 2 };
inline const char* state_name(State s) { return s == State::Coil ? "C" : (s == State::R ? "R" : "L"); }
inline bool is_helix(State s) { return s != State::Coil; }
// Ising-like spin: L = -1, C = 0, R = +1
inline int   spin(State s) { return s == State::R ? 1 : (s == State::L ? -1 : 0); }
inline State from_spin(int q) { return q > 0 ? State::R : (q < 0 ? State::L : State::Coil); }

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
    std::string init       = "rod";   // initial conformation: "rod" (straight) or "walk" (random walk)
    std::string init_state = "coil";  // initial states: "coil" (all C) or "random" (R or L, equal probability)
    int n_states = 3;                 // 3: C/R/L with +-1 spin steps (R<->L only via C)
                                      // 2: R/L only, the state move flips R<->L directly (requires init_state = random)
    // bend_key: which states select the bend parameters at residue j
    //   "triple": triple class of (j-1, j, j+1)                           [default]
    //   "pair"  : pair class of (j-1, j), i.e. of the bond entering j; the pair classes are mapped onto
    //             the bend table as CC->CCC, CH->CHC, HH->HHH, RL->RLR (freely-jointed-with-hinges models)
    std::string bend_key   = "triple";
    BondParams bond[N_PAIRS];         // indexed by Pair
    BendParams bend[N_TRIPLES];       // indexed by Triple

    // nearest-neighbour state (Ising-like) energies, keyed by pair class:
    //   E_HH = -J0   (helix propagation)      E_CH = +J1   (helix/coil boundary)
    //   E_CC =  0                              E_RL = +J2   (R/L junction)
    double J0 = 1.0;
    double J1 = 2.0;
    double J2 = 2.0;

    // thermodynamics / MC
    double kT          = 1.0;
    long   n_equil     = 1000;
    long   n_sweeps    = 10000;
    double max_disp    = 0.1;         // trial displacement amplitude
    int    n_pivot     = 0;           // pivot moves per sweep (0 = off): rotate the chain beyond a random bead by a
                                      // random angle in [-max_rot, max_rot] about a random axis (full energy change)
    double max_rot     = 3.14159265358979;
    int    n_hinge     = 0;           // hinge moves per sweep (0 = off). Flips the spin of a residue and redraws the
                                      // directions of the affected bonds from their new hinge distributions, rotating
                                      // the downstream chain rigidly. Exact for chains WITHOUT non-bonded interactions.
    unsigned long seed = 12345;

    // non-bonded (steric) interaction between all bead pairs except bonded neighbours |i-j| = 1
    //   nb_type = none | gauss | wca
    //   gauss : U(r) = nb_A * exp(-(r/nb_sigma)^2)                       (soft isotropic core)
    //   wca   : U(r) = 4 nb_A [(nb_sigma/r)^12 - (nb_sigma/r)^6] + nb_A   for r < 2^(1/6) nb_sigma
    //   nb_rcut: cutoff (0 = automatic: 3.5 sigma for gauss, 2^(1/6) sigma for wca)
    // (isotropic and state-independent for now; the pair energy is routed through the two states so a
    //  state-dependent table can be added without touching the moves)
    std::string nb_type  = "none";
    double      nb_A     = 9.01;
    double      nb_sigma = 1.5;
    double      nb_rcut  = 0.0;
    // helix-helix pairs: nb_hh = same (the isotropic potential above) | gb (purely repulsive, WCA-truncated
    // Gay-Berne between rods of length rod_L and diameter D = 2 rod_r, axes = local chain tangents):
    //   U = 4 eps [(s0/rho)^12 - (s0/rho)^6] + eps  for rho < 2^(1/6) s0,  rho = r - sigma(u1,u2,r^) + s0,  s0 = D,
    //   sigma = s0 [1 - chi/2 ( (a+b)^2/(1+chi c) + (a-b)^2/(1-chi c) )]^(-1/2),  a = r^.u1, b = r^.u2, c = u1.u2,
    //   chi = (k^2-1)/(k^2+1), k = L/D;  eps = eps0 (gb_aniso_eps = 0) or the standard anisotropic form with
    //   mu, nu, kappa' (gb_aniso_eps = 1).  Helix-coil pairs use the isotropic potential for now.
    std::string nb_hh        = "same";
    double      gb_eps0      = 1.0;
    int         gb_aniso_eps = 1;
    double      gb_mu        = 2.0;
    double      gb_nu        = 1.0;
    double      gb_kappa_p   = 5.0;

    // rod geometry written to the trajectory for helical residues (coil residues are spheres of radius rod_r)
    double rod_L = 2.5;               // rod length
    double rod_r = 1.0;               // rod radius (column rod_D = 2 rod_r)

    // output
    long        log_every  = 100;
    long        dump_every = 0;       // 0 = never dump configurations
    std::string out_prefix = "run";

    const BondParams& p(Pair c)   const { return bond[static_cast<int>(c)]; }
    double pair_energy(Pair c) const {
        switch (c) {
        case Pair::HH: return -J0;
        case Pair::CH: return  J1;
        case Pair::RL: return  J2;
        default:       return 0.0;
        }
    }
    const BendParams& p(Triple t) const { return bend[static_cast<int>(t)]; }
};

bool read_input(const std::string& filename, Input& in);
void print_input(const Input& in);
