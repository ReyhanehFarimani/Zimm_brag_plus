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
    std::string state_pattern;        // explicit initial states, one character per residue (C/R/L), e.g. LLLLCCCC;
                                      // overrides init_state when given (must have exactly N characters)
    int freeze_states = 0;            // 1: no state (or hinge) moves, the states stay as initialised
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
    // on-site ("ground") energy of a helical residue, the same for R and L; coil = 0
    double E_helix = 0.0;

    // thermodynamics / MC
    double kT          = 1.0;
    long   n_equil     = 1000;
    long   n_sweeps    = 10000;
    double max_disp    = 0.1;         // trial displacement amplitude
    int    n_flip      = 0;           // domain sense-flip (R <-> L cluster) moves per sweep (0 = off)
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
    //   mu, nu, kappa' (gb_aniso_eps = 1).  Helix-coil pairs use the same isotropic potential as coil-coil.
    std::string nb_hh        = "same";
    double      gb_eps0      = 1.0;
    int         gb_aniso_eps = 1;
    double      gb_mu        = 2.0;
    double      gb_nu        = 1.0;
    double      gb_kappa_p   = 5.0;
    // nb_hh = fit: the measured helix-helix radial pair potential of the
    // umbrella campaign (see pair_potential.cpp): erfc radial sectors x
    // quadratic angular maps, handedness-resolved (R.R / L.L / R.L).
    //   hf_theta0 : which campaign, 100 or 45 (deg)
    //   hf_scale  : overall amplitude factor
    //   hf_eps_s  : chiral amplification (multiplies the sin sector only;
    //               0 = racemic reference, 1 = as measured)
    //   hf_len    : how many code length units make one campaign length unit a. The parameter files
    //               built from the measured fits (bonds, coil core) are ALREADY in a, so hf_len = 1.
    //               [Until 2026-09-17 the t100 scans used 1.2311 = rod_L / L_block, which evaluated the
    //               helix-helix potential on a scale stretched by 23 % relative to the chain.]
    //   hf_clamp  : 1 (default) = the CHIRAL sector is evaluated only inside the fitted domain: its
    //               angular maps at e1, e2 clamped to [-0.5, 0.5], its radial factor at r clamped to
    //               [1.5, 2.5] a and shifted to vanish at 2.5 a (continuous: constant inside 1.5 a, zero
    //               beyond 2.5 a).  The narrow chiral radial factor otherwise outgrows the achiral one
    //               when extrapolated inward.  0 = extrapolate the fitted form everywhere (legacy).
    //   hf_sig_c  : radial width (RANGE) of the chiral sector in units of a; 0 (default) = the campaign
    //               family median (0.70 a for theta0 = 100, 0.51 a for 45). The chiral force is normalised
    //               to eps_s C at r0 = 1.7 a for every width; a larger width lets it decay more slowly
    //               outward. Weakly constrained by the data (t100: 5 families, median 0.70, weighted
    //               mean 0.46, spread 0.26 a). The eps_s at which attraction first appears depends on it.
    int    hf_theta0 = 100;
    double hf_scale  = 1.0;
    double hf_eps_s  = 1.0;
    double hf_len    = 1.0;
    int    hf_clamp  = 1;
    double hf_sig_c  = 0.0;
    //   hf_eps_rl : chiral amplitude of OPPOSITE-handed (R.L) pairs; < 0 (default) = the same as hf_eps_s.
    //               hf_eps_s then acts on same-handed pairs only: e.g. hf_eps_s = 12.2, hf_eps_rl = 1 gives a
    //               same-handed well of -3 kT while R.L pairs keep the measured, purely repulsive potential.
    double hf_eps_rl = -1.0;
    // nb_hh = db: tabulated helix-helix potential (prelim/cg/helix_pair_db_full.py, HELIX_PAIR_DB v2)
    // read from db_file, nearest-node lookup in (r, betaA, betaB, psi, alphaA, alphaB); the registries
    // alpha come from the per-residue registry directions (registry.h) and the twist move
    // (n_twist attempts per sweep, step +- twist_step degrees).  hf_theta0 / hf_len select the
    // table's theta0 entry and the length unit as for nb_hh = fit.
    std::string db_file  = "";
    int    n_twist    = 0;
    double twist_step = 30.0;          // degrees
    // harmonic registry-twist term between sequential same-handed residues (registry.h):
    //   E = 1/2 twist_kappa (dalpha - s twist_alpha0)^2, s = +1 R.R / -1 L.L, dalpha in (-pi, pi].
    // twist_alpha0 in degrees (R.R sign), twist_kappa in kT/rad^2; a negative value means "the
    // measured default of hf_theta0": 100 -> -135.3 deg, kappa = 1/sigma^2 with sigma = 42.5 deg
    // (1.82 kT/rad^2); 45 -> +25.3 deg, sigma = 22.3 deg (6.60 kT/rad^2).  twist_kappa = 0 switches
    // the term off (the registries are then coupled by the pair table only).
    double twist_alpha0 = -999.0;
    double twist_kappa  = -1.0;
    // the same term for OPPOSITE-handed neighbours (R.L), measured 2026-09-23 on the graded runs: the
    // registries across an R|L junction are anti-aligned, dalpha_0 = 180 deg; harmonic fits (grant repo,
    // cg_bonded_grd_folded_fits.csv) K = 2.47 kT/rad^2 at theta0 = 45 (sigma 26 deg, as sharp as R.R) and
    // 0.11 at theta0 = 100 (broad, ~1 kT deep), kappa = 2K; leaving them free (kappa = 0) gives every
    // R|L junction ~1.7 kT of spurious registry entropy against R.R.  Coil-helix junctions ARE free
    // (measured uniform, R = 0.06) and carry no term.  Negative = the measured default of hf_theta0.
    double twist_alpha0_rl = -999.0;
    double twist_kappa_rl  = -1.0;
    double twist_alpha0_rad = 0.0;     // derived
    double twist_alpha0_rl_rad = 0.0;  // derived
    bool   twist_on = false;           // derived: nb_hh = db and twist_kappa > 0
    int    debug_registry = 0;         // 1: check the registry invariants every sweep (slow)

    // Verlet neighbour list for the non-bonded energy. nl_skin = total skin width (code units):
    //   < 0 = automatic (4 sqrt(3) max_disp: a single displacement never exceeds half of it), 0 = off
    //   (all-pairs loops). The list changes the cost only, the energies are identical.
    double nl_skin   = -1.0;
    // Minimum residue separation |i - j| of a NON-BONDED pair (2026-09-23): the bonded bend potentials
    // are Boltzmann inversions of the monomer-resolved chains and already contain every 1-3
    // interaction, so the coil core / helix pair potential must start at 1-4 (nb_min_sep = 3).
    // With 2 (the old behaviour) the 1-3 core is counted twice and the chain loses its large-bend
    // tail (P(theta > 90 deg) 0.13 vs 0.17 measured; Rg2 +10 % for the all-coil 15-block chain).
    int    nb_min_sep = 3;

    // STAR architecture (user 2026-09-23, "star with 20 arms"): n_arms > 0 makes the chain n_arms arms of
    // N residues each, grafted to a fixed spherical core at the origin (radius core_radius, bead units a).
    // Bonds, bends, state couplings and registry twists exist only within an arm; non-bonded pairs act
    // between all residues (within an arm from separation nb_min_sep, across arms always).  Every residue
    // feels a harmonic wall at core_radius + 1/2 (core_k), and the first residue of each arm is tethered
    // (graft_k) to its graft site, one of n_arms Fibonacci points on the sphere of radius core_radius + 1/2.
    int    n_arms      = 0;
    // PERIODIC BOX (user 2026-09-24, "a box of density 0.7 of 500 50-segment chains"): box > 0 = cubic periodic
    // box of side box (bead units a).  With n_arms > 0 the arms are n_arms FREE chains of N residues (no core,
    // no graft), placed as randomly oriented rods at random positions; every non-bonded distance (pair energies,
    // neighbour list) uses the minimum image; positions stay unwrapped (bonded terms use direct differences).
    // Ree and Rg2 in the obs file are then per-chain means.  Needs box > 2 (cutoff + nl_skin).
    // Segment moves for dense systems (2026-09-24): n_seg crankshaft rotations per sweep of 1..seg_len_max beads
    // about the axis through their anchors by up to seg_rot rad (update_segment.cpp); n_segflip handedness flips
    // of a random helical segment of 1..segflip_len_max residues (update_flipdomain.cpp); n_torsion tail
    // rotations about the bond (i, i+1) by up to torsion_rot rad (a pivot with the bond as axis).
    int    n_seg         = 0;
    int    seg_len_max   = 3;
    double seg_rot       = 0.5;
    int    n_segflip     = 0;
    int    segflip_len_max = 8;
    int    n_torsion     = 0;
    double torsion_rot   = 0.5;
    double box         = 0.0;
    // RESTART (2026-09-24, user: "restart using the last point"): read states, positions and registries of one
    // frame of a trajectory written by this code (<out_prefix>_conf.xyz, extended xyz) instead of the fresh
    // initial condition; restart_frame = -1 is the last frame. Everything else (J, moves, box, ...) comes from
    // this input as usual; a new seed continues the Markov chain from that configuration.
    std::string restart_file;
    int    restart_frame = -1;
    // pivot energies through a cell grid (O(tail) instead of O(tail x N) all-pairs): -1 = automatic (systems of
    // 5000+ residues), 0 = never, 1 = always (used by the harness); the energies agree to round-off
    int    pivot_grid  = -1;
    double core_radius = 4.0;
    double core_k      = 50.0;
    double graft_k     = 10.0;

    // rod geometry written to the trajectory for helical residues (coil residues are spheres of radius rod_r)
    double rod_L = 2.5;               // rod length
    double rod_r = 1.0;               // rod radius (column rod_D = 2 rod_r)

    // output
    long        log_every  = 100;
    long        dump_every = 0;       // 0 = never dump configurations
    std::string out_prefix = "run";

    const BondParams& p(Pair c)   const { return bond[static_cast<int>(c)]; }
    double site_energy(State s) const { return is_helix(s) ? E_helix : 0.0; }
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
