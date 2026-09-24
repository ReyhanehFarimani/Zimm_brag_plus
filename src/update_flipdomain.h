#pragma once
#include <random>
#include "chain.h"

// Sense-flip moves (geometry untouched; only local state/bond/bend energies change). Two symmetric proposals,
// chosen with probability 1/2 each:
//   (a) single-site: flip residue i from R to L or vice versa (moves, creates and annihilates R|L walls);
//   (b) whole-domain: flip the maximal same-sense domain containing i, REJECTED if a neighbouring residue is
//       helical -- flipping would merge domains irreversibly and break detailed balance, so only domains
//       flanked by coil or the chain ends may flip (global sense decorrelation).
// NOTE: assumes the non-bonded interaction is invariant under R <-> L (it depends only on is_helix), which
// holds for the isotropic cores and the Gay-Berne rods. Coil residue: no move (counted as rejected).
bool try_domain_flip(Chain& chain, int i, const Input& in, std::mt19937_64& rng);
// segment flip: residues i .. i + len - 1 (len uniform in 1..segflip_len_max), all helical, flip their handedness
bool try_segment_flip(Chain& chain, int i, const Input& in, std::mt19937_64& rng);
