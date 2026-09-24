#pragma once
#include <random>
#include "chain.h"
// Segment (crankshaft) move (2026-09-24, user: "segment moves" for the dense box where whole-tail pivots are
// accepted at 2 %): the k = 1 .. seg_len_max beads between two anchors i and j = i + k + 1 of one arm rotate
// rigidly about the axis through the anchors by a random angle in [-seg_rot, seg_rot]. Bond lengths and every
// bend inside the segment are preserved (rigid rotation about an axis through both anchors); only the bends at
// i and j, the twist terms of the junctions i-1..j (the anchors' tangents turn, the segment's registries rotate
// with it), the non-bonded pairs of the beads i..j and the core terms of the segment change. The axis is fixed
// by the configuration and the angle is symmetric, so the proposal is symmetric. With the neighbour list on,
// the moving beads are re-referenced first so that the list covers old and new positions whenever the
// displacement stays below skin/2 (seg_len_max 3, seg_rot 0.5, skin 4 in the box generator); otherwise the
// all-pairs path is taken. Returns true if accepted.
bool try_segment_move(Chain& chain, int i, const Input& in, std::mt19937_64& rng);
