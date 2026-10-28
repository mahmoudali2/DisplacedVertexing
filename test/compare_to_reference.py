#!/usr/bin/env python3
"""Compare the vertex collections of two k4run outputs event by event.

usage: compare_to_reference.py <new output.root> <reference output.root>

For each vertex collection, count the vertices in both files and find the
largest difference in vertex position. Both numbers should be zero when the
standalone build reproduces the k4RecTracker-graph build.
"""
import sys
from podio.root_io import Reader

new_file = sys.argv[1]
ref_file = sys.argv[2]
collections = ["MuonSystemHelixVertices", "MuonSystemVertices"]

# one counter set per collection
n_new = {}             # vertices in the new file
n_ref = {}             # vertices in the reference file
n_count_mismatch = {}  # events where the number of vertices differs
max_shift = {}         # largest position difference [mm]
for name in collections:
    n_new[name] = 0
    n_ref[name] = 0
    n_count_mismatch[name] = 0
    max_shift[name] = 0.0

# the event lists can be read only once, so all collections are compared in one pass
new_events = Reader(new_file).get("events")
ref_events = Reader(ref_file).get("events")
n_events = 0
for ev_new, ev_ref in zip(new_events, ref_events):
    n_events += 1
    for name in collections:
        vtx_new = list(ev_new.get(name))
        vtx_ref = list(ev_ref.get(name))
        n_new[name] += len(vtx_new)
        n_ref[name] += len(vtx_ref)

        if len(vtx_new) != len(vtx_ref):
            n_count_mismatch[name] += 1
            continue

        # same number of vertices: compare them in order
        for a, b in zip(vtx_new, vtx_ref):
            pa = a.getPosition()
            pb = b.getPosition()
            shift = max(abs(pa.x - pb.x), abs(pa.y - pb.y), abs(pa.z - pb.z))
            if shift > max_shift[name]:
                max_shift[name] = shift

print(f"events compared: {n_events}")
for name in collections:
    print(f"{name}: new {n_new[name]}  ref {n_ref[name]}  events with different count {n_count_mismatch[name]}  "
          f"max |position shift| {max_shift[name]:.3e} mm")
