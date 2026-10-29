#!/usr/bin/env python3
"""
Compare the C++ HelixPairVertexFinder (muon-system steering with --loose) with the Python
reference (columns h_* of MS_full_outside/events/<id>.npz) on events with exactly two tracks.

usage: compare_helix_cpp_with_python.py <Gaudi output.root> <registry id>
"""
import sys

import numpy as np
from podio.root_io import Reader

gaudi_output_file = sys.argv[1]
registry_id = sys.argv[2]
python_reference = np.load(f"/eos/user/m/maali/HNL/MS_full_outside/events/{registry_id}.npz")

# one row per two-track event: event number, differences C++ - Python, number of C++ vertices
rows = []
for event_number, event in enumerate(Reader(gaudi_output_file).get("events")):
    has_two_tracks = python_reference["ntrk"][event_number] == 2
    has_python_vertex = np.isfinite(python_reference["h_vr"][event_number])
    if not has_two_tracks or not has_python_vertex:
        continue
    cpp_vertices = list(event.get("MuonSystemHelixVertices"))
    if len(cpp_vertices) != 1:
        rows.append((event_number, np.nan, np.nan, np.nan, np.nan, len(cpp_vertices)))
        continue
    vertex = cpp_vertices[0]
    position = vertex.getPosition()
    parameters = list(vertex.getParameters())
    radius_difference = np.hypot(position.x, position.y) - python_reference["h_vr"][event_number]
    z_difference = position.z - python_reference["h_vz"][event_number]
    dca_difference = parameters[3] - python_reference["h_dca"][event_number]
    if np.isfinite(python_reference["h_mass"][event_number]):
        mass_difference = parameters[10] - python_reference["h_mass"][event_number]
    else:
        mass_difference = 0.0
    rows.append((event_number, radius_difference, z_difference, dca_difference, mass_difference, 1))

table = np.array(rows)
with_one_vertex = table[:, 5] == 1
print(f"{registry_id}: two-track events {len(table)}, with exactly one C++ vertex {with_one_vertex.sum()}")
for column, name in ((1, "d r_xy [mm]"), (2, "d z [mm]"), (3, "d DCA [mm]"), (4, "d mass [GeV]")):
    absolute_differences = np.abs(table[with_one_vertex, column])
    print(f"  |{name}|: median {np.median(absolute_differences):.2e}  "
          f"99% below {np.percentile(absolute_differences, 99):.2e}  max {absolute_differences.max():.2e}  "
          f"number above 1e-2: {(absolute_differences > 1e-2).sum()}")
