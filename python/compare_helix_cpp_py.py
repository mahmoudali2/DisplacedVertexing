#!/usr/bin/env python3
"""C++ MuonSystemHelixVertexFinder (--loose) vs the Python reference (MS_full_outside/events/<id>.npz, h_* columns) on 2-track events.
usage: compare_helix_cpp_py.py <gaudi output.root> <registry id>"""
import sys
import numpy as np
from podio.root_io import Reader
out, pid = sys.argv[1], sys.argv[2]
ref = np.load(f"/eos/user/m/maali/HNL/MS_full_outside/events/{pid}.npz")
rows = []
for iev, ev in enumerate(Reader(out).get("events")):
    if ref["ntrk"][iev] != 2 or not np.isfinite(ref["h_vr"][iev]):
        continue
    vs = list(ev.get("MuonSystemHelixVertices"))
    if len(vs) != 1:
        rows.append((iev, np.nan, np.nan, np.nan, np.nan, len(vs))); continue
    v = vs[0]; p = v.getPosition(); par = list(v.getParameters())
    rows.append((iev, np.hypot(p.x, p.y) - ref["h_vr"][iev], p.z - ref["h_vz"][iev], par[3] - ref["h_dca"][iev],
                 (par[10] - ref["h_mass"][iev]) if np.isfinite(ref["h_mass"][iev]) else 0.0, 1))
a = np.array(rows)
ok = a[:, 5] == 1
print(f"{pid}: 2-track events {len(a)}, with exactly one C++ vertex {ok.sum()}")
for k, name in ((1, "d r_xy [mm]"), (2, "d z [mm]"), (3, "d DCA [mm]"), (4, "d mass [GeV]")):
    d = np.abs(a[ok, k])
    print(f"  |{name}|: median {np.median(d):.2e}  p99 {np.percentile(d, 99):.2e}  max {d.max():.2e}  n(>1e-2) {(d > 1e-2).sum()}")
