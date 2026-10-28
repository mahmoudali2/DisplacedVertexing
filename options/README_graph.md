# Graph-based displaced vertex finder — test steerings

Algorithm: J. Kriewald, [arXiv:2510.00856](https://arxiv.org/abs/2510.00856),
[github.com/jkriewald/delphes-LLP](https://github.com/jkriewald/delphes-LLP).
Key4hep port: Mahmoud Althakeel.

## Files

| File | What it does |
|---|---|
| `runGraphVertexFinder.py` | vertex finder **only**, on an existing `FittedTracks` collection |
| `runFullChainGraph.py` | full chain: digi → GGTF → GenFit (DAF) → graph vertexing |

## Build

The login environment on this account injects an LCG view and a local GenFit
that break the build silently. Use the helper at the repository root, which
sanitises both:

```bash
bash build_graph.sh Vertexing     # or: bash build_graph.sh all
```

## Run — vertexer only (fast iteration)

```bash
source <key4hep stack>
B=<repo>/build_nightly
export LD_LIBRARY_PATH=/afs/cern.ch/work/m/maali/public/GenFit/install/lib64:$B:$B/Vertexing:$LD_LIBRARY_PATH
export GAUDI_PLUGIN_PATH=$B/Vertexing:$B/Tracking:$B/.plugins:$GAUDI_PLUGIN_PATH
export PYTHONPATH=$B:$B/Vertexing/genConfDir:$PYTHONPATH

k4run runGraphVertexFinder.py \
    --inputFile  reco_with_fitted_tracks.root \
    --outputFile graph_vertices.root
```

### Diagnostics

```bash
# see the whole chi2 distribution, nothing rejected on chi2
k4run runGraphVertexFinder.py --openChi2 --inputFile in.root --outputFile out.root

# print the finder's internal decisions: prompt/displaced split, and for each
# track pair either ACCEPT (with chi2 and vertex) or the guard that rejected it
k4run runGraphVertexFinder.py --verboseKernel --inputFile in.root --outputFile out.root
```

The verbose output ends each event with a summary line:

```
[event] tracks=2  displaced=2  prompt=0  soft=0  | PV fitted
[pairs] tried=1  ok=0  rejected: guard=1 mahalanobis=0 fitfail=0
[event] clusters formed=0
```

which says directly why no vertex was formed.

## Run — full chain

```bash
k4run runFullChainGraph.py \
    --inputFile  sim.edm4hep.root \
    --outputFile reco.root \
    --modelPath  SimpleGatrIDEAv3o1.onnx \
    --fitter DAF
```

## Tracks beyond the solenoid

For decays outside the solenoid (r > ~2.3 m) there are no inner-tracker tracks;
the muons are only seen by the muon system. Point the finder at that collection
and read its `AtLastHit` state instead:

```bash
k4run runGraphVertexFinder.py \
    --inputTracks StandaloneMuonTracks \
    --trackState 3 \
    --inputFile muon_tracks.root --outputFile muon_vertices.root
```

## Muon-system arc vertexer (decays behind the solenoid)

`MuonSystemArcVertexFinder` vertexes LLP decays **outside the solenoid** with
muon-system tracks only (`StandaloneMuonTracks`). It does not use the tracks'
perigee parameters: each track is refitted locally from its muon-chamber hits
(circle in xy + line in z, generalised least squares with the correlated
multiple scattering of the yoke, 40 X0) and extrapolated back only to the decay,
adding the scattering in the dual-readout calorimeter (90 X0) to its error.
Track pairs are vertexed with a chi2 fit (ndf = 1 for two tracks).

    k4run runMuonSystemArcVertexFinder.py --inputFile tracks.root --outputFile vtx.root
    # --loose : no DCA/upstream/charge cut, all pairs kept (for cut scans)

Default cuts (HNL mN = 20 GeV working point, < 0.1% fake vertices per event on prompt decays):
`MaxPairDCA = 100 mm`, `RequireOppositeCharge = true`, `RequireBetweenIPAndHits = true`
(vertex closer to the origin than both first hits, same side), vertex upstream of both
first hits and outside the solenoid.

Denominator: truth r_xy > 2.3 m with >= 2 muon tracks (>= 3 hits).

| | efficiency | fakes / event | 68% half-width x / y / z |
|---|---|---|---|
| default (DCA 100, OS, IP-to-hits) | 67% | 0.09% | 24 / 22 / 15 mm |
| tighter (DCA 50, OS, IP-to-hits) | 56% | 0.06% | 20 / 18 / 12 mm |
| looser (DCA 300, OS, IP-to-hits) | 77% | 0.19% | 28 / 28 / 19 mm |

Vertex parameters: `[2, index1, index2, pairDCA, s1, s2, upstream, behindSolenoid, chargeProduct, betweenIPAndHits]`.

## Muon-system 3D helix vertexer (decays outside the solenoid, barrel + endcap)

`MuonSystemHelixVertexFinder` is the geometry-free companion of the arc vertexer
(same input, output collection `MuonSystemHelixVertices`, algorithmType 4). Each
track is a helix fitted to its own muon hits (circle in xy, z linear in the
turning angle; 3D line if the hits show no xy curvature) and is **stepped inward
in 3D** through a step-function field map (+2 T inside r < 2.3 m, |z| < 2.62 m,
-1.7 T elsewhere) while the energy loss of the calorimeter and iron is added back
(dE/dx 1.36 / 1.97 MeV/mm, calibrated on truth). Pairs are vertexed at the 3D point
of closest approach. Nothing assumes a barrel geometry, so tracks that cross the
endcap discs are vertexed like barrel tracks; the pair mass at the vertex is
unbiased (p_reco/p_true = 1.00).

    k4run runMuonSystemHelixVertexFinder.py --inputFile tracks.root --outputFile vtx.root
    # runs the arc vertexer in the same job (--noArc to skip it)
    # --region outside (default) | barrel | endcap | any
    # --loose  : no selection, all pairs kept (validation, cut scans)

Default cuts (main selection of the HNL stand-alone study, HNL/MS_full_outside):
`MaxPairDCA = 300 mm`, opposite rotation sense, vertex upstream of and between
the IP and the first hits, `Region = outside` (r_xy > 2300 mm or |z| > 2180 mm).
Mass and opening-angle cuts are hypothesis dependent and off by default
(`MinMass`, `MaxMass`, `MinCosOpening`).

Truth comparison (HNL, exactly two tracks; median vertex distance, fraction with DCA < 300 mm):

| true decay region | point | helix | arc |
|---|---|---|---|
| barrel r_xy > 2.3 m | 20 GeV, 1e-10 | 48 mm, 0.94 | 43 mm, 0.88 |
| forward r_xy < 2.3 m, \|z\| > 2.18 m | 20 GeV, 1e-10 | 72 mm, 0.93 | 76 mm, 0.81 |
| forward r_xy < 2.3 m, \|z\| > 2.18 m | 5 GeV, 10^-6.5 | 161 mm, 0.99 | 214 mm, 0.85 |

m_reco / m_true (16-84 %) at 20 GeV: 0.88-1.25 (helix) vs 0.66-1.41 (arc).
C++ = Python reference (`/eos/user/m/maali/HNL/MS_full_outside/code/helix3d.py`, `pair_vertex_stepped`)
on 2246 two-track events: vertex and DCA agree to < 1 um (float precision).

Vertex parameters: `[2, index1, index2, pairDCA, t1, t2, upstream, inRegion, chargeProduct, betweenIPAndHits, mass, cosOpening, p1, p2]`
(t = path length from the first hit, negative = upstream; chargeProduct of the rotation senses, 0 = straight track).
No covariance is computed: covMatrix = 0, chi2 = pairDCA [mm], ndf = 0.
