# DisplacedVertexing (standalone)

Gaudi/Key4hep algorithms that find displaced (long-lived-particle) decay vertices from
EDM4hep tracks. They were written for HNL → μμ at FCC-ee with the IDEA muon system, where
the decay can happen metres away from the IP.

This is the standalone build, outside k4RecTracker. The same sources live in
`k4RecTracker-graph/Vertexing/` (untracked package there). Building only needs Gaudi,
k4FWCore, EDM4hep and ROOT. DD4hep and GenFit are not needed.

## Algorithms

| Gaudi algorithm | Kernel | What it does |
|---|---|---|
| `ArcPairVertexFinder` | `ArcPairVertexing` | Two-track vertexing with **local** track fits. Each track is refitted from its own hits (circle in xy, line in z), with multiple scattering in the material in its error. It is extrapolated back only as far as the decay, and pairs are vertexed with a χ² fit (gives a covariance). |
| `HelixPairVertexFinder` | `HelixPairVertexing` | 3D version with no barrel/endcap assumption. A helix is fitted per track in the field at its innermost hit and stepped inward through a two-zone field map, adding back the energy lost in a material map. Pairs are vertexed at the 3D point of closest approach (also gives the pair mass). |
| `GraphDisplacedVertexFinder` | `GraphVertexKernel` | Graph-based multi-track vertex finding on fitted inner-tracker tracks (uses the track states, not the hits). Key4hep port of J. Kriewald, arXiv:2510.00856. |

The two pair finders were first written for the IDEA muon system (formerly
`MuonSystemArcVertexFinder` / `MuonSystemHelixVertexFinder`). They contain no detector
numbers now: the field, material, region and resolutions are all properties, and the
defaults describe no detector (`InnerField` / `BFieldForMomentum` must be set). The
detector is described in the steering:

| Steering | Detector |
|---|---|
| `options/runMuonSystemVertexing.py` | IDEA muon system: +2 T / −1.7 T field map, calorimeter + yoke material, region outside the solenoid. Both finders. |
| `options/runInnerTrackerVertexing.py` | silicon inner tracker (CLD in IDEA_w_CLD_IT): uniform 2 T, no material, 8 innermost hits, DCA < 2 mm, r_xy > 10 mm. Helix finder. |
| `options/runGraphVertexFinder.py`, `runFullChainGraph.py` | graph finder (see `options/README_graph.md`) |

`Region` (both pair finders) is `outside | barrel | endcap | inside | any` relative to the
cylinder `RegionR`, `RegionHalfZ` (`include/VertexRegion.h`).

## Layout

```
components/   Gaudi algorithms (EDM4hep in/out, properties)
include/, src/  framework-free kernels
options/      steering examples (k4run)
python/       Python reference implementations the C++ was validated against
  arc_vertexer.py          reference for ArcPairVertexing (matches C++ exactly)
  helix3d.py               reference for HelixPairVertexing (C++ == Python to < 1 um)
  compare_helix_cpp_py.py  the C++ vs Python comparison
test/         run_regression.sh: IDEA muon-system steering vs the k4RecTracker-graph build
              (bit-identical, both finders); run_inner_tracker_smoke.sh: CLD tracker check
```

## Build

On fcc-ironic:

```bash
bash build.sh          # configure + build + install into ./install (stable 2026-04-08)
```

`build.sh` removes the LCG paths from the login environment and sets a private `TMPDIR`
for podio.

## Run

```bash
source /cvmfs/sw.hsf.org/key4hep/releases/2026-04-08/x86_64-almalinux9-gcc14.2.0-opt/key4hep-stack/2026-04-08-i6h4f2/setup.sh
source setup.sh
k4run options/runMuonSystemVertexing.py --inputFile tracks.root --outputFile vertices.root
```

`options/runFullChainGraph.py` also runs digitisation, GGTF and GenFit. Those come from
k4RecTracker, so that steering needs k4RecTracker on the path too.
