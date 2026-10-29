# DisplacedVertexing

Gaudi algorithms for Key4hep that find displaced decay vertices (for example of long-lived
particles) from EDM4hep tracks. Nothing is specific to a detector: the field, material and
cuts are set in the steering file.

## Algorithms

| Algorithm | Input | What it does |
|---|---|---|
| `HelixPairVertexFinder` | track state | Describes each track as a helix from its fitted track state, starting at its innermost hit. Steps it inward in 3D through a field map, adding back the energy lost in the material, and takes each pair's point of closest approach as the vertex. Works in the barrel and the endcaps. Also gives the pair mass. |
| `ArcPairVertexFinder` | track hits | Refits each track locally from its own hits (circle in xy, line in z), with multiple scattering included in the errors. Extrapolates back only to the decay and fits each pair with a χ² vertex fit, which also gives a covariance. |
| `GraphDisplacedVertexFinder` | track states | Graph-based multi-track vertex finding (J. Kriewald, arXiv:2510.00856). A state given at a hit (e.g. `AtFirstHit`) is first moved to the origin. See `options/README_graph.md`. |

## Build

```bash
source /cvmfs/sw.hsf.org/key4hep/setup.sh
git clone <this repository> DisplacedVertexing
cd DisplacedVertexing
cmake -B build -S . -DCMAKE_INSTALL_PREFIX=$PWD/install
cmake --build build -j8 --target install
```

## Set up (in every new shell)

```bash
source /cvmfs/sw.hsf.org/key4hep/setup.sh
export LD_LIBRARY_PATH=<path to>/DisplacedVertexing/install/lib64:$LD_LIBRARY_PATH
export PYTHONPATH=<path to>/DisplacedVertexing/install/python:$PYTHONPATH
```

## Run

```bash
k4run options/runHelixPairVertexing.py --inputFile tracks.root --outputFile vertices.root \
      --inputTracks <track collection> --field 2.0 --stateLocation 2
```

| Steering | Finder |
|---|---|
| `runHelixPairVertexing.py` | helix finder |
| `runArcPairVertexing.py` | arc finder |
| `runGraphVertexing.py` | graph finder |
| `runMuonSystemVertexing.py` | helix and arc finders with the IDEA muon-system field map and material (decays outside the solenoid) |

Each steering lists every property of its finder with the default value, grouped by use:
track input, detector (field, resolutions, material), physics, algorithm, vertex selection
and output.
