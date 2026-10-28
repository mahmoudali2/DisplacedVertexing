# Graph-based displaced vertex finder

Algorithm: J. Kriewald, [arXiv:2510.00856](https://arxiv.org/abs/2510.00856),
[github.com/jkriewald/delphes-LLP](https://github.com/jkriewald/delphes-LLP).
Key4hep port: Mahmoud Althakeel.

| File | What it does |
|---|---|
| `runGraphVertexFinder.py` | vertex finder only, on an existing fitted track collection |
| `runFullChainGraph.py` | digitisation, GGTF, GenFit and graph vertexing (needs k4RecTracker) |

Build and set up the package as in the main `README.md`, then:

```bash
k4run runGraphVertexFinder.py --inputFile tracks.root --outputFile vertices.root
```

Options:

- `--inputTracks <name>`: track collection (default `FittedTracks`)
- `--trackState <n>`: state location used (1 = AtIP, 3 = AtLastHit)
- `--openChi2`: keep every chi2, to see the full distribution
- `--verboseKernel`: print the finder's decisions per event (prompt/displaced split,
  accepted or rejected pairs and why)
