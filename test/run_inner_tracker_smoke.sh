#!/bin/bash
# Smoke test of HelixPairVertexFinder on silicon tracker tracks (CLD, X -> mu mu decays up to
# 2.3 m, tracks from CombinatorialTrackFinder + GenFit): prints the distance of the vertex to
# the true decay.
# usage: run_inner_tracker_smoke.sh <setup script that makes DisplacedVertexing available>
SETUP=$1
DIR=$(cd "$(dirname "$0")/.." && pwd)
IN=/eos/user/m/maali/Displaced_performance/data/ctf/CLD/llp_c0.root
source $SETUP
OUT=$TMPDIR/it_smoke.root
k4run $DIR/options/runHelixPairVertexing.py --inputFile $IN --outputFile $OUT \
      --inputTracks CTFGenfitTracks --stateLocation 2 > $TMPDIR/it_smoke.log 2>&1 || { tail -5 $TMPDIR/it_smoke.log; exit 1; }
python3 $DIR/test/check_truth_distance.py $OUT HelixPairVertices
