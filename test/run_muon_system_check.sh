#!/bin/bash
# Muon-system check: tracks one HNL sample (mN = 15 GeV) with CombinatorialTrackFinder, runs both
# pair finders with the muon-system steering and prints the distance of the vertex to the truth.
# usage: run_muon_system_check.sh <setup script that makes CombinatorialTrackFinder and DisplacedVertexing available>
SETUP=$1
DIR=$(cd "$(dirname "$0")/.." && pwd)
source $SETUP
OUT=$TMPDIR/ms_check; mkdir -p $OUT
export MS_IN=/eos/user/m/maali/HNL/lowcoupling_full/digi/lc_mN15_um9p00.root
export MS_OUT=$OUT/tracks.root
cd $OUT
k4run $DIR/test/run_muon_system_ctf.py > $OUT/ctf.log 2>&1 || { echo "tracking failed"; grep -i "error\|fatal" $OUT/ctf.log | head; exit 1; }
k4run $DIR/options/runMuonSystemVertexing.py -n -1 --inputFile $MS_OUT --outputFile $OUT/vertices.root \
      > $OUT/vertexing.log 2>&1 || { echo "vertexing failed"; tail -5 $OUT/vertexing.log; exit 1; }
for COLLECTION in MuonSystemHelixVertices MuonSystemVertices; do
  echo "== $COLLECTION"
  python3 $DIR/test/check_truth_distance.py $OUT/vertices.root $COLLECTION
done
