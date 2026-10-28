#!/bin/bash
# Smoke test of HelixPairVertexFinder on inner-tracker tracks (CLD tracker in
# IDEA_w_CLD_IT, HNL mN = 20 GeV, |V|^2 = 1e-9, 500 events): runs the
# inner-tracker steering and prints the distance of the vertex to the true decay.
# Run on fcc-ironic. Result 2026-10-27: 435/500 events with a vertex, median 0.15 mm.
REPO=/afs/cern.ch/work/m/maali/public/DisplacedVertexing
IN=/eos/user/m/maali/HNL/IDEA_w_CLD_IT/tracks_it/20_1e-9_c0.root

# private podio scratch area (shared /tmp/podio on ironic belongs to another user)
export TMPDIR=/tmp/dv_$(id -u)
mkdir -p "$TMPDIR"
OUT=$TMPDIR/it_smoke.root

source /cvmfs/sw.hsf.org/key4hep/releases/2026-04-08/x86_64-almalinux9-gcc14.2.0-opt/key4hep-stack/2026-04-08-i6h4f2/setup.sh >/dev/null 2>&1 < /dev/null
source $REPO/setup.sh

k4run $REPO/options/runInnerTrackerVertexing.py -n 500 --inputFile $IN --outputFile $OUT || exit 1
python3 $REPO/test/check_truth_distance.py $OUT InnerTrackerHelixVertices
