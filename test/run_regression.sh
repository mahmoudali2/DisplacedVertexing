#!/bin/bash
# Regression test of the standalone build against the k4RecTracker-graph build.
# Runs HelixPairVertexFinder + ArcPairVertexFinder with the IDEA muon-system steering (--loose) on one HNL
# file and compares every vertex with the reference output made on 2026-10-25 by
# test_vertexing/run_helix_validation.sh. Run on fcc-ironic.
REPO=/afs/cern.ch/work/m/maali/public/DisplacedVertexing
IN=/eos/user/m/maali/HNL/lowcoupling_full/tracks/lc_mN15_um9p00.root
REF=/afs/cern.ch/work/m/maali/public/test_vertexing/helix_gaudi_val_lc_mN15_um9p00.root
OUT=$REPO/test/regression_lc_mN15_um9p00.root

# private podio scratch area (shared /tmp/podio on ironic belongs to another user)
export TMPDIR=/tmp/dv_$(id -u)
mkdir -p "$TMPDIR"

source /cvmfs/sw.hsf.org/key4hep/releases/2026-04-08/x86_64-almalinux9-gcc14.2.0-opt/key4hep-stack/2026-04-08-i6h4f2/setup.sh >/dev/null 2>&1 < /dev/null
source $REPO/setup.sh

k4run $REPO/options/runMuonSystemVertexing.py -n -1 --loose \
      --inputFile $IN --outputFile $OUT || exit 1

python3 $REPO/test/compare_to_reference.py $OUT $REF
