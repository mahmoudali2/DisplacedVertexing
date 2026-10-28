#!/bin/bash
# =============================================================================
# Build and install DisplacedVertexing (standalone, outside k4RecTracker).
#
# Usage:  bash build.sh                       (stable 2026-04-08 release)
#         STACK=<other setup.sh> bash build.sh
#
# Run on fcc-ironic. podio writes into TMPDIR, and the shared /tmp/podio there
# belongs to another user, so a private TMPDIR is set below.
# =============================================================================
# (no set -u: the key4hep setup script references unset variables)
REPO=/afs/cern.ch/work/m/maali/public/DisplacedVertexing
STACK=${STACK:-/cvmfs/sw.hsf.org/key4hep/releases/2026-04-08/x86_64-almalinux9-gcc14.2.0-opt/key4hep-stack/2026-04-08-i6h4f2/setup.sh}
BUILD=${BUILD:-$REPO/build}
INSTALL=${INSTALL:-$REPO/install}

# private podio scratch area
export TMPDIR=/tmp/dv_$(id -u)
mkdir -p "$TMPDIR"

cd "$REPO" || exit 1
source "$STACK" >/dev/null 2>&1 < /dev/null

# --- strip the LCG view that this account's login environment injects ---------
# (its Python 3.13.11 clashes with the stack's 3.13.8 inside CMake's FindPython)
strip_lcg () { printenv "$1" 2>/dev/null | tr ':' '\n' | grep -v "sft.cern.ch/lcg" | paste -sd: ; }
export CMAKE_PREFIX_PATH=$(strip_lcg CMAKE_PREFIX_PATH)
export C_INCLUDE_PATH=$(strip_lcg C_INCLUDE_PATH)
export CPLUS_INCLUDE_PATH=$(strip_lcg CPLUS_INCLUDE_PATH)
export PATH=$(strip_lcg PATH)
unset PYTHONHOME

echo ">>> configuring"
cmake -B "$BUILD" -S . -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$INSTALL" || exit 1

echo ">>> building"
cmake --build "$BUILD" -j8 || exit 1

echo ">>> installing into $INSTALL"
cmake --install "$BUILD" || exit 1

echo ">>> OK"
