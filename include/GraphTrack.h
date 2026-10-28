/*
 * Part of the Key4hep port of the graph-based displaced vertex finder.
 * Algorithm: J. Kriewald, arXiv:2510.00856 (github.com/jkriewald/delphes-LLP).
 * Port author: Mahmoud Althakeel
 */
/*
 * Key4hep-native stand-in for the Delphes `Candidate` used by the graph-based
 * displaced-vertex finder.
 *
 * The ported algorithm (GraphVertexKernel) reads only a handful of quantities
 * from a track. Reproducing exactly those members here lets the ~4500 lines of
 * ported algorithm compile and run unchanged - which is what makes the port
 * verifiable against the published Delphes results.
 *
 * Everything is stored in EDM4hep units: mm, 1/mm, GeV, ns.
 *
 * Origin of the algorithm: J. Kriewald, "No Track left behind: Graph-based
 * Vertexing for long-lived Particle Reconstruction", arXiv:2510.00856,
 * https://github.com/jkriewald/delphes-LLP
 */
#pragma once

#include <TLorentzVector.h>
#include <TMatrixDSym.h>

namespace k4vertexing {

/// Minimal track view: the only track information the graph vertexer consumes.
struct GraphTrack {
  // --- perigee parameters, Delphes/Bedeschi convention, in mm and 1/mm ---
  //     D = d0, Phi = phi0, C = -omega/2, DZ = z0, CtgTheta = tanLambda
  double D0 = 0.0;
  double Phi = 0.0;
  double C = 0.0;
  double DZ = 0.0;
  double CtgTheta = 0.0;

  /// 5x5 perigee covariance in the same convention and units.
  TMatrixDSym TrackCovariance{5};

  // --- kinematics ---
  TLorentzVector Momentum;         ///< used via M(), P(), Beta(), Phi(), Pt()
  TLorentzVector Position;         ///< used via T() for the hit time [ns]
  TLorentzVector InitialPosition;  ///< production point, when available
  double PT = 0.0;
  double Mass = 0.0;
  double L = 0.0;      ///< path length
  double ErrorT = 0.0; ///< time resolution [ns]; <= 0 disables timing
  int Charge = 0;
  int PID = 0;

  /// Index of the originating track in the input EDM4hep collection, so the
  /// Gaudi component can attach the right tracks to each reconstructed vertex.
  int sourceIndex = -1;

  GraphTrack() { TrackCovariance.ResizeTo(5, 5); TrackCovariance.Zero(); }
};

} // namespace k4vertexing

/// The ported algorithm refers to tracks as `Candidate`; keep that spelling so
/// the ~4500 lines of original code need no edits.
using Candidate = k4vertexing::GraphTrack;
