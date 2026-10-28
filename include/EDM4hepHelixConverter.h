/*
 * Part of the Key4hep port of the graph-based displaced vertex finder.
 * Algorithm: J. Kriewald, arXiv:2510.00856 (github.com/jkriewald/delphes-LLP).
 * Port author: Mahmoud Althakeel
 */
/*
 * Shared EDM4hep <-> helix-convention conversion for the Vertexing package.
 *
 * Every vertex algorithm in this package (the graph-based displaced finder now,
 * a Bedeschi/Delphes least-squares fitter and an LCFIPlus-style finder later)
 * needs the same thing: the five perigee parameters and their covariance taken
 * out of an edm4hep::TrackState and expressed in the convention that algorithm
 * was written in. Keeping that in one place means the sign and unit conventions
 * are defined - and validated - exactly once.
 *
 * EDM4hep perigee (mm, 1/mm):   (d0, phi, omega, z0, tanLambda)
 * Delphes/Bedeschi   (mm, 1/mm): (D,  phi0, C,   z0, ct)
 *
 *     D  = d0        phi0 = phi        C = -omega/2
 *     z0 = z0        ct   = tanLambda
 *
 * Only the curvature differs, so the Jacobian is diag(1, 1, -1/2, 1, 1): every
 * covariance element involving omega picks up a factor -1/2, and Var(omega)
 * a factor 1/4.
 */
#pragma once

#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>

#include <edm4hep/TrackState.h>

#include <cmath>

namespace k4vertexing {

/// Index of the EDM4hep perigee parameters, for readability.
enum EDM4hepPerigee { kD0 = 0, kPhi = 1, kOmega = 2, kZ0 = 3, kTanLambda = 4, kNPerigee = 5 };

/// EDM4hep track parameters -> Delphes/Bedeschi convention (D, phi0, C, z0, ct).
inline TVectorD toDelphesPerigee(const edm4hep::TrackState& state) {
  TVectorD parameters(kNPerigee);
  parameters(0) = state.D0;
  parameters(1) = state.phi;
  parameters(2) = -0.5 * state.omega; // C = -omega/2
  parameters(3) = state.Z0;
  parameters(4) = state.tanLambda;
  return parameters;
}

/// Unpack the packed lower-triangular EDM4hep covariance into a full 5x5 matrix.
inline TMatrixDSym unpackCovariance(const edm4hep::TrackState& state) {
  TMatrixDSym covariance(kNPerigee);
  covariance.Zero();
  int packed = 0;
  for (int row = 0; row < kNPerigee; ++row)
    for (int column = 0; column <= row; ++column, ++packed) {
      covariance(row, column) = state.covMatrix[packed];
      covariance(column, row) = covariance(row, column);
    }
  return covariance;
}

/// EDM4hep covariance -> Delphes/Bedeschi convention, via A = diag(1,1,-1/2,1,1).
inline TMatrixDSym toDelphesPerigeeCovariance(const edm4hep::TrackState& state) {
  TMatrixDSym covariance = unpackCovariance(state);
  TMatrixD jacobian(kNPerigee, kNPerigee);
  jacobian.UnitMatrix();
  jacobian(2, 2) = -0.5;
  TMatrixD jacobianTransposed(kNPerigee, kNPerigee);
  jacobianTransposed.Transpose(jacobian);
  covariance.Similarity(jacobianTransposed);
  return covariance;
}

/// True if the covariance is usable (finite, positive diagonal).
inline bool hasUsableCovariance(const TMatrixDSym& covariance) {
  for (int i = 0; i < covariance.GetNrows(); ++i) {
    if (!std::isfinite(covariance(i, i)) || covariance(i, i) <= 0.0)
      return false;
  }
  return true;
}

/// Transverse momentum [GeV] from the EDM4hep curvature and the solenoid field.
/// pT = 0.299792458e-3 * B[T] / |omega|[1/mm]
inline double transverseMomentum(double omega, double magneticFieldZ) {
  constexpr double kPtConversion = 0.000299792458;
  return (omega == 0.0) ? 0.0 : kPtConversion * magneticFieldZ / std::abs(omega);
}

} // namespace k4vertexing
