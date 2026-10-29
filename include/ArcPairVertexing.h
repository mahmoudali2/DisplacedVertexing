/*
 * Two-track vertexing with LOCAL track fits: each track is described from its
 * own hits and extrapolated back only as far as the decay, instead of using a
 * track state at a perigee that may be metres away from the decay.
 *
 * Nothing in the code is specific to a detector: the field, the material and
 * the region are configuration (see Config).
 *
 *   track  : circle in xy + straight line z(s), generalised least-squares fits
 *            in a frame attached to the innermost hit; the hit covariance
 *            includes the correlated multiple scattering between hits
 *   errors : fit covariance propagated to the extrapolation point, plus the
 *            multiple scattering in the material between that point and the
 *            first hit, plus an optional field-mismatch term
 *   vertex : Gauss-Newton fit with two residuals per track (xy distance
 *            perpendicular to the track, z at the same s) -> position,
 *            3x3 covariance, chi2 with ndf = 2N - 3
 *
 * Material model (radiationLengthsPerMm): a "shell" between an inner box (r_xy < shellInnerRadius,
 * |z| < shellInnerHalfLength) and an outer box (shellOuterRadius, shellOuterHalfLength) holding shellRadiationLengths
 * radiation lengths across its radial depth, plus a uniform density
 * outsideRadiationLengths / outsideThickness beyond the outer box. No material inside the
 * inner box. For example shell = calorimeter, outside = iron yoke.
 *
 * Units: mm, GeV, Tesla.
 * Author: Mahmoud Althakeel
 */
#pragma once

#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TVector3.h>

#include <optional>
#include <vector>

namespace k4vertexing::arcpair {

struct Config {
  // ---- detector: magnetic field --------------------------------------------
  double fieldForMomentum = 0.0;    // T, converts the input track omega [1/mm] into pT (must be set)
  // field mismatch between where the circle is fitted and the region it is
  // extrapolated through; adds 0.5*0.3*dB/pT*s^2 to the bending-plane error
  double fieldUncertainty = 0.0;    // T
  // ---- detector: resolutions --------------------------------------------------
  double hitResolution = 1.0;       // mm per coordinate
  double momentumResolution = 0.0;  // sigma(p)/p of the input tracks
  // ---- detector: material (multiple scattering in the track errors) -----------
  bool useMaterial = true;
  // material shell (between inner and outer box) and the region beyond it
  double shellRadiationLengths = 0.0;
  double shellInnerRadius = 0.0, shellOuterRadius = 0.0;
  double shellInnerHalfLength = 0.0, shellOuterHalfLength = 0.0;
  double outsideRadiationLengths = 0.0;
  double outsideThickness = 1.0;    // mm, depth over which outsideRadiationLengths is spread
  // ---- physics: momentum used for the multiple-scattering estimate -------------
  double defaultTransverseMomentum = 10.0; // GeV, used when the input track has omega = 0
  double minimumMomentum = 2.0;     // GeV, floor
  // ---- algorithm: track model and vertex fit -----------------------------------
  bool fitCurvature = true;
  unsigned minimumHits = 3;
  double upstreamTolerance = 100.0; // mm, vertex may lie this far past a first hit
  int integrationSteps = 40;        // steps of the material integrals
  int maximumIterations = 12;       // Gauss-Newton iterations of the vertex fit
  double jacobianStep = 1.0;        // mm, step of the numerical Jacobian
  double convergenceDistance = 0.1; // mm, stop when the vertex moves less than this
};

/// Radiation lengths per mm of path at a point (shell + outside material).
double radiationLengthsPerMm(const TVector3& position, const Config& config);

class ArcTrack {
public:
  /// hits in mm; omega [1/mm] only sets the momentum for the scattering model
  bool build(std::vector<TVector3> hits, double omega, const Config& config);

  /// point of the track at arc position s
  TVector3 point(double arcPosition) const;
  /// local coordinates (arc position s along the track, offset u perpendicular to it)
  /// of a point, in the frame at the first hit
  void localCoordinates(const TVector3& position, double& arcPosition, double& normalOffset) const;
  /// offset u(s) of the fitted circle and its slope du/ds
  void offsetAndSlope(double arcPosition, double& offset, double& slope) const;
  double zAt(double arcPosition) const { return m_firstHit.Z() + m_parametersZ[0] + m_parametersZ[1] * arcPosition; }
  /// uncertainty of the trajectory perpendicular to it in xy, and in z, at arc position s
  void uncertainties(double arcPosition, const Config& config, double& sigmaXY, double& sigmaZ) const;

  const TVector3& firstHit() const { return m_firstHit; }
  double outwardDirectionX() const { return m_outwardDirection[0]; }
  double outwardDirectionY() const { return m_outwardDirection[1]; }
  double zSlope() const { return m_parametersZ[1]; }
  double momentum() const { return m_momentum; }
  double fitChi2() const { return m_chi2; }
  int numberOfHits() const { return m_numberOfHits; }

private:
  TVector3 m_firstHit;
  // local frame at the first hit: outward direction and normal in xy
  double m_outwardDirection[2] = {1, 0}, m_normalDirection[2] = {0, 1};
  std::vector<double> m_parametersXY; // offset a, slope b, (curvature)
  TMatrixDSym m_parameterCovarianceXY;
  double m_parametersZ[2] = {0, 0};   // z offset c, slope d
  TMatrixDSym m_parameterCovarianceZ{2};
  double m_cosLambda = 1, m_momentum = 0, m_transverseMomentum = 0, m_scatteringFactor = 0, m_chi2 = 0;
  int m_numberOfHits = 0;
};

struct VertexResult {
  TVector3 position;
  TMatrixDSym covariance{3};
  double chi2 = 0;
  int degreesOfFreedom = 0;
  double pairDistanceOfClosestApproach = 0; // distance between the two trajectories at the vertex
  std::vector<double> arcPositionOfVertex;  // per track: arc position s of the vertex (negative = upstream)
  bool upstream = false;          // vertex before every track's first hit
  bool betweenIPAndHits = false;  // closer to the origin than every first hit, same side
};

std::optional<VertexResult> fitVertex(const std::vector<const ArcTrack*>& tracks, const Config& config);

} // namespace k4vertexing::arcpair
