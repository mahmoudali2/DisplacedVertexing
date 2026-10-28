/*
 * Two-track vertexing with LOCAL track fits: each track is described from its
 * own hits and extrapolated back only as far as the decay, instead of using a
 * track state at a perigee that may be metres away from the decay.
 *
 * Written for long-lived-particle decays behind the IDEA solenoid with muon-
 * system tracks, but nothing in the code is specific to a detector: the field,
 * the material and the region are configuration (see Config).
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
 * Material model (x0PerMM): a "shell" between an inner box (r_xy < shellInnerR,
 * |z| < shellInnerZ) and an outer box (shellOuterR, shellOuterZ) holding shellX0
 * radiation lengths across its radial depth, plus a uniform density
 * outsideX0 / outsideThickness beyond the outer box. No material inside the
 * inner box. For IDEA: shell = dual-readout calorimeter, outside = iron yoke.
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
  // momentum used for the multiple-scattering estimate
  double bFieldForMomentum = 0.0;   // T, converts the input track omega [1/mm] into pT (must be set)
  double defaultPt = 10.0;          // GeV, used when the input track has omega = 0
  double momentumResolution = 0.0;  // sigma(p)/p of the input tracks
  double minMomentum = 2.0;         // GeV, floor
  double hitResolution = 1.0;       // mm per coordinate
  // field mismatch between where the circle is fitted and the region it is
  // extrapolated through; adds 0.5*0.3*dB/pT*s^2 to the bending-plane error
  double fieldUncertainty = 0.0;    // T
  bool useMaterial = true;          // multiple scattering in the track errors
  // material shell (between inner and outer box) and the region beyond it
  double shellX0 = 0.0;
  double shellInnerR = 0.0, shellOuterR = 0.0;
  double shellInnerZ = 0.0, shellOuterZ = 0.0;
  double outsideX0 = 0.0;
  double outsideThickness = 1.0;    // mm, depth over which outsideX0 is spread
  bool fitCurvature = true;
  unsigned minHits = 3;
  // topology
  double upstreamTolerance = 100.0; // mm, vertex may lie this far past a first hit
  // numerics
  int nSteps = 40;                  // steps of the material integrals
  int maxIterations = 12;           // Gauss-Newton iterations of the vertex fit
  double jacobianStep = 1.0;        // mm, step of the numerical Jacobian
  double convergenceDistance = 0.1; // mm, stop when the vertex moves less than this
};

/// Radiation lengths per mm of path at a point (shell + outside material).
double x0PerMM(const TVector3& point, const Config& cfg);

class ArcTrack {
public:
  /// hits in mm; omega [1/mm] only sets the momentum for the scattering model
  bool build(std::vector<TVector3> hits, double omega, const Config& cfg);

  TVector3 point(double s) const;
  /// local coordinates (s, u) of a point in the frame at the first hit
  void local(const TVector3& v, double& s, double& u) const;
  /// u(s) of the fitted circle and its slope du/ds
  void uOfS(double s, double& u, double& slope) const;
  double zOfS(double s) const { return m_h0.Z() + m_pz[0] + m_pz[1] * s; }
  /// sigma of the trajectory perpendicular to it in xy, and in z, at local s
  void sigmas(double s, const Config& cfg, double& sigmaU, double& sigmaZ) const;

  const TVector3& firstHit() const { return m_h0; }
  double dirX() const { return m_t0[0]; }
  double dirY() const { return m_t0[1]; }
  double slopeZ() const { return m_pz[1]; }
  double momentum() const { return m_p; }
  double fitChi2() const { return m_chi2; }
  int nHits() const { return m_nHits; }

private:
  TVector3 m_h0;
  double m_t0[2] = {1, 0}, m_n0[2] = {0, 1};
  std::vector<double> m_pu; // a, b, (kappa)
  TMatrixDSym m_cu;
  double m_pz[2] = {0, 0};
  TMatrixDSym m_cz{2};
  double m_cosLambda = 1, m_p = 0, m_pt = 0, m_msFactor = 0, m_chi2 = 0;
  int m_nHits = 0;
};

struct VertexResult {
  TVector3 position;
  TMatrixDSym covariance{3};
  double chi2 = 0;
  int ndf = 0;
  double pairDCA = 0;            // distance between the two trajectories at the vertex
  std::vector<double> sFromFirst; // per track: local s of the vertex (negative = upstream)
  bool upstream = false;          // vertex before every track's first hit
  bool betweenIPAndHits = false;  // closer to the origin than every first hit, same side
};

std::optional<VertexResult> fitVertex(const std::vector<const ArcTrack*>& tracks, const Config& cfg);

} // namespace k4vertexing::arcpair
