/*
 * Geometry-free 3D two-track vertexing with stepped helix propagation.
 *
 * Every track is treated as a trajectory in 3D, so tracks crossing endcap discs
 * (small xy motion, large z motion) are handled exactly like barrel tracks:
 *
 *   track  : helix from a fitted track state (d0, phi, omega, z0, tanLambda),
 *            starting at the innermost track hit. A state without curvature,
 *            or a start point in zero field, gives a 3D line.
 *   propagation : stepped INWARD from the start point through a two-zone
 *            field map (innerField inside a cylinder, outerField outside) while
 *            the ionisation loss of a material map is added back, so the
 *            curvature follows the local field and momentum.
 *   vertex : point of closest approach of the two trajectories (polyline
 *            search + segment refinement), DCA = distance; momenta and
 *            directions there give the pair mass. When both states carry a
 *            covariance the vertex is then fitted (see fitVertex): each track's
 *            state covariance is carried to the trajectory point closest to the
 *            vertex with a numerical Jacobian through the same propagation, and
 *            the two tracks are combined in the planes perpendicular to them.
 *            This gives the vertex position, its 3x3 covariance and a chi2 with
 *            ndf = 2 x 2 - 3 = 1.
 *
 * Material map (energyLossPerMmAt): a "shell" between an inner box (r_xy < shellInner,
 * |z| < shellInner) and an outer box (shellOuter), outside a beam hole
 * r_xy < shellBeamHole, with energy loss energyLossInShell; energyLossOutside beyond
 * the outer box; nothing inside the inner box. For example shell = calorimeter,
 * outside = iron yoke. A uniform field and no energy loss describe a plain tracker.
 *
 * No region is assumed: the caller decides where vertices are accepted.
 *
 * Units: mm, GeV, Tesla.
 * Author: Mahmoud Althakeel
 */
#pragma once

#include <TVector3.h>

#include <optional>
#include <vector>

namespace k4vertexing::helixpair {

struct Config {
  // ---- detector: magnetic field --------------------------------------------
  // innerField for r_xy < fieldRadius and |z| < fieldHalfLength, outerField elsewhere
  double innerField = 0.0;         // T, signed
  double outerField = 0.0;         // T, signed
  double fieldRadius = 1.0e12;     // mm (default: innerField everywhere)
  double fieldHalfLength = 1.0e12; // mm
  // ---- detector: material (energy loss, see above) --------------------------
  double shellInner = 0.0;         // mm
  double shellOuter = 0.0;         // mm
  double shellBeamHole = 0.0;      // mm
  double energyLossInShell = 0.0;  // GeV/mm
  double energyLossOutside = 0.0;  // GeV/mm
  // ---- physics: decay products ------------------------------------------------
  double particleMass = 0.1056584; // GeV, mass hypothesis of both tracks for the pair mass
  // ---- algorithm: propagation and vertex search --------------------------------
  double stepLength = 20.0;        // mm
  double backwardFactor = 1.3;     // propagate inward up to backwardFactor x |start point|
  double upstreamTolerance = 100.0; // mm, vertex may lie this far past a start point
};

/// The five helix parameters of a fitted track state (EDM4hep / LCIO convention)
/// and the point they refer to. Units: mm, 1/mm, rad.
struct TrackStateParameters {
  double d0 = 0, phi = 0, omega = 0, z0 = 0, tanLambda = 0;
  TVector3 referencePoint; // parameters describe the point of closest approach to this point
  // covariance of (d0, phi, omega, z0, tanLambda); hasCovariance = false if not given
  double covariance[5][5] = {};
  bool hasCovariance = false;
};

class HelixTrack {
public:
  /// helix from a fitted track state. The innermost hit (if any) only fixes where
  /// the trajectory starts (path length 0); without hits the point of closest approach
  /// to the reference point is used (for a state at the first hit that is the first hit).
  bool build(const TrackStateParameters& state, const std::vector<TVector3>& hits, const Config& config);

  /// helix (or line) at 3D path length from the start point (negative: inward)
  TVector3 point(double pathLength) const;
  /// unit vector of the outward motion at a path length
  TVector3 tangent(double pathLength) const;

  /// stepped inward propagation; path lengths <= upstreamTolerance, ordered from the innermost point outward
  struct Trajectory {
    std::vector<double> pathLength, momentum;
    std::vector<TVector3> position, direction;
  };
  Trajectory propagate(const Config& config) const;

  bool isStraightLine() const { return m_isStraightLine; }
  /// rotation sense of the outward motion in xy (+1 counter-clockwise, 0 for a line);
  /// equal for equal charges in the same field, so only the product is meaningful
  int rotationSense() const { return m_isStraightLine ? 0 : (m_rotationSense > 0 ? 1 : -1); }
  double radius() const { return m_radius; }
  /// momentum at the start point [GeV] (NaN for a line)
  double measuredMomentum() const { return m_measuredMomentum; }
  const TVector3& startPoint() const { return m_startPoint; }
  int numberOfHits() const { return m_numberOfHits; }
  /// the track state this helix was built from
  const TrackStateParameters& state() const { return m_state; }

private:
  TrackStateParameters m_state;
  TVector3 m_startPoint;
  bool m_isStraightLine = true;
  // helix: circle centre and radius, z change per radian of turning, angle and z at the
  // start point, 3D path length per radian, rotation sense (+1 / -1)
  double m_centreX = 0, m_centreY = 0, m_radius = 0, m_zPerRadian = 0;
  double m_angleAtStart = 0, m_zAtStart = 0, m_pathLengthPerRadian = 1, m_rotationSense = 0;
  // field at the start point (gives the momentum from the curvature)
  double m_fieldAtStart = 0;
  // straight line
  TVector3 m_lineStart, m_lineDirection;
  double m_measuredMomentum = 0;
  int m_numberOfHits = 0;
};

double fieldZAt(const TVector3& position, const Config& config);
double energyLossPerMmAt(const TVector3& position, const Config& config);

struct VertexResult {
  TVector3 position;
  double pairDistanceOfClosestApproach = 0;  // distance between the two trajectories at the vertex
  double pathLengthFromStart[2] = {0, 0};    // path length of the vertex from each start point (negative = upstream)
  bool upstream = false;          // vertex before (within the tolerance) every track's start point
  bool betweenIPAndHits = false;  // closer to the origin than every start point, same side
  double mass = 0;                // pair mass at the vertex (NaN if a track is a line)
  double cosOpening = 0;          // cosine of the opening angle at the vertex
  double momentumAtVertex[2] = {0, 0};  // momenta at the vertex (energy loss added back)
  TVector3 closestApproachMidpoint;     // midpoint of the closest approach of the two trajectories
  // vertex fit (only when both states have a covariance)
  bool isFitted = false;
  double chi2 = -1;               // ndf = 1
  double covariance[3][3] = {};   // covariance of the position [mm^2]
};

std::optional<VertexResult> fitVertex(const HelixTrack& firstTrack, const HelixTrack& secondTrack,
                                      const Config& config);

} // namespace k4vertexing::helixpair
