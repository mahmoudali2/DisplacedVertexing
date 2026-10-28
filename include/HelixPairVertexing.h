/*
 * Geometry-free 3D two-track vertexing with stepped helix propagation.
 *
 * Unlike the arc vertexer (ArcPairVertexing.h), which describes a track by its
 * transverse distance s from the first hit with z linear in s, every track is
 * treated as a trajectory in 3D, so tracks crossing endcap discs (small xy
 * motion, large z motion) are handled exactly like barrel tracks:
 *
 *   track  : helix fitted to its own hits in the field at its innermost hit
 *            (circle in xy: algebraic fit + geometric Gauss-Newton refinement;
 *            z linear in the turning angle). Tracks without measurable xy
 *            curvature, or with hits in zero field, fall back to a 3D line.
 *   propagation : stepped INWARD from the innermost hit through a two-zone
 *            field map (innerField inside a cylinder, outerField outside) while
 *            the ionisation loss of a material map is added back, so the
 *            curvature follows the local field and momentum.
 *   vertex : point of closest approach of the two trajectories (polyline
 *            search + segment refinement), vertex = midpoint, DCA = distance;
 *            momenta and directions at the vertex give the pair mass.
 *
 * Material map (dEdxAt): a "shell" between an inner box (r_xy < shellInner,
 * |z| < shellInner) and an outer box (shellOuter), outside a beam hole
 * r_xy < shellHole, with energy loss dEdxShell; dEdxOutside beyond the outer
 * box; nothing inside the inner box. For IDEA: shell = dual-readout
 * calorimeter, outside = iron yoke. A uniform field and zero dE/dx describe a
 * plain tracker.
 *
 * No region is assumed: the caller decides where vertices are accepted.
 * Validated (IDEA muon system) against the Python reference
 * (HNL/MS_full_outside/code/helix3d.py, pair_vertex_stepped) on HNL signal: at
 * m_N = 20 GeV the vertex is within 48 mm (barrel) / 72 mm (forward) of the
 * true decay (median), p_reco/p_true = 1.00, m_reco/m_true 0.88-1.25 (16-84 %).
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
  // field map: innerField for r_xy < fieldR and |z| < fieldHalfZ, outerField elsewhere
  double innerField = 0.0;     // T, signed
  double outerField = 0.0;     // T, signed
  double fieldR = 1.0e12;      // mm (default: innerField everywhere)
  double fieldHalfZ = 1.0e12;  // mm
  // material map for the energy loss (see above)
  double shellInner = 0.0, shellOuter = 0.0, shellHole = 0.0; // mm
  double dEdxShell = 0.0;      // GeV/mm
  double dEdxOutside = 0.0;    // GeV/mm
  // decay products
  double particleMass = 0.1056584; // GeV, mass hypothesis for the pair mass (default muon)
  // track model
  unsigned minHits = 3;
  unsigned maxHitsInFit = 0;   // use only the N innermost hits (0 = all)
  double minXYSpan = 30.0;     // mm, below: straight line
  double maxRadius = 2.0e6;    // mm, above: straight line
  // propagation and vertex search
  double step = 20.0;          // mm
  double backFactor = 1.3;     // propagate inward up to backFactor x |first hit|
  double upstreamTolerance = 100.0; // mm, vertex may lie this far past a first hit
};

class HelixTrack {
public:
  /// hits in mm
  bool build(std::vector<TVector3> hits, const Config& cfg);

  /// fitted helix (or line) at 3D path length t from the innermost hit (t < 0: inward)
  TVector3 point(double t) const;
  /// unit tangent of the outward motion at t
  TVector3 tangent(double t) const;

  /// stepped inward propagation; t <= upstreamTolerance, ordered from the innermost point outward
  struct Trajectory {
    std::vector<double> t, p;
    std::vector<TVector3> x, dir;
  };
  Trajectory propagate(const Config& cfg) const;

  bool isLine() const { return m_line; }
  /// rotation sense of the outward motion in xy (+1 counter-clockwise, 0 for a line);
  /// equal for equal charges in the same field, so only the product is meaningful
  int rotation() const { return m_line ? 0 : (m_h > 0 ? 1 : -1); }
  double radius() const { return m_R; }
  /// momentum at the innermost hit [GeV] (NaN for a line)
  double momentum() const { return m_pMeas; }
  const TVector3& firstHit() const { return m_h0; }
  int nHits() const { return m_nHits; }

private:
  TVector3 m_h0;
  bool m_line = true;
  // helix
  double m_xc = 0, m_yc = 0, m_R = 0, m_k = 0, m_phi0 = 0, m_z0 = 0, m_ds = 1, m_h = 0;
  // field the helix was fitted in (field at the innermost hit)
  double m_fitField = 0;
  // line
  TVector3 m_base, m_dir;
  double m_pMeas = 0;
  int m_nHits = 0;
};

double fieldAt(const TVector3& p, const Config& cfg);
double dEdxAt(const TVector3& p, const Config& cfg);

struct VertexResult {
  TVector3 position;
  double pairDCA = 0;             // distance between the two trajectories at the vertex
  double t[2] = {0, 0};           // path length of the vertex from each first hit (negative = upstream)
  bool upstream = false;          // vertex before (within the tolerance) every track's first hit
  bool betweenIPAndHits = false;  // closer to the origin than every first hit, same side
  double mass = 0;                // pair mass at the vertex (NaN if a track is a line)
  double cosOpening = 0;          // cosine of the opening angle at the vertex
  double p[2] = {0, 0};           // momenta at the vertex (energy loss added back)
};

std::optional<VertexResult> fitVertex(const HelixTrack& a, const HelixTrack& b, const Config& cfg);

} // namespace k4vertexing::helixpair
