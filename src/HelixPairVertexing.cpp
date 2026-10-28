/*
 * Geometry-free 3D two-track vertexing, see HelixPairVertexing.h.
 * Author: Mahmoud Althakeel
 */
#include "HelixPairVertexing.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace k4vertexing::helixpair {

namespace {
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

/// closest approach of the segments p0-p1 and q0-q1: fractions s, t and the distance
double segmentPOCA(const TVector3& p0, const TVector3& p1, const TVector3& q0, const TVector3& q1, double& s,
                   double& t) {
  const TVector3 u = p1 - p0, v = q1 - q0, w = p0 - q0;
  const double a = u.Dot(u), b = u.Dot(v), c = v.Dot(v), d = u.Dot(w), e = v.Dot(w);
  const double den = a * c - b * b;
  auto clip = [](double x) { return std::clamp(x, 0.0, 1.0); };
  t = den > 1e-12 * a * c ? clip((a * e - b * d) / den) : (c > 0 ? clip(e / c) : 0.0);
  s = a > 0 ? clip((b * t - d) / a) : 0.0;
  t = c > 0 ? clip((b * s + e) / c) : 0.0;
  return (p0 + s * u - (q0 + t * v)).Mag();
}
} // namespace

double fieldAt(const TVector3& p, const Config& cfg) {
  return (p.Perp() < cfg.fieldR && std::fabs(p.Z()) < cfg.fieldHalfZ) ? cfg.innerField : cfg.outerField;
}

double dEdxAt(const TVector3& p, const Config& cfg) {
  const double r = p.Perp(), az = std::fabs(p.Z());
  const bool inOuter = r < cfg.shellOuter && az < cfg.shellOuter;
  const bool inInner = r < cfg.shellInner && az < cfg.shellInner;
  // in the shell (outside the beam hole)
  if (inOuter && !inInner && r > cfg.shellHole)
    return cfg.dEdxShell;
  // beyond the shell
  if (!inOuter)
    return cfg.dEdxOutside;
  // inside the inner box or in the beam hole: no energy loss
  return 0.0;
}

bool HelixTrack::build(const StateParams& st, const std::vector<TVector3>& hits, const Config& cfg) {
  // the state describes the track; the hits (possibly none, e.g. a drift-chamber
  // track without silicon hits) only fix where the trajectory starts
  m_nHits = static_cast<int>(hits.size());

  // point of closest approach (PCA) to the reference point, EDM4hep convention
  const double sinPhi = std::sin(st.phi);
  const double cosPhi = std::cos(st.phi);
  const TVector3 pca(st.reference.X() - st.d0 * sinPhi, st.reference.Y() + st.d0 * cosPhi,
                     st.reference.Z() + st.z0);

  // the trajectory starts at the innermost hit (closest to the origin)
  TVector3 start = pca;
  if (!hits.empty()) {
    start = hits.front();
    for (const auto& h : hits)
      if (h.Mag() < start.Mag())
        start = h;
  }
  m_h0 = start;
  m_fitField = fieldAt(m_h0, cfg);

  // ---- straight line (no curvature or no field) ------------------------------
  if (st.omega == 0 || m_fitField == 0) {
    m_line = true;
    m_dir = TVector3(cosPhi, sinPhi, st.tanLambda).Unit();
    // base = point of the line closest to the start point
    m_base = pca + (m_h0 - pca).Dot(m_dir) * m_dir;
    m_R = std::numeric_limits<double>::infinity();
    m_k = 0;
    m_h = 0;
    m_pMeas = kNaN;
    return true;
  }

  // ---- helix ------------------------------------------------------------------
  // omega > 0: the circle centre lies to the right of the direction of motion,
  // i.e. the track turns clockwise (rotation sense -1)
  const double signedRadius = 1.0 / st.omega;
  m_R = std::fabs(signedRadius);
  m_xc = pca.X() + signedRadius * sinPhi;
  m_yc = pca.Y() - signedRadius * cosPhi;
  m_h = st.omega > 0 ? -1.0 : 1.0;
  // z changes by tanLambda * R per radian of turning
  m_k = st.tanLambda * m_R;
  m_ds = std::sqrt(m_R * m_R + m_k * m_k);

  // phase of the start point on the circle, and z there
  const double phasePCA = std::atan2(pca.Y() - m_yc, pca.X() - m_xc);
  const double phaseStart = std::atan2(m_h0.Y() - m_yc, m_h0.X() - m_xc);
  // turning angle from the PCA to the start point, in the direction of motion
  double turn = m_h * (phaseStart - phasePCA);
  while (turn > M_PI)
    turn -= 2 * M_PI;
  while (turn < -M_PI)
    turn += 2 * M_PI;
  m_phi0 = phaseStart;
  m_z0 = pca.Z() + m_k * turn;
  m_line = false;

  // momentum from the curvature in the field at the start point
  const double pt = 0.3 * std::fabs(m_fitField) * 1e-3 / std::fabs(st.omega);
  m_pMeas = pt * std::sqrt(1.0 + st.tanLambda * st.tanLambda);
  return true;
}

TVector3 HelixTrack::point(double t) const {
  if (m_line)
    return m_base + t * m_dir;
  const double u = t / m_ds, ph = m_phi0 + m_h * u;
  return TVector3(m_xc + m_R * std::cos(ph), m_yc + m_R * std::sin(ph), m_z0 + m_k * u);
}

TVector3 HelixTrack::tangent(double t) const {
  if (m_line)
    return m_dir;
  const double ph = m_phi0 + m_h * t / m_ds;
  return TVector3(-m_R * std::sin(ph) * m_h, m_R * std::cos(ph) * m_h, m_k).Unit();
}

HelixTrack::Trajectory HelixTrack::propagate(const Config& cfg) const {
  Trajectory tr;
  // fitted helix up to the tolerance outward (points at tol, tol - step, ... > 0, in increasing order)
  std::vector<double> tf;
  for (double v = cfg.upstreamTolerance; v > 0; v -= cfg.step)
    tf.push_back(v);
  std::reverse(tf.begin(), tf.end());
  // inward stepping
  const int n = static_cast<int>(cfg.backFactor * m_h0.Mag() / cfg.step);
  TVector3 pos = point(0.0), d = tangent(0.0);
  const double cl = std::hypot(d.X(), d.Y()), tl = d.Z() / std::max(cl, 1e-12);
  double phi = std::atan2(d.Y(), d.X());
  double p = m_pMeas;
  const double fitFieldSign = m_fitField < 0 ? -1.0 : 1.0;
  std::vector<TVector3> X{pos}, D{d};
  std::vector<double> P{p};
  for (int i = 0; i < n; ++i) {
    if (!m_line) {
      const double pt = p * cl;
      // signed xy curvature of the outward motion: m_h / R where the hits are, scales with B / pT
      const double kappa = m_h * fitFieldSign * fieldAt(pos, cfg) * 0.3e-3 / std::max(pt, 1e-3);
      pos -= cfg.step * TVector3(cl * std::cos(phi), cl * std::sin(phi), cl * tl);
      phi -= kappa * cfg.step * cl;
      p += dEdxAt(pos, cfg) * cfg.step;
      d = TVector3(cl * std::cos(phi), cl * std::sin(phi), cl * tl);
    } else {
      pos -= cfg.step * d;
    }
    X.push_back(pos);
    D.push_back(d);
    P.push_back(p);
  }
  const std::size_t nb = X.size();
  tr.t.reserve(nb + tf.size());
  for (std::size_t i = 0; i < nb; ++i) {
    const std::size_t j = nb - 1 - i;
    tr.t.push_back(-cfg.step * static_cast<double>(j));
    tr.x.push_back(X[j]);
    tr.dir.push_back(D[j]);
    tr.p.push_back(P[j]);
  }
  for (double v : tf) {
    tr.t.push_back(v);
    tr.x.push_back(point(v));
    tr.dir.push_back(tangent(v));
    tr.p.push_back(m_pMeas);
  }
  return tr;
}

std::optional<VertexResult> fitVertex(const HelixTrack& a, const HelixTrack& b, const Config& cfg) {
  const auto ta = a.propagate(cfg), tb = b.propagate(cfg);
  const std::size_t na = ta.x.size(), nb = tb.x.size();
  if (na < 2 || nb < 2)
    return std::nullopt;
  // coarse: closest pair of trajectory points
  std::size_t i = 0, j = 0;
  double best2 = std::numeric_limits<double>::infinity();
  for (std::size_t ia = 0; ia < na; ++ia)
    for (std::size_t ib = 0; ib < nb; ++ib) {
      const double d2 = (ta.x[ia] - tb.x[ib]).Mag2();
      if (d2 < best2) {
        best2 = d2;
        i = ia;
        j = ib;
      }
    }
  // fine: closest approach of the neighbouring segments
  double dca = std::numeric_limits<double>::infinity(), s = 0, t = 0;
  std::size_t ii = 0, jj = 0;
  for (long ci : {static_cast<long>(i) - 1, static_cast<long>(i)})
    for (long cj : {static_cast<long>(j) - 1, static_cast<long>(j)}) {
      if (ci < 0 || cj < 0 || ci >= static_cast<long>(na) - 1 || cj >= static_cast<long>(nb) - 1)
        continue;
      double ss, tt;
      const double dd = segmentPOCA(ta.x[ci], ta.x[ci + 1], tb.x[cj], tb.x[cj + 1], ss, tt);
      if (dd < dca) {
        dca = dd;
        ii = ci;
        jj = cj;
        s = ss;
        t = tt;
      }
    }
  if (!std::isfinite(dca)) {
    dca = std::sqrt(best2);
    ii = std::min(i, na - 2);
    jj = std::min(j, nb - 2);
    s = t = 0;
  }
  const TVector3 xa = ta.x[ii] + s * (ta.x[ii + 1] - ta.x[ii]);
  const TVector3 xb = tb.x[jj] + t * (tb.x[jj + 1] - tb.x[jj]);
  VertexResult out;
  out.position = 0.5 * (xa + xb);
  out.pairDCA = dca;
  out.t[0] = ta.t[ii] + s * (ta.t[ii + 1] - ta.t[ii]);
  out.t[1] = tb.t[jj] + t * (tb.t[jj + 1] - tb.t[jj]);
  out.upstream = out.t[0] < cfg.upstreamTolerance && out.t[1] < cfg.upstreamTolerance;
  out.betweenIPAndHits = true;
  for (const HelixTrack* trk : {&a, &b})
    if (!(out.position.Mag() < trk->firstHit().Mag() && out.position.Dot(trk->firstHit()) > 0))
      out.betweenIPAndHits = false;
  const TVector3 &da = ta.dir[ii], &db = tb.dir[jj];
  if (a.isLine() || b.isLine()) {
    out.mass = kNaN;
    out.cosOpening = da.Dot(db);
    out.p[0] = out.p[1] = kNaN;
    return out;
  }
  const TVector3 pa = ta.p[ii] * da, pb = tb.p[jj] * db;
  const double m2 = cfg.particleMass * cfg.particleMass;
  const double E = std::sqrt(pa.Mag2() + m2) + std::sqrt(pb.Mag2() + m2);
  out.mass = std::sqrt(std::max(E * E - (pa + pb).Mag2(), 0.0));
  out.cosOpening = pa.Dot(pb) / (pa.Mag() * pb.Mag());
  out.p[0] = ta.p[ii];
  out.p[1] = tb.p[jj];
  return out;
}

} // namespace k4vertexing::helixpair
