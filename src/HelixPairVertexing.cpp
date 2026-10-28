/*
 * Geometry-free 3D two-track vertexing, see HelixPairVertexing.h.
 * Author: Mahmoud Althakeel
 */
#include "HelixPairVertexing.h"

#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>

namespace k4vertexing::helixpair {

namespace {
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

/// least squares A x = y for three parameters, through the normal equations
bool solve3(const std::vector<std::array<double, 3>>& A, const std::vector<double>& y, double x[3]) {
  double M[3][4] = {};
  for (std::size_t r = 0; r < A.size(); ++r)
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j)
        M[i][j] += A[r][i] * A[r][j];
      M[i][3] += A[r][i] * y[r];
    }
  for (int c = 0; c < 3; ++c) { // Gaussian elimination with partial pivoting
    int piv = c;
    for (int r = c + 1; r < 3; ++r)
      if (std::fabs(M[r][c]) > std::fabs(M[piv][c]))
        piv = r;
    if (std::fabs(M[piv][c]) < 1e-300)
      return false;
    std::swap(M[c], M[piv]);
    for (int r = 0; r < 3; ++r) {
      if (r == c)
        continue;
      const double f = M[r][c] / M[c][c];
      for (int k = c; k < 4; ++k)
        M[r][k] -= f * M[c][k];
    }
  }
  for (int i = 0; i < 3; ++i)
    x[i] = M[i][3] / M[i][i];
  return true;
}

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

bool HelixTrack::build(std::vector<TVector3> hits, const Config& cfg) {
  if (hits.size() < cfg.minHits)
    return false;
  std::sort(hits.begin(), hits.end(), [](const TVector3& a, const TVector3& b) { return a.Mag() < b.Mag(); });
  // optionally keep only the innermost hits (closest to the decay)
  if (cfg.maxHitsInFit > 0 && hits.size() > cfg.maxHitsInFit)
    hits.resize(cfg.maxHitsInFit);
  const int n = static_cast<int>(hits.size());
  m_nHits = n;
  m_h0 = hits.front();
  // the helix is fitted in the field at the innermost hit
  m_fitField = fieldAt(m_h0, cfg);

  double span = 0;
  for (const auto& h : hits)
    span = std::max(span, std::hypot(h.X() - m_h0.X(), h.Y() - m_h0.Y()));
  m_line = true;
  // no field at the hits: no curvature to fit, straight line
  if (span >= cfg.minXYSpan && m_fitField != 0) {
    // algebraic circle fit (coordinates centred for conditioning)
    double mx = 0, my = 0;
    for (const auto& h : hits) {
      mx += h.X() / n;
      my += h.Y() / n;
    }
    std::vector<std::array<double, 3>> A(n);
    std::vector<double> y(n);
    for (int i = 0; i < n; ++i) {
      const double x = hits[i].X() - mx, yy = hits[i].Y() - my;
      A[i] = {2 * x, 2 * yy, 1.0};
      y[i] = x * x + yy * yy;
    }
    double sol[3];
    bool ok = solve3(A, y, sol);
    double xc = sol[0] + mx, yc = sol[1] + my;
    double R = ok ? std::sqrt(std::max(sol[2] + sol[0] * sol[0] + sol[1] * sol[1], 1e-9)) : kNaN;
    // geometric refinement (Gauss-Newton on the distances to the circle)
    for (int it = 0; ok && it < 8; ++it) {
      for (int i = 0; i < n; ++i) {
        const double dx = hits[i].X() - xc, dy = hits[i].Y() - yc, d = std::hypot(dx, dy);
        A[i] = {-dx / d, -dy / d, -1.0};
        y[i] = -(d - R);
      }
      double st[3];
      if (!solve3(A, y, st))
        break;
      xc += st[0];
      yc += st[1];
      R += st[2];
      if (std::max({std::fabs(st[0]), std::fabs(st[1]), std::fabs(st[2])}) < 1e-3)
        break;
    }
    if (ok && std::isfinite(R) && R > 0 && R < cfg.maxRadius) {
      std::vector<double> ph(n);
      for (int i = 0; i < n; ++i) {
        ph[i] = std::atan2(hits[i].Y() - yc, hits[i].X() - xc);
        if (i > 0) { // unwrap
          while (ph[i] - ph[i - 1] > M_PI)
            ph[i] -= 2 * M_PI;
          while (ph[i] - ph[i - 1] < -M_PI)
            ph[i] += 2 * M_PI;
        }
      }
      m_h = ph.back() > ph.front() ? 1.0 : -1.0; // rotation sense of the outward motion
      std::vector<double> u(n);
      for (int i = 0; i < n; ++i)
        u[i] = m_h * (ph[i] - ph[0]);
      const auto [umin, umax] = std::minmax_element(u.begin(), u.end());
      double k = 0;
      if (*umax - *umin > 0) { // slope of z(u)
        const double ub = std::accumulate(u.begin(), u.end(), 0.0) / n;
        double zb = 0;
        for (const auto& h : hits)
          zb += h.Z() / n;
        double num = 0, den = 0;
        for (int i = 0; i < n; ++i) {
          num += (u[i] - ub) * (hits[i].Z() - zb);
          den += (u[i] - ub) * (u[i] - ub);
        }
        k = num / den;
      }
      double z0 = 0;
      for (int i = 0; i < n; ++i)
        z0 += (hits[i].Z() - k * u[i]) / n;
      m_xc = xc;
      m_yc = yc;
      m_R = R;
      m_k = k;
      m_phi0 = ph[0];
      m_z0 = z0;
      m_ds = std::sqrt(R * R + k * k);
      m_line = false;
    }
  }
  if (m_line) { // principal axis of the hits, oriented outward
    TVector3 c(0, 0, 0);
    for (const auto& h : hits)
      c += h;
    c *= 1.0 / n;
    TMatrixDSym scatter(3);
    for (const auto& h : hits) {
      const TVector3 d = h - c;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          scatter(i, j) += d[i] * d[j];
    }
    const TMatrixD vectors = TMatrixDSymEigen(scatter).GetEigenVectors();
    m_dir = TVector3(vectors(0, 0), vectors(1, 0), vectors(2, 0)).Unit();
    if ((hits.back() - hits.front()).Dot(m_dir) < 0)
      m_dir = -m_dir;
    m_base = c + (m_h0 - c).Dot(m_dir) * m_dir;
    m_R = std::numeric_limits<double>::infinity();
    m_k = 0;
    m_h = 0;
    m_pMeas = kNaN;
  } else {
    const double pt = 0.3 * std::fabs(m_fitField) * m_R * 1e-3;
    m_pMeas = pt * std::sqrt(1.0 + (m_k / m_R) * (m_k / m_R));
  }
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
