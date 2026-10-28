/*
 * Two-track vertexing with local arc fits, see ArcPairVertexing.h.
 * Author: Mahmoud Althakeel
 */
#include "ArcPairVertexing.h"

#include <TMatrixDSymEigen.h>
#include <TVectorD.h>

#include <algorithm>
#include <array>
#include <cmath>

namespace k4vertexing::arcpair {

namespace {
constexpr double kHighland = 13.6e-3; // GeV

double logTerm(double x0) { return x0 > 0 ? std::pow(1.0 + 0.038 * std::log(std::max(x0, 1e-3)), 2) : 1.0; }

/// generalised least squares: parameters and their covariance
bool gls(const TMatrixD& design, const TMatrixDSym& cov, const TVectorD& y, TVectorD& par, TMatrixDSym& parCov) {
  TMatrixDSym weight(cov);
  double det = 0;
  weight.Invert(&det);
  if (det == 0)
    return false;
  TMatrixDSym normal(weight);
  normal.SimilarityT(design); // A^T W A
  normal.Invert(&det);
  if (det == 0)
    return false;
  parCov.ResizeTo(normal);
  parCov = normal;
  TMatrixD designT(TMatrixD::kTransposed, design);
  par.ResizeTo(design.GetNcols());
  par = TMatrixD(parCov, TMatrixD::kMult, TMatrixD(designT, TMatrixD::kMult, weight)) * y;
  return true;
}
} // namespace

double x0PerMM(const TVector3& p, const Config& cfg) {
  const double r = p.Perp(), az = std::fabs(p.Z());
  const bool insideShellOuter = r < cfg.shellOuterR && az < cfg.shellOuterZ;
  const bool insideShellInner = r < cfg.shellInnerR && az < cfg.shellInnerZ;
  // in the shell: its X0 spread over its radial depth
  if (insideShellOuter && !insideShellInner)
    return cfg.shellX0 > 0 ? cfg.shellX0 / (cfg.shellOuterR - cfg.shellInnerR) : 0.0;
  // beyond the shell: uniform density
  if (!insideShellOuter)
    return cfg.outsideX0 > 0 ? cfg.outsideX0 / cfg.outsideThickness : 0.0;
  // inside the inner box: no material
  return 0.0;
}

bool ArcTrack::build(std::vector<TVector3> hits, double omega, const Config& cfg) {
  if (hits.size() < cfg.minHits)
    return false;
  std::sort(hits.begin(), hits.end(), [](const TVector3& a, const TVector3& b) { return a.Mag() < b.Mag(); });
  const int n = static_cast<int>(hits.size());
  m_nHits = n;
  m_h0 = hits.front();

  // preliminary straight line: principal axis of the hits, oriented outward
  TVector3 centroid(0, 0, 0);
  for (const auto& h : hits)
    centroid += h;
  centroid *= 1.0 / n;
  TMatrixDSym scatter(3);
  for (const auto& h : hits) {
    const TVector3 d = h - centroid;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        scatter(i, j) += d[i] * d[j];
  }
  const TMatrixD vectors = TMatrixDSymEigen(scatter).GetEigenVectors();
  TVector3 axis(vectors(0, 0), vectors(1, 0), vectors(2, 0));
  if ((hits.back() - hits.front()).Dot(axis) < 0)
    axis = -axis;
  const double axisPerp = std::hypot(axis.X(), axis.Y());
  if (axisPerp < 1e-6)
    return false;
  m_t0[0] = axis.X() / axisPerp;
  m_t0[1] = axis.Y() / axisPerp;
  m_n0[0] = -m_t0[1];
  m_n0[1] = m_t0[0];
  const double tanLambda = axis.Z() / axisPerp;
  m_cosLambda = 1.0 / std::sqrt(1.0 + tanLambda * tanLambda);

  // momentum for the scattering model; E[1/p^2] for a Gaussian resolution
  const double pt = omega != 0 ? 0.3 * cfg.bFieldForMomentum * 1e-3 / std::fabs(omega) : cfg.defaultPt;
  m_p = std::max(pt / m_cosLambda, cfg.minMomentum);
  m_pt = m_p * m_cosLambda;
  m_msFactor = std::pow(kHighland / m_p, 2) * (1.0 + 3.0 * cfg.momentumResolution * cfg.momentumResolution);

  // local coordinates
  std::vector<double> s(n), u(n), dz(n), l(n);
  for (int k = 0; k < n; ++k) {
    const TVector3 d = hits[k] - m_h0;
    s[k] = d.X() * m_t0[0] + d.Y() * m_t0[1];
    u[k] = d.X() * m_n0[0] + d.Y() * m_n0[1];
    dz[k] = d.Z();
    l[k] = s[k] / m_cosLambda;
  }

  // hit covariance: resolution + correlated scattering between the hits
  const double lLast = std::max(l[n - 1], 1e-9);
  const TVector3 chord = hits.back() - m_h0;
  double meanRho = 0;
  for (int i = 0; i <= 10; ++i)
    meanRho += x0PerMM(m_h0 + chord * (i / 10.0), cfg);
  const double logFactor = logTerm(meanRho / 11.0 * std::max(l[n - 1], 1.0));
  TMatrixDSym vms(n);
  for (int j = 0; j < n; ++j)
    for (int k = j; k < n; ++k) {
      const double m = std::min(l[j], l[k]);
      if (m <= 0)
        continue;
      const double dt = m / cfg.nSteps;
      double sum = 0;
      for (int i = 0; i < cfg.nSteps; ++i) {
        const double t = (i + 0.5) * dt;
        const double rho = x0PerMM(m_h0 + chord * (t / lLast), cfg);
        sum += rho * (l[j] - t) * (l[k] - t);
      }
      vms(j, k) = vms(k, j) = cfg.useMaterial ? m_msFactor * logFactor * sum * dt : 0.0;
    }
  const double hit2 = cfg.hitResolution * cfg.hitResolution;
  TMatrixDSym vu(vms), vz(n);
  for (int j = 0; j < n; ++j)
    for (int k = 0; k < n; ++k)
      vz(j, k) = vms(j, k) / (m_cosLambda * m_cosLambda);
  for (int j = 0; j < n; ++j) {
    vu(j, j) += hit2;
    vz(j, j) += hit2;
  }

  // xy: u(s) = a + b s (+ kappa s^2 / 2);  z: dz(s) = c + d s
  const int nPar = (cfg.fitCurvature && n >= 3) ? 3 : 2;
  TMatrixD designU(n, nPar), designZ(n, 2);
  TVectorD yu(n), yz(n);
  for (int k = 0; k < n; ++k) {
    designU(k, 0) = 1;
    designU(k, 1) = s[k];
    if (nPar == 3)
      designU(k, 2) = 0.5 * s[k] * s[k];
    designZ(k, 0) = 1;
    designZ(k, 1) = s[k];
    yu[k] = u[k];
    yz[k] = dz[k];
  }
  TVectorD pu, pz;
  m_cu.ResizeTo(nPar, nPar);
  if (!gls(designU, vu, yu, pu, m_cu) || !gls(designZ, vz, yz, pz, m_cz))
    return false;
  m_pu.assign(pu.GetMatrixArray(), pu.GetMatrixArray() + nPar);
  m_pz[0] = pz[0];
  m_pz[1] = pz[1];

  TVectorD resid = yu - designU * pu;
  TMatrixDSym wu(vu);
  wu.Invert();
  m_chi2 = wu.Similarity(resid);
  return true;
}

void ArcTrack::uOfS(double s, double& u, double& slope) const {
  const double a = m_pu[0], b = m_pu[1];
  const double k = m_pu.size() > 2 ? m_pu[2] : 0.0;
  if (std::fabs(k) < 1e-9) {
    u = a + b * s;
    slope = b;
    return;
  }
  // exact circle with tangent slope b and curvature kappa at (0, a)
  const double norm = std::sqrt(1.0 + b * b);
  const double radius = 1.0 / k;
  const double xc = -b / norm * radius, yc = a + radius / norm;
  const double arg = radius * radius - (s - xc) * (s - xc);
  if (arg <= 0) {
    u = a + b * s + 0.5 * k * s * s;
    slope = b + k * s;
    return;
  }
  const double root = std::sqrt(arg);
  const double sign = radius > 0 ? 1.0 : -1.0;
  u = yc - sign * root;
  slope = sign * (s - xc) / root;
}

TVector3 ArcTrack::point(double s) const {
  double u, slope;
  uOfS(s, u, slope);
  return TVector3(m_h0.X() + s * m_t0[0] + u * m_n0[0], m_h0.Y() + s * m_t0[1] + u * m_n0[1], zOfS(s));
}

void ArcTrack::local(const TVector3& v, double& s, double& u) const {
  const double dx = v.X() - m_h0.X(), dy = v.Y() - m_h0.Y();
  s = dx * m_t0[0] + dy * m_t0[1];
  u = dx * m_n0[0] + dy * m_n0[1];
}

void ArcTrack::sigmas(double s, const Config& cfg, double& sigmaU, double& sigmaZ) const {
  const int nPar = static_cast<int>(m_pu.size());
  const double jac[3] = {1.0, s, 0.5 * s * s};
  double vu = 0;
  for (int i = 0; i < nPar; ++i)
    for (int j = 0; j < nPar; ++j)
      vu += jac[i] * m_cu(i, j) * jac[j];
  double vz = m_cz(0, 0) + 2 * s * m_cz(0, 1) + s * s * m_cz(1, 1);
  if (s < 0) {
    // scattering between the extrapolation point and the first hit
    const double length = -s / m_cosLambda;
    const TVector3 start = point(s);
    const TVector3 toFirst = m_h0 - start;
    const double dt = length / cfg.nSteps;
    double x0 = 0, sum = 0;
    for (int i = 0; i < cfg.nSteps; ++i) {
      const double t = (i + 0.5) * dt; // distance of the scatterer from the point
      const double rho = x0PerMM(start + toFirst * (t / length), cfg);
      x0 += rho * dt;
      sum += rho * t * t * dt;
    }
    const double ms = cfg.useMaterial ? m_msFactor * logTerm(x0) * sum : 0.0;
    vu += ms;
    vz += ms / (m_cosLambda * m_cosLambda);
    // curvature error from the field mismatch: sagitta 0.5 * dkappa * s^2
    const double sagitta = 0.5 * 0.3e-3 * cfg.fieldUncertainty / m_pt * s * s;
    vu += sagitta * sagitta;
  }
  sigmaU = std::sqrt(std::max(vu, 0.0));
  sigmaZ = std::sqrt(std::max(vz, 0.0));
}

namespace {
struct Residuals {
  std::vector<double> r, sigU, sigZ, s;
};

void computeResiduals(const std::vector<const ArcTrack*>& tracks, const TVector3& v, const Config& cfg,
                      Residuals& out, bool updateSigmas) {
  const std::size_t n = tracks.size();
  out.r.resize(2 * n);
  out.s.resize(n);
  if (updateSigmas) {
    out.sigU.resize(n);
    out.sigZ.resize(n);
  }
  for (std::size_t i = 0; i < n; ++i) {
    double s, uv, u, slope;
    tracks[i]->local(v, s, uv);
    tracks[i]->uOfS(s, u, slope);
    if (updateSigmas)
      tracks[i]->sigmas(s, cfg, out.sigU[i], out.sigZ[i]);
    out.r[2 * i] = (uv - u) / std::sqrt(1.0 + slope * slope) / out.sigU[i];
    out.r[2 * i + 1] = (v.Z() - tracks[i]->zOfS(s)) / out.sigZ[i];
    out.s[i] = s;
  }
}

std::optional<TVector3> seedFromLines(const ArcTrack& a, const ArcTrack& b) {
  auto line = [](const ArcTrack& t, TVector3& c, TVector3& d) {
    c = t.point(0.0);
    d = TVector3(t.dirX(), t.dirY(), t.slopeZ()).Unit();
  };
  TVector3 c1, u1, c2, u2;
  line(a, c1, u1);
  line(b, c2, u2);
  const TVector3 w = c1 - c2;
  const double bb = u1.Dot(u2), den = 1.0 - bb * bb;
  if (den < 1e-8)
    return std::nullopt;
  const double d = u1.Dot(w), e = u2.Dot(w);
  const double s1 = (bb * e - d) / den, s2 = (e - bb * d) / den;
  return 0.5 * (c1 + s1 * u1 + c2 + s2 * u2);
}

/// J^T J and J^T r from a numerical Jacobian at fixed sigmas
void normalEquations(const std::vector<const ArcTrack*>& tracks, const TVector3& v, const Config& cfg,
                     const Residuals& base, TMatrixDSym& hessian, TVectorD& gradient) {
  const double eps = cfg.jacobianStep; // mm
  const std::size_t m = base.r.size();
  std::vector<std::array<double, 3>> jac(m);
  Residuals shifted = base;
  for (int a = 0; a < 3; ++a) {
    TVector3 dv(0, 0, 0);
    dv[a] = eps;
    computeResiduals(tracks, v + dv, cfg, shifted, false);
    for (std::size_t k = 0; k < m; ++k)
      jac[k][a] = (shifted.r[k] - base.r[k]) / eps;
  }
  hessian.ResizeTo(3, 3);
  hessian.Zero();
  gradient.ResizeTo(3);
  gradient.Zero();
  for (std::size_t k = 0; k < m; ++k)
    for (int a = 0; a < 3; ++a) {
      gradient[a] += jac[k][a] * base.r[k];
      for (int b = 0; b < 3; ++b)
        hessian(a, b) += jac[k][a] * jac[k][b];
    }
}
} // namespace

std::optional<VertexResult> fitVertex(const std::vector<const ArcTrack*>& tracks, const Config& cfg) {
  if (tracks.size() < 2)
    return std::nullopt;
  auto seed = seedFromLines(*tracks[0], *tracks[1]);
  if (!seed)
    return std::nullopt;
  TVector3 v = *seed;
  Residuals res;
  TMatrixDSym hessian(3);
  TVectorD gradient(3);
  for (int it = 0; it < cfg.maxIterations; ++it) {
    computeResiduals(tracks, v, cfg, res, true);
    normalEquations(tracks, v, cfg, res, hessian, gradient);
    double det = 0;
    hessian.Invert(&det);
    if (det == 0 || !std::isfinite(det))
      return std::nullopt;
    const TVectorD step = -1.0 * (hessian * gradient);
    v += TVector3(step[0], step[1], step[2]);
    if (!std::isfinite(v.X()))
      return std::nullopt;
    if (std::sqrt(step.Norm2Sqr()) < cfg.convergenceDistance)
      break;
  }
  computeResiduals(tracks, v, cfg, res, true);
  normalEquations(tracks, v, cfg, res, hessian, gradient);
  double det = 0;
  hessian.Invert(&det);
  if (det == 0 || !std::isfinite(det))
    return std::nullopt;

  VertexResult out;
  out.position = v;
  out.covariance = hessian;
  for (double r : res.r)
    out.chi2 += r * r;
  out.ndf = 2 * static_cast<int>(tracks.size()) - 3;
  out.sFromFirst = res.s;
  out.upstream = std::all_of(res.s.begin(), res.s.end(), [&](double s) { return s < cfg.upstreamTolerance; });
  for (std::size_t i = 0; i < tracks.size(); ++i)
    for (std::size_t j = i + 1; j < tracks.size(); ++j)
      out.pairDCA = std::max(out.pairDCA, (tracks[i]->point(res.s[i]) - tracks[j]->point(res.s[j])).Mag());
  // the decay must lie between the IP and the first hits; an arc extrapolated far
  // enough can otherwise cross another track on the far side of the IP
  out.betweenIPAndHits = std::all_of(tracks.begin(), tracks.end(), [&](const ArcTrack* t) {
    return v.Mag() < t->firstHit().Mag() && v.Dot(t->firstHit()) > 0;
  });
  return out;
}

} // namespace k4vertexing::arcpair
