#!/usr/bin/env python3
"""
Python reference of the arc vertexer (mirrors the Gaudi algorithm
ArcPairVertexFinder). Units: mm, GeV, Tesla.

Per track
  - hits ordered from the IP outward; local frame at the innermost hit h0:
      t0 = outward direction in xy (from a preliminary straight-line fit),
      n0 = z x t0, s = local coordinate along t0, u = coordinate along n0
  - xy: circle through the hits, fitted as u(s) = a + b s + kappa s^2/2
        (generalised least squares, full hit covariance), then evaluated as the
        exact circle with that tangent and curvature
  - z : straight line z(s) = c + d s (generalised least squares)
  - hit covariance = hit resolution (diagonal) + correlated multiple scattering
    from the material between the first hit and each hit (muon-system yoke)

Extrapolation error at arc position s < 0 (upstream of the first hit)
  = fit covariance propagated with J = (1, s, s^2/2) / (1, s)
  + multiple scattering in the material between the point and the first hit
    (dual-readout calorimeter), integrated along the path.

Vertex fit (per track pair)
  Gauss-Newton on 2 residuals per track: the xy distance perpendicular to the
  track and the z difference at the same s. chi2 with ndf = 2N - 3.
"""
import numpy as np

HIGHLAND = 13.6e-3  # GeV


class Config:
    bfield_for_p = 2.0          # T, converts the track omega into pT
    momentum_resolution = 0.20  # standalone muon system sigma(p)/p
    min_momentum = 2.0          # GeV, floor for the scattering estimate
    hit_resolution = 1.0        # mm, per coordinate
    # field mismatch: the circle is fitted in the chambers' field and extrapolated
    # through the calorimeter region, whose field differs. Calibrated on truth
    # (per-track pulls flat vs extrapolation distance): 0.5 T
    field_uncertainty = 0.5     # T
    use_material = True         # multiple scattering in yoke + calorimeter
    # material: dual-readout calorimeter between solenoid and muon system
    calo_x0 = 90.0
    calo_inner_r, calo_outer_r = 2400.0, 4500.0
    calo_inner_z, calo_outer_z = 2400.0, 4300.0
    # muon-system yoke between stations
    yoke_x0 = 40.0
    yoke_thickness = 2000.0
    fit_curvature = True
    min_hits = 3
    # region: only vertices behind the solenoid are kept
    solenoid_r, solenoid_z = 2300.0, 2180.0
    upstream_tolerance = 100.0  # mm, vertex may sit this far past the first hit
    n_steps = 40


def _calo_density(cfg):
    return cfg.calo_x0 / (cfg.calo_outer_r - cfg.calo_inner_r)


def x0_per_mm(p, cfg):
    """Radiation lengths per mm of path at point p."""
    r, az = np.hypot(p[0], p[1]), abs(p[2])
    inside_calo_outer = r < cfg.calo_outer_r and az < cfg.calo_outer_z
    inside_calo_inner = r < cfg.calo_inner_r and az < cfg.calo_inner_z
    if inside_calo_outer and not inside_calo_inner:
        return _calo_density(cfg)
    if not inside_calo_outer:
        return cfg.yoke_x0 / cfg.yoke_thickness
    return 0.0


class ArcTrack:
    def __init__(self, hits, omega, cfg):
        self.ok = False
        hits = hits[np.argsort(np.linalg.norm(hits, axis=1))]
        self.hits = hits
        if len(hits) < cfg.min_hits:
            return
        h0 = hits[0]
        # preliminary 3D line -> outward xy direction and tanLambda
        c = hits.mean(axis=0)
        u = np.linalg.svd(hits - c)[2][0]
        if np.dot(hits[-1] - hits[0], u) < 0:
            u = -u
        txy = u[:2] / np.linalg.norm(u[:2])
        self.t0 = txy
        self.n0 = np.array([-txy[1], txy[0]])
        d = hits - h0
        s = d[:, :2] @ self.t0
        uu = d[:, :2] @ self.n0
        tanl0 = u[2] / np.linalg.norm(u[:2])
        cosl = 1.0 / np.sqrt(1.0 + tanl0 ** 2)

        pt = 0.3 * cfg.bfield_for_p * 1e-3 / abs(omega) if omega else 10.0
        p = max(pt / cosl, cfg.min_momentum)
        self.p = p
        self.pt = p * cosl
        self.cosl = cosl
        # E[1/p^2] with a Gaussian momentum resolution
        self.msfac = (HIGHLAND / p) ** 2 * (1.0 + 3.0 * cfg.momentum_resolution ** 2)

        # ---- hit covariance with correlated multiple scattering ----
        l = s / cosl  # 3D path length from the first hit
        n = len(hits)
        Vms = np.zeros((n, n))
        xtot = sum(x0_per_mm(h0 + (hits[-1] - h0) * f, cfg) for f in np.linspace(0, 1, 11)) / 11 * max(l[-1], 1.0)
        logf = (1.0 + 0.038 * np.log(max(xtot, 1e-3))) ** 2
        for j in range(n):
            for k in range(j, n):
                m = min(l[j], l[k])
                if m <= 0:
                    continue
                ts = np.linspace(0.0, m, cfg.n_steps + 1)
                tm = 0.5 * (ts[1:] + ts[:-1]); dt = ts[1] - ts[0]
                rho = np.array([x0_per_mm(h0 + (hits[-1] - h0) * (t / max(l[-1], 1e-9)), cfg) for t in tm])
                val = np.sum(self.msfac * logf * rho * (l[j] - tm) * (l[k] - tm)) * dt
                Vms[j, k] = Vms[k, j] = val if cfg.use_material else 0.0
        Vu = Vms + np.eye(n) * cfg.hit_resolution ** 2
        Vz = Vms / cosl ** 2 + np.eye(n) * cfg.hit_resolution ** 2

        # ---- GLS fits ----
        npar = 3 if (cfg.fit_curvature and n >= 3) else 2
        A = np.stack([np.ones(n), s, 0.5 * s ** 2][:npar], axis=1)
        Wu = np.linalg.inv(Vu)
        Cu = np.linalg.inv(A.T @ Wu @ A)
        self.pu = Cu @ A.T @ Wu @ uu
        self.Cu = Cu
        Az = np.stack([np.ones(n), s], axis=1)
        Wz = np.linalg.inv(Vz)
        Cz = np.linalg.inv(Az.T @ Wz @ Az)
        self.pz = Cz @ Az.T @ Wz @ d[:, 2]
        self.Cz = Cz
        self.h0 = h0
        self.chi2 = float((uu - A @ self.pu) @ Wu @ (uu - A @ self.pu))
        self.ok = True

    # exact circle with tangent slope b and curvature kappa at (0, a)
    def u_of_s(self, s):
        a, b = self.pu[0], self.pu[1]
        k = self.pu[2] if len(self.pu) > 2 else 0.0
        if abs(k) < 1e-9:
            return a + b * s, b
        norm = np.sqrt(1.0 + b * b)
        R = 1.0 / k
        xc, yc = -b / norm * R, a + R / norm
        arg = R * R - (s - xc) ** 2
        if arg <= 0:
            return a + b * s + 0.5 * k * s * s, b + k * s
        root = np.sqrt(arg)
        uval = yc - np.sign(R) * root
        slope = np.sign(R) * (s - xc) / root
        return uval, slope

    def point(self, s):
        u, _ = self.u_of_s(s)
        xy = self.h0[:2] + s * self.t0 + u * self.n0
        z = self.h0[2] + self.pz[0] + self.pz[1] * s
        return np.array([xy[0], xy[1], z])

    def local(self, v):
        d = v[:2] - self.h0[:2]
        return float(d @ self.t0), float(d @ self.n0)

    def sigmas(self, s, cfg):
        """(sigma_perp_xy, sigma_z) of the trajectory at local s."""
        J = np.array([1.0, s, 0.5 * s * s])[: len(self.pu)]
        vu = float(J @ self.Cu @ J)
        Jz = np.array([1.0, s])
        vz = float(Jz @ self.Cz @ Jz)
        if s < 0:
            # scattering between the point and the first hit (calorimeter)
            L = -s / self.cosl
            start = self.point(s)
            ts = np.linspace(0.0, L, cfg.n_steps + 1)
            tm = 0.5 * (ts[1:] + ts[:-1]); dt = ts[1] - ts[0]
            rho = np.array([x0_per_mm(start + (self.h0 - start) * (t / L), cfg) for t in tm])
            xtot = float(np.sum(rho) * dt)
            logf = (1.0 + 0.038 * np.log(max(xtot, 1e-3))) ** 2 if xtot > 0 else 1.0
            # tm = distance of the scatterer from the point
            ms = float(np.sum(self.msfac * logf * rho * tm ** 2) * dt) if cfg.use_material else 0.0
            vu += ms
            vz += ms / self.cosl ** 2
            # curvature error from the field mismatch: sagitta 0.5 * dkappa * s^2
            vu += (0.5 * 0.3e-3 * cfg.field_uncertainty / self.pt * s * s) ** 2
        return np.sqrt(vu), np.sqrt(vz)


def residuals(tracks, v, cfg, sig=None):
    res, sg, ss = [], [], []
    for i, t in enumerate(tracks):
        s, uv = t.local(v)
        u, slope = t.u_of_s(s)
        zt = t.h0[2] + t.pz[0] + t.pz[1] * s
        r1 = (uv - u) / np.sqrt(1.0 + slope * slope)
        r2 = v[2] - zt
        if sig is None:
            s1, s2 = t.sigmas(s, cfg)
        else:
            s1, s2 = sig[i]
        res += [r1 / s1, r2 / s2]
        sg.append((s1, s2)); ss.append(s)
    return np.array(res), sg, ss


def line_poca(t1, t2):
    """Seed: POCA of the tangent lines at the first hits."""
    def line(t):
        d = np.array([t.t0[0], t.t0[1], t.pz[1]])
        d = d / np.linalg.norm(d)
        return t.point(0.0), d
    (c1, u1), (c2, u2) = line(t1), line(t2)
    w = c1 - c2; b = u1 @ u2; den = 1 - b * b
    if den < 1e-8:
        return None
    dd, e = u1 @ w, u2 @ w
    s1 = (b * e - dd) / den; s2 = (e - b * dd) / den
    return 0.5 * (c1 + s1 * u1 + c2 + s2 * u2)


def fit_vertex(tracks, cfg, n_iter=12):
    v = line_poca(tracks[0], tracks[1])
    if v is None:
        return None
    eps = 1.0
    for _ in range(n_iter):
        r, sig, _ = residuals(tracks, v, cfg)
        J = np.zeros((len(r), 3))
        for a in range(3):
            dv = np.zeros(3); dv[a] = eps
            J[:, a] = (residuals(tracks, v + dv, cfg, sig)[0] - r) / eps
        H = J.T @ J
        try:
            step = -np.linalg.solve(H, J.T @ r)
        except np.linalg.LinAlgError:
            return None
        v = v + step
        if np.linalg.norm(step) < 0.1:
            break
    r, sig, ss = residuals(tracks, v, cfg)
    J = np.zeros((len(r), 3))
    for a in range(3):
        dv = np.zeros(3); dv[a] = eps
        J[:, a] = (residuals(tracks, v + dv, cfg, sig)[0] - r) / eps
    try:
        cov = np.linalg.inv(J.T @ J)
    except np.linalg.LinAlgError:
        return None
    chi2 = float(r @ r)
    ndf = 2 * len(tracks) - 3
    pts = [t.point(s) for t, s in zip(tracks, ss)]
    dca = float(np.linalg.norm(pts[0] - pts[1]))
    upstream = all(s < cfg.upstream_tolerance for s in ss)
    rxy = np.hypot(v[0], v[1])
    behind = rxy > cfg.solenoid_r or abs(v[2]) > cfg.solenoid_z
    # decay must lie between the IP and the chambers: closer to the origin than
    # every first hit and on the same side (rejects arcs extrapolated past the IP)
    between = all(np.linalg.norm(v) < np.linalg.norm(t.h0) and np.dot(v, t.h0) > 0 for t in tracks)
    return dict(v=v, cov=cov, chi2=chi2, ndf=ndf, dca=dca, upstream=upstream,
                behind=behind, s=ss, between=between)
