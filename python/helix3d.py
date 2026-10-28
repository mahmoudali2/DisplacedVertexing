"""Geometry-free 3D pair vertexer for stand-alone muon-system tracks (barrel and endcap alike).

Per track: helix fitted to its own hits in the field where the chambers sit (Bz = -1.7 T everywhere outside the coil in IDEA_o1_v05,
also in the forward volume |z| > 2.62 m):  circle in xy (algebraic fit + geometric refinement), z linear in the turning angle.
The track is parametrised by its 3D path length t from the innermost hit (t < 0 = towards the decay), so a track crossing endcap discs
(small xy motion, large z motion) is described as well as a barrel track.  Tracks without a measurable xy curvature (xy span < XY_MIN
or radius > R_MAX) fall back to a 3D straight line through the hits.
Pair vertex: point of closest approach of the two extrapolated trajectories in 3D (coarse grid in (t1, t2) + local minimisation),
vertex = midpoint, DCA = distance between the two trajectory points.  No region assumption: the caller decides barrel / endcap / outside.
Momentum: pT = 0.3 |B| R with |B| = B_CH (chambers), p = pT sqrt(1 + tan^2 lambda), plus the ionisation loss between vertex and first hit
(DEDX per mm in the dual-readout calorimeter shell and in iron, from the material map below), calibrated on truth (see test/).  Units: mm, GeV, T."""
import numpy as np
from scipy.optimize import minimize

B_CH = 1.7           # |Bz| in the muon chambers and in the calorimeter / forward volume (sim field map)
XY_MIN = 30.0        # mm, minimum xy span for a circle fit
R_MAX = 2.0e6        # mm, larger radius = straight
MMU = 0.1056584
# material map (IDEA_o1_v05): DR calorimeter = box shell between (r<2475, |z|<2475) and (r<4525, |z|<4525), beam hole r<252; iron outside it
CALO_IN, CALO_OUT, CALO_HOLE = 2475.0, 4525.0, 252.0
# GeV/mm.  Calibrated on truth (test/val_*.npz, 3 points, 27k matched muons): (p_true - p_meas) = 1.36 x the loss of the nominal values
# (1.0 calo, 1.45 iron) for all three points (1.355, 1.396, 1.354); with it p_reco/p_true = 1.00 (median), 0.83-1.24 (p16-p84), within +-6% for 3-60 GeV.
DEDX = {"calo": 1.36e-3, "iron": 1.97e-3, "none": 0.0}
# sim field map (IDEA_o1_v05 compact): Bz = +2 T for r_xy < 2300 and |z| < 1.2 x 2180 mm, -1.7 T everywhere else (r_xy < 10 m)
B_IN, R_COIL, Z_FIELD = 2.0, 2300.0, 2616.0


def bz(p):
    return B_IN if (np.hypot(p[0], p[1]) < R_COIL and abs(p[2]) < Z_FIELD) else -B_CH


def dedx_at(p):
    calo, iron = region(np.asarray(p)[None, :])
    return DEDX["calo"] * calo[0] + DEDX["iron"] * iron[0]


def region(p):
    r, az = np.hypot(p[..., 0], p[..., 1]), np.abs(p[..., 2])
    in_out = (r < CALO_OUT) & (az < CALO_OUT)
    in_in = (r < CALO_IN) & (az < CALO_IN)
    calo = in_out & ~in_in & (r > CALO_HOLE)
    iron = ~in_out
    return calo, iron


def eloss(a, b, n=60):
    """ionisation loss [GeV] on the straight segment a -> b"""
    f = (np.arange(n) + 0.5) / n
    pts = a[None, :] + f[:, None] * (b - a)[None, :]
    calo, iron = region(pts)
    dl = np.linalg.norm(b - a) / n
    return float(dl * (DEDX["calo"] * calo.sum() + DEDX["iron"] * iron.sum()))


class HelixTrack:
    def __init__(self, hits):
        self.ok = False
        hits = np.asarray(hits, float)
        hits = hits[np.argsort(np.linalg.norm(hits, axis=1))]
        self.hits, self.n = hits, len(hits)
        if self.n < 3:
            return
        self.h0 = hits[0]
        x, y, z = hits[:, 0], hits[:, 1], hits[:, 2]
        span = np.hypot(x - x[0], y - y[0]).max()
        self.line = True
        if span >= XY_MIN:
            sol, *_ = np.linalg.lstsq(np.c_[2 * x, 2 * y, np.ones(self.n)], x ** 2 + y ** 2, rcond=None)
            xc, yc = sol[0], sol[1]; R = np.sqrt(max(sol[2] + xc * xc + yc * yc, 1e-9))
            for _ in range(8):                                     # geometric refinement (Gauss-Newton)
                dx, dy = x - xc, y - yc; d = np.hypot(dx, dy)
                J = np.c_[-dx / d, -dy / d, -np.ones(self.n)]; res = d - R
                try:
                    st = np.linalg.lstsq(J, -res, rcond=None)[0]
                except np.linalg.LinAlgError:
                    break
                xc, yc, R = xc + st[0], yc + st[1], R + st[2]
                if np.abs(st).max() < 1e-3:
                    break
            if np.isfinite(R) and 0 < R < R_MAX:
                ph = np.unwrap(np.arctan2(y - yc, x - xc))
                self.h = 1.0 if ph[-1] > ph[0] else -1.0           # rotation sense of the outward motion
                u = self.h * (ph - ph[0])                         # turning angle along the motion, 0 at the innermost hit
                k = np.polyfit(u, z, 1)[0] if np.ptp(u) > 0 else 0.0
                self.xc, self.yc, self.R, self.k, self.ph0 = xc, yc, R, k, ph[0]
                self.z0 = float(np.mean(z - k * u))
                self.ds = np.sqrt(R * R + k * k)                  # 3D path length per radian
                self.line = False
                self.chi2xy = float(np.sum((np.hypot(x - xc, y - yc) - R) ** 2))
        if self.line:
            c = hits.mean(axis=0); d = np.linalg.svd(hits - c)[2][0]
            if np.dot(hits[-1] - hits[0], d) < 0:
                d = -d
            self.c, self.d = c, d / np.linalg.norm(d)
            self.R, self.k, self.h = np.inf, 0.0, 0.0
        self.ok = True
        # measured momentum (at the first hit)
        if self.line:
            self.pmeas = np.nan
        else:
            pt = 0.3 * B_CH * self.R * 1e-3
            self.pmeas = pt * np.sqrt(1.0 + (self.k / self.R) ** 2)
        self.tmax = float(np.linalg.norm(self.hits[-1] - self.h0))

    def point(self, t):
        t = np.asarray(t, float)
        if self.line:
            base = self.c + np.dot(self.h0 - self.c, self.d) * self.d       # projection of the innermost hit on the line
            return base + t[..., None] * self.d
        u = t / self.ds; ph = self.ph0 + self.h * u
        return np.stack([self.xc + self.R * np.cos(ph), self.yc + self.R * np.sin(ph), self.z0 + self.k * u], axis=-1)

    def tangent(self, t):
        if self.line:
            return self.d.copy()
        ph = self.ph0 + self.h * t / self.ds
        v = np.array([-self.R * np.sin(ph) * self.h, self.R * np.cos(ph) * self.h, self.k])
        return v / np.linalg.norm(v)


def trajectory(tr, step=20.0, nback=1.3, tol=100.0):
    """Polyline of the track propagated INWARD from the innermost hit through the sim field map, gaining the ionisation loss of the material map
    (curvature follows the local field and momentum); plus the fitted helix up to `tol` mm outward.  Returns t (path length, <= tol), points, p, outward tangents."""
    tf = np.arange(tol, 0.0, -step)[::-1]
    fw_pts = tr.point(tf) if len(tf) else np.zeros((0, 3))
    fw_dir = np.array([tr.tangent(t) for t in tf]).reshape(-1, 3)
    n = int(nback * np.linalg.norm(tr.h0) / step)
    pos = tr.point(0.0).copy(); d = tr.tangent(0.0); cl = np.hypot(d[0], d[1]); tl = d[2] / max(cl, 1e-12); phi = np.arctan2(d[1], d[0])
    p = tr.pmeas if not tr.line else np.nan
    P, D, PP = [pos.copy()], [d.copy()], [p]
    for _ in range(n):
        if not tr.line:
            pt = p * cl
            kap = -tr.h * bz(pos) * 0.3e-3 / max(pt, 1e-3)        # signed xy curvature of the outward motion [1/mm]
            pos = pos - step * np.array([cl * np.cos(phi), cl * np.sin(phi), cl * tl])
            phi = phi - kap * step * cl
            p = p + dedx_at(pos) * step
            d = np.array([cl * np.cos(phi), cl * np.sin(phi), cl * tl])
        else:
            pos = pos - step * d
        P.append(pos.copy()); D.append(d.copy()); PP.append(p)
    P, D, PP = np.array(P[::-1]), np.array(D[::-1]), np.array(PP[::-1])
    t = np.r_[-step * np.arange(n, -1, -1), tf]
    return t, np.r_[P, fw_pts], np.r_[PP, np.full(len(tf), tr.pmeas)], np.r_[D, fw_dir]


def _seg_poca(p0, p1, q0, q1):
    u, v, w = p1 - p0, q1 - q0, p0 - q0
    a, b, c, d, e = u @ u, u @ v, v @ v, u @ w, v @ w
    den = a * c - b * b
    t = np.clip((a * e - b * d) / den, 0.0, 1.0) if den > 1e-12 * a * c else (np.clip(e / c, 0.0, 1.0) if c > 0 else 0.0)
    s = np.clip((b * t - d) / a, 0.0, 1.0) if a > 0 else 0.0
    t = np.clip((b * s + e) / c, 0.0, 1.0) if c > 0 else 0.0
    return s, t, float(np.linalg.norm(p0 + s * u - (q0 + t * v)))


def pair_vertex_stepped(a, b, tol=100.0, step=20.0):
    """3D POCA of the two stepped trajectories; vertex = midpoint; momenta (with the energy gained back to the vertex) and directions at the vertex."""
    if not (a.ok and b.ok):
        return None
    ta, Pa, pa, Da = trajectory(a, step=step, tol=tol); tb, Pb, pb, Db = trajectory(b, step=step, tol=tol)
    Dm = np.linalg.norm(Pa[:, None, :] - Pb[None, :, :], axis=2)
    i, j = np.unravel_index(np.argmin(Dm), Dm.shape)
    best = (np.inf, 0, 0, 0.0, 0.0)
    for ii in (i - 1, i):
        for jj in (j - 1, j):
            if 0 <= ii < len(Pa) - 1 and 0 <= jj < len(Pb) - 1:
                s, t, dd = _seg_poca(Pa[ii], Pa[ii + 1], Pb[jj], Pb[jj + 1])
                if dd < best[0]: best = (dd, ii, jj, s, t)
    dd, ii, jj, s, t = best
    if not np.isfinite(dd):
        dd, ii, jj, s, t = float(Dm[i, j]), min(i, len(Pa) - 2), min(j, len(Pb) - 2), 0.0, 0.0
    xa = Pa[ii] + s * (Pa[ii + 1] - Pa[ii]); xb = Pb[jj] + t * (Pb[jj + 1] - Pb[jj])
    t1 = ta[ii] + s * (ta[ii + 1] - ta[ii]); t2 = tb[jj] + t * (tb[jj + 1] - tb[jj])
    v = 0.5 * (xa + xb)
    between = all(np.linalg.norm(v) < np.linalg.norm(tr.h0) and np.dot(v, tr.h0) > 0 for tr in (a, b))
    out = dict(v=v, dca=float(dd), t=(float(t1), float(t2)), upstream=bool(t1 < tol and t2 < tol), between=bool(between), rxy=float(np.hypot(v[0], v[1])), vz=float(v[2]))
    da, db = Da[ii], Db[jj]
    if a.line or b.line:
        out.update(mass=np.nan, cos_open=float(da @ db), p=(np.nan, np.nan)); return out
    pv = [pa[ii] * da, pb[jj] * db]
    E = sum(np.sqrt(q @ q + MMU ** 2) for q in pv); P = pv[0] + pv[1]
    out.update(mass=float(np.sqrt(max(E * E - P @ P, 0.0))), cos_open=float(pv[0] @ pv[1] / (np.linalg.norm(pv[0]) * np.linalg.norm(pv[1]))), p=(float(pa[ii]), float(pb[jj])))
    return out


def pair_vertex(a, b, tol=100.0, nback=1.3, ngrid=90):
    """3D POCA of two HelixTracks extrapolated inward.  Returns dict or None."""
    if not (a.ok and b.ok):
        return None
    lo_a, lo_b = -nback * np.linalg.norm(a.h0), -nback * np.linalg.norm(b.h0)
    ta, tb = np.linspace(lo_a, tol, ngrid), np.linspace(lo_b, tol, ngrid)
    Pa, Pb = a.point(ta), b.point(tb)
    D = np.linalg.norm(Pa[:, None, :] - Pb[None, :, :], axis=2)
    i, j = np.unravel_index(np.argmin(D), D.shape)
    f = lambda q: float(np.linalg.norm(a.point(q[0]) - b.point(q[1])))
    st = (ta[1] - ta[0], tb[1] - tb[0])
    o = minimize(f, x0=[ta[i], tb[j]], method="Nelder-Mead", options=dict(xatol=0.5, fatol=0.05, maxiter=400, initial_simplex=[[ta[i], tb[j]], [ta[i] + st[0], tb[j]], [ta[i], tb[j] + st[1]]]))
    t1, t2 = o.x
    p1, p2 = a.point(t1), b.point(t2)
    v = 0.5 * (p1 + p2)
    upstream = t1 < tol and t2 < tol
    between = all(np.linalg.norm(v) < np.linalg.norm(t.h0) and np.dot(v, t.h0) > 0 for t in (a, b))
    out = dict(v=v, dca=float(np.linalg.norm(p1 - p2)), t=(float(t1), float(t2)), upstream=bool(upstream), between=bool(between),
               rxy=float(np.hypot(v[0], v[1])), vz=float(v[2]))
    # kinematics at the vertex
    if a.line or b.line:
        out.update(mass=np.nan, cos_open=float(np.dot(a.tangent(t1), b.tangent(t2))), p=(a.pmeas, b.pmeas))
        return out
    pv = []
    for tr, tt in ((a, t1), (b, t2)):
        pm = tr.pmeas + eloss(v, tr.h0)
        pv.append(pm * tr.tangent(tt))
    E = sum(np.sqrt(q @ q + MMU ** 2) for q in pv); P = pv[0] + pv[1]
    out.update(mass=float(np.sqrt(max(E * E - P @ P, 0.0))), cos_open=float(pv[0] @ pv[1] / (np.linalg.norm(pv[0]) * np.linalg.norm(pv[1]))),
               p=(float(np.linalg.norm(pv[0])), float(np.linalg.norm(pv[1]))))
    return out
