"""Simulation drivers.

``simulate_original`` reproduces, line by line, the continuous-time MATLAB
simulations of the submitted manuscript (controller and plant integrated
jointly with fixed-step RK4, identical nominal model in plant and controller,
learning signal evaluated with the exact error acceleration).

``simulate_sampled`` is a sampled-data simulator used for the additional
studies requested in the review: the controller runs at a fixed sampling
period Ts with zero-order hold, the plant may differ from the nominal model,
and optional measurement noise/quantisation, joint friction and actuator
saturation can be included.
"""
from __future__ import annotations

import numpy as np

from .model import NOMINAL_PARAMS, derived, dynamics_terms


# ---------------------------------------------------------------------------
# Reference trajectories (identical to the paper)
# ---------------------------------------------------------------------------
def ref_constant(t):
    r = np.array([np.pi / 4, np.pi / 2, -np.pi / 2, np.pi / 6])
    return r, np.zeros(4), np.zeros(4)


def ref_sinusoidal(t):
    a = np.deg2rad([45.0, 10.0, 25.0, 60.0])
    r = np.array([a[0] * np.cos(t / 2), a[1] * np.sin(t / 4), a[2] * np.cos(t / 4), a[3] * np.sin(t)])
    rd = np.array([-a[0] / 2 * np.sin(t / 2), a[1] / 4 * np.cos(t / 4), -a[2] / 4 * np.sin(t / 4), a[3] * np.cos(t)])
    rdd = np.array([-a[0] / 4 * np.cos(t / 2), -a[1] / 16 * np.sin(t / 4), -a[2] / 16 * np.cos(t / 4), -a[3] * np.sin(t)])
    return r, rd, rdd


def dist_constant(t, d0=0.5):
    return np.full(4, d0)


def dist_pulses(t, d0=0.5, d1=2.5, times=(20.0, 40.0, 60.0, 80.0), width=1.0):
    d = np.full(4, d0)
    for i, ti in enumerate(times):
        if t >= ti and t <= ti + width:   # same (inclusive) test as the MATLAB code
            d[i] = d1
    return d


# ---------------------------------------------------------------------------
# Original continuous-time simulation (port of sim_qarm.m / f_qarm.m)
# ---------------------------------------------------------------------------
def _f_original(t, s, gains, ref_fun, dist_fun, p):
    a1, a2, ad, k1d, k2d = gains
    q, qd = s[0:4], s[4:8]
    k1, k2 = s[8], s[9]
    dh = s[10:14]
    r, rd, rdd = ref_fun(t)
    d = dist_fun(t)
    M, h, _ = dynamics_terms(q, qd, p)
    e1 = r - q
    e2 = rd - qd
    e2d = -k1 * e1 - k2 * e2 - d + dh           # exact error acceleration
    c = e2d + k2d * e2 + k1d * e1
    dk1 = a1 * c @ e1
    dk2 = a2 * c @ e2
    ddh = -ad * c
    u = M @ (k1 * e1 + k2 * e2 + rdd - dh) + h
    qdd = np.linalg.solve(M, u - h) + d
    return np.concatenate([qd, qdd, [dk1, dk2], ddh]), u


def simulate_original(ref_fun, dist_fun, gains=(10.5, 0.2, 5.0, 1.0, 2.0),
                      T=120.0, h=0.05, k0=(1.0, 1.0), p=None, gain_floor=None):
    """gain_floor: optional projection k1, k2 >= gain_floor applied after
    every integration step (not used in the original paper)."""
    p = derived(NOMINAL_PARAMS if p is None else p)
    s = np.zeros(14)
    s[8], s[9] = k0
    # time vector built by repeated addition exactly as in sim_qarm.m
    # (while t(end) < T: t = [t; t(end)+h]) so that sample selection by
    # t >= t_ss is bit-for-bit identical to the original MATLAB scripts
    tl = [0.0]
    while tl[-1] < T:
        tl.append(tl[-1] + h)
    t = np.array(tl)
    n = len(t) - 1
    S = np.zeros((n + 1, 14))
    U = np.zeros((n + 1, 4))
    S[0] = s
    f = lambda tt, ss: _f_original(tt, ss, gains, ref_fun, dist_fun, p)
    for k in range(n):
        tk = t[k]
        f0, u0 = f(tk, s)
        U[k] = u0
        f1, _ = f(tk + h / 2, s + h * f0 / 2)
        f2, _ = f(tk + h / 2, s + h * f1 / 2)
        f3, _ = f(tk + h, s + h * f2)
        s = s + h * (f0 + 2 * f1 + 2 * f2 + f3) / 6
        if gain_floor is not None:
            s[8] = max(s[8], gain_floor)
            s[9] = max(s[9], gain_floor)
        S[k + 1] = s
    U[n] = f(t[n], s)[1]
    return dict(t=t, q=S[:, 0:4], qd=S[:, 4:8], k1=S[:, 8], k2=S[:, 9],
                dhat=S[:, 10:14], tau=U,
                ref=np.array([ref_fun(tt)[0] for tt in t]))


# ---------------------------------------------------------------------------
# Sampled-data simulation with plant/model mismatch
# ---------------------------------------------------------------------------
class AccelFilter:
    """Critically damped 2nd-order low-pass filter applied to the measured
    velocity; its derivative output is used as the acceleration estimate
    (same structure as the QUARC 'Second-Order Low-Pass Filter' block)."""

    def __init__(self, n, wc, Ts):
        self.wc, self.Ts = wc, Ts
        self.y = np.zeros(n)
        self.yd = np.zeros(n)

    def reset(self, v0):
        self.y = v0.copy()
        self.yd[:] = 0.0

    def update(self, v):
        w = self.wc
        ydd = w * w * (v - self.y) - 2 * w * self.yd
        self.y = self.y + self.Ts * self.yd
        self.yd = self.yd + self.Ts * ydd
        return self.y, self.yd


def friction_torque(qd, Fc, Fv, vs=1e-3):
    """Smooth Coulomb + viscous joint friction."""
    return Fc * np.tanh(qd / vs) + Fv * qd


def simulate_sampled(ctrl, ref_fun, dist_fun, p_true=None, T=60.0, Ts=0.002,
                     nsub=2, q0=None, qd0=None, noise=None, rng=None,
                     tau_max=None, friction=None, accel_wc=2 * np.pi * 20.0,
                     log_every=5):
    """Simulate ``ctrl`` on the (possibly perturbed) QArm model.

    noise    : dict(q_quant=rad, q_std=rad, qd_std=rad/s) or None
    tau_max  : array of symmetric torque limits (N m) or None
    friction : dict(Fc=array, Fv=array) or None
    """
    p_true = derived(NOMINAL_PARAMS if p_true is None else p_true)
    rng = np.random.default_rng(0) if rng is None else rng
    q = np.zeros(4) if q0 is None else np.array(q0, float)
    qd = np.zeros(4) if qd0 is None else np.array(qd0, float)
    ctrl.reset()
    filt = AccelFilter(4, accel_wc, Ts)
    filt.reset(qd)
    n = int(round(T / Ts))
    h = Ts / nsub
    nlog = n // log_every + 1
    log = dict(t=np.zeros(nlog), q=np.zeros((nlog, 4)), ref=np.zeros((nlog, 4)),
               tau=np.zeros((nlog, 4)), k1=np.zeros(nlog), k2=np.zeros(nlog),
               dhat=np.zeros((nlog, 4)), sat=np.zeros((nlog, 4), bool))
    Fc = Fv = None
    if friction is not None:
        Fc, Fv = np.asarray(friction["Fc"], float), np.asarray(friction["Fv"], float)

    def accel(tt, qq, vv, tau):
        M, hh, _ = dynamics_terms(qq, vv, p_true)
        rhs = tau - hh
        if Fc is not None:
            rhs = rhs - friction_torque(vv, Fc, Fv)
        return np.linalg.solve(M, rhs) + dist_fun(tt)

    il = 0
    for k in range(n):
        t = k * Ts
        # --- measurement
        qm, qdm = q.copy(), qd.copy()
        if noise is not None:
            if noise.get("q_std", 0) > 0:
                qm = qm + rng.normal(0, noise["q_std"], 4)
            if noise.get("q_quant", 0) > 0:
                qm = np.round(qm / noise["q_quant"]) * noise["q_quant"]
            if noise.get("qd_std", 0) > 0:
                qdm = qdm + rng.normal(0, noise["qd_std"], 4)
        _, qdd_est = filt.update(qdm)
        ref = ref_fun(t)
        tau_cmd = ctrl.step(t, ref, (qm, qdm, qdd_est), Ts)
        tau = tau_cmd if tau_max is None else np.clip(tau_cmd, -tau_max, tau_max)
        if k % log_every == 0:
            log["t"][il] = t
            log["q"][il] = q
            log["ref"][il] = ref[0]
            log["tau"][il] = tau
            if hasattr(ctrl, "k1"):
                log["k1"][il], log["k2"][il] = ctrl.k1, ctrl.k2
                log["dhat"][il] = ctrl.dhat
            if tau_max is not None:
                log["sat"][il] = np.abs(tau_cmd) > tau_max
            il += 1
        # --- plant integration over one sample (ZOH)
        for j in range(nsub):
            tt = t + j * h
            a0 = accel(tt, q, qd, tau)
            q1, v1 = q + h / 2 * qd, qd + h / 2 * a0
            a1 = accel(tt + h / 2, q1, v1, tau)
            q2, v2 = q + h / 2 * v1, qd + h / 2 * a1
            a2 = accel(tt + h / 2, q2, v2, tau)
            q3, v3 = q + h * v2, qd + h * a2
            a3 = accel(tt + h, q3, v3, tau)
            q = q + h / 6 * (qd + 2 * v1 + 2 * v2 + v3)
            qd = qd + h / 6 * (a0 + 2 * a1 + 2 * a2 + a3)
        if not np.all(np.isfinite(q)) or np.max(np.abs(q)) > 1e3:
            log["diverged"] = True
            for key in ("t", "q", "ref", "tau", "k1", "k2", "dhat", "sat"):
                log[key] = log[key][:il]
            return log
    for key in ("t", "q", "ref", "tau", "k1", "k2", "dhat", "sat"):
        log[key] = log[key][:il]
    log["diverged"] = False
    return log


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
def tracking_metrics(t, ref, q, t_ss=5.0, deg=True):
    """Per-joint and pooled RMSE / NRMSE, overall and steady-state."""
    k = 180 / np.pi if deg else 1.0
    e = (ref - q) * k
    r = ref * k
    m = t >= t_ss
    out = {}
    for tag, mask in (("all", np.ones_like(t, bool)), ("ss", m)):
        ee, rr = e[mask], r[mask]
        rmse = np.sqrt(np.mean(ee ** 2, axis=0))
        rms_ref = np.sqrt(np.mean(rr ** 2, axis=0))
        out[f"rmse_{tag}"] = rmse
        out[f"nrmse_{tag}"] = 100 * rmse / rms_ref
        out[f"rmse_{tag}_agg"] = np.sqrt(np.mean(ee ** 2))
        out[f"nrmse_{tag}_agg"] = 100 * out[f"rmse_{tag}_agg"] / np.sqrt(np.mean(rr ** 2))
        out[f"maxe_{tag}"] = np.max(np.abs(ee), axis=0)
    return out
