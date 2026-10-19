"""Joint-space tracking controllers for the 4-DOF QArm model.

All controllers share the interface

    tau = ctrl.step(t, ref, y, Ts)

where ``ref = (r, rd, rdd)`` are the desired joint position, velocity and
acceleration (arrays of length n), ``y = (q, qd, qdd_est)`` are the measured
joint position, velocity and (filtered) acceleration estimate, and ``Ts`` is
the sampling period used to integrate the internal adaptive states with a
forward-Euler rule (as a fixed-step Simulink/QUARC implementation would).

Controllers
-----------
SLController(variant="accel")   proposed law with c = e2dot + k2d e2 + k1d e1,
                                e2dot obtained from the measured acceleration
SLController(variant="model")   acceleration-free (implementable) variant used
                                on the hardware: e2dot replaced by its nominal
                                prediction -k1 e1 - k2 e2 + dhat
FLController                    feedback linearisation with fixed gains,
                                optionally with the disturbance estimator only
PIDGController                  independent-joint PID + nominal gravity
                                compensation, integrator clamping anti-windup
"""
from __future__ import annotations

import numpy as np

from .model import dynamics_terms


class SLController:
    def __init__(self, p_nom, alpha1, alpha2, alpha_d, k1d, k2d,
                 k10=1.0, k20=1.0, n=4, variant="accel", adapt_gains=True,
                 adapt_dhat=True, active=None, gain_floor=None):
        self.p = p_nom
        self.a1, self.a2, self.ad = alpha1, alpha2, alpha_d
        self.k1d, self.k2d = k1d, k2d
        self.k10, self.k20 = k10, k20
        self.n = n
        self.variant = variant
        self.adapt_gains = adapt_gains
        self.adapt_dhat = adapt_dhat
        # mask of joints whose errors enter the learning signal (hardware: wrist only)
        self.active = np.ones(n) if active is None else np.asarray(active, float)
        self.gain_floor = gain_floor  # optional projection k_i >= floor
        self.reset()

    def reset(self):
        self.k1 = float(self.k10)
        self.k2 = float(self.k20)
        self.dhat = np.zeros(self.n)
        self.last_c = np.zeros(self.n)

    def accel_command(self, ref, y):
        r, rd, rdd = ref
        q, qd, _ = y
        e1 = r - q
        e2 = rd - qd
        return self.k1 * e1 + self.k2 * e2 + rdd - self.dhat, e1, e2

    def step(self, t, ref, y, Ts):
        r, rd, rdd = ref
        q, qd, qdd_est = y
        v, e1, e2 = self.accel_command(ref, y)
        M, h, _ = dynamics_terms(q, qd, self.p)
        tau = M @ v + h
        if self.variant == "accel":
            e2dot = rdd - qdd_est
        else:  # nominal prediction of the error derivative (d assumed zero)
            e2dot = -self.k1 * e1 - self.k2 * e2 + self.dhat
        c = (e2dot + self.k2d * e2 + self.k1d * e1) * self.active
        e1a, e2a = e1 * self.active, e2 * self.active
        self.last_c = c
        if self.adapt_gains:
            self.k1 += Ts * self.a1 * float(c @ e1a)
            self.k2 += Ts * self.a2 * float(c @ e2a)
            if self.gain_floor is not None:
                self.k1 = max(self.k1, self.gain_floor)
                self.k2 = max(self.k2, self.gain_floor)
        if self.adapt_dhat:
            self.dhat = self.dhat - Ts * self.ad * c
        return tau

    def state(self):
        return np.concatenate([[self.k1, self.k2], self.dhat])


class FLController(SLController):
    """Feedback linearisation with fixed gains (optionally + disturbance estimate)."""

    def __init__(self, p_nom, k1, k2, n=4, with_dhat=False, alpha_d=0.0,
                 k1d=1.0, k2d=2.0, variant="accel"):
        super().__init__(p_nom, 0.0, 0.0, alpha_d, k1d, k2d, k10=k1, k20=k2, n=n,
                         variant=variant, adapt_gains=False, adapt_dhat=with_dhat)


class PIDGController:
    """tau = Kp e + Ki int(e) + Kd edot + G_nom(q), per joint, with clamping."""

    def __init__(self, p_nom, Kp, Ki, Kd, n=4, int_limit=None, gravity=True):
        self.p = p_nom
        self.Kp, self.Ki, self.Kd = (np.asarray(x, float) for x in (Kp, Ki, Kd))
        self.n = n
        self.int_limit = int_limit
        self.gravity = gravity
        self.reset()

    def reset(self):
        self.xi = np.zeros(self.n)

    def step(self, t, ref, y, Ts):
        r, rd, rdd = ref
        q, qd, _ = y
        e1 = r - q
        e2 = rd - qd
        tau = self.Kp * e1 + self.Ki * self.xi + self.Kd * e2
        if self.gravity:
            _, _, G = dynamics_terms(q, qd, self.p)
            tau = tau + G
        self.xi = self.xi + Ts * e1
        if self.int_limit is not None:
            self.xi = np.clip(self.xi, -self.int_limit, self.int_limit)
        return tau

    def state(self):
        return self.xi.copy()


def pole_matched_pid(p_nom, alpha_d, k1d, k2d, k10=1.0, k20=1.0, q0=None):
    """PID gains whose nominal decoupled closed loop reproduces the linearised
    characteristic polynomial of the proposed controller at its initial gains:

        s^3 + (k2+alpha_d) s^2 + (k1+alpha_d k2d) s + alpha_d k1d.

    For joint i with effective inertia m_i = M_ii(q0):
        Kd_i = m_i (k2+alpha_d), Kp_i = m_i (k1+alpha_d k2d), Ki_i = m_i alpha_d k1d.
    """
    q0 = np.zeros(4) if q0 is None else q0
    M, _, _ = dynamics_terms(q0, np.zeros(4), p_nom)
    m = np.diag(M)
    return m * (k10 + alpha_d * k2d), m * (alpha_d * k1d), m * (k20 + alpha_d)
