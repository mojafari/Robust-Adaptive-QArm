"""Nominal rigid-body dynamics of the Quanser QArm (4-DOF).

Direct port of the MATLAB function ``f_qarm.m`` used for the original
simulations of the paper.  The equations follow Quanser's ancillary material
("Dynamics - Torque Equations", "Dynamic Parameters").

    M(q) qdd + B(q) nu_cr + C_cf(q) nu_sq + g*Gbar(q) = tau + d_a

with nu_cr = [q1d q2d, q1d q3d, q1d q4d, q2d q3d, q2d q4d, q3d q4d]^T and
nu_sq = [q1d^2, q2d^2, q3d^2, q4d^2]^T.
"""
from __future__ import annotations

import numpy as np

NOMINAL_PARAMS = dict(
    g=9.81,
    m2=0.4591, m3=0.269, m4=0.257, mL=0.1,
    L1=0.14, L2=0.35, L3=0.05, L4=0.25, L5=0.15,
    lamc2=0.1071, lamc3=0.1561, lamc4=0.0998,
    I1A=1.489e-3, I2A=1.922e-4, I2L=1.610e-3, I3A=2.679e-4,
    I3L=2.069e-3, I4A=6.528e-4, I4L=1.120e-3,
)


def derived(p: dict) -> dict:
    p = dict(p)
    p["lam2"] = np.sqrt(p["L2"] ** 2 + p["L3"] ** 2)
    p["lam3"] = p["L4"] + p["L5"]
    return p


def dynamics_terms(q, qd, p):
    """Return M (4x4), h = B nu_cr + C nu_sq + g G (4,) and G-vector (4,)."""
    p = p if "lam2" in p else derived(p)
    g = p["g"]; m2, m3, m4, mL = p["m2"], p["m3"], p["m4"], p["mL"]
    lam2, lam3 = p["lam2"], p["lam3"]
    lamc2, lamc3, lamc4 = p["lamc2"], p["lamc3"], p["lamc4"]
    I1A, I2A, I2L, I3A, I3L, I4A, I4L = (p[k] for k in
                                         ("I1A", "I2A", "I2L", "I3A", "I3L", "I4A", "I4L"))
    s2, s3 = np.sin(q[1]), np.sin(q[2])
    c2, c3 = np.cos(q[1]), np.cos(q[2])
    s23, c23 = np.sin(q[1] + q[2]), np.cos(q[1] + q[2])

    M11 = (I1A + m2 * (lam2 - lamc2) ** 2 * c2 ** 2 + I2A * s2 ** 2 + I2L * c2 ** 2
           + m3 * lam2 ** 2 * c2 ** 2 + m3 * lamc3 ** 2 * s23 ** 2
           - 2 * m3 * lam2 * lamc3 * c2 * s23 + I3L * s23 ** 2 + I3A * c23 ** 2
           + m4 * (lam3 - lamc4) ** 2 * s23 ** 2 + m4 * lam2 ** 2 * c2 ** 2
           - 2 * m4 * (lam3 - lamc4) * lam2 * c2 * s23 + I4L * s23 ** 2 + I4A * c23 ** 2
           + mL * lam3 ** 2 * s23 ** 2 + mL * lam2 ** 2 * c2 ** 2
           - 2 * mL * lam2 * lam3 * c2 * s23)
    M14 = -I4A * c23
    M22 = (m2 * (lam2 - lamc2) ** 2 + I2L + m3 * lam2 ** 2 + m3 * lamc3 ** 2
           - 2 * m3 * lam2 * lamc3 * s3 + I3L + m4 * (lam2 + lam3 - lamc4) ** 2 + I4L
           + mL * lam2 ** 2 + mL * lam3 ** 2 - 2 * mL * lam2 * lam3 * s3)
    M23 = (m3 * lamc3 ** 2 - m3 * lamc3 * lam2 * s3 + I3L
           + m4 * (lam3 - lamc4 - lam2 * s3) * (lam3 - lamc4) + I4L
           - mL * (lam3 ** 2 - lam2 * lam3 * s3))
    M33 = m3 * lamc3 ** 2 + I3L + m4 * (lamc4 - lam3) ** 2 + I4L + mL * lam3 ** 2
    M44 = I4A
    M = np.array([[M11, 0, 0, M14],
                  [0, M22, M23, 0],
                  [0, M23, M33, 0],
                  [M14, 0, 0, M44]])

    B11 = (-2 * m2 * (lam2 - lamc2) ** 2 * s2 * c2 + 2 * I2A * s2 * c2 - 2 * I2L * s2 * c2
           - 2 * m3 * lam2 ** 2 * s2 * c2 + 2 * m3 * lamc3 ** 2 * s23 * c23
           + 2 * m3 * lam2 * lamc3 * s2 * s23 - 2 * m3 * lam2 * lamc3 * c2 * c23
           + 2 * I3L * s23 * c23 - 2 * I3A * s23 * c23
           + 2 * m4 * (lam3 - lamc4) ** 2 * s23 * c23 - 2 * m4 * lam2 ** 2 * s2 * c2
           + 2 * m4 * (lam3 - lamc4) * lam2 * s2 * s23
           - 2 * m4 * (lam3 - lamc4) * lam2 * c2 * c23 + 2 * I4L * s23 * c23
           - 2 * I4A * s23 * c23 + 2 * mL * lam3 ** 2 * s23 * c23
           - 2 * mL * lam2 ** 2 * s2 * c2 + 2 * mL * lam2 * lam3 * s2 * s23
           - 2 * mL * lam2 * lam3 * c2 * c23)
    B12 = (2 * m3 * lamc3 ** 2 * s23 * c23 - 2 * m3 * lam2 * lamc3 * c2 * c23
           + 2 * I3L * s23 * c23 - 2 * I3A * s23 * c23
           + 2 * m4 * (lam3 - lamc4) ** 2 * s23 * c23
           - 2 * m4 * (lam3 - lamc4) * lam2 * c2 * c23 + 2 * I4L * s23 * c23
           - 2 * I4A * s23 * c23 + 2 * mL * lam3 ** 2 * s23 * c23
           - 2 * mL * lam2 * lam3 * c2 * c23)
    B15 = I4A * s23
    B16 = I4A * s23
    B23 = -I4A * s23
    B24 = -(2 * m3 * lam2 * lamc3 * c3 + 2 * mL * lam2 * lam3 * c3)
    B33 = -I4A * s23
    B41 = I4A * s23
    B42 = I4A * s23
    B = np.array([[B11, B12, 0, 0, B15, B16],
                  [0, 0, B23, B24, 0, 0],
                  [0, 0, B33, 0, 0, 0],
                  [B41, B42, 0, 0, 0, 0]])

    C = np.zeros((4, 4))
    C[1, 0] = -(-m2 * (lam2 - lamc2) ** 2 * s2 * c2 + I2A * s2 * c2 - I2L * s2 * c2
                - m3 * lam2 ** 2 * s2 * c2 + m3 * lamc3 ** 2 * s23 * c23
                + m3 * lam2 * lamc3 * s2 * s23 - m3 * lam2 * lamc3 * c2 * c23
                + I3L * s23 * c23 - I3A * s23 * c23
                + m4 * (lam3 - lamc4) ** 2 * s23 * c23 - m4 * lam2 ** 2 * s2 * c2
                + m4 * (lam3 - lamc4) * lam2 * s2 * s23
                - m4 * (lam3 - lamc4) * lam2 * c2 * c23 + I4L * s23 * c23
                - I4A * s23 * c23 + mL * lam3 ** 2 * s23 * c23)
    C[2, 0] = -(m3 * lamc3 ** 2 * s23 * c23 - m3 * lam2 * lamc3 * c2 * c23
                + I3L * s23 * c23 - I3A * s23 * c23
                + m4 * (lam3 - lamc4) ** 2 * s23 * c23
                - m4 * (lam3 - lamc4) * lam2 * c2 * c23 + I4L * s23 * c23
                - I4A * s23 * c23 + mL * lam3 ** 2 * s23 * c23
                - mL * lam2 * lam3 * c2 * c23)
    C[2, 1] = m3 * lam2 * lamc3 * c3 + mL * lam2 * lam3 * c3

    G = np.array([
        0.0,
        -(m2 * (lam2 - lamc2) * c2 + m3 * (lam2 * c2 - lamc3 * s23)
          + m4 * (lam2 * c2 - (lam3 - lamc4) * s23) + mL * (lam2 * c2 - lam3 * s23)),
        m3 * lamc3 * s23 + m4 * (lam3 - lamc4) * s23 + mL * lam3 * s23,
        0.0,
    ])

    nu_cr = np.array([qd[0] * qd[1], qd[0] * qd[2], qd[0] * qd[3],
                      qd[1] * qd[2], qd[1] * qd[3], qd[2] * qd[3]])
    nu_sq = qd ** 2
    h = B @ nu_cr + C @ nu_sq + g * G
    return M, h, g * G


def perturbed_params(rng: np.random.Generator, rel: float = 0.2,
                     payload_range=(0.0, 0.3), base=None) -> dict:
    """Randomly perturb masses, inertias and centre-of-mass locations.

    Each mass, inertia and CoM distance is multiplied by an independent factor
    drawn uniformly from [1-rel, 1+rel]; the payload mass m_L is drawn
    uniformly from ``payload_range`` (kg).  Link lengths are kept nominal.
    """
    p = dict(NOMINAL_PARAMS if base is None else base)
    for k in ("m2", "m3", "m4", "lamc2", "lamc3", "lamc4",
              "I1A", "I2A", "I2L", "I3A", "I3L", "I4A", "I4L"):
        p[k] = p[k] * rng.uniform(1 - rel, 1 + rel)
    p["mL"] = rng.uniform(*payload_range)
    return p
