"""Hardware-representative single-joint model of the QArm wrist experiments.

This module reproduces, as closely as the available information allows, the
signal chain used in the hardware experiments of the paper:

* Quanser low-level wrist PID (``LowLevelControlWristCurrent.slx``):
  command rate limiter (pi/4 rad/s), error and derivative paths enabled only
  for t > 5 s, two gain sets selected by the sign of the command,
  integrator with a 0.005 rad error dead-band and +-0.5 rad.s clamping,
  integrator reset whenever the (gated) error crosses zero (Simulink
  'either' external reset driven by the error signal),
  derivative on a 5 Hz second-order filtered velocity, output rate limiter
  (+-1 A/s) and asymmetric current saturation [-0.05, 0.015] A.
* Proposed controller as implemented on the hardware (Nonlinear_Control_Test_Program.slx,
  wrist only): acceleration-free learning signal, joints 1-3 fed to the
  controller as zero, Simulink Derivative blocks for the reference
  derivatives, torque-to-current conversion I = tau / K_tau with
  K_tau = 0.005 N m/A x 353.5 = 1.7675 N m/A, output rate limiter +-1.5 A/s
  and current saturation.
* Wrist plant identified from the logged data (pooled least squares over the
  ten experiments): effective inertia, Coulomb/viscous friction with
  stick-slip, constant offset torque and payload gravity torque, plus
  encoder, velocity and current quantisation of the Dynamixel XM430-W350.

All quantities are SI (rad, rad/s, A, N m).
"""
from __future__ import annotations

import numpy as np

I4A = 6.528e-4                 # nominal wrist inertia used by the controller
K_TAU = 1.7675                 # N m / A used for the torque-to-current conversion
CURRENT_GAIN = I4A / K_TAU     # = 3.6934e-4 A per (rad/s^2)

# identified wrist plant (see scripts/identify_wrist.py)
WRIST_PLANT = dict(J=1.97e-3, dJ_payload=1.30e-4, Fc=0.020, Fs=0.024, Fv=6.5e-3,
                   tau_off=-2.1e-3, g_sin=1.85e-2, g_cos=-6.1e-3,
                   q_quant=2 * np.pi / 4096,          # 12-bit encoder
                   qd_quant=0.229 * 2 * np.pi / 60,   # 0.229 rpm velocity unit
                   i_quant=2.69e-3)                  # 2.69 mA current unit

QUANSER_LIMITS = (-0.05, 0.015)


class RateLimiter:
    def __init__(self, rise, fall, Ts, y0=0.0):
        self.rise, self.fall, self.Ts, self.y = rise, fall, Ts, y0
        self.first = True

    def __call__(self, u):
        if self.first:           # Simulink: first output equals the input
            self.first = False
            self.y = u
            return u
        dy = np.clip(u - self.y, self.fall * self.Ts, self.rise * self.Ts)
        self.y = self.y + dy
        return self.y


class QuanserWristPID:
    def __init__(self, Ts=0.002, gate_time=5.0, limits=QUANSER_LIMITS,
                 gains_pos=(0.01, 0.03, 0.007), gains_neg=(0.007, 0.09, 0.012),
                 deadband=0.005, int_limit=0.5, out_rate=1.0, lpf_wn=2 * np.pi * 5,
                 reset_on_crossing=True):
        self.Ts, self.gate_time, self.limits = Ts, gate_time, limits
        self.reset_on_crossing = reset_on_crossing
        self.gp, self.gn = gains_pos, gains_neg
        self.deadband, self.int_limit = deadband, int_limit
        self.out_rate, self.wn = out_rate, lpf_wn
        self.reset()

    def reset(self):
        self.xi = 0.0
        self.e_prev = 0.0
        self.vf, self.vfd = 0.0, 0.0
        self.cmd_rl = RateLimiter(np.pi / 4, -np.pi / 4, self.Ts)
        self.out_rl = RateLimiter(self.out_rate, -self.out_rate, self.Ts)

    def step(self, t, r, q, qd):
        Ts = self.Ts
        rr = self.cmd_rl(r)
        gate = 1.0 if t > self.gate_time else 0.0
        e = (rr - q) * gate
        # 2nd-order LPF on the measured speed (zeta = 1)
        ydd = self.wn ** 2 * (qd - self.vf) - 2 * self.wn * self.vfd
        vf = self.vf
        self.vf += Ts * self.vfd
        self.vfd += Ts * ydd
        ed = (0.0 - vf) * gate
        Kp, Ki, Kd = self.gp if np.sign(r) > 0 else self.gn
        ep = self.e_prev
        crossed = (ep < 0 <= e) or (ep == 0 < e) or (ep > 0 >= e) or (ep == 0 > e)
        if crossed and self.reset_on_crossing:
            self.xi = 0.0
        self.e_prev = e
        u = Kp * e + Ki * self.xi + Kd * ed
        if abs(e) >= self.deadband:
            self.xi = float(np.clip(self.xi + Ts * e, -self.int_limit, self.int_limit))
        u = self.out_rl(u)
        return float(np.clip(u, *self.limits))


class ProposedWristHW:
    """Acceleration-free learning controller, wrist channel only."""

    def __init__(self, alpha1, alpha2, alpha_d, k1d, k2d, Ts=0.002,
                 limits=QUANSER_LIMITS, k10=1.0, k20=1.0, out_rate=1.5,
                 variant="model", accel_wn=2 * np.pi * 5):
        self.a1, self.a2, self.ad = alpha1, alpha2, alpha_d
        self.k1d, self.k2d = k1d, k2d
        self.Ts, self.limits, self.k10, self.k20 = Ts, limits, k10, k20
        self.out_rate = out_rate
        self.variant = variant
        self.wn = accel_wn
        self.reset()

    def reset(self):
        self.k1, self.k2, self.dh = self.k10, self.k20, 0.0
        self.r_prev, self.r2_prev, self.first = 0.0, 0.0, True
        self.vf, self.vfd = 0.0, 0.0
        self.out_rl = RateLimiter(self.out_rate, -self.out_rate, self.Ts)

    def step(self, t, r, q, qd):
        Ts = self.Ts
        if self.first:
            r2, r2d = 0.0, 0.0
            self.first = False
        else:
            r2 = (r - self.r_prev) / Ts
            r2d = (r2 - self.r2_prev) / Ts
        self.r_prev, self.r2_prev = r, r2
        e1, e2 = r - q, r2 - qd
        if self.variant == "model":
            c = (self.k1d - self.k1) * e1 + (self.k2d - self.k2) * e2 + self.dh
        else:  # filtered-acceleration version of the ideal law
            ydd = self.wn ** 2 * (qd - self.vf) - 2 * self.wn * self.vfd
            self.vf += Ts * self.vfd
            self.vfd += Ts * ydd
            c = (r2d - self.vfd) + self.k2d * e2 + self.k1d * e1
        v = self.k1 * e1 + self.k2 * e2 + r2d - self.dh
        i_cmd = CURRENT_GAIN * v
        self.k1 += Ts * self.a1 * c * e1
        self.k2 += Ts * self.a2 * c * e2
        self.dh -= Ts * self.ad * c
        i_cmd = self.out_rl(i_cmd)
        return float(np.clip(i_cmd, *self.limits))


def simulate_wrist(ctrl, amp_deg, freq_hz, payload=False, T=120.0, Ts=0.002,
                   nsub=10, plant=None, q0=0.0, quantize=True, friction=True,
                   gravity=True, true_inertia=True):
    """Simulate one wrist experiment. Returns dict of logged signals."""
    P = dict(WRIST_PLANT if plant is None else plant)
    J = (P["J"] if true_inertia else I4A) + (P["dJ_payload"] if payload else 0.0)
    Fc = P["Fc"] if friction else 0.0
    Fs = P["Fs"] if friction else 0.0
    Fv = P["Fv"] if friction else 0.0
    ctrl.reset()
    n = int(round(T / Ts))
    h = Ts / nsub
    A = np.deg2rad(amp_deg)
    w = 2 * np.pi * freq_hz
    q, v = q0, 0.0
    tl = np.arange(n + 1) * Ts
    Q = np.zeros(n + 1); R = np.zeros(n + 1); I = np.zeros(n + 1); V = np.zeros(n + 1)
    K1 = np.full(n + 1, np.nan); K2 = np.full(n + 1, np.nan); DH = np.full(n + 1, np.nan)
    for k in range(n + 1):
        t = tl[k]
        r = A * np.sin(w * t)
        qm, vm = q, v
        if quantize:
            qm = np.round(q / P["q_quant"]) * P["q_quant"]
            vm = np.round(v / P["qd_quant"]) * P["qd_quant"]
        i_cmd = ctrl.step(t, r, qm, vm)
        i_app = np.round(i_cmd / P["i_quant"]) * P["i_quant"] if quantize else i_cmd
        Q[k], R[k], I[k], V[k] = q, r, i_cmd, v
        if hasattr(ctrl, "k1"):
            K1[k], K2[k], DH[k] = ctrl.k1, ctrl.k2, ctrl.dh
        if k == n:
            break
        for _ in range(nsub):
            tau_g = (P["g_sin"] * np.sin(q) + P["g_cos"] * np.cos(q)) if (payload and gravity) else 0.0
            drive = K_TAU * i_app - P["tau_off"] - tau_g
            if v == 0.0:
                if abs(drive) <= Fs:
                    continue          # stuck: static friction holds the joint
                a = (drive - Fc * np.sign(drive)) / J
            else:
                a = (drive - Fc * np.sign(v) - Fv * v) / J
            v_new = v + h * a
            if v != 0.0 and np.sign(v_new) != np.sign(v) and Fc > 0:
                v_new = 0.0           # velocity reversal -> re-check sticking
            q += h * 0.5 * (v + v_new)
            v = v_new
    return dict(t=tl, q=Q, ref=R, i=I, qd=V, k1=K1, k2=K2, dhat=DH)


def window_metrics(t, ref, q, i=None, t0=10.0):
    m = t >= t0
    e = np.rad2deg(ref[m] - q[m])
    out = dict(rmse=float(np.sqrt(np.mean(e ** 2))), mae=float(np.mean(np.abs(e))),
               maxe=float(np.max(np.abs(e))),
               nrmse=float(100 * np.sqrt(np.mean(e ** 2)) / np.sqrt(np.mean(np.rad2deg(ref[m]) ** 2))))
    if i is not None:
        out["i_rms"] = float(np.sqrt(np.mean(i[m] ** 2)))
        out["i_peak"] = float(np.max(np.abs(i[m])))
    return out
