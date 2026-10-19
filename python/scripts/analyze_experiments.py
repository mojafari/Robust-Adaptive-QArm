"""Quantitative analysis of the QArm wrist hardware experiments.

Reads data/experimental/{PID,Nonlinear}_Test{1..5}.mat (Simulink 'To File'
format, 19 x N for the proposed controller, 13 x N for PID) and produces
  results/experiments/metrics.json   all metrics, several evaluation windows
  results/experiments/Test{j}.pdf    angle / error / current figure per test
  results/experiments/Exp_Adapt.pdf  adaptive variables of all tests
  results/experiments/Exp_StickSlip.pdf  zoom illustrating stick-slip
  results/experiments/hw_params.json  controller parameters recovered from logs

Row layout of the logged matrices (see data/experimental/README.md):
  1 time, 2 wrist command [rad], 3-6 joint angles [rad], 7-10 joint speeds,
  11 filtered wrist speed, 12 filtered wrist acceleration, 13 wrist current
  command [A], 14-17 dhat_1..4, 18 k1, 19 k2 (rows 14-19 proposed only).
"""
import json
import os
import sys

import numpy as np
import scipy.io as sio

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from qarm_sim import plotting as P  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
DATA = os.path.join(ROOT, "data", "experimental")
OUT = os.path.join(ROOT, "results", "experiments")
os.makedirs(OUT, exist_ok=True)
P.setup()

TESTS = {1: (45, 1 / 25, 0), 2: (75, 1 / 25, 0), 3: (75, 1 / 100, 0),
         4: (45, 1 / 25, 32), 5: (75, 1 / 25, 32)}
WINDOWS = {"full": 0.0, "t5": 5.0, "t10": 10.0, "t25": 25.0}
I_UP, I_LO = 0.015, -0.05


def load(ctrl, j):
    d = sio.loadmat(os.path.join(DATA, f"{ctrl}_Test{j}.mat"))["data"]
    out = dict(t=d[0], ref=d[1], q=d[5], qd=d[9], qd_f=d[10], qdd_f=d[11], i=d[12],
               q_other=d[2:5])
    if d.shape[0] >= 19:
        out.update(dhat=d[13:17], k1=d[17], k2=d[18])
    return out


def metrics(D, t0):
    m = D["t"] >= t0
    e = np.rad2deg(D["ref"][m] - D["q"][m])
    r = np.rad2deg(D["ref"][m])
    i = D["i"][m]
    return dict(rmse=float(np.sqrt(np.mean(e ** 2))),
                nrmse=float(100 * np.sqrt(np.mean(e ** 2)) / np.sqrt(np.mean(r ** 2))),
                mae=float(np.mean(np.abs(e))), maxe=float(np.max(np.abs(e))),
                t_maxe=float(D["t"][m][np.argmax(np.abs(e))]),
                i_rms=float(np.sqrt(np.mean(i ** 2))), i_peak=float(np.max(np.abs(i))),
                frac_at_upper=float(np.mean(np.abs(i - I_UP) < 1e-6)),
                frac_above_upper=float(np.mean(i > I_UP + 1e-6)),
                frac_stuck=float(np.mean(np.abs(D["qd"][m]) < np.deg2rad(0.5))))


def recover_parameters(D):
    """Least-squares identification of (alpha1, alpha2, alpha_d, k1d, k2d) of
    the acceleration-free law  c = (k1d-k1) e1 + (k2d-k2) e2 + dhat,
    k1' = a1 c e1, k2' = a2 c e2, dhat' = -ad c, from the logged signals."""
    t, r, q, qd = D["t"], D["ref"], D["q"], D["qd"]
    k1, k2, dh = D["k1"], D["k2"], D["dhat"][3]
    dt = np.gradient(t)
    r2 = np.gradient(r, t)
    e1, e2 = r - q, r2 - qd
    ddh = np.gradient(dh, t)
    m = t > 1
    X = np.column_stack([e1, e2, -k1 * e1 - k2 * e2 + dh])[m]
    co = np.linalg.lstsq(X, -ddh[m], rcond=None)[0]
    ad = co[2]
    k1d, k2d = co[0] / ad, co[1] / ad
    c = -ddh / ad
    mm = m & (np.abs(e1) > 0.02) & (np.abs(e2) > 0.02)
    a1 = np.linalg.lstsq((c * e1)[mm][:, None], np.gradient(k1, t)[mm], rcond=None)[0][0]
    a2 = np.linalg.lstsq((c * e2)[mm][:, None], np.gradient(k2, t)[mm], rcond=None)[0][0]
    # torque-to-current gain: I = kappa * (k1 e1 + k2 e2 + r2dot - dhat)
    r2d = np.gradient(r2, t)
    v = k1 * e1 + k2 * e2 + r2d - dh
    sel = m & (D["i"] < I_UP - 1e-4) & (D["i"] > -0.0499)
    kappa = np.linalg.lstsq(v[sel][:, None], D["i"][sel], rcond=None)[0][0]
    del dt
    return dict(alpha1=float(a1), alpha2=float(a2), alpha_d=float(ad), k1d=float(k1d),
                k2d=float(k2d), current_gain=float(kappa))


def replay(D, prm):
    """Re-integrate the adaptive laws with the rounded parameters and compare."""
    t, r, q, qd = D["t"], D["ref"], D["q"], D["qd"]
    r2 = np.concatenate([[0], np.diff(r) / np.diff(t)])
    e1, e2 = r - q, r2 - qd
    k1 = np.empty_like(t); k2 = np.empty_like(t); dh = np.empty_like(t)
    k1[0], k2[0], dh[0] = D["k1"][0], D["k2"][0], D["dhat"][3][0]
    for i in range(len(t) - 1):
        c = (prm["k1d"] - k1[i]) * e1[i] + (prm["k2d"] - k2[i]) * e2[i] + dh[i]
        h = t[i + 1] - t[i]
        k1[i + 1] = k1[i] + h * prm["alpha1"] * c * e1[i]
        k2[i + 1] = k2[i] + h * prm["alpha2"] * c * e2[i]
        dh[i + 1] = dh[i] - h * prm["alpha_d"] * c
    rel = lambda a, b: float(np.max(np.abs(a - b)) / np.ptp(b))
    return dict(k1=rel(k1, D["k1"]), k2=rel(k2, D["k2"]), dhat=rel(dh, D["dhat"][3]))


def fig_test(j, Pd, Nd, limits_prop):
    amp, f, pl = TESTS[j]
    fig, ax = plt.subplots(3, 1, figsize=(9, 7.4), sharex=True)
    for a in ax:
        a.axvspan(0, 5, color=P.GREY, alpha=0.15, lw=0)
    ax[0].plot(Pd["t"], np.rad2deg(Pd["q"]), color=P.BLUE, lw=1.2, label="PID (Quanser)")
    ax[0].plot(Nd["t"], np.rad2deg(Nd["q"]), color=P.GREEN, lw=1.2, label="Proposed")
    ax[0].plot(Nd["t"], np.rad2deg(Nd["ref"]), "--", color=P.RED, lw=1.2, label="Desired")
    ax[0].set_ylabel("Wrist angle [deg]")
    ax[0].legend(ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0), frameon=False)
    ax[1].plot(Pd["t"], np.rad2deg(Pd["ref"] - Pd["q"]), color=P.BLUE, lw=1.0)
    ax[1].plot(Nd["t"], np.rad2deg(Nd["ref"] - Nd["q"]), color=P.GREEN, lw=1.0)
    ax[1].set_ylabel("Tracking error [deg]")
    ax[2].plot(Pd["t"], 1e3 * Pd["i"], color=P.BLUE, lw=1.0, label="PID")
    ax[2].plot(Nd["t"], 1e3 * Nd["i"], color=P.GREEN, lw=1.0, label="Proposed")
    ax[2].axhline(1e3 * I_UP, color="k", ls=":", lw=1.0)
    ax[2].text(119, 1e3 * I_UP + 1.0, "+15 mA (Quanser upper limit)", ha="right", va="bottom", fontsize=8.5)
    if limits_prop[1] > I_UP:
        ax[2].text(1, 1e3 * I_UP + 1.0, "upper limit relaxed for the proposed controller",
                   ha="left", va="bottom", fontsize=8.5, color=P.GREEN)
    ax[2].set_ylabel("Current command [mA]")
    ax[2].set_xlabel("Time [s]")
    ax[0].set_xlim(0, 120)
    fig.align_ylabels(ax)
    P.save(fig, os.path.join(OUT, f"Test{j}.pdf"))


def fig_adapt(allN):
    fig, ax = plt.subplots(3, 1, figsize=(9, 7.0), sharex=True)
    cols = [P.GREEN, P.BLUE, P.ORANGE, P.MAGENTA, P.CYAN]
    for j, D in allN.items():
        ax[0].plot(D["t"], D["k1"], color=cols[j - 1], lw=1.0, label=f"Test {j}")
        ax[1].plot(D["t"], D["k2"], color=cols[j - 1], lw=1.0)
        ax[2].plot(D["t"], D["dhat"][3], color=cols[j - 1], lw=1.0)
    ax[0].set_ylabel(r"$k_1$ [s$^{-2}$]")
    ax[1].set_ylabel(r"$k_2$ [s$^{-1}$]")
    ax[2].set_ylabel(r"$\hat d_4$ [rad/s$^2$]")
    ax[2].set_xlabel("Time [s]")
    ax[0].legend(ncol=5, loc="upper left")
    ax[0].set_xlim(0, 120)
    fig.align_ylabels(ax)
    P.save(fig, os.path.join(OUT, "Exp_Adapt.pdf"))


def fig_stickslip(Pd, Nd, t0=30, t1=60):
    fig, ax = plt.subplots(3, 1, figsize=(9, 6.6), sharex=True)
    for D, col, lab in ((Pd, P.BLUE, "PID (Quanser)"), (Nd, P.GREEN, "Proposed")):
        m = (D["t"] >= t0) & (D["t"] <= t1)
        ax[0].plot(D["t"][m], np.rad2deg(D["q"][m]), color=col, label=lab)
        ax[1].plot(D["t"][m], np.rad2deg(D["qd"][m]), color=col, lw=1.0)
        ax[2].plot(D["t"][m], 1e3 * D["i"][m], color=col, lw=1.0)
    m = (Nd["t"] >= t0) & (Nd["t"] <= t1)
    ax[0].plot(Nd["t"][m], np.rad2deg(Nd["ref"][m]), "--", color=P.RED, label="Desired")
    ax[2].axhline(15, color="k", ls=":", lw=1.0)
    ax[0].set_ylabel("Wrist angle [deg]")
    ax[1].set_ylabel("Wrist speed [deg/s]")
    ax[2].set_ylabel("Current command [mA]")
    ax[2].set_xlabel("Time [s]")
    ax[0].legend(ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0), frameon=False)
    ax[0].set_xlim(t0, t1)
    fig.align_ylabels(ax)
    P.save(fig, os.path.join(OUT, "Exp_StickSlip.pdf"))


def main():
    res = {"windows": WINDOWS, "tests": {}}
    allN, params = {}, {}
    for j in TESTS:
        Pd, Nd = load("PID", j), load("Nonlinear", j)
        allN[j] = Nd
        prm = recover_parameters(Nd)
        prm_r = dict(alpha1=round(prm["alpha1"]), alpha2=round(prm["alpha2"]),
                     alpha_d=round(prm["alpha_d"]), k1d=round(prm["k1d"] / 10) * 10,
                     k2d=round(prm["k2d"] / 10) * 10)
        rep = replay(Nd, prm_r)
        limits_prop = (I_LO, 10.0) if Nd["i"].max() > I_UP + 1e-6 else (I_LO, I_UP)
        params[j] = dict(fitted=prm, rounded=prm_r, replay_rel_error=rep,
                         upper_current_limit_exceeded=bool(Nd["i"].max() > I_UP + 1e-6))
        T = {}
        for w, t0 in WINDOWS.items():
            mp, mn = metrics(Pd, t0), metrics(Nd, t0)
            T[w] = dict(PID=mp, Proposed=mn,
                        rmse_red=100 * (mp["rmse"] - mn["rmse"]) / mp["rmse"],
                        mae_red=100 * (mp["mae"] - mn["mae"]) / mp["mae"],
                        maxe_red=100 * (mp["maxe"] - mn["maxe"]) / mp["maxe"])
        T["pid_current_zero_before_5s"] = bool(np.all(Pd["i"][Pd["t"] <= 5.0] == 0))
        T["proposed_current_nonzero_before_5s"] = bool(np.any(Nd["i"][Nd["t"] <= 5.0] != 0))
        T["k1_range"] = [float(Nd["k1"].min()), float(Nd["k1"].max())]
        T["k2_range"] = [float(Nd["k2"].min()), float(Nd["k2"].max())]
        T["dhat4_range"] = [float(Nd["dhat"][3].min()), float(Nd["dhat"][3].max())]
        T["dhat123_maxabs"] = float(np.abs(Nd["dhat"][:3]).max())
        T["other_joints_maxdev_deg"] = float(np.rad2deg(np.abs(Nd["q_other"] - Nd["q_other"][:, :1]).max()))
        T["Ts"] = float(np.median(np.diff(Nd["t"])))
        # within-trial variability: RMSE of each full reference cycle in [20, 120) s
        # (Test 3: two half-cycles of 50 s, since its period is 100 s)
        per = 25.0 if j != 3 else 50.0
        cyc = {}
        for lab, D in (("PID", Pd), ("Proposed", Nd)):
            vals = []
            for k0 in np.arange(20.0, 120.0 - 1e-9, per):
                m = (D["t"] >= k0) & (D["t"] < k0 + per)
                e = np.rad2deg(D["ref"][m] - D["q"][m])
                vals.append(float(np.sqrt(np.mean(e ** 2))))
            cyc[lab] = dict(values=vals, mean=float(np.mean(vals)), std=float(np.std(vals, ddof=1)))
        T["cycle_rmse"] = cyc
        res["tests"][j] = T
        fig_test(j, Pd, Nd, limits_prop)
        print(f"Test {j}: t>=10 RMSE PID {T['t10']['PID']['rmse']:.2f} Prop {T['t10']['Proposed']['rmse']:.2f}"
              f" | params {prm_r} replay {rep}")
    fig_adapt(allN)
    fig_stickslip(load("PID", 1), load("Nonlinear", 1))
    with open(os.path.join(OUT, "metrics.json"), "w") as f:
        json.dump(res, f, indent=1)
    with open(os.path.join(OUT, "hw_params.json"), "w") as f:
        json.dump(params, f, indent=1)


if __name__ == "__main__":
    main()
