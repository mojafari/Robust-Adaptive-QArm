"""Reproduce the three simulation scenarios of the paper (Section 4.1).

Outputs (in ../../results/original):
  Sim_<case>_Traj.pdf   joint angles (0-20 s window) and tracking errors
  Sim_<case>_U.pdf      joint torques
  Sim_<case>_Adapt.pdf  adaptive gains k1, k2 and disturbance estimates
  metrics.json          RMSE / NRMSE values (Table 2 of the manuscript)
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from qarm_sim.simulate import (dist_constant, dist_pulses, ref_constant,  # noqa: E402
                               ref_sinusoidal, simulate_original, tracking_metrics)
from qarm_sim import plotting as P  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "..", "..", "results", "original")
os.makedirs(OUT, exist_ok=True)
P.setup()

CASES = {
    "Const": (ref_constant, dist_constant, 20.0),
    "Sine": (ref_sinusoidal, dist_constant, 20.0),
    "Const_Ex": (ref_constant, dist_pulses, 100.0),
}


def plot_case(name, s, tw):
    t = s["t"]
    q = np.rad2deg(s["q"])
    r = np.rad2deg(s["ref"])
    e = r - q
    m = t <= tw
    # --- trajectories + errors
    fig, ax = plt.subplots(4, 2, figsize=(10, 7.2), sharex="col")
    for i in range(4):
        a = ax[i, 0]
        a.plot(t[m], q[m, i], color=P.GREEN, label="Proposed")
        a.plot(t[m], r[m, i], "--", color=P.RED, label="Desired")
        a.set_ylabel(f"{P.JOINT_NAMES[i]} angle\n[deg]")
        b = ax[i, 1]
        b.plot(t[m], e[m, i], color=P.GREEN)
        b.set_ylabel(f"{P.JOINT_NAMES[i]} error\n[deg]")
        if name == "Const_Ex":
            mm = (t >= 10) & m
            lim = 1.15 * np.abs(e[mm, i]).max()
            b.set_ylim(-lim, lim)
            for tk in (20, 40, 60, 80):
                b.axvspan(tk, tk + 1, color=P.ORANGE, alpha=0.25, lw=0)
    ax[0, 0].legend(ncol=2, loc="best")
    ax[3, 0].set_xlabel("Time [s]")
    ax[3, 1].set_xlabel("Time [s]")
    for a in ax.ravel():
        a.set_xlim(0, tw)
    fig.align_ylabels(ax[:, 0]); fig.align_ylabels(ax[:, 1])
    P.save(fig, os.path.join(OUT, f"Sim_{name}_Traj.pdf"))
    # --- torques
    fig, ax = plt.subplots(4, 1, figsize=(8, 6.4), sharex=True)
    for i in range(4):
        ax[i].plot(t[m], s["tau"][m, i], color=P.BLUE)
        ax[i].set_ylabel(rf"$\tau_{i+1}$ [N$\cdot$m]")
        ax[i].set_xlim(0, tw)
    ax[3].set_xlabel("Time [s]")
    fig.align_ylabels(ax)
    P.save(fig, os.path.join(OUT, f"Sim_{name}_U.pdf"))
    # --- adaptive variables
    fig, ax = plt.subplots(3, 1, figsize=(8, 6.0), sharex=True)
    ax[0].plot(t, s["k1"], color=P.MAGENTA)
    ax[0].set_ylabel(r"$k_1$ [s$^{-2}$]")
    ax[1].plot(t, s["k2"], color=P.PURPLE)
    ax[1].set_ylabel(r"$k_2$ [s$^{-1}$]")
    cols = [P.GREEN, P.BLUE, P.ORANGE, P.CYAN]
    for i in range(4):
        ax[2].plot(t, s["dhat"][:, i], color=cols[i], label=rf"$\hat d_{i+1}$")
    ax[2].axhline(0.5, color=P.RED, ls="--", lw=1.0, label=r"$d_0=0.5$")
    ax[2].set_ylabel(r"$\hat d_i$ [rad/s$^2$]")
    ax[2].legend(ncol=5, loc="lower right")
    ax[2].set_xlabel("Time [s]")
    for a in ax:
        a.set_xlim(0, tw)
    fig.align_ylabels(ax)
    P.save(fig, os.path.join(OUT, f"Sim_{name}_Adapt.pdf"))


def main():
    res = {}
    for name, (rf, df, tw) in CASES.items():
        s = simulate_original(rf, df, T=120.0, h=0.05)
        m = tracking_metrics(s["t"], s["ref"], s["q"], t_ss=5.0)
        res[name] = {k: (v.tolist() if hasattr(v, "tolist") else float(v)) for k, v in m.items()}
        res[name]["k1_range"] = [float(s["k1"].min()), float(s["k1"].max())]
        res[name]["k2_range"] = [float(s["k2"].min()), float(s["k2"].max())]
        res[name]["tau_peak"] = np.abs(s["tau"]).max(0).tolist()
        res[name]["dhat_final"] = s["dhat"][-1].tolist()
        plot_case(name, s, tw)
        print(name, "ss RMSE", np.round(m["rmse_ss"], 4), "agg", round(m["rmse_ss_agg"], 4),
              "k1", res[name]["k1_range"], "k2", res[name]["k2_range"])
    with open(os.path.join(OUT, "metrics.json"), "w") as f:
        json.dump(res, f, indent=1)


if __name__ == "__main__":
    main()
