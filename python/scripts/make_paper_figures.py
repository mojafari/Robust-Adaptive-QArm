"""Compact figures used in the revised manuscript (Revision 1).

Run after run_original.py, run_comparison.py, analyze_experiments.py and
run_wrist_hw.py.  Writes PDF (and PNG previews) to ../../results/paper_figures.
"""
import json
import os
import sys

import numpy as np
import scipy.io as sio

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from qarm_sim.simulate import (dist_constant, dist_pulses, ref_constant,  # noqa: E402
                               ref_sinusoidal, simulate_original)
from qarm_sim import plotting as P  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
RES = os.path.join(ROOT, "results")
DATA = os.path.join(ROOT, "data", "experimental")
OUT = os.path.join(RES, "paper_figures")
os.makedirs(OUT, exist_ok=True)

W = 7.0          # figure width [in]; included at \textwidth
P.setup()
plt.rcParams.update({"font.size": 8.5, "axes.labelsize": 8.5, "axes.titlesize": 8.5,
                     "legend.fontsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
                     "lines.linewidth": 1.1, "figure.constrained_layout.use": True,
                     "figure.constrained_layout.w_pad": 0.06})
JN = P.JOINT_NAMES
CTRL_COL = {"PID+G": P.BLUE, "FL": P.GREY, "FL+DE": P.ORANGE, "Proposed": P.GREEN}
CTRL_LAB = {"PID+G": "PID+G", "FL": "FL", "FL+DE": "FL+DE", "Proposed": "Designed"}


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=170)
    plt.close(fig)


# ---------------------------------------------------------------------------
def simulation_figures():
    cases = {"Const": (ref_constant, dist_constant, 20.0, "Constant reference"),
             "Sine": (ref_sinusoidal, dist_constant, 20.0, "Sinusoidal reference"),
             "Pulse": (ref_constant, dist_pulses, 100.0, "Disturbance pulses")}
    sims = {k: simulate_original(v[0], v[1], T=120.0, h=0.05) for k, v in cases.items()}
    for key, (_, _, tw, _) in cases.items():
        s = sims[key]
        t = s["t"]; m = t <= tw
        q = np.rad2deg(s["q"]); r = np.rad2deg(s["ref"]); e = r - q
        fig, ax = plt.subplots(4, 2, figsize=(W, 3.9), sharex=True)
        for i in range(4):
            ax[i, 0].plot(t[m], q[m, i], color=P.GREEN, label="Designed")
            ax[i, 0].plot(t[m], r[m, i], "--", color=P.RED, label="Desired")
            ax[i, 0].set_ylabel(rf"$q_{i+1}$ [deg]")
            ax[i, 1].plot(t[m], e[m, i], color=P.GREEN)
            ax[i, 1].set_ylabel(rf"$e_{{1,{i+1}}}$ [deg]")
            if key == "Pulse":
                mm = (t >= 10) & m
                lim = 1.15 * np.abs(e[mm, i]).max()
                ax[i, 1].set_ylim(-lim, lim)
                for tk in (20, 40, 60, 80):
                    for a in ax[i]:
                        a.axvspan(tk, tk + 1, color=P.ORANGE, alpha=0.3, lw=0)
            else:
                ax[i, 1].set_ylim(-3, 3)
        ax[3, 0].set_xlabel("Time [s]"); ax[3, 1].set_xlabel("Time [s]")
        ax[0, 0].set_xlim(0, tw)
        h, l = ax[0, 0].get_legend_handles_labels()
        fig.legend(h, l, ncol=2, loc="outside upper center", frameon=False)
        fig.align_ylabels(ax[:, 0]); fig.align_ylabels(ax[:, 1])
        save(fig, f"Sim_{key}_Traj")

    # torques: 4 joints x 3 scenarios
    fig, ax = plt.subplots(4, 3, figsize=(W, 3.9), sharex="col")
    for c, (key, (_, _, tw, title)) in enumerate(cases.items()):
        s = sims[key]; t = s["t"]; m = t <= tw
        for i in range(4):
            ax[i, c].plot(t[m], s["tau"][m, i], color=P.BLUE)
            if c == 0:
                ax[i, c].set_ylabel(rf"$\tau_{i+1}$ [N$\cdot$m]")
        ax[0, c].set_title(title)
        ax[3, c].set_xlabel("Time [s]")
        ax[3, c].set_xlim(0, tw)
    fig.align_ylabels(ax[:, 0])
    save(fig, "Sim_Torques")

    # adaptive variables: k1, k2, dhat x 3 scenarios
    fig, ax = plt.subplots(3, 3, figsize=(W, 3.5), sharex="col")
    cols = [P.GREEN, P.BLUE, P.ORANGE, P.CYAN]
    for c, (key, (_, _, tw, title)) in enumerate(cases.items()):
        s = sims[key]; t = s["t"]; m = t <= tw
        ax[0, c].plot(t[m], s["k1"][m], color=P.MAGENTA)
        ax[1, c].plot(t[m], s["k2"][m], color=P.PURPLE)
        for i in range(4):
            ax[2, c].plot(t[m], s["dhat"][m, i], color=cols[i], label=rf"$\hat d_{i+1}$")
        ax[2, c].axhline(0.5, color=P.RED, ls="--", lw=0.9, label=r"$d_0$")
        ax[0, c].set_title(title)
        ax[2, c].set_xlabel("Time [s]")
        ax[2, c].set_xlim(0, tw)
    ax[0, 0].set_ylabel(r"$k_1$"); ax[1, 0].set_ylabel(r"$k_2$")
    ax[2, 0].set_ylabel(r"$\hat d_i$ [rad/s$^2$]")
    h, l = ax[2, 0].get_legend_handles_labels()
    fig.legend(h, l, ncol=5, loc="outside lower center", frameon=False)
    fig.align_ylabels(ax[:, 0])
    save(fig, "Sim_Adapt")


# ---------------------------------------------------------------------------
def comparison_figure():
    R = json.load(open(os.path.join(RES, "comparison", "comparison_raw.json")))
    fig, ax = plt.subplots(1, 3, figsize=(W, 2.5), width_ratios=[1, 1, 1.15])
    for c in ("PID+G", "FL", "FL+DE", "Proposed"):
        L = np.load(os.path.join(RES, "comparison", f"log_A_{c.replace('+', 'p')}.npz"))
        e = np.rad2deg(L["ref"] - L["q"])
        ax[0].plot(L["t"], e[:, 1], color=CTRL_COL[c], label=CTRL_LAB[c])
        ax[1].plot(L["t"], e[:, 3], color=CTRL_COL[c])
    ax[0].set_ylim(-1.5, 1.5); ax[1].set_ylim(-8, 8)
    ax[0].set_title("(a) Shoulder error, scenario A"); ax[1].set_title("(b) Wrist error, scenario A")
    for a in ax[:2]:
        a.set_xlim(0, 40); a.set_xlabel("Time [s]")
    ax[0].set_ylabel("Error [deg]")
    order = ["PID+G", "FL", "FL+DE", "Proposed"]
    data = [[r["rmse_ss_agg"] for r in R if r["scen"] == "B" and r["ctrl"] == c and not r["diverged"]]
            for c in order]
    bp = ax[2].boxplot(data, tick_labels=[CTRL_LAB[c] for c in order], patch_artist=True, widths=0.6,
                       medianprops=dict(color="k", lw=1.2), flierprops=dict(markersize=3))
    for patch, c in zip(bp["boxes"], order):
        patch.set_facecolor(CTRL_COL[c]); patch.set_alpha(0.5)
    ax[2].set_yscale("log")
    ax[2].set_ylabel("Steady-state RMSE [deg]")
    ax[2].set_title("(c) Scenario B, 30 plants")
    ax[2].tick_params(axis="x", labelsize=7)
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, ncol=4, loc="outside upper center", frameon=False)
    save(fig, "Sim_Comparison")


# ---------------------------------------------------------------------------
def load(ctrl, j):
    return sio.loadmat(os.path.join(DATA, f"{ctrl}_Test{j}.mat"))["data"]


def experiment_figures():
    titles = {1: r"Test 1: $45^\circ$, 1/25 Hz", 2: r"Test 2: $75^\circ$, 1/25 Hz",
              3: r"Test 3: $75^\circ$, 1/100 Hz", 4: r"Test 4: $45^\circ$, 1/25 Hz, 32 g",
              5: r"Test 5: $75^\circ$, 1/25 Hz, 32 g"}
    fig, ax = plt.subplots(5, 2, figsize=(W, 7.2), sharex=True, width_ratios=[1.25, 1])
    for k, j in enumerate(range(1, 6)):
        Pd, Nd = load("PID", j), load("Nonlinear", j)
        for a in ax[k]:
            a.axvspan(0, 5, color=P.GREY, alpha=0.18, lw=0)
        ax[k, 0].plot(Pd[0], np.rad2deg(Pd[5]), color=P.BLUE, lw=0.9, label="PID")
        ax[k, 0].plot(Nd[0], np.rad2deg(Nd[5]), color=P.GREEN, lw=0.9, label="Designed")
        ax[k, 0].plot(Nd[0], np.rad2deg(Nd[1]), "--", color=P.RED, lw=0.9, label="Desired")
        ax[k, 0].set_ylabel("Wrist angle [deg]")
        ax[k, 0].set_title(titles[j], loc="left")
        ax[k, 1].plot(Pd[0], 1e3 * Pd[12], color=P.BLUE, lw=0.7)
        ax[k, 1].plot(Nd[0], 1e3 * Nd[12], color=P.GREEN, lw=0.7)
        ax[k, 1].axhline(15, color="k", ls=":", lw=0.9, label="+15 mA limit")
        ax[k, 1].set_ylabel("Current [mA]")
        lo = 1e3 * min(Pd[12].min(), Nd[12].min()); hi = 1e3 * max(Pd[12].max(), Nd[12].max())
        ax[k, 1].set_ylim(lo - 3, max(hi, 15) + 5)
    ax[4, 0].set_xlabel("Time [s]"); ax[4, 1].set_xlabel("Time [s]")
    ax[0, 0].set_xlim(0, 120)
    h, l = ax[0, 0].get_legend_handles_labels()
    h2, l2 = ax[0, 1].get_legend_handles_labels()
    fig.legend(h + h2, l + l2, ncol=4, loc="outside upper center", frameon=False)
    fig.align_ylabels(ax[:, 0]); fig.align_ylabels(ax[:, 1])
    save(fig, "Exp_Tests")

    # stick-slip zoom, Test 1
    Pd, Nd = load("PID", 1), load("Nonlinear", 1)
    fig, ax = plt.subplots(3, 1, figsize=(W, 3.3), sharex=True)
    for D, col, lab in ((Pd, P.BLUE, "PID"), (Nd, P.GREEN, "Designed")):
        m = (D[0] >= 30) & (D[0] <= 60)
        ax[0].plot(D[0][m], np.rad2deg(D[5][m]), color=col, label=lab)
        ax[1].plot(D[0][m], np.rad2deg(D[9][m]), color=col, lw=0.8)
        ax[2].plot(D[0][m], 1e3 * D[12][m], color=col, lw=0.8)
    m = (Nd[0] >= 30) & (Nd[0] <= 60)
    ax[0].plot(Nd[0][m], np.rad2deg(Nd[1][m]), "--", color=P.RED, label="Desired")
    ax[2].axhline(15, color="k", ls=":", lw=0.9)
    ax[2].set_ylim(-24, 20)
    ax[0].set_ylabel("Angle [deg]"); ax[1].set_ylabel("Speed [deg/s]"); ax[2].set_ylabel("Current [mA]")
    ax[2].set_xlabel("Time [s]"); ax[0].set_xlim(30, 60)
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, ncol=3, loc="outside upper center", frameon=False)
    fig.align_ylabels(ax)
    save(fig, "Exp_StickSlip")

    # adaptive variables of the five tests
    fig, ax = plt.subplots(3, 1, figsize=(W, 3.3), sharex=True)
    cols = [P.GREEN, P.BLUE, P.ORANGE, P.MAGENTA, P.CYAN]
    for j in range(1, 6):
        D = load("Nonlinear", j)
        ax[0].plot(D[0], D[17], color=cols[j - 1], lw=0.9, label=f"Test {j}")
        ax[1].plot(D[0], D[18], color=cols[j - 1], lw=0.9)
        ax[2].plot(D[0], D[16], color=cols[j - 1], lw=0.9)
    ax[0].set_ylabel(r"$k_1$"); ax[1].set_ylabel(r"$k_2$"); ax[2].set_ylabel(r"$\hat d_4$ [rad/s$^2$]")
    ax[2].set_xlabel("Time [s]"); ax[0].set_xlim(0, 120)
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, ncol=5, loc="outside upper center", frameon=False)
    fig.align_ylabels(ax)
    save(fig, "Exp_Adapt")


# ---------------------------------------------------------------------------
def wrist_model_figure():
    fig, ax = plt.subplots(2, 2, figsize=(W, 3.6), sharex="col")
    for col, j in enumerate((1, 3)):
        t1 = 60 if j == 1 else 120
        for row, (ctrl, tag, colr, name) in enumerate((("PID", "pid_hw", P.BLUE, "PID"),
                                                       ("Nonlinear", "prop_hw", P.GREEN, "Designed"))):
            d = load(ctrl, j)
            S = np.load(os.path.join(RES, "wrist_hw", f"sim_{tag}_T{j}.npz"))
            a = ax[row, col]
            m = d[0] <= t1; ms = S["t"] <= t1
            a.plot(d[0][m], np.rad2deg(d[5][m]), color=colr, lw=1.1, label="Experiment")
            a.plot(S["t"][ms], np.rad2deg(S["q"][ms]), color="k", lw=0.9, ls="-.", label="Wrist model")
            a.plot(d[0][m], np.rad2deg(d[1][m]), "--", color=P.RED, lw=0.9, label="Desired")
            a.set_xlim(0, t1)
            a.set_title(f"Test {j}, {name}", loc="left")
            if col == 0:
                a.set_ylabel("Wrist angle [deg]")
            if row == 1:
                a.set_xlabel("Time [s]")
    h, l = ax[0, 0].get_legend_handles_labels()
    fig.legend(h, l, ncol=3, loc="outside upper center", frameon=False)
    save(fig, "Wrist_Model")


if __name__ == "__main__":
    simulation_figures()
    comparison_figure()
    experiment_figures()
    wrist_model_figure()
    print("figures written to", OUT)
