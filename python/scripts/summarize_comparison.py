"""Tables and figures for the comparative simulation study."""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from qarm_sim import plotting as P  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "..", "..", "results", "comparison")
P.setup()
CTRLS = ["PID+G", "FL", "FL+DE", "Proposed"]
COLS = {"PID+G": P.BLUE, "FL": P.GREY, "FL+DE": P.ORANGE, "Proposed": P.GREEN}


def main():
    R = json.load(open(os.path.join(OUT, "comparison_raw.json")))
    summ = {}
    for scen in "ABCD":
        summ[scen] = {}
        for c in CTRLS:
            rows = [r for r in R if r["scen"] == scen and r["ctrl"] == c]
            ok = [r for r in rows if not r["diverged"]]
            agg = np.array([r["rmse_ss_agg"] for r in ok])
            pj = np.array([r["rmse_ss"] for r in ok])
            tr = np.array([r["tau_rms"] for r in ok])
            tp = np.array([r["tau_peak"] for r in ok])
            s = dict(n=len(rows), n_diverged=len(rows) - len(ok),
                     agg_mean=float(agg.mean()), agg_std=float(agg.std(ddof=1)) if len(agg) > 1 else 0.0,
                     agg_median=float(np.median(agg)), agg_max=float(agg.max()),
                     joint_mean=pj.mean(0).tolist(),
                     tau_rms_mean=tr.mean(0).tolist(), tau_peak_max=tp.max(0).tolist(),
                     sat_frac_max=np.max([r["sat_frac"] for r in ok], axis=0).tolist())
            if c == "Proposed":
                k1 = np.array([r["k1_range"] for r in ok]); k2 = np.array([r["k2_range"] for r in ok])
                s["k1_range"] = [float(k1[:, 0].min()), float(k1[:, 1].max())]
                s["k2_range"] = [float(k2[:, 0].min()), float(k2[:, 1].max())]
            summ[scen][c] = s
    # paired comparison Proposed vs others in MC scenarios
    for scen in "BD":
        for c in ("PID+G", "FL+DE"):
            a = {r["seed"]: r["rmse_ss_agg"] for r in R if r["scen"] == scen and r["ctrl"] == "Proposed" and not r["diverged"]}
            b = {r["seed"]: r["rmse_ss_agg"] for r in R if r["scen"] == scen and r["ctrl"] == c and not r["diverged"]}
            common = sorted(set(a) & set(b))
            wins = sum(a[s] < b[s] for s in common)
            summ[scen][f"Proposed_better_than_{c}"] = f"{wins}/{len(common)}"
    json.dump(summ, open(os.path.join(OUT, "comparison_summary.json"), "w"), indent=1)
    for scen in "ABCD":
        print(scen, {c: (round(summ[scen][c]["agg_mean"], 3), round(summ[scen][c]["agg_std"], 3),
                         np.round(summ[scen][c]["joint_mean"], 3).tolist()) for c in CTRLS})
    print({s: {k: v for k, v in summ[s].items() if k.startswith("Proposed_better")} for s in "BD"})

    # --- figure 1: time histories, scenario A (shoulder and wrist errors)
    fig, ax = plt.subplots(2, 2, figsize=(10.5, 5.6), sharex=True)
    for c in CTRLS:
        L = np.load(os.path.join(OUT, f"log_A_{c.replace('+', 'p')}.npz"))
        e = np.rad2deg(L["ref"] - L["q"])
        for k, j in enumerate((1, 3)):
            ax[0, k].plot(L["t"], e[:, j], color=COLS[c], lw=1.2, label=c)
            ax[1, k].plot(L["t"], e[:, j], color=COLS[c], lw=1.2)
    for k, j in enumerate((1, 3)):
        ax[0, k].set_title(f"{P.JOINT_NAMES[j]} joint")
        ax[0, k].set_ylabel("Error [deg]")
        ax[1, k].set_ylabel("Error [deg] (zoom)")
        ax[1, k].set_ylim(-1.5, 1.5) if j == 1 else ax[1, k].set_ylim(-6, 6)
        ax[1, k].set_xlabel("Time [s]")
    h, l = ax[0, 0].get_legend_handles_labels()
    fig.legend(h, l, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.03), frameon=False)
    ax[0, 0].set_xlim(0, 40)
    P.save(fig, os.path.join(OUT, "Sim_Comparison_A.pdf"))

    # --- figure 2: Monte Carlo distributions (B and D)
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.8))
    for k, scen in enumerate("BD"):
        data = [[r["rmse_ss_agg"] for r in R if r["scen"] == scen and r["ctrl"] == c and not r["diverged"]]
                for c in CTRLS]
        bp = ax[k].boxplot(data, tick_labels=CTRLS, patch_artist=True, widths=0.6,
                           medianprops=dict(color="k", lw=1.5))
        for patch, c in zip(bp["boxes"], CTRLS):
            patch.set_facecolor(COLS[c]); patch.set_alpha(0.45)
        ax[k].set_yscale("log")
        ax[k].set_ylabel("Pooled steady-state RMSE [deg]")
        ax[k].set_title("B: parametric uncertainty" if scen == "B" else
                        "D: uncertainty + noise + saturation")
    P.save(fig, os.path.join(OUT, "Sim_Comparison_MC.pdf"))


if __name__ == "__main__":
    main()
