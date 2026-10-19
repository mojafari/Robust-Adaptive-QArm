"""Hardware-representative wrist simulation (new Section 4.4 of the revision).

1. Validation   : the five experimental conditions with the exact controller
                  settings recovered from the logs, compared with the measured
                  RMSE (window t >= 10 s).
2. Counterfactual comparisons (not possible on the hardware any more):
     - PID without the 5 s start-up gate,
     - proposed controller with one common parameter set and the Quanser
       current limits in all five tests,
3. Sim-to-real gap decomposition: starting from the nominal model used for
   design, the identified non-idealities are added one at a time.
"""
import json
import os
import sys
from multiprocessing import Pool

import numpy as np
import scipy.io as sio

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from qarm_sim.wrist_hw import (QUANSER_LIMITS, WRIST_PLANT, ProposedWristHW,  # noqa: E402
                               QuanserWristPID, simulate_wrist, window_metrics)
from qarm_sim import plotting as P  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
OUT = os.path.join(ROOT, "results", "wrist_hw")
DATA = os.path.join(ROOT, "data", "experimental")
os.makedirs(OUT, exist_ok=True)
P.setup()

TESTS = {1: (45, 1 / 25, False), 2: (75, 1 / 25, False), 3: (75, 1 / 100, False),
         4: (45, 1 / 25, True), 5: (75, 1 / 25, True)}
RELAXED = (-15.0, 10.0)
HW = {1: dict(alpha1=2, alpha2=1, alpha_d=1, k1d=200, k2d=100, limits=QUANSER_LIMITS),
      2: dict(alpha1=2, alpha2=1, alpha_d=1, k1d=200, k2d=100, limits=QUANSER_LIMITS),
      3: dict(alpha1=2, alpha2=1, alpha_d=1, k1d=200, k2d=100, limits=QUANSER_LIMITS),
      4: dict(alpha1=4, alpha2=1, alpha_d=1, k1d=400, k2d=200, limits=RELAXED),
      5: dict(alpha1=2, alpha2=1, alpha_d=1, k1d=400, k2d=210, limits=RELAXED)}
COMMON_A = dict(alpha1=2, alpha2=1, alpha_d=1, k1d=200, k2d=100, limits=QUANSER_LIMITS)
COMMON_B = dict(alpha1=2, alpha2=1, alpha_d=1, k1d=400, k2d=200, limits=QUANSER_LIMITS)


def prop(cfg, k0=(1.0, 1.0), variant="model"):
    c = dict(cfg)
    lim = c.pop("limits")
    return ProposedWristHW(c["alpha1"], c["alpha2"], c["alpha_d"], c["k1d"], c["k2d"], limits=lim,
                           k10=k0[0], k20=k0[1], variant=variant)


def job(args):
    tag, j, kind, opts = args
    amp, f, pl = TESTS[j]
    if kind == "pid":
        ctrl = QuanserWristPID(gate_time=opts.get("gate", 5.0))
    else:
        ctrl = prop(opts["cfg"], k0=opts.get("k0", (1.0, 1.0)), variant=opts.get("variant", "model"))
    simkw = {k: opts[k] for k in ("quantize", "friction", "gravity", "true_inertia") if k in opts}
    if "plant" in opts:
        simkw["plant"] = opts["plant"]
    S = simulate_wrist(ctrl, amp, f, payload=pl, **simkw)
    m = window_metrics(S["t"], S["ref"], S["q"], S["i"], t0=10.0)
    if opts.get("save"):
        np.savez_compressed(os.path.join(OUT, f"sim_{tag}_T{j}.npz"), **S)
    return tag, j, m


def exp_metrics():
    out = {}
    for j in TESTS:
        for ctrl in ("PID", "Nonlinear"):
            d = sio.loadmat(os.path.join(DATA, f"{ctrl}_Test{j}.mat"))["data"]
            out[(ctrl, j)] = window_metrics(d[0], d[1], d[5], d[12], t0=10.0)
    return out


def main():
    jobs = []
    for j in TESTS:
        jobs.append(("pid_hw", j, "pid", dict(save=True)))
        jobs.append(("prop_hw", j, "prop", dict(cfg=HW[j], save=True)))
        jobs.append(("pid_nogate", j, "pid", dict(gate=0.0)))
        jobs.append(("prop_commonA", j, "prop", dict(cfg=COMMON_A)))
        jobs.append(("prop_commonB", j, "prop", dict(cfg=COMMON_B)))
        # ablation: fixed-gain feedback linearisation with the target gains
        jobs.append(("fl_target", j, "prop", dict(cfg=dict(COMMON_A, alpha1=0, alpha2=0, alpha_d=0),
                                                  k0=(200.0, 100.0))))
        # acceleration-based learning signal (filtered acceleration), common settings
        jobs.append(("prop_accel", j, "prop", dict(cfg=COMMON_A, variant="accel")))
    # sim-to-real decomposition on Tests 1 and 3 (both controllers; PID without gate)
    NO_OFFSET = dict(WRIST_PLANT, tau_off=0.0)
    steps = [("ideal", dict(quantize=False, friction=False, gravity=False, true_inertia=False,
                            plant=NO_OFFSET)),
             ("+inertia", dict(quantize=False, friction=False, gravity=False, true_inertia=True,
                               plant=NO_OFFSET)),
             ("+friction", dict(quantize=False, friction=True, gravity=False, true_inertia=True,
                                plant=NO_OFFSET)),
             ("+offset", dict(quantize=False, friction=True, gravity=False, true_inertia=True)),
             ("+quantisation", dict(quantize=True, friction=True, gravity=False, true_inertia=True))]
    for j in (1, 3):
        for name, kw in steps:
            jobs.append((f"gap_pid_{name}", j, "pid", dict(kw, gate=0.0)))
            jobs.append((f"gap_prop_{name}", j, "prop", dict(kw, cfg=HW[j])))
    with Pool(2) as pool:
        res = pool.map(job, jobs, chunksize=1)
    R = {}
    for tag, j, m in res:
        R.setdefault(tag, {})[j] = m
    E = exp_metrics()
    R["experiment_PID"] = {j: E[("PID", j)] for j in TESTS}
    R["experiment_Proposed"] = {j: E[("Nonlinear", j)] for j in TESTS}
    with open(os.path.join(OUT, "wrist_hw_results.json"), "w") as f:
        json.dump(R, f, indent=1)
    for tag in ("experiment_PID", "pid_hw", "pid_nogate", "experiment_Proposed", "prop_hw",
                "prop_commonA", "prop_commonB", "fl_target", "prop_accel"):
        print(f"{tag:22s}", " ".join(f"{R[tag][j]['rmse']:6.2f}" for j in TESTS))
    for j in (1, 3):
        print("gap test", j, {n: (round(R[f"gap_pid_{n}"][j]["rmse"], 2), round(R[f"gap_prop_{n}"][j]["rmse"], 2))
                              for n, _ in steps})
    figure_validation()


def figure_validation():
    fig, ax = plt.subplots(2, 2, figsize=(11, 6.4), sharex="col")
    for col, j in enumerate((1, 3)):
        t1 = 60 if j == 1 else 120
        for row, (ctrl, tag, colr) in enumerate((("PID", "pid_hw", P.BLUE),
                                                 ("Nonlinear", "prop_hw", P.GREEN))):
            d = sio.loadmat(os.path.join(DATA, f"{ctrl}_Test{j}.mat"))["data"]
            S = np.load(os.path.join(OUT, f"sim_{tag}_T{j}.npz"))
            a = ax[row, col]
            m = d[0] <= t1
            a.plot(d[0][m], np.rad2deg(d[5][m]), color=colr, lw=1.3, label="experiment")
            ms = S["t"] <= t1
            a.plot(S["t"][ms], np.rad2deg(S["q"][ms]), color="k", lw=1.0, ls="-.", label="simulation")
            a.plot(d[0][m], np.rad2deg(d[1][m]), "--", color=P.RED, lw=1.0, label="desired")
            a.set_xlim(0, t1)
            name = "PID (Quanser)" if ctrl == "PID" else "Proposed"
            a.set_title(f"Test {j}, {name}")
            if col == 0:
                a.set_ylabel("Wrist angle [deg]")
            if row == 1:
                a.set_xlabel("Time [s]")
    ax[0, 0].legend(ncol=3, loc="lower left", fontsize=8.5)
    P.save(fig, os.path.join(OUT, "Wrist_Sim_vs_Exp.pdf"))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "fig":
        figure_validation()
    else:
        main()
