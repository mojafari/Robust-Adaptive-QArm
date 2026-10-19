"""Comparative simulation study requested in the review (new Section 4.2).

Controllers (all with the Table-1 design parameters where applicable):
  PID+G    independent-joint PID + nominal gravity compensation, gains chosen
           so that its nominal decoupled poles equal the linearised poles of
           the proposed controller at k1(0)=k2(0)=1 (equal-bandwidth tuning)
  FL       feedback linearisation, fixed gains k1=k2=1, no disturbance estimate
  FL+DE    FL with the disturbance estimator only (gain adaptation disabled)
  Proposed full controller (gain adaptation + disturbance estimate)

Scenarios (sinusoidal reference of the paper, constant matched disturbance
d0=0.5 rad/s^2 on every joint, controller sampled at 500 Hz, acceleration
estimate obtained by a 20 Hz second-order filter on the measured velocity):
  A  nominal plant
  B  parametric uncertainty, Monte Carlo (masses/inertias/CoM +-20 %,
     payload m_L ~ U[0, 0.3] kg)
  C  nominal plant + measurement noise/quantisation
  D  as B + measurement noise + torque saturation at the continuous ratings
     (joint friction is studied separately on the wrist, Section 4.4, with
     friction parameters identified from the hardware logs)
"""
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from qarm_sim.controllers import (FLController, PIDGController,  # noqa: E402
                                  SLController, pole_matched_pid)
from qarm_sim.model import NOMINAL_PARAMS, derived, perturbed_params  # noqa: E402
from qarm_sim.simulate import (dist_constant, ref_sinusoidal,  # noqa: E402
                               simulate_sampled, tracking_metrics)

OUT = os.path.join(os.path.dirname(__file__), "..", "..", "results", "comparison")
os.makedirs(OUT, exist_ok=True)

P_NOM = derived(NOMINAL_PARAMS)
A1, A2, AD, K1D, K2D = 10.5, 0.2, 5.0, 1.0, 2.0
T, TS, T_SS = 40.0, 0.002, 10.0
NOISE = dict(q_quant=2 * np.pi / 4096, q_std=0.0, qd_std=0.01)
TAU_MAX = np.array([2.65, 5.3, 2.65, 1.0])          # continuous ratings
N_MC = 30
CTRLS = ["PID+G", "FL", "FL+DE", "Proposed"]


def make_ctrl(name):
    if name == "PID+G":
        Kp, Ki, Kd = pole_matched_pid(P_NOM, AD, K1D, K2D)
        return PIDGController(P_NOM, Kp, Ki, Kd, int_limit=1.0)
    if name == "FL":
        return FLController(P_NOM, 1.0, 1.0)
    if name == "FL+DE":
        return FLController(P_NOM, 1.0, 1.0, with_dhat=True, alpha_d=AD, k1d=K1D, k2d=K2D)
    if name == "Proposed":
        return SLController(P_NOM, A1, A2, AD, K1D, K2D)
    raise ValueError(name)


def run_one(args):
    scen, ctrl_name, seed = args
    rng = np.random.default_rng(1000 + seed)
    p_true = P_NOM
    kw = {}
    if scen in ("B", "D"):
        p_true = derived(perturbed_params(np.random.default_rng(seed), rel=0.2,
                                          payload_range=(0.0, 0.3)))
    if scen in ("C", "D"):
        kw["noise"] = NOISE
    if scen == "D":
        kw["tau_max"] = TAU_MAX
    c = make_ctrl(ctrl_name)
    L = simulate_sampled(c, ref_sinusoidal, dist_constant, p_true=p_true, T=T, Ts=TS,
                         nsub=1, rng=rng, log_every=5, **kw)
    if L["diverged"]:
        return dict(scen=scen, ctrl=ctrl_name, seed=seed, diverged=True)
    m = tracking_metrics(L["t"], L["ref"], L["q"], t_ss=T_SS)
    ss = L["t"] >= T_SS
    out = dict(scen=scen, ctrl=ctrl_name, seed=seed, diverged=False,
               rmse_ss=m["rmse_ss"].tolist(), rmse_ss_agg=float(m["rmse_ss_agg"]),
               nrmse_ss_agg=float(m["nrmse_ss_agg"]),
               maxe_ss=m["maxe_ss"].tolist(),
               tau_rms=np.sqrt(np.mean(L["tau"][ss] ** 2, axis=0)).tolist(),
               tau_peak=np.abs(L["tau"]).max(0).tolist(),
               sat_frac=L["sat"].mean(0).tolist())
    if hasattr(c, "k1") and ctrl_name == "Proposed":
        out["k1_range"] = [float(L["k1"].min()), float(L["k1"].max())]
        out["k2_range"] = [float(L["k2"].min()), float(L["k2"].max())]
    if seed == 0:
        np.savez_compressed(os.path.join(OUT, f"log_{scen}_{ctrl_name.replace('+', 'p')}.npz"),
                            t=L["t"], q=L["q"], ref=L["ref"], tau=L["tau"], k1=L["k1"],
                            k2=L["k2"], dhat=L["dhat"])
    return out


def main():
    jobs = []
    for scen in ("A", "C"):
        jobs += [(scen, c, 0) for c in CTRLS]
    for scen in ("B", "D"):
        jobs += [(scen, c, s) for c in CTRLS for s in range(N_MC)]
    t0 = time.time()
    with Pool(2) as pool:
        res = pool.map(run_one, jobs, chunksize=1)
    print(f"{len(jobs)} runs in {time.time() - t0:.0f} s")
    with open(os.path.join(OUT, "comparison_raw.json"), "w") as f:
        json.dump(res, f)
    Kp, Ki, Kd = pole_matched_pid(P_NOM, AD, K1D, K2D)
    with open(os.path.join(OUT, "pid_gains.json"), "w") as f:
        json.dump(dict(Kp=Kp.tolist(), Ki=Ki.tolist(), Kd=Kd.tolist()), f, indent=1)


if __name__ == "__main__":
    main()
