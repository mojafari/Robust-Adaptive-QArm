"""Parameter sensitivity, gain projection, acceleration-free variant and
computation time (new Sections 4.3-4.4 / response to Reviewer 1).
"""
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from qarm_sim.controllers import FLController, PIDGController, SLController, pole_matched_pid  # noqa: E402
from qarm_sim.model import NOMINAL_PARAMS, derived, dynamics_terms  # noqa: E402
from qarm_sim.simulate import (dist_constant, dist_pulses, ref_constant, ref_sinusoidal,  # noqa: E402
                               simulate_original, simulate_sampled, tracking_metrics)

OUT = os.path.join(os.path.dirname(__file__), "..", "..", "results", "sensitivity")
os.makedirs(OUT, exist_ok=True)
P_NOM = derived(NOMINAL_PARAMS)
BASE = dict(alpha1=10.5, alpha2=0.2, alpha_d=5.0, k1d=1.0, k2d=2.0)
SCALES = [0.25, 0.5, 1.0, 2.0, 4.0]


def rh_margin(prm, k1=1.0, k2=1.0):
    a2 = k2 + prm["alpha_d"]
    a1 = k1 + prm["alpha_d"] * prm["k2d"]
    a0 = prm["alpha_d"] * prm["k1d"]
    return a2 * a1 - a0


def sens_job(args):
    name, scale = args
    prm = dict(BASE)
    if name != "base":
        prm[name] = BASE[name] * scale
    c = SLController(P_NOM, prm["alpha1"], prm["alpha2"], prm["alpha_d"], prm["k1d"], prm["k2d"])
    L = simulate_sampled(c, ref_sinusoidal, dist_constant, T=40.0, Ts=0.002, nsub=1, log_every=5)
    if L["diverged"]:
        return dict(param=name, scale=scale, diverged=True)
    m = tracking_metrics(L["t"], L["ref"], L["q"], t_ss=10.0)
    ss = L["t"] >= 10.0
    return dict(param=name, scale=scale, diverged=False, rmse_ss=m["rmse_ss"].tolist(),
                rmse_ss_agg=float(m["rmse_ss_agg"]),
                tau_rms=float(np.sqrt(np.mean(L["tau"][ss] ** 2))),
                k1_range=[float(L["k1"].min()), float(L["k1"].max())],
                k2_range=[float(L["k2"].min()), float(L["k2"].max())],
                rh_margin=rh_margin(prm))


def accel_free_job(k1d):
    k2d = 2.0 * np.sqrt(k1d)
    c = SLController(P_NOM, BASE["alpha1"], BASE["alpha2"], BASE["alpha_d"], k1d, k2d, variant="model")
    L = simulate_sampled(c, ref_constant, dist_constant, T=60.0, Ts=0.002, nsub=1, log_every=5)
    e_end = (L["ref"][-1] - L["q"][-1])
    return dict(k1d=k1d, k2d=k2d, diverged=bool(L["diverged"]), e_final=e_end.tolist(),
                predicted=float(-0.5 / k1d))


def timing():
    """Mean wall-clock time per control-law evaluation (Python/NumPy, one core)."""
    rng = np.random.default_rng(0)
    q, qd = rng.normal(size=4), rng.normal(size=4)
    ref = (rng.normal(size=4), rng.normal(size=4), rng.normal(size=4))
    y = (q, qd, rng.normal(size=4))
    Kp, Ki, Kd = pole_matched_pid(P_NOM, 5.0, 1.0, 2.0)
    ctrls = {"PID+G": PIDGController(P_NOM, Kp, Ki, Kd),
             "FL": FLController(P_NOM, 1.0, 1.0),
             "Proposed": SLController(P_NOM, 10.5, 0.2, 5.0, 1.0, 2.0)}

    class RBFNN:
        """Representative RBF-NN adaptive term (Gaussian nodes on z=[e1,e2,q,qd],
        weight update W' = Gamma*phi*s^T) added to a PD law, N nodes."""

        def __init__(self, N):
            self.C = rng.normal(size=(N, 16))
            self.W = np.zeros((N, 4))

        def step(self, t, ref, y, Ts):
            e1, e2 = ref[0] - y[0], ref[1] - y[1]
            s = e2 + e1
            z = np.concatenate([e1, e2, y[0], y[1]])
            phi = np.exp(-np.sum((self.C - z) ** 2, axis=1) / 2.0)
            tau = 5 * s + self.W.T @ phi
            self.W += Ts * 10.0 * np.outer(phi, s)
            return tau

    for N in (64, 256, 1024):
        ctrls[f"RBFNN-{N}"] = RBFNN(N)
    def bench(fn, n=20000, reps=7):
        for _ in range(500):
            fn()
        best = np.inf
        for _ in range(reps):
            t0 = time.perf_counter()
            for _ in range(n):
                fn()
            best = min(best, (time.perf_counter() - t0) / n)
        return 1e6 * best

    out = {}
    for name, c in ctrls.items():
        out[name] = bench(lambda c=c: c.step(0.0, ref, y, 0.002))
        c.reset() if hasattr(c, "reset") else None
    out["model_terms_only"] = bench(lambda: dynamics_terms(q, qd, P_NOM))
    out["n_adaptive_states"] = {"PID+G": 4, "FL": 0, "Proposed": 6,
                                "RBFNN-64": 256, "RBFNN-256": 1024, "RBFNN-1024": 4096}
    return out


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "timing":
        res = json.load(open(os.path.join(OUT, "sensitivity.json")))
        res["timing_us"] = timing()
        json.dump(res, open(os.path.join(OUT, "sensitivity.json"), "w"), indent=1)
        print(res["timing_us"])
        return
    jobs = [("base", 1.0)] + [(k, s) for k in BASE for s in SCALES if s != 1.0]
    with Pool(2) as pool:
        sens = pool.map(sens_job, jobs, chunksize=1)
        af = pool.map(accel_free_job, [1.0, 10.0, 100.0], chunksize=1)
    proj = {}
    for floor in (None, 0.1):
        s = simulate_original(ref_constant, dist_pulses, gain_floor=floor)
        m = tracking_metrics(s["t"], s["ref"], s["q"], t_ss=5.0)
        proj[str(floor)] = dict(rmse_ss=m["rmse_ss"].tolist(), rmse_ss_agg=float(m["rmse_ss_agg"]),
                                k1_range=[float(s["k1"].min()), float(s["k1"].max())])
    res = dict(sensitivity=sens, accel_free=af, projection=proj, timing_us=timing())
    with open(os.path.join(OUT, "sensitivity.json"), "w") as f:
        json.dump(res, f, indent=1)
    for r in sens:
        print(r["param"], r["scale"], r.get("rmse_ss_agg"), r.get("k1_range"), r.get("rh_margin"))
    print(af)
    print(proj)
    print(res["timing_us"])


if __name__ == "__main__":
    main()
