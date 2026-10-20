# Learning-Based Adaptive Robust Control of a 4-DOF Robot Manipulator

Code and data for the paper

> M. Jafari, D. Cadigan, M. Reyhanoglu, *Learning-Based Adaptive Robust Control of a
> 4-DOF Robot Manipulator*, Robotics (MDPI), under review (manuscript robotics-4620744, revision 1).

The controller combines feedback linearization of the Quanser QArm dynamics with
gradient-based online adaptation of two feedback gains `$k_1$`, `k2` and a lumped
disturbance estimate `dhat`. All three are driven by the learning signal
`c = e2_dot + k2d*e2 + k1d*e1`, which measures the mismatch from a target error dynamics.
The repository covers the full 4-DOF simulations and the wrist-joint hardware experiments.

## Contents

The original simulations of the paper were run in MATLAB (`matlab/original_simulations/`); all
simulations added in the revision were run in Python (`python/`). The experiments were run in
MATLAB/Simulink/QUARC. The QArm parameters used by the model are listed in Table 1 of the paper
and in `python/qarm_sim/model.py`.

```
data/experimental/         raw logs of the 10 hardware runs (+ README with the row layout)
matlab/original_simulations/  MATLAB scripts of the original simulations (3 scenarios)
matlab/hardware/           plotting scripts of the experiments and a reference
                           implementation of the controller block used on the QArm
python/qarm_sim/           Python package: QArm model, controllers, simulators,
                           hardware-representative wrist model
python/scripts/            scripts that produce every number and figure of the paper
results/                   outputs of the scripts (JSON + PDF/PNG figures)
```

## Reproducing the paper (Python 3.10+, `pip install -r python/requirements.txt`)

| Paper item | Script | Runtime (1 core) |
|---|---|---|
| Sec. 4.1, Table 3 (original scenarios; reproduces the MATLAB results to all digits) | `python scripts/run_original.py` | ~5 s |
| Sec. 4.2, Table 4 (PID+G / FL / FL+DE / designed; Monte Carlo, noise, limits) | `python scripts/run_comparison.py` then `python scripts/summarize_comparison.py` | ~25 min (2 cores) |
| Sec. 4.3, Table 5 (sensitivity, projection, timing) | `python scripts/run_sensitivity.py` | ~3 min |
| Sec. 4.4, Tables 6-8 (hardware parameters; metrics over 0-120 s, t>=5 s, t>=10 s) | `python scripts/analyze_experiments.py` | ~20 s |
| Sec. 4.5, Table 9 (identified wrist model) | `python scripts/run_wrist_hw.py` | ~3 min |
| All figures of the revised manuscript | `python scripts/make_paper_figures.py` | ~30 s |

Run the scripts from the `python/` folder. Results are written to `results/`.

The MATLAB scripts in `matlab/original_simulations/<case>/sim_qarm.m` are the original
simulation code (fixed: they now delete the torque log of a previous run before
starting, and the "Wrist" axis labels are spelled correctly).

## Hardware set-up (summary)

* Quanser QArm, QUARC/Simulink, 500 Hz, fixed-step `ode4`. Base/shoulder/elbow in
  position mode at home; wrist (Dynamixel XM430-W350-R) in **current mode**.
* **Baseline:** vendor model `LowLevelControlWristCurrent.slx` (Quanser, not
  redistributed here). Its configuration: parallel PID in amperes, direction-dependent
  gains (r>0: Kp=0.010, Ki=0.030, Kd=0.007; r<=0: 0.007, 0.090, 0.012), derivative on a 5 Hz
  filtered speed, integrator dead band 0.005 rad, clamp +-0.5 rad s, reset at error zero crossings,
  output enabled only for t > 5 s, rate limiter +-1 A/s, saturation [-0.05, +0.015] A.
* **Proposed controller (code that ran on the QArm):** `matlab/hardware/f_qarm_hardware.m`
  and `matlab/hardware/torque_to_current.m`, copied verbatim from the authors' real-time
  model (only the authors' blocks are distributed, not the Quanser blocks). Joints 1-3 are
  fed to the controller as zero; reference derivatives come from Simulink Derivative blocks;
  the adaptive states are integrated with Integrator blocks (k1(0)=k2(0)=1, dhat(0)=0);
  I = tau / (0.005*353.5) A; output rate limiter +-1.5 A/s, then saturation.
  `matlab/hardware/proposed_wrist_controller.m` is a compact wrist-only equivalent
  (Law B as used; Law A and gain projection selectable).
  Per-test parameters are given in the paper and in `results/experiments/hw_params.json`
  (the saved model holds the Test-5 values).

## Citation

```bibtex
@article{jafari2026learning,
  title   = {Learning-Based Adaptive Robust Control of a 4-{DOF} Robot Manipulator},
  author  = {Jafari, Mohammad and Cadigan, Dylan and Reyhanoglu, Mahmut},
  journal = {Robotics},
  year    = {2026},
  note    = {Under review}
}
```

## License

MIT (see `LICENSE`). Quanser software and documentation are the property of Quanser Inc.
