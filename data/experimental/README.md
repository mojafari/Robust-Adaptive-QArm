# Hardware experiment logs (QArm wrist, 500 Hz, 120 s each)

Files: `PID_Test{1..5}.mat` (vendor low-level wrist PID) and
`Nonlinear_Test{1..5}.mat` (proposed controller). Each file contains a matrix
`data` written by a Simulink *To File* block; one column per sample.

| Row | Signal | Unit |
|----|--------|------|
| 1 | time | s |
| 2 | wrist command r(t) = A sin(2 pi f t) | rad |
| 3-6 | measured joint angles q1..q4 (q4 = wrist) | rad |
| 7-10 | measured joint speeds | rad/s |
| 11 | wrist speed, 5 Hz 2nd-order filtered | rad/s |
| 12 | wrist acceleration, 5 Hz 2nd-order filtered | rad/s^2 |
| 13 | wrist current command (after rate limiter and saturation) | A |
| 14-17 | dhat_1..dhat_4 (proposed only; dhat_1..3 disabled = 0) | rad/s^2 |
| 18 | k1 (proposed only) | 1/s^2 |
| 19 | k2 (proposed only) | 1/s |

Test conditions: 1) 45 deg, 1/25 Hz; 2) 75 deg, 1/25 Hz; 3) 75 deg, 1/100 Hz;
4) 45 deg, 1/25 Hz, +32 g; 5) 75 deg, 1/25 Hz, +32 g (payload 7.6 cm from the axis).

Important notes for anyone re-analysing the data
* The vendor PID output is **zero for t <= 5 s** (its error and derivative paths
  are gated by `t > 5`). Compare controllers on t >= 10 s.
* Vendor current limits are [-0.05, +0.015] A. The proposed controller used the
  same limits in Tests 1-3 but a relaxed upper limit in Tests 4-5.
* Proposed-controller parameters differ between tests; they are recovered from
  the logs by `python/scripts/analyze_experiments.py` (`results/experiments/hw_params.json`).
* One trial per controller and condition.
