function [i_cmd, k1_dot, k2_dot, dhat_dot] = proposed_wrist_controller( ...
    r, r_dot, r_ddot, q, q_dot, q_ddot_f, k1, k2, dhat, prm)
%PROPOSED_WRIST_CONTROLLER  Learning-based controller, QArm wrist channel.
%
%   Reference implementation of the controller block used in the hardware
%   experiments of the paper (wrist joint, current mode). It is intended to be
%   placed in a MATLAB Function block of the QUARC/Simulink model, with three
%   Integrator (1/s) blocks closing the loops k1_dot -> k1, k2_dot -> k2 and
%   dhat_dot -> dhat (initial conditions k1(0)=k2(0)=1, dhat(0)=0). The
%   current command must then pass through the rate limiter and saturation of
%   the vendor model before the HIL Write block.
%
%   Compact, wrist-only equivalent of f_qarm_hardware.m + torque_to_current.m
%   (the code that ran on the QArm). It was checked against the logged signals (k1, k2, dhat,
%   current), which it reproduces to within 0.2 % (see
%   python/scripts/analyze_experiments.py). It replaces the 4-DOF f_qarm
%   function used on the hardware, which, with joints 1-3 fed as zero and
%   their learning signals therefore zero, reduces exactly to the wrist expressions
%   below (I = I4A*v/K_tau, R^2 = 1.000 against the logs).
%
%   Inputs
%     r, r_dot, r_ddot : wrist reference and its discrete derivatives (rad, rad/s, rad/s^2)
%     q, q_dot         : measured wrist angle and speed (rad, rad/s)
%     q_ddot_f         : filtered wrist acceleration (rad/s^2), used only by Law A
%                        (e.g. 'yd' output of the 2nd-order LPF of the vendor model)
%     k1, k2, dhat     : current values of the adaptive states
%     prm              : struct with fields
%                        alpha1, alpha2, alpha_d, k1d, k2d, law ('B' or 'A'),
%                        I4A (6.528e-4 kg m^2), K_tau (1.7675 N m/A),
%                        k_min (optional gain projection, [] to disable)
%   Hardware settings used in the paper (table of experimental conditions):
%     Tests 1-3: alpha1=2, alpha2=1, alpha_d=1, k1d=200, k2d=100, law 'B'
%     Test 4   : alpha1=4, alpha2=1, alpha_d=1, k1d=400, k2d=200, law 'B'
%     Test 5   : alpha1=2, alpha2=1, alpha_d=1, k1d=400, k2d=210, law 'B'

e1 = r - q;
e2 = r_dot - q_dot;

% acceleration command and torque (wrist row of M(q)[...] + C qdot + G;
% with joints 1-3 at rest M44 = I4A and the C/G terms of the wrist vanish)
v   = k1*e1 + k2*e2 + r_ddot - dhat;
tau = prm.I4A * v;
i_cmd = tau / prm.K_tau;

% learning signal
if prm.law == 'A'          % acceleration-based (requires filtered acceleration)
    c = (r_ddot - q_ddot_f) + prm.k2d*e2 + prm.k1d*e1;
else                       % acceleration-free (used in the experiments)
    c = (prm.k1d - k1)*e1 + (prm.k2d - k2)*e2 + dhat;
end

k1_dot   = prm.alpha1 * c * e1;
k2_dot   = prm.alpha2 * c * e2;
dhat_dot = -prm.alpha_d * c;

% optional projection (recommended): stop decreasing a gain at k_min
if isfield(prm,'k_min') && ~isempty(prm.k_min)
    if k1 <= prm.k_min && k1_dot < 0, k1_dot = 0; end
    if k2 <= prm.k_min && k2_dot < 0, k2_dot = 0; end
end
end
