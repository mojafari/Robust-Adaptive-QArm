% MATLAB Function block "Wrist Controller/MATLAB Function" of
% Nonlinear_Control_Test_Program.slx, copied verbatim. Its output passes through
% a Rate Limiter (+-1.5 A/s) and a Saturation Dynamic block ([-15, 10] A in the
% saved Test-5 state; [-0.05, 0.015] A in Tests 1-3) before the HIL write.

function I = fcn(u)

% Torque = kt * current
kt = 0.005*353.5;       % Motor torque constant (N*m/A) Value from dynamic model simulink file

%kt = 0.3*272.5;       % Motor torque constant (N*m/A) Value from dynamic model simulink file

%kt = 0.3*353.5;       % Motor torque constant (N*m/A) Value from dynamic model simulink file - change gear ratio




I = u/kt;

