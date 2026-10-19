%% ============================================================
%  Quantitative Experimental Assessment
%  ============================================================

clear; clc; close all;

%% Experimental test files

testFiles.proposed = { ...
    'Nonlinear_Test1.mat', ...
    'Nonlinear_Test2.mat', ...
    'Nonlinear_Test3.mat', ...
    'Nonlinear_Test4.mat', ...
    'Nonlinear_Test5.mat'};

testFiles.pid = { ...
    'PID_Test1.mat', ...
    'PID_Test2.mat', ...
    'PID_Test3.mat', ...
    'PID_Test4.mat', ...
    'PID_Test5.mat'};

nTests = numel(testFiles.proposed);

%% Analysis settings

t_ss = 5;       % beginning of steady-state evaluation [s]

% Storage
results = struct();

%% ============================================================
%  Process each experiment
%  ============================================================

for j = 1:nTests

    % ----------------------------------------------------------
    % Load proposed-controller data
    % ----------------------------------------------------------

    S = load(testFiles.proposed{j});
    data = S.data;

    t_prop = data(1,:);
    ref_prop = rad2deg(data(2,:));
    pos_prop = rad2deg(data(6,:));
    current_prop = data(13,:);

    % ----------------------------------------------------------
    % Load PID data
    % ----------------------------------------------------------

    S = load(testFiles.pid{j});
    data = S.data;

    t_pid = data(1,:);
    ref_pid = rad2deg(data(2,:));
    pos_pid = rad2deg(data(6,:));
    current_pid = data(13,:);

    % ----------------------------------------------------------
    % Check reference consistency
    % ----------------------------------------------------------

    % The two experiments should use the same reference trajectory.
    % Interpolate PID reference onto proposed time vector if needed.

    if length(t_prop) ~= length(t_pid) || ...
            max(abs(t_prop(:) - t_pid(:))) > 1e-9

        ref_pid_interp = interp1(t_pid, ref_pid, t_prop, ...
                                 'linear', 'extrap');

        t = t_prop;
        ref = ref_prop;

        % Check difference between references
        ref_difference = max(abs(ref - ref_pid_interp));

    else

        t = t_prop;
        ref = ref_prop;

        ref_difference = max(abs(ref_prop - ref_pid));

    end

    %% ========================================================
    %  Proposed controller
    %  ========================================================

    % Interpolate proposed position/current if necessary
    posP = pos_prop;
    currentP = current_prop;

    % Tracking error
    eP = ref - posP;

    % Overall metrics
    RMSE_P = sqrt(mean(eP.^2));
    MAE_P  = mean(abs(eP));
    MaxE_P = max(abs(eP));

    RMS_ref = sqrt(mean(ref.^2));

    NRMSE_P = 100 * RMSE_P / RMS_ref;

    % ----------------------------------------------------------
    % Steady-state
    % ----------------------------------------------------------

    idx_ss = t >= t_ss;

    eP_ss = eP(idx_ss);
    ref_ss = ref(idx_ss);
    currentP_ss = currentP(idx_ss);

    RMSE_P_ss = sqrt(mean(eP_ss.^2));
    MAE_P_ss  = mean(abs(eP_ss));
    MaxE_P_ss = max(abs(eP_ss));

    RMS_ref_ss = sqrt(mean(ref_ss.^2));

    NRMSE_P_ss = 100 * RMSE_P_ss / RMS_ref_ss;

    % ----------------------------------------------------------
    % Current metrics
    % ----------------------------------------------------------

    RMS_current_P = sqrt(mean(currentP.^2));
    Peak_current_P = max(abs(currentP));

    RMS_current_P_ss = sqrt(mean(currentP_ss.^2));
    Peak_current_P_ss = max(abs(currentP_ss));

    %% ========================================================
    %  PID controller
    %  ========================================================

    % Use PID time vector as the analysis vector if it differs.
    t2 = t_pid;

    % Interpolate proposed reference if necessary
    if length(t_pid) ~= length(t_prop) || ...
            max(abs(t_pid(:) - t_prop(:))) > 1e-9

        ref2 = interp1(t_prop, ref_prop, t_pid, ...
                       'linear', 'extrap');

    else
        ref2 = ref_pid;
    end

    ePID = ref2 - pos_pid;

    % Overall metrics
    RMSE_PID = sqrt(mean(ePID.^2));
    MAE_PID  = mean(abs(ePID));
    MaxE_PID = max(abs(ePID));

    RMS_ref_PID = sqrt(mean(ref2.^2));

    NRMSE_PID = 100 * RMSE_PID / RMS_ref_PID;

    % ----------------------------------------------------------
    % Steady-state
    % ----------------------------------------------------------

    idx_ss_PID = t2 >= t_ss;

    ePID_ss = ePID(idx_ss_PID);
    refPID_ss = ref2(idx_ss_PID);
    currentPID_ss = current_pid(idx_ss_PID);

    RMSE_PID_ss = sqrt(mean(ePID_ss.^2));
    MAE_PID_ss  = mean(abs(ePID_ss));
    MaxE_PID_ss = max(abs(ePID_ss));

    RMS_ref_PID_ss = sqrt(mean(refPID_ss.^2));

    NRMSE_PID_ss = 100 * RMSE_PID_ss / RMS_ref_PID_ss;

    % ----------------------------------------------------------
    % Current metrics
    % ----------------------------------------------------------

    RMS_current_PID = sqrt(mean(current_pid.^2));
    Peak_current_PID = max(abs(current_pid));

    RMS_current_PID_ss = sqrt(mean(currentPID_ss.^2));
    Peak_current_PID_ss = max(abs(currentPID_ss));

    %% ========================================================
    %  Controller-to-controller comparison
    %  ========================================================

    % RMSE reduction
    RMSE_change = 100 * ...
        (RMSE_PID - RMSE_P) / RMSE_PID;

    % Steady-state RMSE reduction
    RMSE_change_ss = 100 * ...
        (RMSE_PID_ss - RMSE_P_ss) / RMSE_PID_ss;

    % MAE reduction
    MAE_change = 100 * ...
        (MAE_PID - MAE_P) / MAE_PID;

    % Maximum-error reduction
    MaxE_change = 100 * ...
        (MaxE_PID - MaxE_P) / MaxE_PID;

    %% ========================================================
    %  Store results
    %  ========================================================

    results(j).test = j;

    results(j).RMSE_PID = RMSE_PID;
    results(j).RMSE_Proposed = RMSE_P;

    results(j).NRMSE_PID = NRMSE_PID;
    results(j).NRMSE_Proposed = NRMSE_P;

    results(j).RMSE_PID_ss = RMSE_PID_ss;
    results(j).RMSE_Proposed_ss = RMSE_P_ss;

    results(j).NRMSE_PID_ss = NRMSE_PID_ss;
    results(j).NRMSE_Proposed_ss = NRMSE_P_ss;

    results(j).MAE_PID = MAE_PID;
    results(j).MAE_Proposed = MAE_P;

    results(j).MAE_PID_ss = MAE_PID_ss;
    results(j).MAE_Proposed_ss = MAE_P_ss;

    results(j).MaxE_PID = MaxE_PID;
    results(j).MaxE_Proposed = MaxE_P;

    results(j).MaxE_PID_ss = MaxE_PID_ss;
    results(j).MaxE_Proposed_ss = MaxE_P_ss;

    results(j).RMS_current_PID = RMS_current_PID;
    results(j).RMS_current_Proposed = RMS_current_P;

    results(j).Peak_current_PID = Peak_current_PID;
    results(j).Peak_current_Proposed = Peak_current_P;

    results(j).RMS_current_PID_ss = RMS_current_PID_ss;
    results(j).RMS_current_Proposed_ss = RMS_current_P_ss;

    results(j).Peak_current_PID_ss = Peak_current_PID_ss;
    results(j).Peak_current_Proposed_ss = Peak_current_P_ss;

    results(j).RMSE_change = RMSE_change;
    results(j).RMSE_change_ss = RMSE_change_ss;

    results(j).MAE_change = MAE_change;
    results(j).MaxE_change = MaxE_change;

    results(j).reference_difference = ref_difference;

end

%% ============================================================
%  Display results
%  ============================================================

fprintf('\n============================================================\n');
fprintf('QUANTITATIVE EXPERIMENTAL ASSESSMENT\n');
fprintf('Steady-state begins at t = %.2f s\n', t_ss);
fprintf('============================================================\n');

for j = 1:nTests

    fprintf('\n-------------------- Test %d --------------------\n', j);

    fprintf('\nOverall RMSE [deg]\n');
    fprintf('PID       = %.6f\n', results(j).RMSE_PID);
    fprintf('Proposed  = %.6f\n', results(j).RMSE_Proposed);

    fprintf('\nOverall NRMSE [%%]\n');
    fprintf('PID       = %.6f\n', results(j).NRMSE_PID);
    fprintf('Proposed  = %.6f\n', results(j).NRMSE_Proposed);

    fprintf('\nSteady-state RMSE [deg]\n');
    fprintf('PID       = %.6f\n', results(j).RMSE_PID_ss);
    fprintf('Proposed  = %.6f\n', results(j).RMSE_Proposed_ss);

    fprintf('\nSteady-state NRMSE [%%]\n');
    fprintf('PID       = %.6f\n', results(j).NRMSE_PID_ss);
    fprintf('Proposed  = %.6f\n', results(j).NRMSE_Proposed_ss);

    fprintf('\nMAE [deg]\n');
    fprintf('PID       = %.6f\n', results(j).MAE_PID);
    fprintf('Proposed  = %.6f\n', results(j).MAE_Proposed);

    fprintf('\nMaximum absolute error [deg]\n');
    fprintf('PID       = %.6f\n', results(j).MaxE_PID);
    fprintf('Proposed  = %.6f\n', results(j).MaxE_Proposed);

    fprintf('\nRMS current [A]\n');
    fprintf('PID       = %.6f\n', results(j).RMS_current_PID);
    fprintf('Proposed  = %.6f\n', results(j).RMS_current_Proposed);

    fprintf('\nPeak absolute current [A]\n');
    fprintf('PID       = %.6f\n', results(j).Peak_current_PID);
    fprintf('Proposed  = %.6f\n', results(j).Peak_current_Proposed);

    fprintf('\nRMSE reduction [%%]\n');
    fprintf('Overall       = %.3f\n', results(j).RMSE_change);
    fprintf('Steady-state  = %.3f\n', results(j).RMSE_change_ss);

    fprintf('\nMAE reduction [%%] = %.3f\n', ...
        results(j).MAE_change);

    fprintf('Maximum-error reduction [%%] = %.3f\n', ...
        results(j).MaxE_change);

    fprintf('--------------------------------------------------\n');
end

%% ============================================================
%  Create summary table
%  ============================================================

Test = (1:nTests)';

PID_RMSE = [results.RMSE_PID]';
Proposed_RMSE = [results.RMSE_Proposed]';

PID_NRMSE = [results.NRMSE_PID]';
Proposed_NRMSE = [results.NRMSE_Proposed]';

PID_RMSE_ss = [results.RMSE_PID_ss]';
Proposed_RMSE_ss = [results.RMSE_Proposed_ss]';

PID_NRMSE_ss = [results.NRMSE_PID_ss]';
Proposed_NRMSE_ss = [results.NRMSE_Proposed_ss]';

RMSE_Reduction = [results.RMSE_change]';
RMSE_Reduction_ss = [results.RMSE_change_ss]';

SummaryTable = table(Test, ...
    PID_RMSE, Proposed_RMSE, ...
    PID_NRMSE, Proposed_NRMSE, ...
    PID_RMSE_ss, Proposed_RMSE_ss, ...
    PID_NRMSE_ss, Proposed_NRMSE_ss, ...
    RMSE_Reduction, RMSE_Reduction_ss);

disp(SummaryTable);

%% Optional:
% writetable(SummaryTable,'experimental_quantitative_results.csv');
