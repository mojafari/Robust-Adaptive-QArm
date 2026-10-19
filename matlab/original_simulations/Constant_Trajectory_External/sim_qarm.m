%% Simulate the learning control for the qarm

clc
clear variables
close all

%% Initial conditions (base, shoulder, elbow and wrist)

th10=0;
th20=0;
th30=0;
th40=0;
th1d0=0;
th2d0=0;
th3d0=0;
th4d0=0;
k10=1;
k20=1;
dh10=0;
dh20=0;
dh30=0;
dh40=0;
z0=[th10 th20 th30 th40 th1d0 th2d0 th3d0 th4d0 k10 k20 dh10 dh20 dh30 dh40];

% delete the torque log of a previous run (fopen(...,'a') appends)
if exist('exp_const_Ext.txt','file'), delete('exp_const_Ext.txt'); end

%% solve ODE using RK4

h = 0.05;      % 50ms
t0 = 0;
t = [t0];

Q = z0;

sim_t = 120;

while t(end)<sim_t
    t0 = t(end);
    x0 = (Q(end,:))';
    X = rk4('f_qarm',t0,x0,h);
    t = [t;t0+h];
    Q = [Q;X'];
end
%% extract states
th1 = rad2deg(Q(:,1));
th2= rad2deg(Q(:,2));
th3 = rad2deg(Q(:,3));
th4 = rad2deg(Q(:,4));
k1 = Q(:,9);
k2 = Q(:,10);
dh1 = Q(:,11);
dh2 = Q(:,12);
dh3 = Q(:,13);
dh4 = Q(:,14);

% Reference trajectory
r11 = 45.*ones(size(t));
r12 = 90.*ones(size(t));
r13 = -90.*ones(size(t));
r14 = 30.*ones(size(t));



fileID = fopen('exp_const_Ext.txt','r');
formatSpec = '%f %f %f %f';
sizeA = [5 Inf];
A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
t_c = A(1,:);
u_1 = A(2,:);
u_2 = A(3,:);
u_3 = A(4,:);
u_4 = A(5,:);

% %% Plotting
close all;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextInterpreter','latex');
green = [0, 150, 0]/255;
blue = [0, 0, 150]/255;
figure(1)
subplot(4,1,1); hold on;
plot(t,th1,'color',green);% actual trajectory
plot(t,r11,'--r'); % desired trajectory
legend({'Designed','Desired'},'Orientation','horizontal')
xlim([0, t(end)]);
% ylim([-1.1, 1.75]);
ylabel({'Base Angle', '[deg]'});
xticks([]);
yticklabels(strrep(yticklabels,'-','$-$'))
Layout; % Fancy layout for figure

subplot(4,1,2); hold on;
plot(t,th2,'color',green);% actual trajectory
plot(t,r12,'--r'); % desired trajectory
% ylim([-1.3, 1.3]);
xlim([0, t(end)]);
ylabel({'Shoulder Angle', '[deg]'});
xticks([]);
yticklabels(strrep(yticklabels,'-','$-$'))
Layout; % Fancy layout for figure

subplot(4,1,3); hold on;
plot(t,th3,'color',green);% actual trajectory
plot(t,r13,'--r'); % desired trajectory
xlim([0, t(end)]);
ylabel({'Elbow Angle', '[deg]'});
xticks([]);
yticklabels(strrep(yticklabels,'-','$-$'))
Layout;

subplot(4,1,4); hold on;
plot(t,th4,'color',green);% actual trajectory
plot(t,r14,'--r'); % desired trajectory
xlim([0, t(end)]);
ylabel({'Wrist Angle', '[deg]'});
xlabel('Time [s]');
yticklabels(strrep(yticklabels,'-','$-$'))
Layout;



figure(2) 
subplot(411); hold on;
plot(t_c,u_1,'color',blue);

xticks([]);
yticklabels(strrep(yticklabels,'-','$-$'))
ylabel('$\tau_1$ [N.m]');
xlim([0, t_c(end)]);

Layout; 
subplot(412); hold on;
plot(t_c,u_2,'color',blue);

xticks([]);
yticklabels(strrep(yticklabels,'-','$-$'))
ylabel('$\tau_2$ [N.m]');
xlim([0, t_c(end)]);

Layout;
subplot(413); hold on;
plot(t_c,u_3,'color',blue);
xticks([]);
yticklabels(strrep(yticklabels,'-','$-$'))
ylabel('$\tau_3$ [N.m]');
xlim([0, t_c(end)]);
Layout;

subplot(414); hold on;
plot(t_c,u_4,'color',blue);

ylabel('$\tau_4$ [N.m]');
xlabel('Time [s]');
xlim([0, t_c(end)]);
yticklabels(strrep(yticklabels,'-','$-$'))
Layout;


figure(3)
subplot(2,1,1); hold on;
plot(t,k1,'m')

xticks([]);
ylabel('$k_{1}$');
xlim([0, t(end)]);
yticklabels(strrep(yticklabels,'-','$-$'))
Layout;


subplot(2,1,2); hold on;
plot(t,k2,'m')
ylabel('$k_{2}$');
xlabel('Time [s]');
xlim([0, t(end)]);
yticklabels(strrep(yticklabels,'-','$-$'))
Layout;


figure(4)
subplot(4,1,1); hold on;
plot(t,dh1,'c');

xticks([]);
ylabel('$\hat{d}_{1}$');
xlim([0, t(end)]);
yticklabels(strrep(yticklabels,'-','$-$'))
Layout;

subplot(4,1,2); hold on;
plot(t,dh2,'c');

xticks([]);
ylabel('$\hat{d}_{2}$');
xlim([0, t(end)]);
yticklabels(strrep(yticklabels,'-','$-$'))
Layout;

subplot(4,1,3); hold on;
plot(t,dh3,'c');
xticks([]);
ylabel('$\hat{d}_{3}$');
xlim([0, t(end)]);
yticklabels(strrep(yticklabels,'-','$-$'))
Layout;

subplot(4,1,4); hold on;
plot(t,dh4,'c');

ylabel('$\hat{d}_{4}$');
xlabel('Time [s]');
xlim([0, t(end)]);
yticklabels(strrep(yticklabels,'-','$-$'))
Layout;




% figure(10)
% subplot(2,2,1); hold on;
% plot(t,th1*180/pi);% actual trajectory
% plot(t,r11,'--'); % desired trajectory
% plot(t,th2*180/pi);% actual trajectory
% plot(t,r12,'--'); % desired trajectory
% plot(t,th3*180/pi);% actual trajectory
% plot(t,r13,'--'); % desired trajectory
% plot(t,th4*180/pi);% actual trajectory
% plot(t,r14,'--'); % desired trajectory
% lgd = legend({'Designed (base)','Desired (base)','Designed (shoulder)','Desired (shoulder)','Designed (elbow)','Desired (elbow)','Designed (wrist)','Desired (wrist)'},'Orientation','horizontal');
% lgd.NumColumns = 2;
% % xlim([0, t(end)]);
% % ylim([3.5, 11]);
% xlabel('Time [s]');
% ylabel('$\theta_{1,2,3,4}$ [deg]');
% % xticks([]);
% yticklabels(strrep(yticklabels,'-','$-$'))
% Layout; % Fancy layout for figure
% 
% subplot(2,2,2); hold on;
% plot(t_c,u_1);
% plot(t_c,u_2);
% plot(t_c,u_3);
% plot(t_c,u_4);
% lgd = legend({'$u_{1}$','$u_{2}$','$u_{3}$','$u_{4}$'},'Orientation','horizontal');
% lgd.NumColumns = 4;
% %xlim([0, t(end)]);
% % ylim([3.5, 11]);
% xlabel('Time [s]');
% ylabel('$\tau_{1,2,3,4}$ [N.m]');
% % xticks([]);
% yticklabels(strrep(yticklabels,'-','$-$'))
% Layout;
% subplot(2,2,3); hold on;
% plot(t,k1,'m')
% plot(t,k2,'b')
% lgd = legend({'$k_1$','$k_2$'},'Orientation','horizontal');
% lgd.NumColumns = 2;
% %xlim([0, t(end)]);
% % ylim([3.5, 11]);
% xlabel('Time [s]');
% ylabel('$k_{1,2}$');
% % xticks([]);
% yticklabels(strrep(yticklabels,'-','$-$'))
% Layout;
% subplot(2,2,4); hold on;
% plot(t,dh1);
% plot(t,dh2);
% plot(t,dh3);
% plot(t,dh4);
% lgd = legend({'$\hat{d}_{1}$','$\hat{d}_{2}$','$\hat{d}_{3}$','$\hat{d}_{4}$'},'Orientation','horizontal');
% lgd.NumColumns = 4;
% %xlim([0, t(end)]);
% % ylim([3.5, 11]);
% xlabel('Time [s]');
% ylabel('$\hat{d}_{1,2,3,4}$');
% % xticks([]);
% yticklabels(strrep(yticklabels,'-','$-$'))
% Layout;


%% Quantitative simulation assessment - Normalized RMSE

% Tracking errors [deg]
e1 = r11 - th1;
e2 = r12 - th2;
e3 = r13 - th3;
e4 = r14 - th4;

% ============================================================
% Overall NRMSE - entire simulation
% ============================================================

RMSE1 = sqrt(mean(e1.^2));
RMSE2 = sqrt(mean(e2.^2));
RMSE3 = sqrt(mean(e3.^2));
RMSE4 = sqrt(mean(e4.^2));

RMS_ref1 = sqrt(mean(r11.^2));
RMS_ref2 = sqrt(mean(r12.^2));
RMS_ref3 = sqrt(mean(r13.^2));
RMS_ref4 = sqrt(mean(r14.^2));

NRMSE1 = 100 * RMSE1 / RMS_ref1;
NRMSE2 = 100 * RMSE2 / RMS_ref2;
NRMSE3 = 100 * RMSE3 / RMS_ref3;
NRMSE4 = 100 * RMSE4 / RMS_ref4;

RMSE_all = sqrt(mean([e1.^2; e2.^2; e3.^2; e4.^2], 'all'));

RMS_ref_all = sqrt(mean([r11.^2; r12.^2; r13.^2; r14.^2], 'all'));

NRMSE_all = 100 * RMSE_all / RMS_ref_all;


% ============================================================
% Steady-state NRMSE
% ============================================================

t_ss = 5;  % beginning of steady-state evaluation [s]

idx_ss = t >= t_ss;

% Steady-state errors
e1_ss = e1(idx_ss);
e2_ss = e2(idx_ss);
e3_ss = e3(idx_ss);
e4_ss = e4(idx_ss);

% Steady-state references
r11_ss = r11(idx_ss);
r12_ss = r12(idx_ss);
r13_ss = r13(idx_ss);
r14_ss = r14(idx_ss);

% Steady-state RMSE
RMSE1_ss = sqrt(mean(e1_ss.^2));
RMSE2_ss = sqrt(mean(e2_ss.^2));
RMSE3_ss = sqrt(mean(e3_ss.^2));
RMSE4_ss = sqrt(mean(e4_ss.^2));

% Steady-state RMS reference
RMS_ref1_ss = sqrt(mean(r11_ss.^2));
RMS_ref2_ss = sqrt(mean(r12_ss.^2));
RMS_ref3_ss = sqrt(mean(r13_ss.^2));
RMS_ref4_ss = sqrt(mean(r14_ss.^2));

% Steady-state NRMSE [%]
NRMSE1_ss = 100 * RMSE1_ss / RMS_ref1_ss;
NRMSE2_ss = 100 * RMSE2_ss / RMS_ref2_ss;
NRMSE3_ss = 100 * RMSE3_ss / RMS_ref3_ss;
NRMSE4_ss = 100 * RMSE4_ss / RMS_ref4_ss;

% Aggregate steady-state NRMSE
RMSE_all_ss = sqrt(mean( ...
    [e1_ss.^2; e2_ss.^2; e3_ss.^2; e4_ss.^2], 'all'));

RMS_ref_all_ss = sqrt(mean( ...
    [r11_ss.^2; r12_ss.^2; r13_ss.^2; r14_ss.^2], 'all'));

NRMSE_all_ss = 100 * RMSE_all_ss / RMS_ref_all_ss;


% ============================================================
% Display results
% ============================================================

fprintf('\n========================================\n');
fprintf('Quantitative Simulation Assessment\n');
fprintf('========================================\n');
fprintf('\nOverall RMSE (0-%g s)\n', t(end));
fprintf('Joint 1 RMSE     = %.7f %%\n', RMSE1);
fprintf('Joint 2 RMSE     = %.7f %%\n', RMSE2);
fprintf('Joint 3 RMSE     = %.7f %%\n', RMSE3);
fprintf('Joint 4 RMSE     = %.7f %%\n', RMSE4);
fprintf('Aggregate RMSE   = %.7f %%\n', RMSE_all);
fprintf('========================================\n\n');
fprintf('\nSteady-State RMSE (%g-%g s)\n', t_ss, t(end));
fprintf('Joint 1 RMSE     = %.7f %%\n', RMSE1_ss);
fprintf('Joint 2 RMSE     = %.7f %%\n', RMSE2_ss);
fprintf('Joint 3 RMSE     = %.7f %%\n', RMSE3_ss);
fprintf('Joint 4 RMSE     = %.7f %%\n', RMSE4_ss);
fprintf('Aggregate RMSE   = %.7f %%\n', RMSE_all_ss);
fprintf('========================================\n\n');
fprintf('\nOverall NRMSE (0-%g s)\n', t(end));
fprintf('Joint 1 NRMSE     = %.7f %%\n', NRMSE1);
fprintf('Joint 2 NRMSE     = %.7f %%\n', NRMSE2);
fprintf('Joint 3 NRMSE     = %.7f %%\n', NRMSE3);
fprintf('Joint 4 NRMSE     = %.7f %%\n', NRMSE4);
fprintf('Aggregate NRMSE   = %.7f %%\n', NRMSE_all);
fprintf('========================================\n\n');
fprintf('\nSteady-State NRMSE (%g-%g s)\n', t_ss, t(end));
fprintf('Joint 1 NRMSE     = %.7f %%\n', NRMSE1_ss);
fprintf('Joint 2 NRMSE     = %.7f %%\n', NRMSE2_ss);
fprintf('Joint 3 NRMSE     = %.7f %%\n', NRMSE3_ss);
fprintf('Joint 4 NRMSE     = %.7f %%\n', NRMSE4_ss);
fprintf('Aggregate NRMSE   = %.7f %%\n', NRMSE_all_ss);

fprintf('========================================\n\n');
