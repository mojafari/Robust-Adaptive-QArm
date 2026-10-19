clear all; clc; close all;
%% Loads and plots results
load Nonlinear_Test1.mat;
t_sl = data(1,:);
wrist_cmd_sl = rad2deg(data(2,:)); % (rad)
phi_meas_sl = rad2deg(data(3:6,:)); % (rad)
phi_dot_sl = rad2deg(data(7:10,:)); % (rad/s)
wrist_spd_sl = rad2deg(data(11,:)); % (rad/s)
wrist_acc_sl = rad2deg(data(12,:)); % (rad/s^2)
wrist_current_sl = data(13,:); % (A)
dhat_sl = data(14:17,:);
k1_sl = data(18,:);
k2_sl = data(19,:);

load PID_Test1.mat;
t_pid = data(1,:);
wrist_cmd_pid = rad2deg(data(2,:)); % (rad)
phi_meas_pid = rad2deg(data(3:6,:)); % (rad)
phi_dot_pid = rad2deg(data(7:10,:)); % (rad/s)
wrist_spd_pid = rad2deg(data(11,:)); % (rad/s)
wrist_acc_pid = rad2deg(data(12,:)); % (rad/s^2)
wrist_current_pid = data(13,:); % (A)


% Figures
close all;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextInterpreter','latex');
green = [0, 150, 0]/255;
blue = [0, 0, 150]/255;
figure;
subplot(311); hold on;
plot(t_sl,phi_meas_sl(4,:),'color',green); % actual trajectory
plot(t_pid,phi_meas_pid(4,:),'color',blue); % actual trajectory
plot(t_sl,wrist_cmd_sl(1,:),'--r'); % desired trajectory
plot(t_pid,wrist_cmd_pid(1,:),':b'); % desired trajectory
legend({'Designed','PID','Desired 1','Desired 2'},'Orientation','horizontal')
% ylim([-1.5, 1.5]);
% xlim([0, t(end)]);
ylabel('Wrist Angle [deg]');
xticks([]);
yticklabels(strrep(yticklabels,'-','$-$'))
Layout; % Fancy layout for figure


subplot(312); hold on;
plot(t_sl,wrist_spd_sl(1,:),'color',green); % actual trajectory
plot(t_pid,wrist_spd_pid(1,:),'color',blue); % actual trajectory
% ylim([-1.5, 1.5]);
% xlim([0, t(end)]);
ylabel('Wrist Angular Speed [deg/s]');
xticks([]);
yticklabels(strrep(yticklabels,'-','$-$'))
Layout; % Fancy layout for figure

subplot(313); hold on;
plot(t_sl,wrist_acc_sl(1,:),'color',green); % actual trajectory
plot(t_pid,wrist_acc_pid(1,:),'color',blue); % actual trajectory
% xlim([0, t(end)]);

ylabel('Wrist Angular Acceleration [deg/$s^{2}$]');
xlabel('Time [s]');
yticklabels(strrep(yticklabels,'-','$-$'))
Layout;


figure;
subplot(211); hold on;
plot(t_sl,phi_meas_sl(4,:),'color',green); % actual trajectory
plot(t_pid,phi_meas_pid(4,:),'color',blue); % actual trajectory
plot(t_sl,wrist_cmd_sl(1,:),'--r'); % desired trajectory
legend({'Designed','PID','Desired'},'Orientation','horizontal')
% ylim([-1.5, 1.5]);
% xlim([0, t(end)]);
ylabel('Wrist Angle [deg]');
xticks([]);
yticklabels(strrep(yticklabels,'-','$-$'))
Layout; % Fancy layout for figure

subplot(212); hold on;
plot(t_sl,wrist_current_sl(1,:),'color',green); % actual trajectory
plot(t_pid,wrist_current_pid(1,:),'color',blue); % actual trajectory
% xlim([0, t(end)]);

ylabel('Applied Current to Wrist [A]');
xlabel('Time [s]');
yticklabels(strrep(yticklabels,'-','$-$'))
Layout;


figure; 
subplot(211); hold on;
plot(t_sl,wrist_current_sl(1,:),'r');

xticks([]);
yticklabels(strrep(yticklabels,'-','$-$'))
ylabel('Applied Current to Wrist [A]');
% xlim([0, t(end)]);

Layout; 
subplot(212); hold on;
plot(t_pid,wrist_current_pid(1,:),'r');

ylabel('Applied Current to Wrist [A]');
xlabel('Time [s]');
% xlim([0, t(end)]);
yticklabels(strrep(yticklabels,'-','$-$'))
Layout;


figure; 
subplot(211); hold on;
plot(t_sl,k1_sl,'b');

xticks([]);
yticklabels(strrep(yticklabels,'-','$-$'))
ylabel('$k_{1}$');
% xlim([0, t(end)]);

Layout; 
subplot(212); hold on;
plot(t_sl,k2_sl,'b');

ylabel('$k_{2}$');
xlabel('Time [s]');
% xlim([0, t(end)]);
yticklabels(strrep(yticklabels,'-','$-$'))
Layout;



figure; 
subplot(221); hold on;
plot(t_sl,dhat_sl(1,:),'b');

xticks([]);
yticklabels(strrep(yticklabels,'-','$-$'))
ylabel('$\hat{d}_{1}$');
% xlim([0, t(end)]);

Layout; 
subplot(222); hold on;
plot(t_sl,dhat_sl(2,:),'b');

xticks([]);
yticklabels(strrep(yticklabels,'-','$-$'))
ylabel('$\hat{d}_{2}$');
% xlim([0, t(end)]);

Layout;
subplot(223); hold on;
plot(t_sl,dhat_sl(3,:),'b');

xticks([]);
yticklabels(strrep(yticklabels,'-','$-$'))
ylabel('$\hat{d}_{3}$');
% xlim([0, t(end)]);

Layout;
subplot(224); hold on;
plot(t_sl,dhat_sl(4,:),'b');

ylabel('$\hat{d}_{4}$');
xlabel('Time [s]');
yticklabels(strrep(yticklabels,'-','$-$'))
% xlim([0, t(end)]);
Layout;




figure(10)
subplot(2,2,1); hold on;
plot(t_sl,phi_meas_sl(1,:));% actual trajectory
plot(t_sl,wrist_cmd_sl(1,:),'--'); % desired trajectory
plot(t_sl,phi_meas_sl(2,:));% actual trajectory
plot(t_sl,wrist_cmd_sl(1,:),'--'); % desired trajectory
plot(t_sl,phi_meas_sl(3,:));% actual trajectory
plot(t_sl,wrist_cmd_sl(1,:),'--'); % desired trajectory
plot(t_sl,phi_meas_sl(4,:));% actual trajectory
plot(t_sl,wrist_cmd_sl(1,:),'--'); % desired trajectory
lgd = legend({'Designed (base)','Desired (base)','Designed (shoulder)','Desired (shoulder)','Designed (elbow)','Desired (elbow)','Designed (wrist)','Desired (wrist)'},'Orientation','horizontal');
lgd.NumColumns = 2;
% xlim([0, t(end)]);
% ylim([3.5, 11]);
xlabel('Time [s]');
ylabel('$x_{1,2,3,4}$');
% xticks([]);
yticklabels(strrep(yticklabels,'-','$-$'))
Layout; % Fancy layout for figure

subplot(2,2,2); hold on;
plot(t_sl,wrist_current_sl(1,:));
plot(t_sl,wrist_current_sl(1,:));
plot(t_sl,wrist_current_sl(1,:));
plot(t_sl,wrist_current_sl(1,:));
lgd = legend({'$u_{1}$','$u_{2}$','$u_{3}$','$u_{4}$'},'Orientation','horizontal');
lgd.NumColumns = 4;
%xlim([0, t(end)]);
% ylim([3.5, 11]);
xlabel('Time [s]');
ylabel('$u_{1,2,3,4}$');
% xticks([]);
yticklabels(strrep(yticklabels,'-','$-$'))
Layout;
subplot(2,2,3); hold on;
plot(t_sl,k1_sl,'m')
plot(t_sl,k2_sl,'b')
lgd = legend({'$k_1$','$k_2$'},'Orientation','horizontal');
lgd.NumColumns = 2;
%xlim([0, t(end)]);
% ylim([3.5, 11]);
xlabel('Time [s]');
ylabel('$k_{1,2}$');
% xticks([]);
yticklabels(strrep(yticklabels,'-','$-$'))
Layout;
subplot(2,2,4); hold on;
plot(t_sl,dhat_sl(1,:));
plot(t_sl,dhat_sl(2,:));
plot(t_sl,dhat_sl(3,:));
plot(t_sl,dhat_sl(4,:));
lgd = legend({'$\hat{d}_{1}$','$\hat{d}_{2}$','$\hat{d}_{3}$','$\hat{d}_{4}$'},'Orientation','horizontal');
lgd.NumColumns = 4;
%xlim([0, t(end)]);
% ylim([3.5, 11]);
xlabel('Time [s]');
ylabel('$\hat{d}_{1,2,3,4}$');
% xticks([]);
yticklabels(strrep(yticklabels,'-','$-$'))
Layout;