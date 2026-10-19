% MATLAB Function block "Wrist Controller/Subsystem/MATLAB Function" of
% Nonlinear_Control_Test_Program.slx, copied verbatim (state saved after Test 5:
% alpha1=2, alpha2=1, alphadh=1, k1d=400, k2d=210; see the paper for Tests 1-4).
% Inputs phi, phi_dot = [0 0 0 q4], [0 0 0 q4dot] (joints 1-3 fed as zero);
% r11=r12=r13=0, r14 = wrist command; r2*, r2*d from Simulink Derivative blocks.
% Outputs dk1, dk2, ddh1..4 are integrated by Integrator blocks (k1(0)=k2(0)=1,
% dhat(0)=0); u(4) goes to torque_to_current.m.
% NOTE: e2d below uses d = 0, i.e. the acceleration-free learning signal (Law B).

function [dk1, dk2, ddh1, ddh2, ddh3, ddh4, u] = f_qarm(phi, phi_dot, k1, k2, d_hat1, d_hat2, d_hat3, d_hat4, r11, r12, r13, r14, r21, r22, r23, r24, r21d, r22d, r23d, r24d)
% Reference trajectory
    % r11 = (1*pi/180)*cos(t);
    % r21 = -(1*pi/180)*sin(t);
    % r21d = -(1*pi/180)*cos(t);
    % 
    % r12 = (1*pi/180)*sin(t/4);
    % r22 = ((1/4)*pi/180)*cos(t/4);
    % r22d = -((1/16)*pi/180)*sin(t/4);
    % 
    % r13 = (45*pi/180)*sin(t/4);
    % r23 = ((45/4)*pi/180)*cos(t/4);
    % r23d = -((45/16)*pi/180)*sin(t/4);
    % 
    % r14 = (30*pi/180)*sin(t/4);
    % r24 = ((30/4)*pi/180)*cos(t/4);
    % r24d = -((30/16)*pi/180)*sin(t/4);


    
    r1=[r11 r12 r13 r14]';
    r2=[r21 r22 r23 r24]';
    r2d=[r21d r22d r23d r24d]';
    
%     d_x = rectangularPulse(0,1,t);    % actual disturbance 
%     d_y = rectangularPulse(0,1,t);    % actual disturbance
    d_1 = 0;    % actual disturbance 
    d_2 = 0;    % actual disturbance
    d_3= 0;
    d_4= 0;
    d=[d_1 d_2 d_3 d_4]';

%% States assignment
    th1 = phi(1);
    th2 = phi(2);
    th3= phi(3);
    th4=phi(4);
    th1d=phi_dot(1);
    th2d=phi_dot(2);
    th3d=phi_dot(3);
    th4d=phi_dot(4);
    % k1=s(9);
    % k2=s(10);
    dh1=d_hat1;
    dh2=d_hat2;
    dh3=d_hat3;
    dh4=d_hat4;
    x1=[th1 th2 th3 th4]';
    x2=[th1d th2d th3d th4d]';
    th=x1;
    thd=x2;
    dh=[dh1 dh2 dh3 dh4]';

    dk1 = k1;
    dk2 = k2;

%% control gains and learning rate of gradient descent
    alpha1 = 2;
    alpha2 = 1;
    %alpha3 = 0.05;
    %alpha4=0.05;
    alphadh= 1;
    
    k1d=400;
    k2d=210;
    %k3d=3;
    %k4d=4;

%% Definition of qarm parameters
g=9.81;
m2=0.4591;
m3=0.269;
m4=0.257;
mL=0.1;
L1=0.14;
L2=0.35;
L3=0.05;
L4=0.25;
L5=0.15;
lam1=L1;
lam2=sqrt(L2*L2+L3*L3);
lam3=L4+L5;
lamc2=0.1071;
lamc3=0.1561;
lamc4=0.0998;
I1A=1.489*10^(-3);
I2A=1.922*10^(-4);
I2L=1.610*10^(-3);
I3A=2.679*10^(-4);
I3L=2.069*10^(-3);
I4A=6.528*10^(-4);
I4L=1.120*10^(-3);

s2=sin(th(2));
s3=sin(th(3));
c2=cos(th(2));
c3=cos(th(3));
s23=sin(th(2)+th(3));
c23=cos(th(2)+th(3));

%-----------------------------------------------
%The Inertia matrix M
%-----------------------------------------------
M11=I1A+m2*(lam2-lamc2)^2*c2^2+I2A*s2^2+I2L*c2^2+m3*lam2^2*c2^2+m3*lamc3^2*s23^2-2*m3*lam2*lamc3*c2*s23+I3L*s23^2+I3A*c23^2+m4*(lam3-lamc4)^2*s23^2+m4*lam2^2*c2^2-2*m4*(lam3-lamc4)*lam2*c2*s23+I4L*s23^2+I4A*c23^2+mL*lam3^2*s23^2+mL*lam2^2*c2^2-2*mL*lam2*lam3*c2*s23;
M14=-I4A*c23;
M22=m2*(lam2-lamc2)^2+I2L+m3*lam2^2+m3*lamc3^2-2*m3*lam2*lamc3*s3+I3L+m4*(lam2+lam3-lamc4)^2+I4L+mL*lam2^2+mL*lam3^2-2*mL*lam2*lam3*s3;
M23=m3*lamc3^2-m3*lamc3*lam2*s3+I3L+m4*(lam3-lamc4-lam2*s3)*(lam3-lamc4)+I4L-mL*(lam3^2-lam2*lam3*s3);
M32=m3*lamc3^2-m3*lamc3*lam2*s3+I3L+m4*(lam3-lamc4-lam2*s3)*(lam3-lamc4)+I4L-mL*(lam3^2-lam2*lam3*s3);
M33=m3*lamc3^2+I3L+m4*(lamc4-lam3)^2+I4L+mL*lam3^2;
M41=-(I4A*c23);
M44=I4A;
M=[M11 0 0 M14;
    0 M22 M23 0;
    0 M32 M33 0;
    M41 0 0 M44];
%-----------------------------------------------


%-----------------------------------------------
%The Coriolis matrix B
%-----------------------------------------------
%Note: Parameter mapping changed from B(x,y) to resolve warnings
B11=-2*m2*(lam2-lamc2)^2*s2*c2+2*I2A*s2*c2-2*I2L*s2*c2-2*m3*lam2^2*s2*c2+2*m3*lamc3^2*s23*c23+2*m3*lam2*lamc3*s2*s23-2*m3*lam2*lamc3*c2*c23+2*I3L*s23*c23-2*I3A*s23*c23+2*m4*(lam3-lamc4)^2*s23*c23-2*m4*lam2^2*s2*c2+2*m4*(lam3-lamc4)*lam2*s2*s23-2*m4*(lam3-lamc4)*lam2*c2*c23+2*I4L*s23*c23-2*I4A*s23*c23+2*mL*lam3^2*s23*c23-2*mL*lam2^2*s2*c2+2*mL*lam2*lam3*s2*s23-2*mL*lam2*lam3*c2*c23;
B12=2*m3*lamc3^2*s23*c23-2*m3*lam2*lamc3*c2*c23+2*I3L*s23*c23-2*I3A*s23*c23+2*m4*(lam3-lamc4)^2*s23*c23-2*m4*(lam3-lamc4)*lam2*c2*c23+2*I4L*s23*c23-2*I4A*s23*c23+2*mL*lam3^2*s23*c23-2*mL*lam2*lam3*c2*c23;
B15=I4A*s23;
B16=I4A*s23;
B23=-(I4A*s23);
B24=-(2*m3*lam2*lamc3*c3+2*mL*lam2*lam3*c3);
B33=-(I4A*s23);
B41=I4A*s23;
B42=I4A*s23;
B = [B11 B12 0 0 B15 B16;
     0 0 B23 B24 0 0;
     0 0 B33 0 0 0;
     B41 B42 0 0 0 0];
%-----------------------------------------------


%-----------------------------------------------
%The Centrifugal matrix C
%-----------------------------------------------
C=zeros(4);
C(2,1)=-(-m2*(lam2-lamc2)^2*s2*c2+I2A*s2*c2-I2L*s2*c2-m3*lam2^2*s2*c2+m3*lamc3^2*s23*c23+m3*lam2*lamc3*s2*s23-m3*lam2*lamc3*c2*c23+I3L*s23*c23-I3A*s23*c23+m4*(lam3-lamc4)^2*s23*c23-m4*lam2^2*s2*c2+m4*(lam3-lamc4)*lam2*s2*s23-m4*(lam3-lamc4)*lam2*c2*c23+I4L*s23*c23-I4A*s23*c23+mL*lam3^2*s23*c23);
C(3,1)=-(m3*lamc3^2*s23*c23-m3*lam2*lamc3*c2*c23+I3L*s23*c23-I3A*s23*c23+m4*(lam3-lamc4)^2*s23*c23-m4*(lam3-lamc4)*lam2*c2*c23+I4L*s23*c23-I4A*s23*c23+mL*lam3^2*s23*c23-mL*lam2*lam3*c2*c23);
C(3,2)=m3*lam2*lamc3*c3+mL*lam2*lam3*c3;
%-----------------------------------------------


%-----------------------------------------------
%The Gravity vector G
%-----------------------------------------------
G1=0;
G2=-(m2*(lam2-lamc2)*c2+m3*(lam2*c2-lamc3*s23)+m4*(lam2*c2-(lam3-lamc4)*s23)+mL*(lam2*c2-lam3*s23));
G3=m3*lamc3*s23+m4*(lam3-lamc4)*s23+mL*lam3*s23;    
G4=0;
G=[G1;G2;G3;G4];
%-----------------------------------------------

%-----------------------------------------------
%The thdij vector for the Coriolis matrix
%-----------------------------------------------
thd12=thd(1)*thd(2);
thd13=thd(1)*thd(3);
thd14=thd(1)*thd(4);
thd23=thd(2)*thd(3);
thd24=thd(2)*thd(4);
thd34=thd(3)*thd(4);
thdij=[thd12;thd13;thd14;thd23;thd24;thd34];
%-----------------------------------------------


%-----------------------------------------------
%The thdis vector for the Centrifugal matrix
%-----------------------------------------------
thd1s=thd(1)*thd(1);
thd2s=thd(2)*thd(2);
thd3s=thd(3)*thd(3);
thd4s=thd(4)*thd(4);
thdis=[thd1s;thd2s;thd3s;thd4s];
%-----------------------------------------------

ff=-inv(M)*(B*thdij+C*thdis+g*G);
gg=inv(M);



%% Error definitions

e1=r1-x1;
e2=r2-x2;
e2d=-k1*e1-k2*e2-d+dh;

c=e2d+k2d*e2+k1d*e1;  
    
%% First order ODEs
dk1=alpha1*c'*e1;
dk2=alpha2*c'*e2;

ddh=-alphadh*c;

ddh1 = ddh(1);
ddh2 = ddh(2);
ddh3 = ddh(3);
ddh4 = ddh(4);

u=inv(gg)*(-ff+k1*e1+k2*e2+r2d-dh);

% pwm_s = pwm(u); 
    
end
 
