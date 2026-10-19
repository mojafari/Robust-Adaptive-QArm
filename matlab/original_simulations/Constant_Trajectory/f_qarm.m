function dx = f_qarm(t,s)
% Reference trajectory
    r11 = pi/4;
    r21 = 0;
    r21d = 0;
    
    r12 = pi/2;
    r22 = 0;
    r22d = 0;
    
    r13 = -pi/2;
    r23 = 0;
    r23d = 0;

    r14 = pi/6;
    r24 = 0;
    r24d = 0;
    
    r1=[r11 r12 r13 r14]';
    r2=[r21 r22 r23 r24]';
    r2d=[r21d r22d r23d r24d]';
    
%     d_x = rectangularPulse(0,1,t);    % actual disturbance 
%     d_y = rectangularPulse(0,1,t);    % actual disturbance
    d_1 = 0.5;    % actual disturbance 
    d_2 = 0.5;    % actual disturbance
    d_3 = 0.5;
    d_4 = 0.5;
    d = [d_1 d_2 d_3 d_4]';

%% States assignment
    th1 = s(1);
    th2 = s(2);
    th3 = s(3);
    th4 = s(4);
    th1d = s(5);
    th2d = s(6);
    th3d = s(7);
    th4d = s(8);
    k1 = s(9);
    k2 = s(10);
    dh1 = s(11);
    dh2 = s(12);
    dh3 = s(13);
    dh4 = s(14);
    x1 = [th1 th2 th3 th4]';
    x2 = [th1d th2d th3d th4d]';
    th = x1;
    thd = x2;
    dh = [dh1 dh2 dh3 dh4]';

%% control gains and learning rate of gradient descent
    alpha1 = 10.5;
    alpha2 = 0.2;
    alphadh = 5;
    
    k1d=1;
    k2d=2;

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

u=inv(gg)*(-ff+k1*e1+k2*e2+r2d-dh);

fileID = fopen('exp_const.txt','a');
fprintf(fileID,'%f %f %f %f\n',[t;u(1);u(2);u(3);u(4)]);
fclose(fileID);

dx=[x2;ff+gg*u+d;dk1;dk2;ddh];
    
end
 