clc
clear all
close all

syms tau1 tau2 tau5 R1 R5 H2 V2 H3 V3 H4 V4 H5 V5 M1 M5 Fx Fy Mf ...
ddxc1 ddxc5 ddq2 ddx2 ddy2 ddPx ddPy ddq3 ddx4 ddy4 ddalpha ddq4 ...
I3 I2 I5 ma m1 d g l dxc1 dxc5 dq2 dq3 dq4 xc1 xc5 q2 q3 q4 q5 Px Py dPx dPy alpha ep real

var1=load('C:\Users\giola\Desktop\tau1.mat');
var2=load('C:\Users\giola\Desktop\tau2.mat');
var5=load('C:\Users\giola\Desktop\tau5.mat');

tau1=var1.tau1;
tau5=var5.tau5;
tau2=var2.tau2;

%Calcolo matrice di massa M
M11=diff(tau1,ddxc1);
M12=diff(tau1,ddxc5);
M13=diff(tau1,ddq2);

M21=diff(tau5,ddxc1);
M22=diff(tau5,ddxc5);
M23=diff(tau5,ddq2);

M31=diff(tau2,ddxc1);
M32=diff(tau2,ddxc5);
M33=diff(tau2,ddq2);

M=[M11 M12 M13; M21 M22 M23; M31 M32 M33];

%Calcolo matrice di gravità G
G1=diff(tau1,g)
G2=diff(tau5,g)
G3=diff(tau2,g)

G=[G1;G2;G3]*g;

%Calcolo matrice di Coriolus V
V=[tau1;tau5;tau2]-M*[ddxc1;ddxc5;ddq2]-G;


save('C:\Users\giola\Desktop\M.mat','M');
save('C:\Users\giola\Desktop\G.mat','G');
save('C:\Users\giola\Desktop\V.mat','V');

