%% DINAMICA
clc
close all
clear all

global ma m1 d g l

syms tau1 tau2 tau5 R1 R5 H2 V2 H3 V3 H4 V4 H5 V5 M1 M5 Fx Fy Mf ...
     ddxc1 ddxc5 ddq2 ddx2 ddy2 ddPx ddPy ddq3 ddx4 ddy4 ddalpha ddq4 ...
     I3 I2 I5 ma m1 d g l dxc1 dxc5 dq2 dq3 dq4 xc1 xc5 q2 q3 q4 q5 Px Py dPx dPy alpha ep real

%Inerzia asta 
% I2 = ma*l^2/4 + 1/12 *ma*l^2;
% I5 = I2;
% I3 = I2 + 1/12 * ma*d^2 + ma*(d^2/4+l^2/4);

% Scrivo le derivate delle dPx dPy per poi risolvere le equazioni
% successive senza problemi :)
dPx = dxc1-dq2*l*sin(q2)-(dq2+dq3)*l/2*sin(q2+q3)-(dq2+dq3)*d*sin(q2+q3-pi/2);
dPy = +dq2*l*cos(q2)+(dq2+dq3)*l/2*cos(q2+q3)+(dq2+dq3)*d*cos(q2+q3-pi/2);

ddPx = diff(dPx,dxc1)*ddxc1+diff(dPx,dq2)*ddq2+diff(dPx,dq3)*ddq3+diff(dPx,q2)*dq2+diff(dPx,q3)*dq3;
ddPy = diff(dPy,dxc1)*ddxc1+diff(dPy,dq2)*ddq2+diff(dPy,dq3)*ddq3+diff(dPy,q2)*dq2+diff(dPy,q3)*dq3;
ddPx = simplify(ddPx);
ddPy = simplify(ddPy);


eq1 = tau1-H2-m1*ddxc1;
eq2 = R1-V2-m1*g;
%eq3 = -tau2-M1;

eq4 = H2+H3-ma*(ddxc1-l/2*(cos(q2)*dq2^2+sin(q2)*ddq2)); % l'ultimo terime è ddx2
eq5 = V2-ma*g-V3-ma*ddy2;
eq6 = tau2-l/2*ma*g*cos(q2)-l*V3*cos(q2)-l*H3*sin(q2)-I2*ddq2;

eq7 = H4-H3+Fx-2*ma*ddPx;
eq8 = V3-V4-2*ma*g+Fy-2*ma*ddPy;
eq9 = Mf+2*ma*g*( (l/2-d/4*tan(ep)) *cos(ep) )+V4*l*cos(ep)-H4*l*sin(ep)-I3*(ddq2+ddq3);

eq10 = -H4-H5-ma*ddx4;
eq11 = V4-V5-ma*g-ma*(ddxc5+(ddq2+ddq3+ddq4)*l/2*cos(q5)-(dq2+dq3+dq4)^2*l/2*sin(q5)); %l'ultimo termine è ddy4
eq12 = -ma*g*l/2*cos(q5)+H4*l*sin(q5)+l*V4*cos(q5)-I5*(ddq2+ddq3+ddq4);

eq13 = R5+H5;
eq14 = tau5-m1*g+V5-m1*ddxc5;
%eq15 = M5;

eq16 = ddalpha-ddq2-ddq3;
eq17 = q5-q4+pi-q2-q3;
eq18 = ep-pi+q2+q3;

sol2 = solve(eq1,eq2,eq4,eq5,eq6,eq7,eq8,eq9,eq10,eq11,eq12,eq13,eq14,eq16,eq17,eq18,...
      tau1,tau2,tau5,H2,V2,H3,V3,H4,V4,H5,V5,R1,R5,ddalpha,ep,q5);

  
% Devo esprimere ddq3 e ddq4 in funzione delle varaibili di giunto!!
tau1 = simplify(sol2.tau1)
tau2 = simplify(sol2.tau2)
tau5 = simplify(sol2.tau5)


% Matrice di massa 
M11 = simplify(diff(tau1,ddxc1));
M12 = simplify(diff(tau1,ddxc5));
M13 = simplify(diff(tau1,ddq2));

M21 = simplify(diff(tau5,ddxc1));
M22 = simplify(diff(tau5,ddxc5));
M23 = simplify(diff(tau5,ddq2));

M31 = simplify(diff(tau2,ddxc1));
M32 = simplify(diff(tau2,ddxc5));
M33 = simplify(diff(tau2,ddq2));

M = [M11 M12 M13; M21 M22 M23; M31 M32 M33];

% Matrice di gravità
G1 = simplify(diff(tau1,g));
G2 = simplify(diff(tau5,g));
G3 = simplify(diff(tau2,g));

G = [G1 G2 G3]';

% Matrice di Coriolis
V = simplify([tau1 tau5 tau2]' - M * [ddxc1 ddxc5 ddq2]' - G);







