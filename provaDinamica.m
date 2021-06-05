%% DINAMICA
clc
close all
clear all

global ma m1 d g l

syms tau1 tau2 tau5 R1 R5 H2 V2 H3 V3 H4 V4 H5 V5 M1 M5 Fx Fy Mf ...
     ddxc1 ddxc5 ddq2 ddx2 ddy2 ddPx ddPy ddx4 ddy4 ddalpha  ...
     I3 I2 I5 ma m1 d g l dxc1 dxc5 dq2 dq3 dq4 xc1 xc5 q2 q3 q4 q5 Px Py dPx dPy alpha ep ddq3 ddq4 real


% Inerzia asta 
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

ddq3_perM= ...
(2*ddxc5*((2*l + (2*sin(q2)*l^3 - 2*sin(2*q2)*l^2*xc1 - 2*cos(2*q2)*l^2*xc5 - 2*sin(q2)*l*xc1^2 + 4*cos(q2)*l*xc1*xc5 - 6*sin(q2)*l*xc5^2 + 2*xc1^2*xc5 + 2*xc5^3)/(3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2)) + ((2*xc5 - 2*l*sin(q2))*(2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2))/((2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5)^2/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2 + 1) - ddq2*(((2*(2*l^2*cos(q2) + (2*l*(- sin(q2)*l^2*xc1 - cos(q2)*l^2*xc5 + sin(2*q2)*l*xc1^2 + 2*cos(2*q2)*l*xc1*xc5 - sin(2*q2)*l*xc5^2 + sin(q2)*xc1^3 + cos(q2)*xc1^2*xc5 + sin(q2)*xc1*xc5^2 + cos(q2)*xc5^3))/(3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2)))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2)) + (4*l*(xc5*cos(q2) - l*sin(q2) + xc1*sin(q2))*(2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2)/((2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5)^2/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2 + 1) + 1) + (2*ddxc1*(((2*xc1 - 2*l + 2*l*cos(q2))*(2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2 + (2*(xc1 + l*cos(q2))*(- l^2 + 2*cos(q2)*l*xc1 - 2*sin(q2)*l*xc5 + xc1^2 + xc5^2))/((l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))*(3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2))))/((2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5)^2/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2 + 1);
 
ddq4_perM=...
- ddxc1*(((2*(2*xc1 - 2*l + 2*l*cos(q2))*((3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l^2*sin(q2) + 2*l*xc5))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2 + (2*(2*xc1 + 2*l*cos(q2))*(- l^2 + 2*cos(q2)*l*xc1 - 2*sin(q2)*l*xc5 + xc1^2 + xc5^2))/((l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))*(3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2)))/(((3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l^2*sin(q2) + 2*l*xc5)^2/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2 + 1) + ((2*(2*xc1 - 2*l + 2*l*cos(q2))*(2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2 + (2*(2*xc1 + 2*l*cos(q2))*(- l^2 + 2*cos(q2)*l*xc1 - 2*sin(q2)*l*xc5 + xc1^2 + xc5^2))/((l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))*(3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2)))/((2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5)^2/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2 + 1)) - ddxc5*(((2*(2*l + (2*sin(q2)*l^3 - 2*sin(2*q2)*l^2*xc1 - 2*cos(2*q2)*l^2*xc5 - 2*sin(q2)*l*xc1^2 + 4*cos(q2)*l*xc1*xc5 - 6*sin(q2)*l*xc5^2 + 2*xc1^2*xc5 + 2*xc5^3)/(3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2)))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2)) + (2*(2*xc5 - 2*l*sin(q2))*(2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2)/((2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5)^2/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2 + 1) - ((2*(2*l - (2*sin(q2)*l^3 - 2*sin(2*q2)*l^2*xc1 - 2*cos(2*q2)*l^2*xc5 - 2*sin(q2)*l*xc1^2 + 4*cos(q2)*l*xc1*xc5 - 6*sin(q2)*l*xc5^2 + 2*xc1^2*xc5 + 2*xc5^3)/(3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2)))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2)) - (2*(2*xc5 - 2*l*sin(q2))*((3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l^2*sin(q2) + 2*l*xc5))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2)/(((3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l^2*sin(q2) + 2*l*xc5)^2/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2 + 1)) - ddq2*(((2*(2*l^2*cos(q2) - (2*l*(- sin(q2)*l^2*xc1 - cos(q2)*l^2*xc5 + sin(2*q2)*l*xc1^2 + 2*cos(2*q2)*l*xc1*xc5 - sin(2*q2)*l*xc5^2 + sin(q2)*xc1^3 + cos(q2)*xc1^2*xc5 + sin(q2)*xc1*xc5^2 + cos(q2)*xc5^3))/(3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2)))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2)) - (4*l*(xc5*cos(q2) - l*sin(q2) + xc1*sin(q2))*((3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l^2*sin(q2) + 2*l*xc5))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2)/(((3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l^2*sin(q2) + 2*l*xc5)^2/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2 + 1) - ((2*(2*l^2*cos(q2) + (2*l*(- sin(q2)*l^2*xc1 - cos(q2)*l^2*xc5 + sin(2*q2)*l*xc1^2 + 2*cos(2*q2)*l*xc1*xc5 - sin(2*q2)*l*xc5^2 + sin(q2)*xc1^3 + cos(q2)*xc1^2*xc5 + sin(q2)*xc1*xc5^2 + cos(q2)*xc5^3))/(3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2)))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2)) + (4*l*(xc5*cos(q2) - l*sin(q2) + xc1*sin(q2))*(2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2)/((2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5)^2/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2 + 1));
 
% Matrice di massa 
M11 = simplify( diff(tau1,ddxc1)+diff(tau1,ddq3)*diff(ddq3_perM,ddxc1)+diff(tau1,ddq4)*diff(ddq4_perM,ddxc1) );
M12 = simplify( diff(tau1,ddxc5)+diff(tau1,ddq3)*diff(ddq3_perM,ddxc5)+diff(tau1,ddq4)*diff(ddq4_perM,ddxc5) );
M13 = simplify( diff(tau1,ddq2)+diff(tau1,ddq3)*diff(ddq3_perM,ddq2)+diff(tau1,ddq4)*diff(ddq4_perM,ddq2) );

M21 = simplify( diff(tau5,ddxc1)+diff(tau5,ddq3)*diff(ddq3_perM,ddxc1)+diff(tau5,ddq4)*diff(ddq4_perM,ddxc1) );
M22 = simplify( diff(tau5,ddxc5)+diff(tau5,ddq3)*diff(ddq3_perM,ddxc5)+diff(tau5,ddq4)*diff(ddq4_perM,ddxc5) );
M23 = simplify( diff(tau5,ddq2)+diff(tau5,ddq3)*diff(ddq3_perM,ddq2)+diff(tau5,ddq4)*diff(ddq4_perM,ddq2) );

M31 = simplify( diff(tau2,ddxc1)+diff(tau2,ddq3)*diff(ddq3_perM,ddxc1)+diff(tau2,ddq4)*diff(ddq4_perM,ddxc1) );
M32 = simplify( diff(tau2,ddxc5)+diff(tau2,ddq3)*diff(ddq3_perM,ddxc5)+diff(tau2,ddq4)*diff(ddq4_perM,ddxc5) );
M33 = simplify( diff(tau2,ddq2)+diff(tau2,ddq3)*diff(ddq3_perM,ddq2)+diff(tau2,ddq4)*diff(ddq4_perM,ddq2) );


M = [M11 M12 M13; M21 M22 M23; M31 M32 M33];

% Matrice di gravità
G1 = simplify(diff(tau1,g));
G2 = simplify(diff(tau5,g));
G3 = simplify(diff(tau2,g));

G = [G1 G2 G3]';
syms ddq4_1 ddq3_1 real

ddq3 = ddq3_1+ddq3_perM;
ddq4 = ddq4_1+ddq4_perM;

 
tau1 =...
-(4*I3*ddq2*cos(q2 + q3 + q4) - 4*Mf*cos(q2 + q3 + q4) + 4*I3*ddq3*cos(q2 + q3 + q4) - 4*I5*ddq2*cos(q2 + q3) - 4*I5*ddq3*cos(q2 + q3) - 4*I5*ddq4*cos(q2 + q3) + 4*Fx*l*sin(q4) + 3*g*l*ma*cos(q4) + 3*g*l*ma*cos(2*q2 + 2*q3 + q4) - d*g*ma*sin(q4) + d*g*ma*sin(2*q2 + 2*q3 + q4) - 4*ddxc1*l*m1*sin(q4) - 12*ddxc1*l*ma*sin(q4) + 10*dq2^2*l^2*ma*cos(q2)*sin(q4) + 4*ddq2*l^2*ma*sin(q2 + q3)*sin(q4) + 4*ddq3*l^2*ma*sin(q2 + q3)*sin(q4) + 10*ddq2*l^2*ma*sin(q2)*sin(q4) + 4*dq2^2*l^2*ma*cos(q2 + q3)*sin(q4) + 4*dq3^2*l^2*ma*cos(q2 + q3)*sin(q4) - 8*d*ddq2*l*ma*cos(q2 + q3)*sin(q4) - 8*d*ddq3*l*ma*cos(q2 + q3)*sin(q4) + 8*dq2*dq3*l^2*ma*cos(q2 + q3)*sin(q4) + 8*d*dq2^2*l*ma*sin(q2 + q3)*sin(q4) + 8*d*dq3^2*l*ma*sin(q2 + q3)*sin(q4) + 16*d*dq2*dq3*l*ma*sin(q2 + q3)*sin(q4))/(4*l*sin(q4))
 
 
tau2 =... 
(2*Mf*cos(q3)*sin(q4) + 2*Mf*cos(q4)*sin(q3) + 2*I2*ddq2*sin(q4) + 2*I5*ddq2*sin(q3) + 2*I5*ddq3*sin(q3) + 2*I5*ddq4*sin(q3) - 2*I3*ddq2*cos(q3)*sin(q4) - 2*I3*ddq2*cos(q4)*sin(q3) - 2*I3*ddq3*cos(q3)*sin(q4) - 2*I3*ddq3*cos(q4)*sin(q3) - 2*Fy*l*cos(q2)*sin(q4) + 2*Fx*l*sin(q2)*sin(q4) + 4*ddq2*l^2*ma*sin(q4) - 2*dq2^2*l^2*ma*sin(q3)*sin(q4) - 2*dq3^2*l^2*ma*sin(q3)*sin(q4) - d*g*ma*cos(q2)*cos(q4) + 6*g*l*ma*cos(q2)*sin(q4) + 3*g*l*ma*cos(q4)*sin(q2) - 4*ddxc1*l*ma*sin(q2)*sin(q4) + 2*ddq2*l^2*ma*cos(q3)*sin(q4) + 2*ddq3*l^2*ma*cos(q3)*sin(q4) + 4*d*dq2^2*l*ma*cos(q3)*sin(q4) + 4*d*dq3^2*l*ma*cos(q3)*sin(q4) - 4*dq2*dq3*l^2*ma*sin(q3)*sin(q4) + d*g*ma*cos(q2)*cos(q3)^2*cos(q4) - 3*g*l*ma*cos(q2)*cos(q3)^2*sin(q4) - 3*g*l*ma*cos(q3)^2*cos(q4)*sin(q2) - d*g*ma*cos(q3)^2*sin(q2)*sin(q4) + 4*d*ddq2*l*ma*sin(q3)*sin(q4) + 4*d*ddq3*l*ma*sin(q3)*sin(q4) - 3*g*l*ma*cos(q2)*cos(q3)*cos(q4)*sin(q3) - d*g*ma*cos(q2)*cos(q3)*sin(q3)*sin(q4) - d*g*ma*cos(q3)*cos(q4)*sin(q2)*sin(q3) + 8*d*dq2*dq3*l*ma*cos(q3)*sin(q4) + 3*g*l*ma*cos(q3)*sin(q2)*sin(q3)*sin(q4))/(2*sin(q4))
 
 
tau5 =... 
ddxc5*m1 + g*m1 + g*ma + ma*(ddxc5 - (l*cos(q2 + q3 + q4)*(ddq2 + ddq3 + ddq4))/2 + (l*sin(q2 + q3 + q4)*(dq2 + dq3 + dq4)^2)/2) - (4*Mf*sin(q2 + q3 + q4) - 4*I3*ddq2*sin(q2 + q3 + q4) - 4*I3*ddq3*sin(q2 + q3 + q4) + 4*I5*ddq2*sin(q2 + q3) + 4*I5*ddq3*sin(q2 + q3) + 4*I5*ddq4*sin(q2 + q3) - d*g*ma*cos(q4) + d*g*ma*cos(2*q2 + 2*q3 + q4) - g*l*ma*sin(q4) - 3*g*l*ma*sin(2*q2 + 2*q3 + q4))/(4*l*sin(q4))
 
 
% Matrice di Coriolis
 V = ([tau1 tau5 tau2]' - M * [ddxc1 ddxc5 ddq2]' - G);
