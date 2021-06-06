close all
clear all
clc

syms tau1 tau2 tau5 R1 R5 H2 V2 H3 V3 H4 V4 H5 V5 M1 M5 Fx Fy Mf ...
ddxc1 ddxc5 ddq2 ddx2 ddy2 ddPx ddPy ddq3 ddx4 ddy4 ddalpha ddq4 ...
I3 I2 I5 ma m1 d g l dxc1 dxc5 dq2 dq3 dq4 xc1 xc5 q2 q3 q4 q5 Px Py dPx dPy alpha ep real


%calcolo tau1 tau2 tau5
eq1 = tau1-H2-m1*ddxc1;
eq2 = R1-V2-m1*g;
%eq3 = -tau2-M1;

eq4 = H2+H3-ma*ddx2; 
eq5 = V2-ma*g-V3-ma*ddy2;
eq6 = tau2-l/2*ma*g*cos(q2)-l*V3*cos(q2)-l*H3*sin(q2)-I2*ddq2;

eq7 = H4-H3+Fx-2*ma*ddPx;
eq8 = V3-V4-2*ma*g+Fy-2*ma*ddPy;
eq9 = Mf+2*ma*g*( (l/2-d/4*tan(ep)) *cos(ep) )+V4*l*cos(ep)-H4*l*sin(ep)-I3*(ddq2+ddq3);

eq10 = -H4-H5+ma*ddx4;
eq11 = V4-V5-ma*g-ma*ddy4; 
eq12 = -ma*g*l/2*cos(q5)+H4*l*sin(q5)+l*V4*cos(q5)-I5*(ddq2+ddq3+ddq4);

eq13 = R5+H5;
eq14 = tau5-m1*g+V5-m1*ddxc5;
%eq15 = M5;

eq16 = ddalpha-ddq2-ddq3;
eq17 = q5-q4+pi-q2-q3;
eq18 = ep-pi+q2+q3;

sol2 = solve(eq1,eq2,eq4,eq5,eq6,eq7,eq8,eq9,eq10,eq11,eq12,eq13,eq14,eq16,eq17,eq18,...
      tau1,tau2,tau5,H2,V2,H3,V3,H4,V4,H5,V5,R1,R5,ddalpha,ep,q5);

tau1 = simplify(sol2.tau1);
tau2 = simplify(sol2.tau2);
tau5 = simplify(sol2.tau5);

% save('C:\Users\giola\Desktop\tau1_s.mat','tau1')
% save('C:\Users\giola\Desktop\tau2_s.mat','tau2')
% save('C:\Users\giola\Desktop\tau5_s.mat','tau5')

%calcolo ddPx ddPy 
% %q3
% q3=-q2 - 2*atan((2*l*xc5 + (-(2*l*xc1 + l^2 + l^2*tan(q2/2)^2 + xc1^2 + xc5^2 + xc1^2*tan(q2/2)^2 + xc5^2*tan(q2/2)^2 - 4*l*xc5*tan(q2/2) - 2*l*xc1*tan(q2/2)^2)*(2*l*xc1 - 3*l^2 - 3*l^2*tan(q2/2)^2 + xc1^2 + xc5^2 + xc1^2*tan(q2/2)^2 + xc5^2*tan(q2/2)^2 - 4*l*xc5*tan(q2/2) - 2*l*xc1*tan(q2/2)^2))^(1/2) - 4*l^2*tan(q2/2) + 2*l*xc5*tan(q2/2)^2)/(3*l^2*tan(q2/2)^2 - l^2 + xc1^2 + xc5^2 + xc1^2*tan(q2/2)^2 + xc5^2*tan(q2/2)^2 - 4*l*xc5*tan(q2/2) - 4*l*xc1*tan(q2/2)^2) - (4*(l*xc5 - 2*l^2*tan(q2/2) + l*xc5*tan(q2/2)^2))/(3*l^2*tan(q2/2)^2 - l^2 + xc1^2 + xc5^2 + xc1^2*tan(q2/2)^2 + xc5^2*tan(q2/2)^2 - 4*l*xc5*tan(q2/2) - 4*l*xc1*tan(q2/2)^2));
% 
% %q4
% q4=2*atan((2*l*xc5 + (-(2*l*xc1 + l^2 + l^2*tan(q2/2)^2 + xc1^2 + xc5^2 + xc1^2*tan(q2/2)^2 + xc5^2*tan(q2/2)^2 - 4*l*xc5*tan(q2/2) - 2*l*xc1*tan(q2/2)^2)*(2*l*xc1 - 3*l^2 - 3*l^2*tan(q2/2)^2 + xc1^2 + xc5^2 + xc1^2*tan(q2/2)^2 + xc5^2*tan(q2/2)^2 - 4*l*xc5*tan(q2/2) - 2*l*xc1*tan(q2/2)^2))^(1/2) - 4*l^2*tan(q2/2) + 2*l*xc5*tan(q2/2)^2)/(3*l^2*tan(q2/2)^2 - l^2 + xc1^2 + xc5^2 + xc1^2*tan(q2/2)^2 + xc5^2*tan(q2/2)^2 - 4*l*xc5*tan(q2/2) - 4*l*xc1*tan(q2/2)^2)) + 2*atan((2*l*xc5 + (-(2*l*xc1 + l^2 + l^2*tan(q2/2)^2 + xc1^2 + xc5^2 + xc1^2*tan(q2/2)^2 + xc5^2*tan(q2/2)^2 - 4*l*xc5*tan(q2/2) - 2*l*xc1*tan(q2/2)^2)*(2*l*xc1 - 3*l^2 - 3*l^2*tan(q2/2)^2 + xc1^2 + xc5^2 + xc1^2*tan(q2/2)^2 + xc5^2*tan(q2/2)^2 - 4*l*xc5*tan(q2/2) - 2*l*xc1*tan(q2/2)^2))^(1/2) - 4*l^2*tan(q2/2) + 2*l*xc5*tan(q2/2)^2)/(3*l^2*tan(q2/2)^2 - l^2 + xc1^2 + xc5^2 + xc1^2*tan(q2/2)^2 + xc5^2*tan(q2/2)^2 - 4*l*xc5*tan(q2/2) - 4*l*xc1*tan(q2/2)^2) - (4*(l*xc5 - 2*l^2*tan(q2/2) + l*xc5*tan(q2/2)^2))/(3*l^2*tan(q2/2)^2 - l^2 + xc1^2 + xc5^2 + xc1^2*tan(q2/2)^2 + xc5^2*tan(q2/2)^2 - 4*l*xc5*tan(q2/2) - 4*l*xc1*tan(q2/2)^2));
% 
% %dq3
% dq3=(2*dxc5*((2*l + (4*sin(q2)*l^3 - 4*sin(2*q2)*l^2*xc1 - 4*cos(2*q2)*l^2*xc5 - 4*sin(q2)*l*xc1^2 + 8*cos(q2)*l*xc1*xc5 - 12*sin(q2)*l*xc5^2 + 4*xc1^2*xc5 + 4*xc5^3)/(2*(3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2)))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2)) + ((2*xc5 - 2*l*sin(q2))*(2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2))/((2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5)^2/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2 + 1) - dq2*((2*(((- 4*sin(q2)*l^3*xc1 - 4*cos(q2)*l^3*xc5 + 4*sin(2*q2)*l^2*xc1^2 + 8*cos(2*q2)*l^2*xc1*xc5 - 4*sin(2*q2)*l^2*xc5^2 + 4*sin(q2)*l*xc1^3 + 4*cos(q2)*l*xc1^2*xc5 + 4*sin(q2)*l*xc1*xc5^2 + 4*cos(q2)*l*xc5^3)/(2*(3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2)) + 2*l^2*cos(q2))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2)) + ((2*l*xc5*cos(q2) - 2*l^2*sin(q2) + 2*l*xc1*sin(q2))*(2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2))/((2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5)^2/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2 + 1) + 1) + (2*dxc1*(((2*xc1 - 2*l + 2*l*cos(q2))*(2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2 + (- 4*cos(q2)*l^3 + 4*cos(2*q2)*l^2*xc1 - 4*sin(2*q2)*l^2*xc5 + 12*cos(q2)*l*xc1^2 - 8*sin(q2)*l*xc1*xc5 + 4*cos(q2)*l*xc5^2 + 4*xc1^3 + 4*xc1*xc5^2)/(2*(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))*(3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2))))/((2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5)^2/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2 + 1);
% 
% %dq4
% dq4=- dxc1*((2*(((2*xc1 - 2*l + 2*l*cos(q2))*((3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l^2*sin(q2) + 2*l*xc5))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2 + (2*(xc1 + l*cos(q2))*(- l^2 + 2*cos(q2)*l*xc1 - 2*sin(q2)*l*xc5 + xc1^2 + xc5^2))/((l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))*(3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2))))/(((3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l^2*sin(q2) + 2*l*xc5)^2/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2 + 1) + (2*(((2*xc1 - 2*l + 2*l*cos(q2))*(2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2 + (2*(xc1 + l*cos(q2))*(- l^2 + 2*cos(q2)*l*xc1 - 2*sin(q2)*l*xc5 + xc1^2 + xc5^2))/((l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))*(3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2))))/((2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5)^2/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2 + 1)) - dxc5*((2*((2*l + (2*sin(q2)*l^3 - 2*sin(2*q2)*l^2*xc1 - 2*cos(2*q2)*l^2*xc5 - 2*sin(q2)*l*xc1^2 + 4*cos(q2)*l*xc1*xc5 - 6*sin(q2)*l*xc5^2 + 2*xc1^2*xc5 + 2*xc5^3)/(3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2)) + ((2*xc5 - 2*l*sin(q2))*(2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2))/((2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5)^2/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2 + 1) - (2*((2*l - (2*sin(q2)*l^3 - 2*sin(2*q2)*l^2*xc1 - 2*cos(2*q2)*l^2*xc5 - 2*sin(q2)*l*xc1^2 + 4*cos(q2)*l*xc1*xc5 - 6*sin(q2)*l*xc5^2 + 2*xc1^2*xc5 + 2*xc5^3)/(3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2)) - ((2*xc5 - 2*l*sin(q2))*((3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l^2*sin(q2) + 2*l*xc5))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2))/(((3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l^2*sin(q2) + 2*l*xc5)^2/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2 + 1)) - dq2*((2*((2*l^2*cos(q2) - (2*l*(- sin(q2)*l^2*xc1 - cos(q2)*l^2*xc5 + sin(2*q2)*l*xc1^2 + 2*cos(2*q2)*l*xc1*xc5 - sin(2*q2)*l*xc5^2 + sin(q2)*xc1^3 + cos(q2)*xc1^2*xc5 + sin(q2)*xc1*xc5^2 + cos(q2)*xc5^3))/(3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2)) - (2*l*(xc5*cos(q2) - l*sin(q2) + xc1*sin(q2))*((3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l^2*sin(q2) + 2*l*xc5))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2))/(((3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l^2*sin(q2) + 2*l*xc5)^2/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2 + 1) - (2*((2*l^2*cos(q2) + (2*l*(- sin(q2)*l^2*xc1 - cos(q2)*l^2*xc5 + sin(2*q2)*l*xc1^2 + 2*cos(2*q2)*l*xc1*xc5 - sin(2*q2)*l*xc5^2 + sin(q2)*xc1^3 + cos(q2)*xc1^2*xc5 + sin(q2)*xc1*xc5^2 + cos(q2)*xc5^3))/(3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2)) + (2*l*(xc5*cos(q2) - l*sin(q2) + xc1*sin(q2))*(2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5))/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2))/((2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5)^2/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))^2 + 1));
% 

var1=load('C:\Users\giola\Desktop\q3');
var2=load('C:\Users\giola\Desktop\q4');
var3=load('C:\Users\giola\Desktop\dq3');
var4=load('C:\Users\giola\Desktop\dq4');

q3=var1.q3;
q4=var2.q4;
dq3=var3.dq3;
dq4=var4.dq4;

dPx = dxc1-dq2*l*sin(q2)-(dq2+dq3)*l/2*sin(q2+q3)-(dq2+dq3)*d*sin(q2+q3-pi/2);
dPy = +dq2*l*cos(q2)+(dq2+dq3)*l/2*cos(q2+q3)+(dq2+dq3)*d*cos(q2+q3-pi/2);

ddPx = diff(dPx,dxc1)*ddxc1+diff(dPx,xc1)*dxc1+diff(dPx,dq2)*ddq2+diff(dPx,q2)*dq2+diff(dPx,dxc5)*ddxc5+diff(dPx,xc5)*dxc5
ddPy = diff(dPy,dxc1)*ddxc1+diff(dPy,xc1)*dxc1+diff(dPy,dq2)*ddq2+diff(dPy,q2)*dq2+diff(dPy,dxc5)*ddxc5+diff(dPy,xc5)*dxc5

save('C:\Users\giola\Desktop\ddPx.mat','ddPx')
save('C:\Users\giola\Desktop\ddPy.mat','ddPy')
