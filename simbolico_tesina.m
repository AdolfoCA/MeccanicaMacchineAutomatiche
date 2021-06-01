clc
close all
clear all

syms xc1 xc5 q2 q3 q4 d l Px Py alpha beta gamma dxc1 dxc5 dq2 real

global l d m1 ma g

%PER TROVARE Q3 E Q4
eq1=xc1+l*cos(q2)+l*cos(gamma)+l*cos(beta);
eq2=-xc5+l*sin(q2)+l*sin(gamma)+l*sin(beta);
eq3=gamma-q2-q3;
eq4=beta-q2-q3-q4;
eq5=q3 - (- q2 - 2*atan((2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5)/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))));
%
sol=solve(eq1,eq2,eq3,eq4,gamma,beta,q3,q4);              
q3=simplify(sol.q3(2));
q4=simplify(sol.q4(2));
% %
% dq3=simplify(diff(q3,xc1,1)*dxc1+diff(q3,xc5,1)*dxc5+diff(q3,q2,1)*dq2);
% dq4=simplify(diff(q4,xc1,1)*dxc1+diff(q4,xc5,1)*dxc5+diff(q4,q2,1)*dq2);
% 
% 
% %PER TROVARE XC1 E Q2
% eq1=xc1+l*cos(q2)+l/2*cos(alpha+pi/2)+d*cos(alpha)-Px;
% eq2=l*sin(q2)+l/2*sin(alpha+pi/2)+d*sin(alpha)-Py;
% %
% sol_1=solve(eq1,eq2,xc1,q2);
% xc1_1=simplify(sol_1.xc1(1));
% xc1_2=simplify(sol_1.xc1(2));
% q2_1=simplify(sol_1.q2(1));
% q2_2=simplify(sol_1.q2(2));
% 
% 
% %PER TROVARE XC5
% eq1=xc1+l*cos(q2)+l*cos(q2+q3)+l*cos(q2+q3+q4);
% eq2=l*sin(q2)+l*sin(q2+q3)+l*sin(q2+q3+q4)-xc5;
% % 
% sol_2=solve(eq1,eq2,xc5,q4);
% xc5_1=simplify(sol_2.xc5(1))
% xc5_2=simplify(sol_2.xc5(2))


%PER TROVARE LA JJACOBIANA
% P(1) = xc1+l*cos(q2)+l/2*cos(q2+q3)+d*cos(alpha);
% P(2) = l*sin(q2)+l/2*sin(q2+q3)+d*sin(alpha);
% P(3) = q2 + q3 - pi/2;
% 
% J=simplify(jacobian([P(1) P(2) P(3)],[xc1 xc5 q2]))

%% DINAMICA
syms tau1 tau2 tau5 R1 R5 H2 V2 H3 V3 H4 V4 H5 V5 M1 M5 Fx Fy Mf ...
     ddxc1 ddxc5 ddq2 ddx2 ddy2 ddPx ddPy ddq3 ddx4 ddy4 ddalpha ddq4 ...
     I3 I2 I5 ma m1 d g l q5 ep real

%Inerzia asta 
I2 = ma*l^2/4 + 1/12 *ma*l^2;
I5 = I2;
I3 = I2 + 1/12 * ma*d^2 + ma*(d^2/4+l^2/4);

eq1 = tau1-H2-m1*ddxc1;
%eq2 = R1-V2-m1*g;
eq3 = -tau2-M1;

eq4 = H2+H3-tau1-ma*ddx2;
eq5 = V2-ma*g-V3-ma*ddy2;
eq6 = tau2-l/2*ma*g*cos(q2)-l*V3*cos(q2)-l*H3*sin(q2)-I2*ddq2;

eq7 = H4-H3+Fx-2*ma*ddPx;
eq8 = V3-V4-2*ma*g+Fy-2*ma*ddPy;
eq9 = Mf+2*ma*g*( (l/2-d/4*tan(ep)) *cos(ep) )+V4*l*cos(ep)-H4*l*sin(ep)-I3*(ddq2+ddq3);

eq10 = -H4+H5-ma*ddx4;
eq11 = V4+V5-ma*g-tau5-ma*ddy4;
eq12 = -ma*g*l/2*cos(q5)+H4*l*sin(q5)+l*V4*cos(q5)-I5*(ddq2+ddq3+ddq4);

%eq13 = R5-H5;
eq14 = tau5-m1*g-V5-m1*ddxc5;
%eq15 = M5;
eq16 = ddalpha-ddq2-ddq3;

eq17 = q5-q4+pi-q2-q3;
eq18 = ep-pi+q2+q3;

sol2 = solve(eq1,eq3,eq4,eq5,eq6,eq7,eq8,eq9,eq10,eq11,eq12,eq14,eq16,eq17,eq18,...
      ddPx,ddPy,Fx,Fy,Mf,I3,ddx2,ddy2,ddalpha,H3,H2,ddq4,ddq3,ep,q5);

ddPx = simplify(sol2.ddPx);
ddPy = simplify(sol2.ddPy);
ddalpha = simplify(sol2.ddalpha);







