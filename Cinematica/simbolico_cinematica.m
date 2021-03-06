clc
close all
clear all

syms xc1 xc5 q2 q3 q4 d l Px Py alpha beta gamma dxc1 dxc5 dq2 ddxc1 ddxc5 ddq2 real

global l d m1 ma g

%% PER TROVARE Q3 Q4 dQ3 dQ4 ddQ3 ddQ4
% eq1=xc1+l*cos(q2)+l*cos(gamma)+l*cos(beta);
% eq2=-xc5+l*sin(q2)+l*sin(gamma)+l*sin(beta);
% eq3=gamma-q2-q3;
% eq4=beta-q2-q3-q4;
% eq5=q3 - (- q2 - 2*atan((2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5)/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))));
% 
% sol=solve(eq1,eq2,eq3,eq4,gamma,beta,q3,q4);              
% q3=simplify(sol.q3(2));
% q4=simplify(sol.q4(2));
%
% dq3=simplify(diff(q3,xc1,1)*dxc1+diff(q3,xc5,1)*dxc5+diff(q3,q2,1)*dq2);
% dq4=simplify(diff(q4,xc1,1)*dxc1+diff(q4,xc5,1)*dxc5+diff(q4,q2,1)*dq2);
% 
% ddq3=(diff(dq3,xc1,1)*dxc1+diff(dq3,xc5,1)*dxc5+diff(dq3,q2,1)*dq2+diff(dq3,dxc1,1)*ddxc1+diff(dq3,dxc5,1)*ddxc5+diff(dq3,dq2,1)*ddq2);
% ddq4=(diff(dq4,xc1,1)*dxc1+diff(dq4,xc5,1)*dxc5+diff(dq4,q2,1)*dq2+diff(dq4,dxc1,1)*ddxc1+diff(dq4,dxc5,1)*ddxc5+diff(dq4,dq2,1)*ddq2);


%% PER TROVARE XC1 e Q2
% eq1=xc1+l*cos(q2)+l/2*cos(alpha+pi/2)+d*cos(alpha)-Px;
% eq2=l*sin(q2)+l/2*sin(alpha+pi/2)+d*sin(alpha)-Py;
% 
% sol_1=solve(eq1,eq2,xc1,q2);
% xc1_1=simplify(sol_1.xc1(1));
% xc1_2=simplify(sol_1.xc1(2));
% q2_1=simplify(sol_1.q2(1));
% q2_2=simplify(sol_1.q2(2));
% 
% 
%% PER TROVARE XC5
% eq1=xc1+l*cos(q2)+l*cos(q2+q3)+l*cos(q2+q3+q4);
% eq2=l*sin(q2)+l*sin(q2+q3)+l*sin(q2+q3+q4)-xc5;
% % 
% sol_2=solve(eq1,eq2,xc5,q4);
% xc5_1=simplify(sol_2.xc5(1))
% xc5_2=simplify(sol_2.xc5(2))


%% PER TROVARE LA JACOBIANA
% P(1) = xc1+l*cos(q2)+l/2*cos(q2+q3)+d*cos(alpha);
% P(2) = l*sin(q2)+l/2*sin(q2+q3)+d*sin(alpha);
% P(3) = q2 + q3 - pi/2;
% 
% J=simplify(jacobian([P(1) P(2) P(3)],[xc1 xc5 q2]))
% dJ=diff(J,xc1)*dxc1+diff(J,xc5)*dxc5+diff(J,q2)*dq2

