clc
clear all

syms xc1 xc5 q2 q3 q4 l alpha beta dxc1 dxc5 dq2 real

global l

eq1=xc1+l*cos(q2)+l*cos(alpha)+l*cos(beta);
eq2=-xc5+l*sin(q2)+l*sin(alpha)+l*sin(beta);
eq3=alpha-q2-q3;
eq4=beta-q2-q3-q4;

sol=solve(eq1,eq2,eq3,eq4,alpha,beta,q3,q4);              
q3=simplify(sol.q3(2));
q4=simplify(sol.q4(2));

dq3=simplify(diff(q3,xc1,1)*dxc1+diff(q3,xc5,1)*dxc5+diff(q3,q2,1)*dq2);
dq4=simplify(diff(q4,xc1,1)*dxc1+diff(q4,xc5,1)*dxc5+diff(q4,q2,1)*dq2);


% l=1;
% xc1=1*l;
% xc5=1.5*l;
% q2=pi/4;
% 
% 
