clc
clear all

syms x y q2 q3 q4 l real

eq1=x+l*cos(q2)-l*cos(q2+q3)-l*cos(q4)==0;
eq2=y+l*sin(q4)-l*sin(q2)-l*sin(pi-q2-q3)==0;

sol=solve(eq1,eq2,q3,q4);
q3=simplify(sol.q3)
q4=simplify(sol.q4)
beta=simplify(q4(2)+pi-q2-q3(2))

% q3_1=-q2 - 2*atan((2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*x -4*sin(q2)*l^3*y - 2*cos(2*q2)*l^2*x^2 + 4*sin(2*q2)*l^2*x*y + 2*cos(2*q2)*l^2*y^2 - 4*cos(q2)*l*x^3 + 4*sin(q2)*l*x^2*y - 4*cos(q2)*l*x*y^2 + 4*sin(q2)*l*y^3 - x^4 - 2*x^2*y^2 - y^4)^(1/2) - 2*l*y)/(2*l*x + l^2 + x^2 + y^2 + 2*l^2*cos(q2) + 2*l*x*cos(q2) - 2*l*y*sin(q2)));
% q3_2=2*atan(((3*l^4 + 4*cos(q2)*l^3*x - 4*sin(q2)*l^3*y -2*cos(2*q2)*l^2*x^2 + 4*sin(2*q2)*l^2*x*y + 2*cos(2*q2)*l^2*y^2 - 4*cos(q2)*l*x^3 + 4*sin(q2)*l*x^2*y - 4*cos(q2)*l*x*y^2 + 4*sin(q2)*l*y^3 - x^4 - 2*x^2*y^2 - y^4)^(1/2) - 2*l^2*sin(q2) + 2*l*y)/(2*l*x + l^2 + x^2 + y^2 + 2*l^2*cos(q2) + 2*l*x*cos(q2) - 2*l*y*sin(q2))) - q2;
%  
% q4_1=-2*atan(((3*l^4 + 4*cos(q2)*l^3*x -4*sin(q2)*l^3*y -2*cos(2*q2)*l^2*x^2 + 4*sin(2*q2)*l^2*x*y + 2*cos(2*q2)*l^2*y^2 - 4*cos(q2)*l*x^3 + 4*sin(q2)*l*x^2*y - 4*cos(q2)*l*x*y^2 + 4*sin(q2)*l*y^3 - x^4 - 2*x^2*y^2 - y^4)^(1/2) - 2*l^2*sin(q2) + 2*l*y)/(2*l*x + l^2 + x^2 + y^2 + 2*l^2*cos(q2) + 2*l*x*cos(q2) - 2*l*y*sin(q2)));
% q4_2=2*atan((2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*x - 4*sin(q2)*l^3*y- 2*cos(2*q2)*l^2*x^2 + 4*sin(2*q2)*l^2*x*y + 2*cos(2*q2)*l^2*y^2 - 4*cos(q2)*l*x^3 + 4*sin(q2)*l*x^2*y - 4*cos(q2)*l*x*y^2 + 4*sin(q2)*l*y^3 - x^4 - 2*x^2*y^2 - y^4)^(1/2) - 2*l*y)/(2*l*x + l^2 + x^2 + y^2 + 2*l^2*cos(q2) + 2*l*x*cos(q2) - 2*l*y*sin(q2)));
%  
% beta=q2 + pi/2 - 2*atan(((3*l^4 + 4*cos(q2)*l^3*x - 4*sin(q2)*l^3*y - 2*cos(2*q2)*l^2*x^2 + 4*sin(2*q2)*l^2*x*y + 2*cos(2*q2)*l^2*y^2 - 4*cos(q2)*l*x^3 + 4*sin(q2)*l*x^2*y - 4*cos(q2)*l*x*y^2 + 4*sin(q2)*l*y^3 - x^4 - 2*x^2*y^2 - y^4)^(1/2) - 2*l^2*sin(q2) + 2*l*y)/(2*l*x + l^2 + x^2 + y^2 + 2*l^2*cos(q2) + 2*l*x*cos(q2) - 2*l*y*sin(q2))) + 2*atan((2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*x - 4*sin(q2)*l^3*y - 2*cos(2*q2)*l^2*x^2 + 4*sin(2*q2)*l^2*x*y + 2*cos(2*q2)*l^2*y^2 - 4*cos(q2)*l*x^3 + 4*sin(q2)*l*x^2*y - 4*cos(q2)*l*x*y^2 + 4*sin(q2)*l*y^3 - x^4 - 2*x^2*y^2 - y^4)^(1/2) - 2*l*y)/(2*l*x + l^2 + x^2 + y^2 + 2*l^2*cos(q2) + 2*l*x*cos(q2) - 2*l*y*sin(q2)))
