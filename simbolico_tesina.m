clc
clear all

syms xc1 xc5 q2 q3 q4 l alpha beta real

eq1=xc1+l*cos(q2)+l*cos(alpha)+l*cos(beta);
eq2=-xc5+l*sin(q2)+l*sin(alpha)+l*sin(beta);
eq3=alpha-q2-q3;
eq4=beta-q2-q3-q4;

%%%%%%%%%%%%
% cos_23=cos(q2)*cos(q3)-sin(q2)*sin(q3);
% sin_23=sin(q2)*cos(q3)+cos(q2)*sin(q3);
% cos_234=cos_23*cos(q4)-sin_23*sin(q4);
% sin_234=sin_23*cos(q4)+cos_23*sin(q4);
% 
% eq1=x+l*cos(q2)+l*cos_23+l*cos_234==0;
% eq2=y+l*sin(q2)+l*sin_23+l*sin_234==0;
%%%%%%%%%%%%
% eq1=x+l*cos(q2)+l*(cos(q2)*cos(q3)-sin(q2)*sin(q3))+l*((cos(q2)*cos(q3)-sin(q2)*sin(q3))*cos(q4)-(sin(q2)*cos(q3)+cos(q2)*sin(q3))*sin(q4))==0;
% eq2=y+l*sin(q2)+l*(sin(q2)*cos(q3)+cos(q2)*sin(q3))+l*((sin(q2)*cos(q3)+cos(q2)*sin(q3))*cos(q4)+(cos(q2)*cos(q3)-sin(q2)*sin(q3))*sin(q4))==0;

sol=solve(eq1,eq2,eq3,eq4,alpha,beta,q3,q4);              
q3=simplify(sol.q3)
q4=simplify(sol.q4)


% q3_1=-q2 - 2*atan((2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*x -4*sin(q2)*l^3*y - 2*cos(2*q2)*l^2*x^2 + 4*sin(2*q2)*l^2*x*y + 2*cos(2*q2)*l^2*y^2 - 4*cos(q2)*l*x^3 + 4*sin(q2)*l*x^2*y - 4*cos(q2)*l*x*y^2 + 4*sin(q2)*l*y^3 - x^4 - 2*x^2*y^2 - y^4)^(1/2) - 2*l*y)/(2*l*x + l^2 + x^2 + y^2 + 2*l^2*cos(q2) + 2*l*x*cos(q2) - 2*l*y*sin(q2)));
% q3_2=2*atan(((3*l^4 + 4*cos(q2)*l^3*x - 4*sin(q2)*l^3*y -2*cos(2*q2)*l^2*x^2 + 4*sin(2*q2)*l^2*x*y + 2*cos(2*q2)*l^2*y^2 - 4*cos(q2)*l*x^3 + 4*sin(q2)*l*x^2*y - 4*cos(q2)*l*x*y^2 + 4*sin(q2)*l*y^3 - x^4 - 2*x^2*y^2 - y^4)^(1/2) - 2*l^2*sin(q2) + 2*l*y)/(2*l*x + l^2 + x^2 + y^2 + 2*l^2*cos(q2) + 2*l*x*cos(q2) - 2*l*y*sin(q2))) - q2;
%  
% q4_1=-2*atan(((3*l^4 + 4*cos(q2)*l^3*x -4*sin(q2)*l^3*y -2*cos(2*q2)*l^2*x^2 + 4*sin(2*q2)*l^2*x*y + 2*cos(2*q2)*l^2*y^2 - 4*cos(q2)*l*x^3 + 4*sin(q2)*l*x^2*y - 4*cos(q2)*l*x*y^2 + 4*sin(q2)*l*y^3 - x^4 - 2*x^2*y^2 - y^4)^(1/2) - 2*l^2*sin(q2) + 2*l*y)/(2*l*x + l^2 + x^2 + y^2 + 2*l^2*cos(q2) + 2*l*x*cos(q2) - 2*l*y*sin(q2)));
% q4_2=2*atan((2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*x - 4*sin(q2)*l^3*y- 2*cos(2*q2)*l^2*x^2 + 4*sin(2*q2)*l^2*x*y + 2*cos(2*q2)*l^2*y^2 - 4*cos(q2)*l*x^3 + 4*sin(q2)*l*x^2*y - 4*cos(q2)*l*x*y^2 + 4*sin(q2)*l*y^3 - x^4 - 2*x^2*y^2 - y^4)^(1/2) - 2*l*y)/(2*l*x + l^2 + x^2 + y^2 + 2*l^2*cos(q2) + 2*l*x*cos(q2) - 2*l*y*sin(q2)));
%  
% beta=q2 + pi/2 - 2*atan(((3*l^4 + 4*cos(q2)*l^3*x - 4*sin(q2)*l^3*y - 2*cos(2*q2)*l^2*x^2 + 4*sin(2*q2)*l^2*x*y + 2*cos(2*q2)*l^2*y^2 - 4*cos(q2)*l*x^3 + 4*sin(q2)*l*x^2*y - 4*cos(q2)*l*x*y^2 + 4*sin(q2)*l*y^3 - x^4 - 2*x^2*y^2 - y^4)^(1/2) - 2*l^2*sin(q2) + 2*l*y)/(2*l*x + l^2 + x^2 + y^2 + 2*l^2*cos(q2) + 2*l*x*cos(q2) - 2*l*y*sin(q2))) + 2*atan((2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*x - 4*sin(q2)*l^3*y - 2*cos(2*q2)*l^2*x^2 + 4*sin(2*q2)*l^2*x*y + 2*cos(2*q2)*l^2*y^2 - 4*cos(q2)*l*x^3 + 4*sin(q2)*l*x^2*y - 4*cos(q2)*l*x*y^2 + 4*sin(q2)*l*y^3 - x^4 - 2*x^2*y^2 - y^4)^(1/2) - 2*l*y)/(2*l*x + l^2 + x^2 + y^2 + 2*l^2*cos(q2) + 2*l*x*cos(q2) - 2*l*y*sin(q2)))

% l=1;
% x=0.99*l;
% y=1.99*l;
% q2=pi/2;


