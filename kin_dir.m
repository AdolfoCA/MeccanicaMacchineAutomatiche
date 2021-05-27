function [X,Y,gamma] = kin_dir(xc1,xc5,q2)
% kin_dir: determina la posiozne del terminale dalla configurazione delle
% varaibili di giunto 
global l 


% prendo le prime soluzioni 
q3 = ...
2*atan(((3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l^2*sin(q2) + 2*l*xc5)/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))) - q2;
 
 
q4 = ...
- 2*atan(((3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l^2*sin(q2) + 2*l*xc5)/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2))) - 2*atan((2*l^2*sin(q2) + (3*l^4 + 4*cos(q2)*l^3*xc1 - 4*sin(q2)*l^3*xc5 - 2*cos(2*q2)*l^2*xc1^2 + 4*sin(2*q2)*l^2*xc1*xc5 + 2*cos(2*q2)*l^2*xc5^2 - 4*cos(q2)*l*xc1^3 + 4*sin(q2)*l*xc1^2*xc5 - 4*cos(q2)*l*xc1*xc5^2 + 4*sin(q2)*l*xc5^3 - xc1^4 - 2*xc1^2*xc5^2 - xc5^4)^(1/2) - 2*l*xc5)/(l^2 - 2*l*xc1 + xc1^2 + xc5^2 - 2*l^2*cos(q2) + 2*l*xc1*cos(q2) - 2*l*xc5*sin(q2)));
 

alpha = q2+q3;
beta = q2+q3+q4;


X = xc1+l*cos(q2)+l*cos(alpha)+l*cos(beta)
Y = -xc5+l*sin(q2)+l*sin(alpha)+l*sin(beta)
gamma = q2 + q3 - pi/2;


end

