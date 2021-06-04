function [q3,q4]=clc_q3_q4(Q)
global l

xc1=Q(1);
xc5=Q(2);
q2=Q(3);

q3 = ...
 -q2 - 2*atan((2*l*xc5 + (-(2*l*xc1 + l^2 + l^2*tan(q2/2)^2 + xc1^2 + xc5^2 + xc1^2*tan(q2/2)^2 + xc5^2*tan(q2/2)^2 - 4*l*xc5*tan(q2/2) - 2*l*xc1*tan(q2/2)^2)*(2*l*xc1 - 3*l^2 - 3*l^2*tan(q2/2)^2 + xc1^2 + xc5^2 + xc1^2*tan(q2/2)^2 + xc5^2*tan(q2/2)^2 - 4*l*xc5*tan(q2/2) - 2*l*xc1*tan(q2/2)^2))^(1/2) - 4*l^2*tan(q2/2) + 2*l*xc5*tan(q2/2)^2)/(3*l^2*tan(q2/2)^2 - l^2 + xc1^2 + xc5^2 + xc1^2*tan(q2/2)^2 + xc5^2*tan(q2/2)^2 - 4*l*xc5*tan(q2/2) - 4*l*xc1*tan(q2/2)^2) - (4*(l*xc5 - 2*l^2*tan(q2/2) + l*xc5*tan(q2/2)^2))/(3*l^2*tan(q2/2)^2 - l^2 + xc1^2 + xc5^2 + xc1^2*tan(q2/2)^2 + xc5^2*tan(q2/2)^2 - 4*l*xc5*tan(q2/2) - 4*l*xc1*tan(q2/2)^2));
q4 = ...
2*atan((2*l*xc5 + (-(2*l*xc1 + l^2 + l^2*tan(q2/2)^2 + xc1^2 + xc5^2 + xc1^2*tan(q2/2)^2 + xc5^2*tan(q2/2)^2 - 4*l*xc5*tan(q2/2) - 2*l*xc1*tan(q2/2)^2)*(2*l*xc1 - 3*l^2 - 3*l^2*tan(q2/2)^2 + xc1^2 + xc5^2 + xc1^2*tan(q2/2)^2 + xc5^2*tan(q2/2)^2 - 4*l*xc5*tan(q2/2) - 2*l*xc1*tan(q2/2)^2))^(1/2) - 4*l^2*tan(q2/2) + 2*l*xc5*tan(q2/2)^2)/(3*l^2*tan(q2/2)^2 - l^2 + xc1^2 + xc5^2 + xc1^2*tan(q2/2)^2 + xc5^2*tan(q2/2)^2 - 4*l*xc5*tan(q2/2) - 4*l*xc1*tan(q2/2)^2)) + 2*atan((2*l*xc5 + (-(2*l*xc1 + l^2 + l^2*tan(q2/2)^2 + xc1^2 + xc5^2 + xc1^2*tan(q2/2)^2 + xc5^2*tan(q2/2)^2 - 4*l*xc5*tan(q2/2) - 2*l*xc1*tan(q2/2)^2)*(2*l*xc1 - 3*l^2 - 3*l^2*tan(q2/2)^2 + xc1^2 + xc5^2 + xc1^2*tan(q2/2)^2 + xc5^2*tan(q2/2)^2 - 4*l*xc5*tan(q2/2) - 2*l*xc1*tan(q2/2)^2))^(1/2) - 4*l^2*tan(q2/2) + 2*l*xc5*tan(q2/2)^2)/(3*l^2*tan(q2/2)^2 - l^2 + xc1^2 + xc5^2 + xc1^2*tan(q2/2)^2 + xc5^2*tan(q2/2)^2 - 4*l*xc5*tan(q2/2) - 4*l*xc1*tan(q2/2)^2) - (4*(l*xc5 - 2*l^2*tan(q2/2) + l*xc5*tan(q2/2)^2))/(3*l^2*tan(q2/2)^2 - l^2 + xc1^2 + xc5^2 + xc1^2*tan(q2/2)^2 + xc5^2*tan(q2/2)^2 - 4*l*xc5*tan(q2/2) - 4*l*xc1*tan(q2/2)^2));
end