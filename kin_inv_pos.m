function [Q]=kin_inv_pos(P)
global l d
Px=P(1);
Py=P(2);
alpha=P(3);

% Q(1)=...
%     Px + l*(1 - (l*cos(alpha) - 2*Py + 2*d*sin(alpha))^2/(4*l^2))^(1/2) - d*cos(alpha) + (l*sin(alpha))/2;
Q(1)=...
    Px - l*(1 - (l*cos(alpha) - 2*Py + 2*d*sin(alpha))^2/(4*l^2))^(1/2) - d*cos(alpha) + (l*sin(alpha))/2;

% Q(3)=...
%     pi + asin((l*cos(alpha) - 2*Py + 2*d*sin(alpha))/(2*l));
Q(3)=...
    -asin((l*cos(alpha) - 2*Py + 2*d*sin(alpha))/(2*l));


xc1=Q(1);
q2=Q(3);
q3=alpha-q2+pi/2;

% Q(2)=...
%     l*sin(q2 + q3) + l*sin(q2) - l*(1 - (xc1 + l*cos(q2 + q3) + l*cos(q2))^2/l^2)^(1/2);
Q(2)=...
    l*(sin(q2 + q3) + sin(q2) + (1 - (xc1 + l*cos(q2 + q3) + l*cos(q2))^2/l^2)^(1/2));
  
end