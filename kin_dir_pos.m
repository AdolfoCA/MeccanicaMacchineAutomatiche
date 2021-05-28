function [P] = kin_dir_pos(Q) %P[X Y alpha]
% kin_dir: determina la posiozne del terminale dalla configurazione delle
% varaibili di giunto 
global l d
xc1=Q(1);
%xc5=Q(2);
q2=Q(3);

% prendo le seconde soluzioni 
[q3,~]=clc_q3_q4(Q);

gamma = q2+q3;
%beta = q2+q3+q4;
alpha=gamma-pi/2;

P(1) = xc1+l*cos(q2)+l/2*cos(gamma)+d*cos(alpha);
P(2) = l*sin(q2)+l/2*sin(gamma)+d*sin(alpha);
P(3) = q2 + q3 - pi/2;

end

