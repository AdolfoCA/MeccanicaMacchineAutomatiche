function [dP]=kin_dir_vel(Q,dQ)
global l d
% xc1=Q(1);
% xc5=Q(2);
% q2=Q(3);
% dxc1=dQ(1);
% % dxc5=dQ(2);
% dalpha=dQ(3);
% 
% 
% [q3,~]=clc_q3_q4(Q);
% [dq3,~]=clc_dq3_dq4(Q,dQ);
% 
% dP(1)=dxc1-dq2*l*sin(q2)-(dq2+dq3)*l/2*sin(q2+q3)-dalpha*d*sin(alpha);
% dP(2)=+dq2*l*cos(q2)+(dq2+dq3)*l/2*cos(q2+q3)+dalpha*d*cos(alpha);
% dP(3)=dq2+dq3;

dP=jacobian_tesina(Q)*dQ;

end