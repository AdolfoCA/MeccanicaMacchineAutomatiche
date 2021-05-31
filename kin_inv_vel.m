function [dQ]=kin_inv_vel(dP,Q)
global l d

dQ=jacobian_tesina(Q)\dP;

end

