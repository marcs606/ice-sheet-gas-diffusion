function [deriv] = calc_deriv_ftcs_avgs(C,dzbnds,D,dt)

r=D*dt./dzbnds./dzbnds;

deriv=r.*(C(3:end)-2*C(2:end-1)+C(1:end-2));

deriv=[0;deriv;0];


end

