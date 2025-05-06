function [depth, age, w] = steady_depth_age(BDOT, MDOT, H, p)
% Gives age-depth profile

%BDOT = 0.02 ; This was set up w/o comment. Does this redefine the acc
%              value given in the temperature model script? keep commented
%MDOT = 0;%.002 ; Is this the melt rate? yes, keep commented
HDOT = 0 ;  % can keep as zero
%age = 2e6 ; % purpose of this variable? nothing
dz = 1/1e6 ;
zpos = linspace(0,1,1/dz) ;
%z = linspace(0,1000,1e6) ;
%H = 1000 ; This was set up w/o comment. does this redefine the ice depth
%           variable defined in the temperature model script?
h = p ; % parameter for vertical velocity profile. don't change for now, 
        % keep same as 'p' in temperature_model script

psi = (1-(h+2)/(h+1)*zpos + 1/(h+1)*zpos.^(h+2) ) ;
w = -( (BDOT - MDOT - HDOT) * psi ) - MDOT ; % vertical velocity    
tf = -w/BDOT ; % thinning function

ddepth = 1 ;
depth = 0:ddepth:H ; 
years = zeros(size(depth)) ;
for i = 2:length(depth)
    years(i) = interp1(zpos*H,1./tf,mean(depth(i-1:i)))/BDOT*ddepth ;
end

age = cumsum(years) ;

w = interp1(zpos*H,w,depth);