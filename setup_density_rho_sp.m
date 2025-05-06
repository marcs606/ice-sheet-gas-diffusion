%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function   [rho_ice, rho_s, h_firn, rho_rock, g ] = setup_density_rho_sp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

% I am repurposing the variables rho_s and h_firn
global rho_ice rho_s h_firn
%
%   set parameters to control density of firn
%  ------------------------------------------
     rho_ice   = 917;           % ice density  (kg m^-3)
%     rho_s     = 450;          % surface snow density   (kg m^-3)
%
%  set rock density
%  ----------------

     rho_rock  = 2400.;         % bedrock density   (kg m^-3) 
     %made consistent with setup_thermal_parameter
%

%  acceleration due to gravity (m s^-2)
%  ------------------------------------
     g  = 9.8;           
%
%
%  set thickness of firn column 
%   h_firn is firn compaction e-folding depth (m)
%   The thickness of air is   h_firn * (1 - rho_s/rho_ice)
%
%     h_firn    = 27;
%
%
%
    load sp_firn_density.mat
    h_firn = [0 h h(end)+10 100e3] ; %add deep points to avoid NaN, the rock mask will take care of the proper value later
    rho_s = [1000*rho_h(1) 1000*rho_h 917 917] ;
    
%     %make equivalent to ice
%     
%     rho_s = ones(size(rho_s))*917;
    
    a=1;
