   function [L, c_vals, k_vals, VD_vals, ZERO_C, lapse ] ...
       = setup_thermal_parameters;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function sets up
%   - thermal constitutive parameters for firn, ice, and rock.
%   - freezing-point parameters
%   - geothermal flux
%   - atmospheric lapse rate (needed when surface height changes) 
%========================================================================
% add in global call for depress and L
global depress L

%
%  REFERENCE TEMPERATURE
%  ----------------------
      ZERO_C    = 273.15;   %  conversion between Kelvin and Celsius
%
%
%  SPECIFIC HEAT c(TEMP)
%  ----------------------
%      (Paterson, 1994, p.205)
%        c = 152.5 + 7.122 TEMP
      c_1   = 152.5;
      c_2   = 7.122;
      c_3   = ZERO_C;
      c_vals = [c_1; c_2; c_3;] ;
%
%
%  ICE CONDUCTIVITY k(TEMP)
%  --------------------------
%  (Paterson, 1994, p.205)
%         k = k_1 exp( -k_2 TEMP )
% %  in units of W m^-1 deg^-1.   Note TEMP is in deg Kelvin (hence k_3).
%       k_1   = 9.828;
%       k_2   = 5.7e-3;
%       k_3   = ZERO_C;
%       k_vals = [k_1; k_2; k_3; NaN; k5; k6] ;
% %     
%     
% %        from Yen 1981 fit for T(K) > 195
%       k_1 = 6.727 ;
%       k_2 = 0.0041 ;
%       k_3 = ZERO_C ; %0.5962
%       k_vals = [k_1; k_2; k_3; NaN; k5; k6] ;
%       
      %fitting Waite 2006
      k_1 = 8.985 ;
      k_2 = .005182 ;
      k_3 = ZERO_C ;
      k_vals = [k_1; k_2; k_3; NaN;] ;
%
%  FIRN CONDUCTIVITY k_firn(TEMP)
%  -------------------------------
%    Van Dusen formula  (Paterson, 1994, p.205)
%          k_firn = 2.1e-2 + 4.2e-4 * rho + 2.2e-9 * rho^3
%    in units of W m^-1 deg^-1
%    At the limit rho -> rho_ice, this appears to give conductivity
%    at warm temperatures (-5 C). To get k_firn for other temperatures,
%    I scale it by [ k(TEMP) / k(-5 deg C)]
      VD_1  = 2.1e-2;
      VD_2  = 4.2e-4;
      VD_3  = 2.2e-9;
      VD_vals = [VD_1; VD_2; VD_3];
      % had to change VD_1 to 2.1e-2
%
%
%  SENSIBLE HEAT CONTENT
%  ---------------------
%    function sens_heat_ice(TEMP) = integral{0,TEMP} c(T') dT'
%    latent heat of fusion L = 3.34 10^5 J kg^-1
         L = 3.34e05;
%
%
%  GEOTHERMAL FLUX
%  ----------------
%   Gary recommends 80 mW m^-2 at Taylor Dome, based on transient model
%    that matches borehole measurements by scaling delta 18_O 
%    (as in this code).
     %q_geo = 0.077;       %  W m^-2  %original setting
    % q_geo = 0.100;
%
%
%  LAPSE RATE for surface air temperature
%  -----------------------
        lapse = -0.01;          %  lapse rate (deg/m)
%
%
%  THERMAL PROPERTIES FOR BEDROCK
%  ------------------------------
     k_rock     = 2.8;                      % rock conductivity
     k_vals(4) = k_rock ;
  %   rho_rock   = 2400;                     % rock density
     c_rock     =  750;                     % rock specific heat
     c_vals(4) = c_rock ;

           %add two k vals for adjsuting basal temperature
%       k_vals(5)   = k5;%1.02 ;
%       k_vals(6)   = k6;%.08
      %
%
%  MELTING POINT DEPRESSION WITH PRESSURE
%  --------------------------------------
   depress = -7.42e-08;            % (deg Pa^-1) no dissolved air
   %  depress = -7.42e-08*1.33;       % (deg Pa^-1)
 %  The 1.33 factor if basal water contains dissolved firn air
 %
 
    %depress = -7.42e-08 ;
    %depress = -7.42e-08*1.053 ; %adjusted to match EDC value for dZ of 0.005, but doesn't work for dZ of 0.001
    %depress = -7.42e-08*.9685 ;%*1.06 ; %adjusted to match EDC value for dZ of 0.005, but doesn't work for dZ of 0.001
    %the measurements (extrapolation) gives a bottom valye of -2.11 which
    %is warmer than predicted theoretically, so adjsut the depress to get
    %the same value, now chekc this works with dz spacing
    %this doesn't quite work because it doesn't set the temperature at the
    %rock bed interface but rather the lowest grid point
    %depress = -7.42e-08*.9857 ;%*1.06 ; %adjusted to match EDC value for dZ of 0.005, but doesn't work for dZ of 0.001
    
%
%  If there is more gas in the water than the amount brought with the ice ...
    %%% depress = -7.42e-08*1.66;  % (deg Pa^-1)
     
     % pressure_melting_point = depress * 9.81*3465*917


