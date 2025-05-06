%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function [ ss_TEMP, Q_melt, Z0] = ss_T_init( X0, Z0, TEMP_surf, H_TIME, ...
          k_vals, c_vals, h, BDOT, HDOT, VD_vals, phi_B, q_geo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  calculate steady state 1-D vertical temperature distribution
%  following Firestone et al., 1990, J. Glaciol. 36(123) 163-168.
%
%  Problem is nonlinear because kappa depends on temperature.
%  Program starts with constant thermal parameters 
%  K(TEMP_ref), c(TEMP_ref), kappa(TEMP_ref)
%  determined by initial guess at temperature TEMP_ref.
%  Program uses this kappa in call to ss_T_calc.m to find TEMP_trial,
%  and uses this TEMP_trial to estimate better kappa(TEMP_trial).
%  Loop continues until changes in TEMP_trial are smaller than 
%  CONV_TEST = 0.001.
%
%  Also checks that bed is not above pressure melting point (pmp)
%  If bed is too hot, reduces geothermal flux and tries again.
%
%    Q_melt is the heat flux used up in melting in steady state
%                                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  =============================================================
%
     %global H_NOW %HDOT %BDOT %H_TIME h     %  ice dynamics parameters
     global N_P M_P                        %  mesh parameters
     global rho_s  h_firn   rho_ice        %  firn parameters
     %global TEMP_BAR TEMP_AMPL YEAR        %  seasonality parameters
     %global TEMP_ref                       %  reference temperature
     %global c_1 c_2                        % constants for specific heat c(TEMP)
     %global k_1 k_2                        % constants for conductivity k(TEMP)
     %global K_rock C_rock RHO_rock q_geo   %  rock thermal conditions
     global g depress L                    % basal melting parameters
%
%  =============================================================
% add constants that I can't find other places
%TEMP_ref = -26 ;


%  initialize output matrix ss_TEMP steady state temperature profile
      ss_TEMP = zeros( size(Z0) );

      CONV_TEST = 0.001;       %  convergence test for variable coefficients
      BED_TEST  = 0.01 ;       %  convergence test for melting bed
      Q_base    = q_geo;       % trial geothermal flux reduced if bed melts
      d_Q       = q_geo * 0.005;   % flux perturbation for Newton Raphson
                               % when finding (d T_base/d Q_base )

      T_depress = ( rho_ice * g * H_TIME ) * depress;

%  set up finely spaced depth array in ice for trapezoid rule integration
%    note that Z_vec(1) is at the top surface
      dZ = 0.005;
      Z = [ 1: -dZ: 0 ]' * H_TIME;
      X = repmat( [1:length(H_TIME)], size(Z,1), 1); 

%  get initial constant trial temperature distribution
      TEMP_trial_old = TEMP_surf * ones( size(Z) ); %I switched TEMP_ref to be TEMP_surf

%  get depth-variable thermal diffusivity array
    rho_var = rho(X, Z, H_TIME) ;  
    %kappa = ( k(TEMP_trial_old, k_vals) ./ ( rho_var .* c(TEMP_trial_old) ) );
    kappa = ( K_firn(rho_var,TEMP_trial_old, k_vals, VD_vals) ./ ( rho_var .* c(TEMP_trial_old, c_vals) ) );


%  get first realistic guess at temperature
      TEMP_trial = ss_T_calc( X, Z, TEMP_surf, TEMP_trial_old, kappa, Q_base,...
          H_TIME, k_vals, h, BDOT, HDOT, VD_vals, rho_var, phi_B);


%     loop on trial temperature solutions until maximum change is small
%     -----------------------------------------------------------------
    counter=0 ; 
    while ( max( abs( TEMP_trial - TEMP_trial_old ) ) >= CONV_TEST )
        counter=counter+1 ;
        %disp(counter) ;
        TEMP_trial_old = TEMP_trial;
        kappa = ( K_firn(rho_var, TEMP_trial_old, k_vals, VD_vals) ./...
            ( rho_var .* c(TEMP_trial_old, c_vals) ) );
        TEMP_trial = ss_T_calc( X, Z, TEMP_surf, TEMP_trial_old, kappa, Q_base,...
            H_TIME, k_vals, h, BDOT, HDOT, VD_vals, rho_var, phi_B);
     end    %  while ( max( abs( TEMP_trial ...



%  test whether bed is above pressure melting temperature
%     If yes, adjust heat flux into base to bring bed to pmt
%  ---------------------------------------------------------
     T_base = TEMP_trial(end);
     hot_bed = (T_base - T_depress); 

     if( hot_bed > 0 )
       del_Q = 0;


%    loop on basal flux Q_base until hot_bed -> 0
%    --------------------------------------------
       while ( abs(hot_bed) > BED_TEST )

%     find basal temperature for Q_base
%     ---------------------------------
           TEMP_trial_old = TEMP_trial + 2*CONV_TEST;
          while ( max( abs( TEMP_trial - TEMP_trial_old ) ) >= CONV_TEST )
             TEMP_trial_old = TEMP_trial;
             kappa = ( K_firn(rho_var, TEMP_trial_old, k_vals, VD_vals) ./...
            ( rho_var .* c(TEMP_trial_old, c_vals) ) );
             TEMP_trial =  ...
                   ss_T_calc( X, Z, TEMP_surf, TEMP_trial_old, kappa, Q_base, H_TIME, k_vals,...
                   h, BDOT, HDOT, VD_vals, rho_var, phi_B);
           end    %  while ( max( abs( TEMP_trial ...

           T_base_0 = TEMP_trial(end);


%     find basal temperature for Q_base + d_Q
%     ---------------------------------------

%     loop on trial temperature solutions until maximum change is small
%         (Newton-Raphson)
%     -----------------------------------------------------------------
           TEMP_trial_old = TEMP_trial + 2*CONV_TEST;
          while ( max( abs( TEMP_trial - TEMP_trial_old ) ) >= CONV_TEST )
             TEMP_trial_old = TEMP_trial;
             kappa = ( K_firn(rho_var, TEMP_trial_old, k_vals, VD_vals) ./...
            ( rho_var .* c(TEMP_trial_old, c_vals) ) );
	     TEMP_trial =  ...
                ss_T_calc( X, Z, TEMP_surf, TEMP_trial_old, kappa, Q_base+d_Q, H_TIME,...
                k_vals, h, BDOT, HDOT, VD_vals, rho_var, phi_B);
           end    %  while ( max( abs( TEMP_trial ...

           T_base_1 = TEMP_trial(end);      



%     find new del_Q for Newton Raphson adjustment
%     --------------------------------------------
          hot_bed = T_base_1 - T_depress;

          del_Q = - hot_bed / ( (T_base_1 - T_base_0)/d_Q );  
          Q_base = Q_base + del_Q;


      end      % while ( abs(hotbed) > BED_TEST ) ...


    end         % if( hot_bed >= 0 ) ...



      Q_melt = q_geo - Q_base;




%  Now put temperatures back on Z grid
%  ===================================
       Z_bed = Z(end, :);

%    find depths that are in ice
      index = find( Z0 >= Z_bed);
      Z_ice = Z0(index);

%  interpolate ice temperatures at depths Z_ice
      ss_TEMP(index) = interp1( Z, TEMP_trial, Z_ice);

%  set up depth array in rock
      index = find( Z0 <= Z_bed);
      Z_rock = Z0(index);

%  set constant temperature gradient in rock using Q_geo
      ss_TEMP(index) = TEMP_trial(end) ...
                             + (q_geo/k_vals(4)) * ( Z_bed - Z_rock ); %k_vals(4) is k_rock
      %ss_TEMP(index) = TEMP_trial(end) ...
      %                       + (.0193/2.8) * ( Z_bed - Z_rock );
      ss_TEMP = ss_TEMP';
      save('ss_TEMP.mat','ss_TEMP','Z0') ;
%
%
%
%
%
