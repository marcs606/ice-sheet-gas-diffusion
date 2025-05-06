%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function ss_T_calc = ss_T_calc( X, Z, TEMP_surf, TEMP_trial, kappa, Q_base,...
      H_TIME, k_vals, h, BDOT, HDOT, VD_vals, rho_var, phi_B )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  calculate steady state 1-D vertical temperature distribution
%  following Firestone et al., 1990, J. Glaciol. 36(123) 163-168.
%
%  this version is called repeatedly by ss_T_init.m
%  until kappa and ss_T_calc are mutually consistent and
%  boundary conditions TEMP_surf and Q_base are matched.
% -------------------------------------------------------
%
     %global K_rock C_rock RHO_rock   %  rock thermal conditions

     
     %set MDOT = 0 so that this will initialize
     MDOT = 0 ;

%  initialize output matrix ss_T_calc for steady state temperature profile
      ss_T_calc = zeros( size(Z) );


      dZ = diff(Z);
 
%  get integrating factor F(z)
%   "top" refers to value at top of trapezoidal integration interval
%   "bott" refers to value at bottom of interval
%   because of the first "flipud", intgrand_1(1) is at the bedrock interface
%   the final flipud for F restores F(1) to the ice sheet surface
%   I_dirn = 1 for integration surface to bed
%   I_dirn = -1 for integration bed to surface 
      %v_0 = v(X, Z, H_TIME, h, BDOT, HDOT, MDOT);
      w_0 = w_velo(X, Z, H_TIME, h, BDOT, HDOT, MDOT, phi_B);
      
      w_k_top = w_0(1:end-1,:) ./ kappa(1:end-1,:);
      w_k_bott = w_0(2:end,:) ./ kappa(2:end,:);
      I_dirn = -1;
      intgrand_1 = I_dirn * dZ .* ( w_k_top + w_k_bott )/2;
      intgral_1 = cumsum( flipud(  ...
                         [ intgrand_1' zeros( size(TEMP_surf') ) ]' ) );

      F = flipud( exp( - intgral_1 ) );


%  integrate 1/F(Z) get temperature
%   "flipud" restores the first element of F_top etc 
%   to the ice sheet surface
      F_top = 1 ./ F(1:end-1,:);
      F_bott = 1 ./ F(2:end,:);
      I_dirn = 1;
      intgrand_2 = I_dirn * dZ .* ( F_top + F_bott )/2;
      intgral_2 = cumsum( [  zeros( size(TEMP_surf') )  intgrand_2' ]' );

%  set Z_bed and TEMP_bed (initial guess) conditions and ice-rock interface
      
      %ss_T_calc = repmat(TEMP_surf, size(X,1),1)  ...
      %         - repmat( ( Q_base./k( TEMP_trial(end,:), k_vals ) ), size(X,1), 1 ) ...
      %              .* intgral_2;
%      ss_T_calc = repmat(TEMP_surf, size(X,1),1)  ...
%               - repmat( ( Q_base./K_firn(rho_var(end), TEMP_trial(end,:), k_vals, VD_vals) ), size(X,1), 1 ) ...
%                    .* intgral_2;
      ss_T_calc = repmat(TEMP_surf, size(X,1),1)  ...
               - ( Q_base./K_firn(rho_var, TEMP_trial(end,:),  ...
                                      k_vals, VD_vals) .* intgral_2 );


clear v_0 v_k_top v_k_bott intgrand_1 intgral_1  ...
         F  F_top  F_bott intgrand_2 intgral_2

%
%
%
%
%
