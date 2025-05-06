%% develop inverse model for borehole temeprature reconstructions
function [ss_TEMP, Q_melt, Z] = func_run_steady(H,bdot,T0,p,q_geo) 

%%
%I want to try a recent warming to better approximate the borehole
%temperature profile. Let's jsut make the last three time steps 1.5 degrees
%warmer
%%

% I've jsut changed to basal melt rather than geothermal flux, but I
% haven't tested it

%I will have an isotope profile, depth-age relationship, and borehole
%temperature profile as observables

%I will seek to find:
%1) the geothermal flux
%2) the D-J kink height
%3) the modern surface temperature
%4) the late Holocene isotope gradient
%5) the early holocene isotope gradient
%6) the deglacial isotope gradient
%
% I will assume no thickness change


% so I need a forward model that predicts the borehole temperature from
% these five parameters

%at first, I will not consider the uncertainty of the depth-age
%relationship, so that I can computer the accumulation rate history
%directly from the depth-age relationship for the thinning function implied
%by the D_J kink height

%so I will make initial guesses of these 5 parameters, find the misfit with
%the borehole temperature profile, and then update them


%  =============================================================
%
     global s_per_year                    %  seconds per year
     global rho_s rho_ice h_firn
     global k_rock c_rock rho_rock %q_geo  
     global g depress L                % basal melting parameters
%     global lapse
     global rock_mask_P
%     global HDOT BDOT 
%     global tau_test itmax_tau_2
%
%  =============================================================
     s_per_year = 365.25*24*3600;   % seconds per year
     
HDOT = 0 ;

[rho_ice, rho_s, h_firn, rho_rock, g ] = setup_density_rho_sp;
[L, c_vals, k_vals, VD_vals, ZERO_C, lapse  ] = setup_thermal_herc;
%need to define accumulation and surface temperature


%% define forward model parameters
%these are the items that get created in the setup file
%[calc_times, H_vec, Hdot_vec, h_vec, slope, w_vec, TEMP_climate, BDOT_climate, q_geo, MDOT_vec, particle_times, phi_B ] = herc_setupcalc_times = 

    %bdot = .043/.917 ;
    %T0 = -51 ;
    mdot = 0;%
    %H = 2700 ;
    
    
    %define geothermal flux
    %q_geo = .07 ; %
           
% Define fraction of horizontal motion coming from basal sliding
    phi_B = 0 ;
        %  Get slope of ice-sheet surface
    slope = 0.001;

    
     
    %p = 2 ;%

    %%
    [ X_P, Z_P ] = setup_grid( H );
    
   [ ss_TEMP, Q_melt, Z] = ss_T_init(X_P, Z_P, T0, H, k_vals, c_vals, p, bdot/s_per_year, 0, VD_vals, phi_B, q_geo) ;
  

