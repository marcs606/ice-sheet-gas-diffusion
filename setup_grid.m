%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   function [ X_P, Z_P ] = ...
               setup_grid( H_TIME )
           %%%function [ X_corners, Z_corners, X_P, Z_P ] = ...%%%
               %%%setup_grid( X_vec, H_TIME, H_NOW )%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   set up up 2-D grid for temperature model using time-varying 
%   ice thickness H_TIME
%   In fact I am probably using this only for a 1-D model
%   but being generalized and ready (sort of ...) 
%   for 2-D ain't so bad neither
% ------------------------------------------------------------
%
      global Z_u Z_d X_e X_w dX_P dZ_P dZ_U dZ_D dX_E dX_W
      global N_P M_P M_ice_P M_rock_P   %  mesh parameters
      global rock_mask_P

  %    global H_NOW                % ice sheet thickness at end of run
  %    global x_vec_P z_vec_P      % positions of centers 


%  Define the edges of the Control Volumes
%  ===============================================
     x_vec_corners = H_TIME * [0 1];
     dZ = 0.001;
     %dZ = 0.005;
     %dZ = 0.1 ;
     z_vec_ice = H_TIME * [ 1: -dZ: 0.0 ];
     %z_vec_ice = H_TIME * [ 1: -dZ: 0.0 ]; %original for above
     %z_vec_rock = -H_TIME * cumsum( logspace(log10(dZ), 0, 15) );
     z_vec_rock = -3465 * cumsum( logspace(log10(dZ), 0, 15) ); %tj jan19 using H_TIME causes an interpolation error in the bedrock, set fixed bedrock depth
    
     z_vec_corners = [ z_vec_ice z_vec_rock ];

     [ X_corners, Z_corners ] = meshgrid( x_vec_corners, z_vec_corners );

     N_corners = length(x_vec_corners);
     M_ice_corners = length(z_vec_ice);
     M_rock_corners = length(z_vec_rock);
     M_corners = M_ice_corners + M_rock_corners;



%  Find the Control Volume centers
%  ===============================
     x_vec_P = ( x_vec_corners( 2: N_corners ) ...
                               +  x_vec_corners( 1: N_corners - 1 ) ) /2;
     z_vec_P = ( z_vec_corners( 2: M_corners ) ...
                               +  z_vec_corners( 1: M_corners - 1 ) ) /2;
     N_P = N_corners-1;
     M_P = M_corners-1;
     M_ice_P = M_ice_corners - 1;

     [ X_P, Z_P ] = meshgrid( x_vec_P, z_vec_P );

 %%  X_P = ones(1,M_P)' * x_vec_P;
 %%  Z_P = z_vec_P' * ones(1,N_P);

%  Make mask for the Points in the rock portion of the grid
%  ========================================================
%  (Recall that ice_corners includes ice-rock boundary, rock_corners do not)
     M_ice_P = M_ice_corners - 1;
     M_rock_P = M_rock_corners;

     rock_mask_P = ones( M_P, N_P);
     rock_mask_P(M_ice_P+1: M_P, : ) = zeros( M_rock_P, N_P);


%  Find the centers of the Control Volume faces
%  ============================================
     Z_u = z_vec_corners(1: M_corners -1)' * ones(1,N_corners -1);
     Z_d = z_vec_corners(2: M_corners)' * ones(1,N_corners -1);

     X_e = ones(1,M_P)' * x_vec_corners(2:N_corners);
     X_w = ones(1,M_P)' * x_vec_corners(1:N_corners -1);


%  Find the lengths of the sides of the Control Volumes
%  ====================================================
     dX_P = X_e - X_w;
     dZ_P = Z_u - Z_d;


%  Find the separations of the Control Volume center points
%  ========================================================
%    ( remember that z_vec counts down (-ve z) x_vec counts up (+ve x) )
     dx_P_vec = ( x_vec_P(2:N_P) - x_vec_P(1:N_P -1) );
     dz_P_vec = -( z_vec_P(2:M_P) - z_vec_P(1:M_P -1) );
%
%  fill in intervals at top and bottom by extrapolating to external point
     dz_P_vec = [ -2*(z_vec_P(1) - z_vec_corners(1) )  ...
                dz_P_vec  -2*( z_vec_corners(M_corners) - z_vec_P(M_P) ) ];

     dZ_U = dz_P_vec(1: M_P)' * ones(1,N_P);
     dZ_D = dz_P_vec(2: M_P+1)' * ones(1,N_P);
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Additional form for 2-D
%% =========================
%% fill in intervals at East and West by extrapolating to external point
     dx_P_vec = [ 2*(x_vec_P(1) - x_vec_corners(1) )  ...
              dx_P_vec  2*(  x_vec_corners(N_corners) - x_vec_P(N_P) ) ];
  
     dX_E = ones(1,M_P)' * dx_P_vec(2:N_P +1);
     dX_W = ones(1,M_P)' * dx_P_vec(1:N_P);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%