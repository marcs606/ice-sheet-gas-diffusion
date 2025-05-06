%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   function [w_velo w_velo_no_firn]= w_velo( xpos, zpos, H_TIME, h, BDOT, HDOT, MDOT, phi_B )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%  finds the vertical velocity at position (xpos, zpos) 
%  with density structure rho( xpos, zpos )
%  total velocity = dynamic velocity + compaction velocity
%
%  Dynamic velocity given by Dansgaard Johnsen model (1969)
%  ignoring dynamic effects of expanded firn layer
%  H_TIME is ice sheet thickness at current TIME
%  h is always dimensionless fraction of H_TIME
%
%  Compaction velocity given by vertical continuity i.e.
%          d_(rho w)/dz = 0
%
%     w_firn = BDOT * ( RHO_ice ./ rho ) - 1 );
%
%  =============================================================
%
     %global H_NOW HDOT BDOT %H_TIME h
                                           %  ice dynamics parameters
     global rho_s rho_ice h_firn           %  firn parameters
     %global T_BAR T_AMPL YEAR              %  seasonality parameters
     global N_P M_ice_P M_P                %  mesh parameters
%
     global s_per_year
     %global L_char TAU_char T_char K_char C_char RHO_char
                                           %  characteristic scalings
%
%  =============================================================
%
%
%  phi_B = sliding fraction in Dansgaard-Johnsen model
%  m_dot = melting rate in DJ model
%
%      phi_B = 0.15 ;        % TJ's value of sliding (flux fraction)
   %   phi_B = 0.30 ;    % Jakob's value 
%
%
%
%   put HDOT into m/s
%      HDOT_s = HDOT/s_per_year ;
%
%    dynamic velocity given by Dansgaard-Johnsen model
%    =================================================
      w_velo = zeros( size( xpos ) );
  %    shape = zeros( size( xpos ) );
  %    shape1 = zeros( size( xpos ) );
      psi = zeros( size( xpos ) );
  
%   create thickness matrix - column "j" contains thickness at site "j"
      H_matrix = repmat( H_TIME, size(zpos, 1), 1 );
       
      %define parameters for sliding vertical velocity
      H1 = H_TIME ;
      h1 = h*H1 ;
      %phi_B = 0.125 ;

      if ( BDOT ~= 0 )

%          index1  = find( (zpos./H_matrix < h) & ( zpos./H_matrix > 0) ); 
%             if( ~isempty( index1 ) )
%                %shape(index1) = ( (zpos(index1)./H_matrix(index1) ).^2 )./ ...
%                %          ( 2.0 .*h .* ( 1 - h./2.0 ) );
%                      
%          %%%      shape1(index1) = ( phi_B*zpos(index1) + 0.5*(1-phi_B)*zpos(index1).^2/h1 ) / (H1 -0.5*h1*(1-phi_B)) ;     
%                psi(index1) = ( phi_B*zpos(index1) + 0.5*(1-phi_B)*zpos(index1).^2/h1 ) / (H1 -0.5*h1*(1-phi_B)) ;     
%             end     % if( ~isempty( index1 ) ...
% 
%          index2  = find( zpos./H_matrix >= h );
%             if( ~isempty( index2 ) )
%                %shape(index2) = ( zpos(index2)./H_matrix(index2) - h/2.0 ) ...
%                %            / ( 1 - h/2.0 );
%          %%%      shape1(index2) = ( zpos(index2) - 0.5*h1*(1-phi_B) ) / (H1 -0.5*h1*(1-phi_B)) ;        
%                psi(index2) = ( zpos(index2) - 0.5*h1*(1-phi_B) ) / (H1 -0.5*h1*(1-phi_B)) ;        
%             end     % if( ~isempty( index2 ) ...

            
            psi = (1-(h+2)/(h+1)*(1-zpos./H_matrix) + 1/(h+1)*(1-zpos./H_matrix).^(h+2) ) ;
            psi = phi_B*(zpos./H_matrix) + (1-phi_B)*psi ;

%             %it looks like there is a numerical rounding error for the
%             %bottommost point
%             in = psi < 0 ; 
%             psi(in)=0 ;
            
            %adding this to test the tempeature profile when there is no
            %basal motion
%             in = zpos < 500 ;
%             psi(in)= 0 ;
%             in = zpos >=500 & zpos < 1000 ;
%             psi(in) = psi(in)*2 ;
        end       % if( BDOT ~= 0 ) ...


%    add firn compaction velocity given by density distribution
%    and scale with accumulation rate
%    BDOT and HDOT are row vectors defining dynamics and climate 
%    at each site
%    ==========================================================
    %  v = -(BDOT - HDOT) * ( v + ( RHO_ice ./ rho( xpos, zpos ) - 1 ) );
    %  %old velocity, I think its for a single location
    %   v = -repmat( (BDOT - HDOT_s), size(zpos,1), 1) .* ...
    %         ( v + ( rho_ice ./ rho( xpos, zpos, H_TIME ) - 1 ) );
         % working vertical velocity, does not include melting
       %v = -repmat( (BDOT + MDOT - HDOT_s), size(zpos,1), 1) .* ...
       %      ( v + ( rho_ice ./ rho( xpos, zpos, H_TIME ) - 1 ) ) + MDOT;
         % above doesn't work because the firn velocity scaling is being
         % applied to the melting and elevation change also, but it should only be applied to
         % the BDOT
         
       %v = -( (BDOT + MDOT - HDOT_s) * shape ) - ( ( (rho_ice ./ rho( xpos, zpos, H_TIME)) - 1 ) * BDOT ) + MDOT ;
       
              
       w_velo = -( (BDOT - MDOT - HDOT) * psi ) - MDOT ...
               - ( ( (rho_ice ./ rho( xpos, zpos, H_TIME)) - 1 ) * BDOT );
       w_velo_no_firn = -( (BDOT - MDOT - HDOT) * psi ) - MDOT ;
       %this will only be good for 1D
       
       %v = -( (BDOT + MDOT - HDOT_s) * shape ) + MDOT ;
       %v = v - ( ( (rho_ice ./ rho( xpos, zpos, H_TIME)) - 1 ) * BDOT )
       
       %using carlos' relationship
       %w_carlos = -((BDOT - MDOT - HDOT) * (rho_ice ./ rho( xpos, zpos, H_TIME))) .* psi - MDOT ;
       
       %
       %figure(1)
          %plot(v)
          %hold on
%          figure(2)
%          plot(v)
%          hold on
%          [HDOT BDOT]
       clear index1 index2 H_matrix
%
%
%
%
%
%
