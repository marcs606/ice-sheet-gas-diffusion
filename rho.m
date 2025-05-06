%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   function rho = rho( xpos, zpos, H_TIME )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%  finds the density at position (xpos, zpos) 
%  using simple exponential with e-folding depth H_firn
%
%  =============================================================
%
     global rho_s rho_ice h_firn           %  firn parameters
     
%  =============================================================


%   create thickness matrix - column "j" contains thickness at site "j"
   %  H_matrix = repmat( H_thick, size(zpos, 1), 1 );

 %   rho = RHO_s +  ( RHO_ice -RHO_s ) * ...
 %       ( 1 - exp( -( repmat( H_thick, size(zpos, 1), 1 ) - zpos )/H_firn ) );

% clear H_matrix

     %rho = 450 + ( rho_ice -450 ) * ...
     %                        ( 1 - exp( -( H_TIME - zpos )/(27+1) ) );
                         
     %rho = rho_ice*ones(length(rho),1) ;
   %
%

    rho = interp1(H_TIME-h_firn,rho_s,zpos) ;
    %rho = rho_ice*ones(length(rho),1) ;
    a=1;
    