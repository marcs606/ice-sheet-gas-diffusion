%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   function K_firn = K_firn( RHO, T, k_vals, VD_vals )
%   function [K_firn K_firn_vd K_firn_s K_firn_c] = K_firn( RHO, T, k_vals, VD_vals )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%  finds the thermal conductivity in ice/firn 
%  with density structure RHO( xpos, zpos )
%
%  The routine takes average of Van Dusen and Schwertfeger formulae 
%   (Paterson, 1994, p.205)
%
%    Van Dusen formula  (SI units) 
%        k_firn = 2.1e-3 + 4.2e-4 * RHO + 2.2e-9 * RHO^3
%
%    Schwertfeger formula  (SI units)
%        k_firn = [ 2 k_ice RHO ] / [ 3 RHO_ice - RHO ]
%
%   I scale the Van Dusen formula by the conductivity of solid ice 
%   at the actual temperature
%
%  =============================================================
%
     %global L_char TAU_char T_char K_char C_char RHO_char
                                           %  characteristic scalings
     %global VD_1 VD_2 VD_3                 % constants in Van Dusen formula
     global rho_ice
     %
%  =============================================================
%
%     K_firn = k(T,k_vals) .* ...
%          ( ( VD_vals(1) + VD_vals(2) * RHO + VD_vals(3) * RHO.^3 ) ...
%              ./ k(-5, k_vals) ...
%             + ( 2 * RHO )./( 3 * rho_ice - RHO ) ) /2;
% %         
% %    
% % %   VD only
%    K_firn_vd = k(T,k_vals) .* ...
%           ( VD_vals(1) + VD_vals(2) * RHO + VD_vals(3) * RHO.^3 ) ...
%              ./ k(-5, k_vals) ;%...
% %            + ( 2 * RHO )./( 3 * rho_ice - RHO ) ) /2;
% %         
% %   Schwertfeger  only
K_firn = k(T,k_vals) .* ( 2 * RHO )./( 3 * rho_ice - RHO ) ;
% %    
% %         rho_char seems to be rho_ice so use instead
% %         
% % use approximate of Calonne, 2019

% a=0.02 ;
% rho_transition = 450 ;
% theta = 1./(1+exp(-2*a*(RHO-rho_transition))) ;
% k_ref_firn = 2.107 + 0.003618 * (RHO - rho_ice) ;
% k_ref_snow = 0.024 - 1.23e-4*RHO + 2.5e-6*RHO.^2 ;
% %ka = .02 * ones(size(k(T,k_vals) ;
% ka_ref = .024 ;
% ka = ones(size(RHO))*ka_ref ; %make thermal conducitivity of air constant
% ki_ref = 2.107 ;
% K_firn = (1-theta) .* k(T,k_vals).*ka/(ki_ref*ka_ref).*k_ref_snow + theta.*k(T,k_vals)/ki_ref.*k_ref_firn ;

%aa=1;
% 
%     K_firn2 = k(T,k_vals) ; %ignore firn conductivity effects
%     
%     K_firn=(K_firn1+K_firn2)/2 ;
% % 
%
%
%
%
