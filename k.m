%%%%%%%%%%%%%%%%%%%%%%%%%%  
   function k = k( TEMP, k_vals )
%%%%%%%%%%%%%%%%%%%%%%%%%%  
%  finds the thermal conductivity in ice at temperature TEMP
%
%   k =  9.828 * exp( -5.7e-3 .* ( 273.15 + TEMP ) );
%
%   (Paterson, 1994, p.205)
%
%    This works only for rho and K_ice in SI units !
%    Fortunately this routine should be called only from K_firn.m
%    where T has been dimensionalized.
%
%  =============================================================
%
     %global h H_TIME H_NOW BDOT            %  ice dynamics parameters
     %global RHO_s RHO_ice H_firn           %  firn parameters
     %global TEMP_BAR TEMP_AMPL YEAR        %  seasonality parameters
     %global N_P M_P                        %  mesh parameters

     %global L_char TAU_char TEMP_char K_char C_char RHO_char
                                           %  characteristic scalings
     %global TEMP_ref ZERO_C                %  reference temperatures

     %global k_1 k_2 k_3                    % constants for conductivity k(T)
     %global VD_1 VD_2 VD_3                 % constants in Van Dusen formula
%
%  =============================================================
%
          k =  k_vals(1) * exp( -k_vals(2) .* ( k_vals(3) + TEMP ) );
%          k = ones(size(k))*2.1 ;%remove this for temperature sensitivity
          
% %          %average teh two together
% %          k1 =  k_vals(1) * exp( -k_vals(2) .* ( k_vals(3) + TEMP ) ) * k_vals(8);
% %          k2 =  k_vals(5) * exp( -k_vals(6) .* ( k_vals(7) + TEMP ) );
% %          k = (k1+k2)/2 ;
% %          
%          %k=k*1.02 ; %follow Cuffey 2016 using Waite 2006 
%         k=k*k_vals(5) ;
% %         k(1:end-31) = k(1:end-31).*linspace(1,1.05,length(k(1:end-31)))' ;
%         k(end-30:end) = k(end-30:end)*1.05 ;

%% hard code Ross 1978 Table IIIa values
% a = 335 ;
% b = 3.27 ;
% c = -.00881 ;
% k = a./(TEMP+k_vals(3))+b+c*(TEMP+k_vals(3)) ;
% 
%         T0 = -13 ;
%         in = TEMP > T0 ;
% %        k(in) = k(in) + k(in) * 0.08 .* (TEMP(in) -T0)/-T0 ;
%         k(in) = k(in) + k(in) * k_vals(6) .* (TEMP(in) -T0)/-T0 ;
% %         
%         
%  
%    %     
%          %k = repmat( 2.0418, size( TEMP ) );



aa=1;
         


