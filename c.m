%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   function c = c( TEMP, c_vals )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%  finds the specific heat capacity of ice at temperature TEMP
%
%   c =  152.5 + 7.122 * ( 273.15 + TEMP );
%     =   c_1  +  c_2  * (  c_3   + TEMP );
%
%   (Paterson, 1994, p.205)
%
%    This works only for rho and K_ice in SI units !
%    Fortunately this routine should be called only from *temp*.m
%    where T has not been nondimensionalized.
%  =============================================================
%
  %   global c_1 c_2 c_3      % constants for specific heat c(T)
  %   global C_char

  
        
      c =  c_vals(1) + c_vals(2) * ( c_vals(3) + TEMP );
%      c=ones(size(c))*2097 ;%remove this for temperature sensitive
 
      % c = repmat(2000, size(TEMP) );

%
%
%
%
