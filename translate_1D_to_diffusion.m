function [AGE,T,Z] = translate_1D_to_diffusion(sim_years,age,depth,ss_temp,z)
%translate_1D_to_diffusion Converts results from 1D model to a readable format
%for the diffusion model

AGE = 0:sim_years+3e5;

if max(age) < sim_years+3e5
    AGE(round(max(age))+1:end) = NaN;
% else
%     AGE = 0:round(max(age));
end
Z = interp1(age,depth,AGE,'pchip','extrap');
z = z(1:1000); % removes data points from 'z' that are under ice bed
ss_temp = ss_temp(1:1000); % ditto
T_formatted = interp1(z(end:-1:1),ss_temp,depth,"linear","extrap"); 
% makes temperature array same size as 'depth' so we can relate temperature to age
T = interp1(age,T_formatted+273.15,AGE,"pchip","extrap");

end