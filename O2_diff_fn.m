function [delta_out,O2_ratio,age] = O2_diff_fn(sim_years,phi,save_years,ice_sheet_age,ice_sheet_z,ice_sheet_T)
% computes O2/N2 diffusion
%   Detailed explanation goes here

%% User inputs
% sim_years=1.9e6;%total simulation years
% phi=20e3;%period of oscillation (years)
% site='data_065'; %'OIC' from Bereiter 2009 or 'COLDEX'
% save_years=5e4:5e4:sim_years;%saves xO2_out at initial time and these years;
%py=1e6*[0 .4 .8 1.2 1.6 1.9];%plot years to make a plot like Figure 6 in Bereiter 2014;


% Model
% props=O2_get_ice_sheet_props(site,sim_years);%age, temperature, and depth in props structure
% fname=strcat(site,'_',props.nm(15:end-4),'_',string(datetime,'yyyyMMddHHmmss')); %file name to save outputsave_years=5e4:5e4:sim_years;%xCO2_out saves initial conditions + mixing ratios at these years
% fname = strcat(site,'_',string(datetime,'yyyyMMddHHmmss')); % file name data set + date/time

TAC=.091;%cm^3 STP g^-1, from Raynaud 2007
TAC=TAC*1e-6*1e3;%m^3 kg^-1

rho=917; %kg m^-3;
R=8.314; %gas constant

N_layers=2*phi;%simulate this many annual layers at a time

age_resolution=phi/200; %how far apart are your layers in age

switcher=50e3; %recalculate properties after this many years

sw_l=round(switcher/age_resolution);%average at least switcher years for properties
if switcher>N_layers
N_layers=switcher;
end

z=((ice_sheet_z(1:age_resolution:N_layers+1)+ice_sheet_z(age_resolution:age_resolution:N_layers+age_resolution))/2)';
z_bnds=(ice_sheet_z(1:age_resolution:N_layers+age_resolution+1))';

T=interp1(ice_sheet_z(~isnan(ice_sheet_z)),ice_sheet_T(~isnan(ice_sheet_T)),z);
T=mean(T(1:sw_l));

dzbnds=diff(z_bnds);
dzbnds=mean(dzbnds(1:sw_l));

age=interp1(ice_sheet_z(~isnan(ice_sheet_z)),ice_sheet_age(~isnan(ice_sheet_age)),z);

age_dummy=0:round(age(end)+1);%calculate initial concentrations
xO2_dummy=.21+1.1e-6*sin(2*pi/phi*age_dummy);
xO2=(interp1(age_dummy,xO2_dummy,age'))';

[~,ind]=min(abs(age-2*phi));%simulate only the annual layers you need after getting switcher years averages
xO2=xO2(1:ind);
age=age(1:ind);

mass=dzbnds*rho;%mass of layer, kg
V_ice=dzbnds;%volume of each layer, m^3
n_air=mass*TAC*100000/R/273.15;%moles of air per layer
nCO2=n_air.*xO2; %moles of CO2 per layer;

sec_per_year=365.25*24*3600;
t=0;

D0=3.5e-9;%from Ikeda-Fukuzawa,2005 
Q_D=9.7e3;
X0 = 3.7e-7 * 1e-6; % fast set
Q_X = 7.9e3; % fast set
% P0=7.5e-8/1e-6/sec_per_year;%salamatin 2001
% Q_P=50e3;
a=9.68;%salamatin 2001
b=717;
% a=10.17;%bereiter 2009
% b=850;

Pdiss=10.^(a-b./T);
Pdiss220=10.^(a-b/220);
D=D0*exp(-Q_D/R./T);
X=X0*exp(-Q_X/R./T); % fast set
% P=P0*exp(-Q_P/R./T)./Pdiss220;
% X=P./D;
X=X*917/.018;%convert from mol_gas mol^-1_ice Pa^-1 to Henry's law

save_cntr=1;
layer_cntr=1;
sec_cntr=0;

xO2_out=NaN(length(save_years)+1,length(nCO2));%initialize matrix for saving
xO2_out(1,:)=xO2;%save initial conditions

dt=.5*dzbnds^2/D;%keep timestep stable


while t<sim_years

    xO2=nCO2./(Pdiss.*X.*V_ice+n_air);%mass balance
    C_ice=xO2.*Pdiss.*X;
    
    [deriv]=calc_deriv_ftcs_avgs(C_ice,dzbnds,D,dt);%calculate derivative
    
    C_ice=C_ice+deriv;%timestep
    
    nCO2=C_ice.*V_ice+xO2.*n_air;%recalculate nCO2
    
    t=t+dt/sec_per_year;
    sec_cntr=sec_cntr+dt;
    
    if sec_cntr/sec_per_year>=switcher %recalculate properties ever switcher years

        xO2=nCO2./n_air;
              
        layer_cntr=layer_cntr+switcher;
        
        z=((ice_sheet_z(layer_cntr:age_resolution:layer_cntr+N_layers)+ice_sheet_z(layer_cntr+age_resolution:age_resolution:age_resolution+layer_cntr+N_layers))/2)';
        z_bnds=(ice_sheet_z(layer_cntr:age_resolution:layer_cntr+N_layers+age_resolution))';
        
        T=interp1(ice_sheet_z(~isnan(ice_sheet_z)),ice_sheet_T(~isnan(ice_sheet_T)),z);
        T=mean(T(1:sw_l));
        
        dzbnds=diff(z_bnds);
        dzbnds=mean(dzbnds(1:sw_l));
        
        mass=dzbnds*rho;%mass of layer, kg
        V_ice=dzbnds;%volume of each layer, m^3
        n_air=mass*TAC*100000/R/273.15;%moles of air per layer
        nCO2=n_air.*xO2; %moles of CO2 per layer;
        Pdiss=10.^(a-b./T);
        D=D0*exp(-Q_D/R./T);
        X=X0*exp(-Q_X/R./T); % fast set
        % P=P0*exp(-Q_P/R./T)./Pdiss220;
        % X=P./D;
        X=X*917/.018;%convert from mol_gas mol^-1_ice Pa^-1 to Henry's law
        
        dt=.5*dzbnds^2/D;
        sec_cntr=0;
    end
    
    if t>=save_years(save_cntr)
       save_cntr=save_cntr+1;
       xO2_out(save_cntr,:)=nCO2./n_air;
       disp(strcat(num2str(save_years(save_cntr-1)),' years have passed'))
    end
    
end
  
xO2_out(end,:)=nCO2./n_air;

delta_out=((xO2_out/.78)/(.21/.78)-1)*1e6;

amp_0 = delta_out(1,151);
O2_ratio = 1 - delta_out(:,151)./amp_0;

end