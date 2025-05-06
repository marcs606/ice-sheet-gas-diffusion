function [xCO2_out,CO2_ratio,age] = CO2_diff_fn(sim_years,phi,save_years,AGE,Znew,TEMP)
% perfroms CO2 diffusion model
%   Detailed explanation goes here

%% User inputs
% sim_years=1.5e6;%total simulation years
% phi=40e3;%period of oscillation (years)
% site='site_226_-126'; %'OIC' from Bereiter 2009 or 'COLDEX'
% save_years=5e4:5e4:sim_years;%saves xCO2_out at initial time and these years;
py=1e6*[0 .3 .6 .9 1.2 1.5];%plot years to make a plot like Figure 6 in Bereiter 2014;


%Model
%props=get_ice_sheet_props(site,sim_years);%age, temperature, and depth in props structure
% props=CO2_get_ice_sheet_props(site,sim_years);%age, temperature, and depth in props structure
% fname = strcat(site,'_',string(datetime,'yyyyMMddHHmmss'));
% fname=strcat('040_ky_',site,'_',props.nm(15:end-4),'_',string(datetime,'yyyyMMddHHmmss')); %file name to save outputsave_years=5e4:5e4:sim_years;%xCO2_out saves initial conditions + mixing ratios at these years

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

z=((Znew(1:age_resolution:N_layers+1)+Znew(age_resolution:age_resolution:N_layers+age_resolution))/2)';
z_bnds=(Znew(1:age_resolution:N_layers+age_resolution+1))';

T=interp1(Znew(~isnan(Znew)),TEMP(~isnan(TEMP)),z);
T=mean(T(1:sw_l));

dzbnds=diff(z_bnds);
dzbnds=mean(dzbnds(1:sw_l));

age=interp1(Znew(~isnan(Znew)),AGE(~isnan(AGE)),z);

age_dummy=0:round(age(end)+1);%calculate initial concentrations
xCO2_dummy=230+50*sin(2*pi/phi*age_dummy);
xCO2=(interp1(age_dummy,xCO2_dummy,age'))';

[~,ind]=min(abs(age-2*phi));%simulate only the annual layers you need after getting switcher years averages
xCO2=xCO2(1:ind);
age=age(1:ind);

mass=dzbnds*rho;%mass of layer, kg
V_ice=dzbnds;%volume of each layer, m^3
n_air=mass*TAC*100000/R/273.15;%moles of air per layer
nCO2=n_air.*xCO2; %moles of CO2 per layer;

sec_per_year=365.25*24*3600;
t=0;

% K = 1; %uncertainty ratios
% C = 1;
%all gas from parameters from Bereiter, 2009; a is Bereiter +3 to make it make sense
D0=9.1e-8; %ikeda-fukazawa 2004
Q_D=14.7e3; %ikeda-fukazawa 2004
S0=1.02e-15;
Q_S=0;
a=10.16;
b=1121;

% r0 = 1;
D=D0*exp(-Q_D/R./T);
X=S0*exp(-Q_S/R./T);
X=X*917/.018;%convert from mol_gas mol^-1_ice Pa^-1 to Henry's law
Pdiss=10.^(a-b./T);

save_cntr=1;
layer_cntr=1;
sec_cntr=0;

xCO2_out=NaN(length(save_years)+1,length(nCO2));%initialize matrix for saving
xCO2_out(1,:)=xCO2;%save initial conditions

dt=.3*dzbnds^2/D;%keep timestep stable


while t<sim_years

    xCO2=nCO2./(Pdiss.*X.*V_ice+n_air);%mass balance
    C_ice=xCO2.*Pdiss.*X;
    
    [deriv]=calc_deriv_ftcs_avgs(C_ice,dzbnds,D,dt);%calculate derivative
    
    C_ice=C_ice+deriv;%timestep
    
    nCO2=C_ice.*V_ice+xCO2.*n_air;%recalculate nCO2
    
    t=t+dt/sec_per_year;
    sec_cntr=sec_cntr+dt;
    
    if sec_cntr/sec_per_year>=switcher %recalculate properties ever switcher years

        xCO2=nCO2./n_air;
              
        layer_cntr=layer_cntr+switcher;
        
        z=((Znew(layer_cntr:age_resolution:layer_cntr+N_layers)+Znew(layer_cntr+age_resolution:age_resolution:age_resolution+layer_cntr+N_layers))/2)';
        z_bnds=(Znew(layer_cntr:age_resolution:layer_cntr+N_layers+age_resolution))';
        
        T=interp1(Znew(~isnan(Znew)),TEMP(~isnan(TEMP)),z);
        T=mean(T(1:sw_l));
        
        dzbnds=diff(z_bnds);
        dzbnds=mean(dzbnds(1:sw_l));
        
        mass=dzbnds*rho;%mass of layer, kg
        V_ice=dzbnds;%volume of each layer, m^3
        n_air=mass*TAC*100000/R/273.15;%moles of air per layer
        nCO2=n_air.*xCO2; %moles of CO2 per layer;
        
        D=D0*exp(-Q_D/R./T);
        X=S0*exp(-Q_S/R./T);
        X=X*917/.018;%convert from mol_gas mol^-1_ice Pa^-1 to Henry's law
        Pdiss=10.^(a-b./T);
        
        dt=.3*dzbnds^2/D;
        sec_cntr=0;
    end
    
    if t>=save_years(save_cntr)
       save_cntr=save_cntr+1;
       xCO2_out(save_cntr,:)=nCO2./n_air;
       disp(strcat(num2str(save_years(save_cntr-1)),' years have passed'))
    end
    
end
  
xCO2_out(end,:)=nCO2./n_air;

% for jj=1:length(py)
%     L{jj}=strcat(num2str(py(jj)/1e6),' My');%make the legend
% end
% 
% dir = 'C:\Users\msail\Downloads\ice_core_diffusion_research\temperature_steady_for_marc\data_profiles';
% load(fullfile(dir,props.nm),'acc','T0','q_geo','H')
% 
% figure(1)
% hold on; grid on; box on;
% plot(age(151:250)-.75*phi,xCO2_out(1,151:250),'LineWidth',1.5)
ind_plot = zeros([1 6]);
ind_plot(1) = 1;
for jj=2:length(py)
    ind=find(save_years==py(jj))+1;
    ind_plot(jj) = ind; % shows which indices of xCO2_out are plotted
%     plot(age(151:250)-.75*phi,xCO2_out(ind,151:250),'LineWidth',1.5)
end
 
CO2_int = xCO2_out(1,151);
CO2_mid = (xCO2_out(1,151)+xCO2_out(1,250))/2;

CO2_ratio = 1 - (xCO2_out(:,151)-CO2_mid)./(CO2_int-CO2_mid);
% age of gas = (n-1)*50e3 yrs where n = 1,...,31 is the index


end