% Integrated temp- age-depth model + diffusion model 

%% temp- age-depth model

clear;

% Input parameters
acc = 0.02 ; % [m yr^-1 ice eq.]
T0 = -60 ; % [deg C]
p = 4 ; % Lliboutry shape parameter appropriate for slow flanks
q_geo = .05 ; % [W m^-2]
H = 2500; % [m]

disp('Calculating temperature profile...')
[ss_TEMP, Q_melt, z] = func_run_steady(H,acc,T0,p,q_geo) ;
%Q_melt is the geothermal flux going into melting, so the melt rate is
%approximately 1/10 of this

disp('Calculating age profile...')
[depth, age] = steady_depth_age(acc, Q_melt/10, H, p) ;

%% diffusion model

sim_years = 1.5e6;
save_years=5e4:5e4:sim_years;%saves xCO2_out at initial time and these years;
py=round(1e6*(0:0.3:1.5),-1);

[age_diff,T_diff,z_diff] = translate_1D_to_diffusion(sim_years,age,depth,ss_TEMP,z);

disp('Diffusing CO2')
[xCO2,CO2_SDR,CO2_age_plot] = CO2_diff_fn(sim_years,40e3,save_years,age_diff,z_diff,T_diff);

disp('Diffusing O2')
[xO2,O2_SDR,O2_age_plot] = O2_diff_fn(sim_years,20e3,save_years,age_diff,z_diff,T_diff);

[maxCO2SDR, maxind] = max(CO2_SDR);
maxO2SDR = max(O2_SDR);
n = ceil((maxind-1)/5);

%% Plot

% make the legend
for jj=1:n
    L{jj}=strcat(num2str(py(jj)/1e6),' Ma');
end

% parameter values
txt = {'Forcings:',...
    ['A = ' num2str(round(acc*100,1)) ' cm yr^{-1}'],...
    ['H = ' num2str(H) ' m'],...
    ['Q = ' num2str(round(q_geo*1000),2) ' mW m^{-2}'],...
    ['T_S = ' num2str(T0) ' \circC']}';
dim = [0.68 0.64 0.2 0.25];

co = '#000000';
figure(10)
ax(1) = subplot(1,6,1); % temp-depth
hold on; box on;
title('Temperature Profile','Fontname','SansSerif','Fontsize',14)
plot(ss_TEMP,z,'Color',co,'linewidth',2)
ylim([0 H])
text(T0*15/16,H/16,'(a)',FontSize=24,FontName='SansSerif',FontWeight='bold')
xlabel('\circC','Fontname','SansSerif','FontSize',14)
ylabel('Height above bed (m)','Fontname','SansSerif','FontSize',14)
text(T0+9,H-200,['Melt = ' num2str(round(Q_melt/10*1000,1)) ' mm yr^{-1}'], ...
    'Fontname','SansSerif','Fontsize',14)

ax(2) = subplot(162); % age-depth
hold on; box on;
plot(age/1e6,H-depth,'Color',co,'linewidth',2)
xlim([0 2])
ylim([0 H])
text(0.125,H/16,'(b)',FontSize=24,FontName='SansSerif',FontWeight='bold')
title('Age Profile','Fontname','SansSerif','Fontsize',14)
xlabel('Age (Ma)','Fontname','SansSerif','FontSize',14)
text(0.4,H-600,txt,'Fontname','SansSerif','Fontsize',14,'LineWidth',1)
ax(2).YTickLabel = {[]};

ax(3) = subplot(1,6,[3,4]); % CO2
hold on; grid on; box on;
plot((CO2_age_plot(151:250)-.75*40e3)/1e3,xCO2(1,151:250),'LineWidth',1.5)
ind_plot = zeros([1 length(py)]);
ind_plot(1) = 1;
for jj=2:length(py)
    ind=find(save_years==py(jj))+1;
    ind_plot(jj) = ind; % shows which indices of xCO2_out are plotted
    plot((CO2_age_plot(151:250)-.75*40e3)/1e3,xCO2(ind,151:250),'LineWidth',1.5)
end
title('CO_2 Signal Amplitudes','Fontname','SansSerif','Fontsize',14)
text(2,270,'SDR ('+string(round((maxind-1)*5e-2,3))+' Ma) = '+string(round(maxCO2SDR,2)), ...
    'Fontname','SansSerif','Fontsize',15)
text(1,172.5,'(c)',FontSize=24,FontName='SansSerif',FontWeight='bold')
xlabel('Years (kyr)','Fontname','SansSerif','FontSize',14)
ylabel('CO_2 (ppm)','Fontname','SansSerif','FontSize',14)
legend(L,'Color','#ffffff','Fontname','SansSerif','Fontsize',15,'location','southeast')
ylim([160 300])

ax(4) = subplot(1,6,[5,6]); % O2/N2
hold on; grid on; box on;
plot((O2_age_plot(151:250)-.75*20e3)/1e3,xO2(1,151:250),'LineWidth',1.5)
for jj=2:length(py)
    ind=find(save_years==py(jj))+1;
    plot((O2_age_plot(151:250)-.75*20e3)/1e3,xO2(ind,151:250),'LineWidth',1.5)
end
text(0.5,-5.75,'(d)',FontSize=24,FontName='SansSerif',FontWeight='bold')
title('O_2/N_2 Ratio Amplitudes','Fontname','SansSerif','Fontsize',14)
xlabel('Years (kyr)','Fontname','SansSerif','FontSize',14)
ylabel(['\deltaO_2/N_2 (',char(8240),')'],'Fontname','SansSerif','FontSize',14)
legend(L,'Color','#ffffff','Fontname','SansSerif','Fontsize',14,'Location','southeast')
text(1,4,'SDR ('+string(round((maxind-1)*5e-2,3))+' Ma) = '+string(round(maxO2SDR,2)), ...
    'Fontname','SansSerif','Fontsize',15)
ylim([-7 7])

% set size and postion
set(gcf,'units','normalized','OuterPosition',[0 0.25 1 0.65])

ax(1).Position = [0.10 0.11 0.12 0.815];
ax(2).Position = [0.23 0.11 0.12 0.815];
ax(3).Position = [0.40 0.11 0.24 0.815];
ax(3).PlotBoxAspectRatio = [1 1 1];
ax(4).Position = [0.68 0.11 0.24 0.815];
ax(4).PlotBoxAspectRatio = [1 1 1];
