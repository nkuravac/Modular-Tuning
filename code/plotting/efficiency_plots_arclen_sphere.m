close all
clearvars

set(groot,'DefaultAxesFontSize',14)

%%
plotdir='../../reports/figures/efficiency_plots_test';
if ~exist(plotdir,'dir')
    mkdir(plotdir)
end

%%
full_sweep_data=table2array(readtable('../../data/new_arclen_sweep_sphere/summary/arclen_sweep_table_sphere_new_opt_ineff.txt'));
opt_ineff_struct=jsondecode(fileread('../../data/ineff_optimize_sphere/summary/optimal_parameters_purcell_ineff.json'));
config_parameters=jsondecode(fileread('../../data/ineff_optimize_sphere/summary/config_parameters.json'));

mu=config_parameters.Viscosity;
r=config_parameters.r;

%%%restrict to only "reasonable" arclengths (arclength >= wavelength)
%full_sweep_data=full_sweep_data(full_sweep_data(:,1)>opt_ineff_struct.Wavelength,:);

arclengths=full_sweep_data(:,1);
forces=full_sweep_data(:,4);
torques=full_sweep_data(:,5);
translational_speeds=full_sweep_data(:,6);

stokes_drags=6*pi*mu*r*translational_speeds;
max_power=2*pi*1200*200;

%% Part I: speed at fixed power

speeds_on_power_curve=translational_speeds.*sqrt(max_power/2/pi./torques);
purcell_efficiencies=stokes_drags.*translational_speeds./(torques*2*pi);
miscalc_purcell_efficiencies=forces.*translational_speeds./(torques*2*pi);
[max_eff,opt_eff_index]=max(purcell_efficiencies);
[max_speed,max_speed_index]=max(speeds_on_power_curve);

arclength_axis_label='Arc length ({\mu}m)';
speed_axis_label='Speed ({\mu}m/s)';
efficiency_axis_label='Purcell efficiency';

figure
plot(arclengths,speeds_on_power_curve);
xlabel(arclength_axis_label)
xline(arclengths(opt_eff_index),'Label','Maximum propulsive efficiency','LabelVerticalAlignment','middle',FontSize=14)
ylabel(speed_axis_label)
exportgraphics(gca,[plotdir,'/speed_vs_S.png'])

figure
hold on
plot(arclengths,miscalc_purcell_efficiencies)
plot(arclengths,purcell_efficiencies)
xline(arclengths(max_speed_index),'Label','Maximum speed','LabelVerticalAlignment','bottom',FontSize=14)
xlabel(arclength_axis_label)
ylabel(efficiency_axis_label)
legend('\epsilon (simulated drag)','\epsilon (Stokes drag)')
exportgraphics(gca,[plotdir,'/eff_vs_S.png'])

figure
hold on
plot(purcell_efficiencies,speeds_on_power_curve,'LineStyle','none','Marker','o')
plot(purcell_efficiencies,sqrt(max_power*purcell_efficiencies/(6*pi*config_parameters.r*config_parameters.Viscosity)))
legend('Simulation','Theory','Location','northwest')
xlabel(efficiency_axis_label)
ylabel(speed_axis_label)
hold off
exportgraphics(gca,[plotdir,'/speed_vs_eff.png'],'Width',8,'Height',6,'Units','inches','Resolution',300)

%% Part I.1: efficiency and load line slope

ll_axis_label='Load line slope (pN{\cdot}nm{\cdot}s)';

figure
plot(arclengths,torques)
xlabel(arclength_axis_label)
ylabel(ll_axis_label)
exportgraphics(gca,[plotdir,'/ll_slope_vs_S.png'])

figure
plot(torques,speeds_on_power_curve);
xlabel(ll_axis_label)
ylabel(speed_axis_label)
exportgraphics(gca,[plotdir,'/U_vs_ll_slope.png'])

figure
plot(torques,purcell_efficiencies);
xlabel(ll_axis_label)
ylabel(efficiency_axis_label)
exportgraphics(gca,[plotdir,'/eff_vs_ll_slope.png'])

figure
hold on
freq_array=linspace(10,1350,100);
plot(freq_array,max_power/2/pi./freq_array,'Color','k');
for i=1:length(torques)
    if i==1
        handle_toggle='on';
    else
        handle_toggle='off';
    end
    freq_points=[0,sqrt(max_power/2/pi/torques(i))];
    plot(freq_points,torques(i)*freq_points,'HandleVisibility',handle_toggle,'Color','b')
end
hold off
xlabel('Frequency (Hz)')
ylabel('Torque (pN{\cdot}nm)')
legend('Fixed power curve','Load lines')
ylim([0,2000])
exportgraphics(gca,[plotdir,'/all_load_lines_geometry.png'])

%% Part II: Geometry vs. tuning

Omega_0=2*sqrt(torques*max_power/(2*pi))./(torques(opt_eff_index)+torques);
Omega_m=sqrt(max_power/2/pi./torques);

U_0=translational_speeds(opt_eff_index)*Omega_0;
U_m=translational_speeds.*Omega_m;

Omega_ratio=Omega_m./Omega_0;
m_ratio=torques(opt_eff_index)./torques;
U_ratio=U_m./U_0;

power_function_omega=(2-1./Omega_ratio)./Omega_ratio;   %(2*Omega_ratio./(1+Omega_ratio.^2)).^2;
power_function_m=4*m_ratio./(1+m_ratio).^2;
U_ratio_omega_theory=1./sqrt(power_function_omega*max_eff./purcell_efficiencies);
U_ratio_m_theory=1./sqrt(power_function_m*max_eff./purcell_efficiencies);

for repeat=1:2
    figure
    hold on
    plot(Omega_ratio,U_ratio,'LineStyle','none','Marker','o')
    plot(Omega_ratio,U_ratio_omega_theory)
    xlabel('\Omega_m/\Omega_0')
    ylabel('U_m/U_0')
    legend('Simulation','Theory')
    xline(1,'LineStyle','--','HandleVisibility','off')
    hold off
    if repeat==1
        exportgraphics(gca,[plotdir,'/speed_ratio_vs_freq_ratio_geometry_change.png'])
    else
        xscale log
        exportgraphics(gca,[plotdir,'/logx_speed_ratio_vs_freq_ratio_geometry_change.png'])
    end
end

for repeat=1:2
    figure
    hold on
    plot(m_ratio,U_ratio,'LineStyle','none','Marker','o')
    plot(m_ratio,U_ratio_m_theory)
    xlabel('m_0/m_1')
    ylabel('U_1/U_0')
    legend('Simulation','Theory')
    xline(1,'LineStyle','--','HandleVisibility','off')
    hold off
    if repeat==1
        exportgraphics(gca,[plotdir,'/speed_ratio_vs_ll_slope_ratio_geometry_change.png'])
    else
        xscale log
        exportgraphics(gca,[plotdir,'/logx_speed_ratio_vs_ll_slope_ratio_geometry_change.png'])
    end
end
%%
figure
hold on
plot(m_ratio,purcell_efficiencies/max_eff)
plot(m_ratio,power_function_m)
xlabel('m_0/m_1')
legend('\epsilon_1/\epsilon_0','P_{max}/P')
exportgraphics(gca,[plotdir,'/eff_ratio_and_power_ratio.png'])




%% Part III: Viscosity vs. tuning

mu_array=linspace(0.93,16.5,12);
ll_slopes_mu=torques(opt_eff_index)*mu_array/mu;

Omega_0_mu=2*sqrt(ll_slopes_mu*max_power/(2*pi))./(mu_array(1)/mu_array(1)*torques(opt_eff_index)+ll_slopes_mu); %tuned for other viscosity but using water motor
Omega_m_mu=sqrt(max_power/2/pi./ll_slopes_mu);

figure
hold on
freq_array=linspace(10,400,100);
plot(freq_array,max_power/2/pi./freq_array,'Color','k');
for i=1:length(ll_slopes_mu)
    if i==1
        handle_toggle='on';
    else
        handle_toggle='off';
    end
    freq_points=[0,sqrt(max_power/2/pi/ll_slopes_mu(i))];
    plot(freq_points,ll_slopes_mu(i)*freq_points,'HandleVisibility',handle_toggle,'Color','b')
end
hold off
xlabel('Frequency (Hz)')
ylabel('Torque (pN{\cdot}nm)')
legend('Fixed power curve','Load lines')
ylim([0,5000])
exportgraphics(gca,[plotdir,'/all_load_lines_viscosity.png'])

U_0_mu=translational_speeds(opt_eff_index)*Omega_0_mu;
U_m_mu=translational_speeds(opt_eff_index)*Omega_m_mu;

Omega_ratio_mu=Omega_m_mu./Omega_0_mu;
U_ratio_mu=1./(U_0_mu./U_m_mu);

power_function_mu=(2-1./Omega_ratio_mu)./Omega_ratio_mu;  %(2*Omega_ratio_mu./(1+Omega_ratio_mu.^2)).^2;
U_ratio_theory_mu=1./sqrt(power_function_mu.*mu_array/mu_array(1));

for repeat=1:2
    figure
    hold on
    plot(Omega_ratio_mu,U_ratio_mu,'LineStyle','none','Marker','o')
    plot(Omega_ratio_mu,U_ratio_theory_mu)
    xlabel('\Omega_m/\Omega_0')
    ylabel('U_m/U_0')
    legend('Simulation','Theory')
    hold off
    if repeat==1
        exportgraphics(gca,[plotdir,'/speed_ratio_vs_freq_ratio_viscosity_change.png'])
    else
        xscale log
        exportgraphics(gca,[plotdir,'/logx_speed_ratio_vs_freq_ratio_viscosity_change.png'])
    end
end

%%
figure
plot(Omega_ratio,power_function_omega)

%% Purcell efficiency and viscosity

viscosity_sweep_data=readtable('../../data/viscosity_sweep/summary/viscosity_sweep_table.txt');
stokes_drag_coeffs=6*pi*r*viscosity_sweep_data.Viscosity_mu;
purcell_eff=stokes_drag_coeffs.*viscosity_sweep_data.Speed_U_per_freq.^2./(viscosity_sweep_data.Torque_Tau_per_freq*2*pi);
plot(viscosity_sweep_data.Viscosity_mu,purcell_eff)
%evidently purcell efficiency even with the new definition is independent of
%viscosity to within working precision

%% Speed vs. viscosity ratio

max_power=2*pi*1200*200;

viscosity_sweep_data=table2array(readtable('../../data/viscosity_sweep/summary/viscosity_sweep_table.txt'));
viscosity_col=viscosity_sweep_data(:,7);
speed_col_visc=viscosity_sweep_data(:,6); %should just be the same entry repeated down the column, but don't want to hard code that
torque_col_visc=viscosity_sweep_data(:,5);

water_visc_index=find(isapprox(viscosity_col,0.93));

speed_on_power_curve_visc=speed_col_visc.*sqrt(max_power/2/pi./torque_col_visc);
speed_on_power_curve_ratio_visc_sim=(speed_on_power_curve_visc/speed_on_power_curve_visc(water_visc_index));
speed_on_power_curve_ratio_visc_theory=1./sqrt(viscosity_col/viscosity_col(water_visc_index));

visc_ratio=viscosity_col/viscosity_col(water_visc_index);

%compare speeds at various viscosities on the power curve
figure
hold on
scatter(visc_ratio,speed_on_power_curve_ratio_visc_sim,'filled','diamond')
plot(visc_ratio,speed_on_power_curve_ratio_visc_theory)
xscale log
xlabel('\eta_b/\eta_a')
ylabel('U_b/U_a')
legend('Simulation','Theory')
exportgraphics(gca,[plotdir,'/speed_ratio_vs_visc_ratio_on_power_curve.png'])


speed_off_power_curve_ratio_visc_sim=speed_col_visc/speed_col_visc(water_visc_index)*2.*sqrt(torque_col_visc*torque_col_visc(water_visc_index))./(torque_col_visc+torque_col_visc(water_visc_index));
speed_off_power_curve_ratio_visc_theory=2*sqrt(visc_ratio)./(1+visc_ratio);
visc_ratio_inverse=1./visc_ratio;

%viscosity is fixed at eta_a, tuning not fixed. U_a: speed if tuned to eta_a (on power
%curve). U_c: speed on eta_a load line but if tuned to eta_c.
figure
hold on
scatter(visc_ratio,speed_off_power_curve_ratio_visc_sim,'filled','diamond')
plot(visc_ratio,speed_off_power_curve_ratio_visc_theory)
xscale log
xlabel('\eta_c/\eta_a')
ylabel('U_c/U_a')
legend('Simulation','Theory')
exportgraphics(gca,[plotdir,'/speed_ratio_vs_visc_ratio_tuning_comparison_same_viscosity.png'])

%viscosity not fixed, tuning fixed to best for eta_a. U_a: speed at eta_a with eta_a load line. U_c: speed at
%eta_c with eta_c load line.
figure
hold on
scatter(visc_ratio,speed_off_power_curve_ratio_visc_sim.*sqrt(visc_ratio_inverse),'filled','diamond')
plot(visc_ratio,speed_off_power_curve_ratio_visc_theory.*sqrt(visc_ratio_inverse))
xscale log
xlabel('\eta_c/\eta_a')
ylabel('U_c/U_a')
legend('Simulation','Theory')
exportgraphics(gca,[plotdir,'/speed_ratio_vs_visc_ratio_power_and_viscosity_comparison.png'])
%%

%this compares speed on power curve to speed on ts curve at other tuning,
%at single viscosity
mu_ratio=linspace(0.01,10,100);
speed_ratio=2*sqrt(mu_ratio)./(1+mu_ratio); %this function looks a lot like something lognormal, but the tail ends up slightly heavier

mu_ratio_log=logspace(-3,3,100);
speed_ratio_log=2*sqrt(mu_ratio_log)./(1+mu_ratio_log);

mu_ratio_axis_label='\eta_0/\eta_m';
speed_ratio_axis_label='U_0/U_m';

figure
plot(mu_ratio,speed_ratio)
xlabel(mu_ratio_axis_label)
ylabel(speed_ratio_axis_label)
exportgraphics(gca,[plotdir,'/speed_vs_viscosity_ratio.png'])

figure
plot(mu_ratio_log,speed_ratio_log)
xlabel(mu_ratio_axis_label)
ylabel(speed_ratio_axis_label)
xscale log
exportgraphics(gca,[plotdir,'/logx_speed_vs_viscosity_ratio.png'])


