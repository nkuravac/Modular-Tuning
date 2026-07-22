close all
clearvars

set(groot,'DefaultAxesFontSize',14)

%%
plotdir='../../reports/figures/all_parameter_efficiency_plots';
if ~exist(plotdir,'dir')
    mkdir(plotdir)
end

%%
full_sweep_data=table2array(readtable('../../data/sweep_near_global_max_sphere/summary/3_parameter_1d_sweeps.txt'));
opt_ineff_struct=jsondecode(fileread('../../data/ineff_optimize_sphere/summary/optimal_parameters_purcell_ineff.json'));
opt_config_parameters=jsondecode(fileread('../../data/ineff_optimize_sphere/summary/config_parameters.json'));
new_config_parameters=jsondecode(fileread('../../data/sweep_near_global_max_sphere/summary/config_parameters.json'));
row_counts=jsondecode(fileread('../../data/sweep_near_global_max_sphere/summary/row_counts.json'));

mu_old=opt_config_parameters.Viscosity;
r_old=opt_config_parameters.r;
mu_new=new_config_parameters.mu;
r_new=new_config_parameters.r;

arclen_indices=1:row_counts.num_arclens;
wavelength_indices=(row_counts.num_arclens+1):row_counts.num_arclens+row_counts.num_wl;
hrad_indices=(row_counts.num_arclens+row_counts.num_wl+1):length(full_sweep_data(:,1));

arclengths=full_sweep_data(:,1);
wavelengths=full_sweep_data(:,2);
radii=full_sweep_data(:,3);
forces=full_sweep_data(:,4);
torques=full_sweep_data(:,5);
translational_speeds=full_sweep_data(:,6);

stokes_drags=6*pi*mu_new*r_new*translational_speeds;
max_power=2*pi*1200*200;

speeds_on_power_curve=translational_speeds.*sqrt(max_power/2/pi./torques);
purcell_efficiencies=stokes_drags.*translational_speeds./(torques*2*pi);
miscalc_purcell_efficiencies=forces.*translational_speeds./(torques*2*pi);
[max_eff,opt_eff_index]=max(purcell_efficiencies);
[max_speed,max_speed_index]=max(speeds_on_power_curve);

%TESTING (arclength curve is so flat that a very slight numerical error shifts the index over by one compared to what we know should be true from the optimization script)
max_speed_index=19;
opt_eff_index=19;


%%

arclength_axis_label='Arc length ({\mu}m)';
wavelength_axis_label='Wavelength ({\mu}m)';
hrad_axis_label='Helix radius ({\mu}m)';
geo_parameter_labels={arclength_axis_label,wavelength_axis_label,hrad_axis_label};
geo_parameter_plot_snippets={'S','wl','R'};
speed_axis_label='Speed ({\mu}m/s)';
efficiency_axis_label='Purcell efficiency';
ll_axis_label='Load line slope (pN{\cdot}nm{\cdot}s)';
indices={arclen_indices,wavelength_indices,hrad_indices};

Omega_0=2*sqrt(torques*max_power/(2*pi))./(torques(opt_eff_index)+torques);
Omega_m=sqrt(max_power/2/pi./torques);

U_0=translational_speeds(opt_eff_index)*Omega_0;
U_m=translational_speeds.*Omega_m;

m_ratio=torques(opt_eff_index)./torques;
Omega_ratio=Omega_m./Omega_0;
U_ratio=U_m./U_0;

power_function_omega=(2-1./Omega_ratio)./Omega_ratio;
power_function_m=4*m_ratio./(1+m_ratio).^2;
U_ratio_omega_theory=1./sqrt(power_function_omega*max_eff./purcell_efficiencies);
U_ratio_m_theory=1./sqrt(power_function_m*max_eff./purcell_efficiencies);
%%

figure
hold on
markers={'square','diamond','o'};
colors={'r','b','g'};
for i=1:3
    scatter(m_ratio(indices{i}),U_ratio(indices{i}),'filled',markers{i},'MarkerFaceColor',colors{i})
    plot(m_ratio(indices{i}),U_ratio_m_theory(indices{i}),'Color',colors{i})
end
xlabel('m_0/m_1')
ylabel('U_1/U_0')
legend('Simulation (S)','Theory (S)','Simulation (\lambda)','Theory (\lambda)','Simulation (R)','Theory (R)')
xline(1,'LineStyle','--','HandleVisibility','off')
xscale log
hold off
exportgraphics(gca,[plotdir,'/speed_ratio_vs_m_ratio_all_parameters.png'])
savefig(gcf,[plotdir,'/speed_ratio_vs_m_ratio_all_parameters.fig'])

%%
figure
hold on
line_styles={'-','--','-.'};
for i=1:3
    plot(full_sweep_data(indices{i},i)/full_sweep_data(opt_eff_index,i),purcell_efficiencies(indices{i}),'LineWidth',1,'LineStyle',line_styles{i})
end
xlabel('Fractional deviation from optimal value')
ylabel(efficiency_axis_label)
legend('S','\lambda','R')
xline(1,'LineStyle','--','Color','k','HandleVisibility','off')
exportgraphics(gca,[plotdir,'/eff_vs_fractional_deviation_all_parameters.png'])
savefig(gcf,[plotdir,'/eff_vs_fractional_deviation_all_parameters.fig'])


%% 9-panel plot

figure

for i=1:3
    
    subplot(3,3,i)
    plot(full_sweep_data(indices{i},i),speeds_on_power_curve(indices{i}));
    xlabel(geo_parameter_labels{i})
    xline(full_sweep_data(opt_eff_index,i),'Label','Max \epsilon','Linestyle','--','LabelVerticalAlignment','bottom','HandleVisibility','off')
    if i==1
        ylabel(speed_axis_label)
    end
    ylim([0,40])
    %exportgraphics(gca,[plotdir,'/speed_vs_',geo_parameter_plot_snippets{i},'.png'])

    %figure
    subplot(3,3,i+6)
    hold on
    %plot(full_sweep_data(indices{i},i),miscalc_purcell_efficiencies(indices{i}))
    plot(full_sweep_data(indices{i},i),purcell_efficiencies(indices{i}))
    xline(full_sweep_data(max_speed_index,i),'Label','Max U','Linestyle','--','LabelVerticalAlignment','bottom','HandleVisibility','off')
    xlabel(geo_parameter_labels{i})
    if i==1
        ylabel(efficiency_axis_label)
    end
    ylim([0,0.01])
    %legend('\epsilon (simulated drag)','\epsilon (Stokes drag)')
    %exportgraphics(gca,[plotdir,'/eff_vs_',geo_parameter_plot_snippets{i},'.png'])

    %figure
    subplot(3,3,i+3)
    hold on
    freq_array=linspace(10,1350,100);
    plot(freq_array,max_power/2/pi./freq_array,'Color','k');
    torque_section=torques(indices{i});
    for j=1:length(torque_section)
        if j==1
            handle_toggle='on';
        else
            handle_toggle='off';
        end
        freq_points=[0,sqrt(max_power/2/pi/torque_section(j))];
        plot(freq_points,torque_section(j)*freq_points,'HandleVisibility',handle_toggle,'Color','b')
    end
    hold off
    xlabel('Frequency (Hz)')
    if i==1
        ylabel('Torque (pN{\cdot}nm)')
    end
    %legend('Fixed power curve','Load lines')
    ylim([0,2000])
    %exportgraphics(gca,[plotdir,'/all_load_lines_',geo_parameter_plot_snippets{i},'.png'])
end
exportgraphics(gca,[plotdir,'/9_panel_speed_eff_and_ll_all_parameters.png'])
savefig(gcf,[plotdir,'/9_panel_speed_eff_and_ll_all_parameters.fig'])

%% arc length plot alone for ppt
ppt_plot_font_size=24; %points

figure
plot(full_sweep_data(arclen_indices,1),speeds_on_power_curve(arclen_indices),'LineWidth',1);
fontsize(ppt_plot_font_size,'points')
xlabel(geo_parameter_labels{1})
xline(full_sweep_data(opt_eff_index,1),'Label','Maximum Purcell efficiency','Linestyle','--','LabelVerticalAlignment','middle','HandleVisibility','off','FontSize',get(gca, 'FontSize'))
ylabel(speed_axis_label)
ylim([0,40])
exportgraphics(gca,[plotdir,'/speed_vs_S.png'])
savefig(gcf,[plotdir,'/speed_vs_S.fig'])

figure
hold on
scatter(m_ratio(indices{1}),U_ratio(indices{1}),'filled',markers{1},'MarkerFaceColor','b')
plot(m_ratio(indices{1}),U_ratio_m_theory(indices{1}),'Color','b')
fontsize(ppt_plot_font_size,'points')
xlabel('m_0/m_1')
ylabel('U_1/U_0')
legend('Simulation','Theory')
xline(1,'LineStyle','--','HandleVisibility','off')
xscale log
exportgraphics(gca,[plotdir,'/logx_speed_ratio_vs_m_ratio_',geo_parameter_plot_snippets{1},'_change_big_font.png'])
savefig(gcf,[plotdir,'/logx_speed_ratio_vs_m_ratio_',geo_parameter_plot_snippets{1},'_change_big_font.fig'])

figure
hold on
plot(full_sweep_data(indices{1},1),miscalc_purcell_efficiencies(indices{1}),'LineWidth',1)
plot(full_sweep_data(indices{1},1),purcell_efficiencies(indices{1}),'LineWidth',1)
fontsize(ppt_plot_font_size,'points')
xline(full_sweep_data(max_speed_index,1),'Label','Max U','Linestyle','--','LabelVerticalAlignment','bottom','HandleVisibility','off','FontSize',get(gca,'FontSize'))
xlabel(geo_parameter_labels{1})
ylabel(efficiency_axis_label)
legend('\epsilon (simulated drag)','\epsilon (Stokes drag)')
exportgraphics(gca,[plotdir,'/eff_vs_',geo_parameter_plot_snippets{1},'_with_miscalc_eff.png'])
savefig(gcf,[plotdir,'/eff_vs_',geo_parameter_plot_snippets{1},'_with_miscalc_eff.fig'])



%%
for i=1:3
    
    for repeat=1:2
        figure
        hold on
        plot(Omega_ratio(indices{i}),U_ratio(indices{i}),'LineStyle','none','Marker','+')
        plot(Omega_ratio(indices{i}),U_ratio_omega_theory(indices{i}))
        xlabel('\Omega_m/\Omega_0')
        ylabel('U_m/U_0')
        legend('Simulation','Theory')
        xline(1,'LineStyle','--','HandleVisibility','off')
        hold off
        if repeat==1
            exportgraphics(gca,[plotdir,'/speed_ratio_vs_freq_ratio_',geo_parameter_plot_snippets{i},'_change.png'])
            savefig(gcf,[plotdir,'/speed_ratio_vs_freq_ratio_',geo_parameter_plot_snippets{i},'_change.fig'])
        else
            xscale log
            exportgraphics(gca,[plotdir,'/logx_speed_ratio_vs_freq_ratio_',geo_parameter_plot_snippets{i},'_change.png'])
            savefig(gcf,[plotdir,'/logx_speed_ratio_vs_freq_ratio_',geo_parameter_plot_snippets{i},'_change.fig'])
        end

        figure
        hold on
        plot(m_ratio(indices{i}),U_ratio(indices{i}),'LineStyle','none','Marker','+')
        plot(m_ratio(indices{i}),U_ratio_m_theory(indices{i}))
        xlabel('m_0/m_1')
        ylabel('U_1/U_0')
        legend('Simulation','Theory')
        xline(1,'LineStyle','--','HandleVisibility','off')
        hold off
        if repeat==1
            exportgraphics(gca,[plotdir,'/speed_ratio_vs_m_ratio_',geo_parameter_plot_snippets{i},'_change.png'])
            savefig(gcf,[plotdir,'/speed_ratio_vs_m_ratio_',geo_parameter_plot_snippets{i},'_change.fig'])
        else
            xscale log
            exportgraphics(gca,[plotdir,'/logx_speed_ratio_vs_m_ratio_',geo_parameter_plot_snippets{i},'_change.png'])
            savefig(gcf,[plotdir,'/logx_speed_ratio_vs_m_ratio_',geo_parameter_plot_snippets{i},'_change.fig'])
        end

    end
    

    figure
    plot(full_sweep_data(indices{i},i),torques(indices{i}))
    xlabel(geo_parameter_labels{i})
    ylabel(ll_axis_label)
    exportgraphics(gca,[plotdir,'/ll_slope_vs_',geo_parameter_plot_snippets{i},'.png'])
    savefig(gcf,[plotdir,'/ll_slope_vs_',geo_parameter_plot_snippets{i},'.fig'])

end

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
savefig(gcf,[plotdir,'/speed_ratio_vs_visc_ratio_on_power_curve.fig'])


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
savefig(gcf,[plotdir,'/speed_ratio_vs_visc_ratio_tuning_comparison_same_viscosity.fig'])

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
savefig(gcf,[plotdir,'/speed_ratio_vs_visc_ratio_power_and_viscosity_comparison.fig'])