close all
clearvars


%%
%full_sweep_data=table2array(readtable("../../data/fuller_sweep_table.txt"));
%modified back to PA
%Modified back for checks
%Modified for PA

full_sweep_data=table2array(readtable("../../data/arclen_sweep_PA/summary/arclen_sweep_table_PA.txt"));
%full_sweep_data=table2array(readtable("../../data/fuller_sweep_table.txt"));

forces=full_sweep_data(:,4);
torques=full_sweep_data(:,5);
translational_speeds=full_sweep_data(:,6);

%%
dir='../../reports/figures/PA/ineff_plots';
if ~exist(dir,'dir')
    mkdir(dir)
end

%%
purcell_ineff=2*pi*torques./(forces.*translational_speeds); %note these are all "per frequency" quantities, the 2pi accounts for the spare omega
%[min_ineff,min_ineff_index]=min(purcell_ineff); %this doesn't quite work
%because of very slight numerical issues
opt_ineff_struct=jsondecode(fileread('../../data/optimal_parameters_PA_purcell_ineff.json'));


opt_arclen_purcell=opt_ineff_struct.Arc_length;
opt_wl_purcell=opt_ineff_struct.Wavelength;
opt_hrad_purcell=opt_ineff_struct.Helix_radius;
min_ineff=opt_ineff_struct.Minimum_Purcell_inefficiency;
min_ineff_index=find(isapprox(full_sweep_data(:,1),opt_arclen_purcell,'tight') & isapprox(full_sweep_data(:,2),opt_wl_purcell,'tight') & isapprox(full_sweep_data(:,3),opt_hrad_purcell,'tight'));
%min_ineff_index=16; %16 makes the plot work, but 10 is correct...

%%
radius_indices=find(full_sweep_data(:,1)==full_sweep_data(min_ineff_index,1) & full_sweep_data(:,2)==full_sweep_data(min_ineff_index,2));
radii=full_sweep_data(radius_indices,3);

wavelength_indices=find(full_sweep_data(:,1)==full_sweep_data(min_ineff_index,1) & full_sweep_data(:,3)==full_sweep_data(min_ineff_index,3));
wavelengths=full_sweep_data(wavelength_indices,2);

arclen_indices=find(full_sweep_data(:,2)==full_sweep_data(min_ineff_index,2) & full_sweep_data(:,3)==full_sweep_data(min_ineff_index,3));
arclengths=full_sweep_data(arclen_indices,1);
arclen_torque_slopes=torques(arclen_indices);
arclen_speed_slopes=translational_speeds(arclen_indices);
arclen_load_lines=[arclen_torque_slopes,purcell_ineff(arclen_indices)];
arclen_load_lines=sortrows(arclen_load_lines,1,"ascend");

%%  Don't run this section unless using full sweep data
%{
figure
plot(arclen_load_lines(:,1),arclen_load_lines(:,2))
xlabel('Load line slope (pN{\cdot}nm{\cdot}s)')
ylabel('Purcell inefficiency')
exportgraphics(gca,[dir,'/ineff_vs_load_line_slope_arclen_only.png'])

figure
plot(arclengths,purcell_ineff(arclen_indices))
xlabel('Arc length ({\mu}m)')
ylabel('Purcell inefficiency')
exportgraphics(gca,[dir,'/ineff_vs_arc_length.png'])

figure
plot(wavelengths,purcell_ineff(wavelength_indices))
xlabel('Wavelength ({\mu}m)')
ylabel('Purcell inefficiency')
exportgraphics(gca,[dir,'/ineff_vs_wavelength.png'])

figure
plot(radii,purcell_ineff(radius_indices))
xlabel('Helical radius ({\mu}m)')
ylabel('Purcell inefficiency')
exportgraphics(gca,[dir,'/ineff_vs_radius.png'])

figure
plot(arclengths,arclen_torque_slopes)
xlabel('Arc length ({\mu}m)')
ylabel('Load line slope (pN{\cdot}nm{\cdot}s)')
exportgraphics(gca,[dir,'/ll_slope_vs_arc_length.png'])
%}

%% Ineff on torque bundles (irrelevant if restricting to arc length)
%{
%num_torques=[20,30,50,70,100];
%bundle_width=[0.1,0.25,0.5,0.9];
num_torques=30;
bundle_width=0.5;

torque_bundle_dir=[dir,'/ineff_vs_torque_slope'];
if(~exist(torque_bundle_dir,'dir'))
    mkdir(torque_bundle_dir)
end

for num_torques_index=1:length(num_torques)
for bundle_width_index=1:length(bundle_width)

torque_array=linspace(min(torques),max(torques),num_torques(num_torques_index));
bundle_ineff=zeros(length(torque_array),1);

for i=1:length(torque_array)
    nearby_torque_indices=find(abs(torques-torque_array(i))<min(diff(torque_array))*bundle_width(bundle_width_index)/2);
    %want to take min ineff over "bundle" of torques, but don't want
    %overlap, hence restricting bundle width to 1/2 (that is, 2*1/4) of torque array
    %discretization width.
    if ~isempty(nearby_torque_indices)
        bundle_ineff(i)=min(purcell_ineff(nearby_torque_indices));
    else
        bundle_ineff(i)=NaN;
    end
end

figure
plot(torque_array,bundle_ineff)
xlabel("Load line slope (pN\cdot nm/Hz)")
ylabel("Local minimum of Purcell inefficiency")
exportgraphics(gca,[torque_bundle_dir,'/ineff_vs_load_line_slope',num2str(num_torques(num_torques_index)),'_points_',num2str(bundle_width(bundle_width_index)),'bundle_width.png'])

figure
zoom_cutoff=round(num_torques(num_torques_index)/5);
plot(torque_array(zoom_cutoff:end),bundle_ineff(zoom_cutoff:end))
xlabel("Load line slope (pN\cdot nm/Hz)")
ylabel("Local minimum of Purcell inefficiency")
exportgraphics(gca,[torque_bundle_dir,'/zoomed_ineff_vs_load_line_slope',num2str(num_torques(num_torques_index)),'_points_',num2str(bundle_width(bundle_width_index)),'bundle_width.png'])

end
end
%}

%% Compare effects of geometry off optimal with tuning off optimal

power=2*pi*150*1250; %this is maybe a typical power, but I did check to make sure this does not change the shape of the plot whatsoever

good_geo_speeds=zeros(length(arclengths),1);
good_motor_speeds=good_geo_speeds;

fixed_power_curve_speeds=translational_speeds(arclen_indices).*sqrt(power./(2*pi*torques(arclen_indices)));
[max_speed_on_curve_mu,max_speed_on_curve_index_mu]=max(fixed_power_curve_speeds)
[min_ineff_on_curve,min_ineff_on_curve_index]=min(purcell_ineff(arclen_indices))
%PA: index 16 lines up with plot but not with min_ineff at 11. 
%EC: same but 24, 17
for repeat=1:2

    if repeat==1
        fastest_torque_slope=torques(min_ineff_on_curve_index);
        fastest_speed_slope=translational_speeds(min_ineff_on_curve_index);
        %max speed 'should' be just speed at optimal purcell geo
        plotname=[dir,'geo_vs_tuning_assume_purcell_has_max_speed.png'];
    else
        fastest_torque_slope=torques(max_speed_on_curve_index_mu);
        fastest_speed_slope=translational_speeds(max_speed_on_curve_index_mu);
        %artificially impose that 'best' geometry has the real max speed
        %across geometries even if it isn't optimal for purcell efficiency
        plotname=[dir,'geo_vs_tuning_real_max_speed.png'];
    end


    for i=1:length(arclengths)
        good_motor_speeds(i)=arclen_speed_slopes(i)*sqrt(power/(2*pi*arclen_torque_slopes(i)));
        good_geo_speeds(i)=fastest_speed_slope*2*sqrt(arclen_torque_slopes(i)*power/(2*pi))/(fastest_torque_slope+arclen_torque_slopes(i));
    end
    
    max_speed=fastest_speed_slope*sqrt(power/(2*pi*fastest_torque_slope));

    good_motor_speeds=good_motor_speeds/max_speed;
    good_geo_speeds=good_geo_speeds/max_speed;

    figure
    hold on
    plot(arclengths,good_motor_speeds)
    plot(arclengths,good_geo_speeds)
    xlabel('Arc length for motor tuning ({\mu}m)')
    ylabel('Normalized translational speed')
    legend('Geometry tuned to motor','Optimal Purcell efficiency geometry','Location','southeast')
    yline(1,'LineStyle','--','HandleVisibility','off')
    exportgraphics(gca,plotname)
end

%%check to make sure power is irrelevant
%{

%power=[0.5;1;1.5]*2*pi*150*1250; %I'm not sure if this matters or not so checking - turns out does not matter.
power=2*pi*150*1250;
max_speed=zeros(length(power),1);

good_geo_speeds=zeros(length(arclengths),length(power));
good_motor_speeds=good_geo_speeds;

purcell_opt_torque_slope=torques(min_ineff_index);
purcell_opt_speed_slope=translational_speeds(min_ineff_index);


for j=1:length(power)
    for i=1:length(arclengths)
        good_motor_speeds(i,j)=arclen_speed_slopes(i)*sqrt(power(j)/(2*pi*arclen_torque_slopes(i)));
        good_geo_speeds(i,j)=purcell_opt_speed_slope*2*sqrt(arclen_torque_slopes(i)*power(j)/(2*pi))/(purcell_opt_torque_slope+arclen_torque_slopes(i));
    end
    max_speed(j)=purcell_opt_speed_slope*sqrt(power(j)/(2*pi*purcell_opt_torque_slope));
    good_motor_speeds(:,j)=good_motor_speeds(:,j)/max_speed(j);
    good_geo_speeds(:,j)=good_geo_speeds(:,j)/max_speed(j);

    figure
    hold on
    plot(arclengths,good_motor_speeds(:,j))
    plot(arclengths,good_geo_speeds(:,j))
    xlabel('Arc length for motor tuning ({\mu}m)')
    ylabel('Normalized translational speed')
    legend('Geometry tuned to motor','Optimal Purcell efficiency geometry')
    yline(1,'LineStyle','--','HandleVisibility','off')
    exportgraphics(gca,[dir,'/geo_vs_tuning.png'])
end
%}

%% compare fixed motor but various viscosity to various motor but fixed viscosity

viscosity_sweep_data=table2array(readtable("../../data/viscosity_sweep_PA/summary/viscosity_sweep_table_PA.txt"));
forces_mu=viscosity_sweep_data(:,4);
torques_mu=viscosity_sweep_data(:,5);
translational_speeds_mu=viscosity_sweep_data(:,6);
mu_array=viscosity_sweep_data(:,7);
base_mu_index=ceil(length(mu_array)/2); %'middle' element, round up
base_mu=mu_array(base_mu_index);


power=2*pi*150*1250; %this is maybe a typical power, but I did check to make sure this does not change the shape of the plot whatsoever

base_torque_slope=torques_mu(base_mu_index);
base_speed_slope=translational_speeds_mu(base_mu_index);

fixed_viscosity_speeds=base_speed_slope*2*sqrt(torques_mu*power/2/pi)./(torques_mu+base_torque_slope);
fixed_motor_speeds=translational_speeds_mu.*2*sqrt(base_torque_slope*power/2/pi)./(torques_mu+base_torque_slope);
fixed_motor_speeds_normed=translational_speeds_mu.*2*sqrt(base_torque_slope*power/2/pi)./(torques_mu+base_torque_slope).*sqrt(mu_array/base_mu);
%scale by viscosity

max_speed=base_speed_slope*sqrt(power/(2*pi*base_torque_slope));

fixed_viscosity_speeds=fixed_viscosity_speeds/max_speed;
fixed_motor_speeds=fixed_motor_speeds/max_speed;
fixed_motor_speeds_normed=fixed_motor_speeds_normed/max_speed;

figure
hold on
plot(mu_array,fixed_viscosity_speeds)
plot(mu_array,fixed_motor_speeds)
xlabel('Viscosity (mPa{\cdot}s)')
ylabel('Normalized translational speed')
legend('Tuned motor, constant load line','Constant motor, variable load line','Location','southeast')
yline(1,'LineStyle','--','HandleVisibility','off')
exportgraphics(gca,[dir,'/speed_vs_mu.png'])

figure
hold on
plot(mu_array,fixed_viscosity_speeds)
plot(mu_array,fixed_motor_speeds_normed)
xlabel('Viscosity (mPa{\cdot}s)')
ylabel('Normalized translational speed')
legend('Tuned motor, constant load line','Constant motor, variable load line, normalized by $\sqrt{\mu/\mu_1}$','Interpreter','latex','Location','southeast')
yline(1,'LineStyle','--','HandleVisibility','off')
exportgraphics(gca,[dir,'/speed_vs_mu_normed.png'])

%%





