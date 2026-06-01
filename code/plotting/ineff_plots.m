full_sweep_data=table2array(readtable("../../data/full_sweep_table.txt"));
forces=full_sweep_data(:,4);
torques=full_sweep_data(:,5);
translational_speeds=full_sweep_data(:,6);

dir='../../reports/figures/ineff_plots';
if(~exist(dir,'dir'))
    mkdir(dir)
end

%%
purcell_ineff=2*pi*torques./(forces.*translational_speeds); %note these are all "per frequency" quantities, the 2pi accounts for the spare omega
[min_ineff,min_ineff_index]=min(purcell_ineff);

radius_indices=find(full_sweep_data(:,1)==full_sweep_data(min_ineff_index,1) & full_sweep_data(:,2)==full_sweep_data(min_ineff_index,2));
radii=full_sweep_data(radius_indices,3);

wavelength_indices=find(full_sweep_data(:,1)==full_sweep_data(min_ineff_index,1) & full_sweep_data(:,3)==full_sweep_data(min_ineff_index,3));
wavelengths=full_sweep_data(wavelength_indices,2);

arclen_indices=find(full_sweep_data(:,2)==full_sweep_data(min_ineff_index,2) & full_sweep_data(:,3)==full_sweep_data(min_ineff_index,3));
arclengths=full_sweep_data(arclen_indices,1);

figure
plot(arclengths,purcell_ineff(arclen_indices))
xlabel('Arc length (\mu m)')
ylabel('Purcell inefficiency')
exportgraphics(gca,[dir,'/ineff_vs_arc_length.png'])

figure
plot(wavelengths,purcell_ineff(wavelength_indices))
xlabel('Wavelength (\mu m)')
ylabel('Purcell inefficiency')
exportgraphics(gca,[dir,'/ineff_vs_wavelength.png'])

figure
plot(radii,purcell_ineff(radius_indices))
xlabel('Helical radius (\mu m)')
ylabel('Purcell inefficiency')
exportgraphics(gca,[dir,'/ineff_vs_radius.png'])

%%
torque_array=linspace(min(torques),max(torques),50);
bundle_ineff=zeros(length(torque_array),1);

for i=1:length(torque_array)
    nearby_torque_indices=find(abs(torques-torque_array(i))<min(diff(torque_array))/4);
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
xlabel("Slope of torque-speed curve (pN\cdot nm/Hz)")
ylabel("Local minimum of Purcell inefficiency")
exportgraphics(gca,[dir,'/ineff_vs_torque-speed_slope.png'])

figure
plot(torque_array(4:end),bundle_ineff(4:end))
xlabel("Slope of torque-speed curve (pN\cdot nm/Hz)")
ylabel("Local minimum of Purcell inefficiency")
exportgraphics(gca,[dir,'/ineff_vs_torque-speed_slope_zoomed.png'])