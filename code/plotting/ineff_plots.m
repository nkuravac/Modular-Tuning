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
arclen_torque_slopes=full_sweep_data(arclen_indices,5);
arclen_load_lines=[arclen_torque_slopes,purcell_ineff(arclen_indices)];
arclen_load_lines=sortrows(arclen_load_lines,1,"ascend");

figure
plot(arclen_load_lines(:,1),arclen_load_lines(:,2))
xlabel('Load line slope (pN\cdot nm/Hz)')
ylabel('Purcell inefficiency')
exportgraphics(gca,[dir,'/ineff_vs_load_line_slope_arclen_only.png'])

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
num_torques=[20,30,50,70,100];
bundle_width=[0.1,0.25,0.5,0.9];

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
