full_sweep_data=table2array(readtable("../../data/full_sweep_table.txt"));
forces=full_sweep_data(:,4);
torques=full_sweep_data(:,5);
translational_speeds=full_sweep_data(:,6);

dir='../../reports/figures/energy_per_dist_plots';
if(~exist(dir,'dir'))
    mkdir(dir)
else
    rmdir(dir,'s')
    mkdir(dir)
end

%%
freq_array=linspace(1,500,100);
torque_array=linspace(1,4000,100);

energy_per_dist=zeros(length(freq_array),length(torque_array));

for i=1:length(freq_array)
    for j=1:length(torque_array)
       local_torque_indices=find(abs(torques)-torque_array(j)/freq_array(i)<0.1);
       local_energies=torques(local_torque_indices)*2*pi*freq_array(i)./translational_speeds(local_torque_indices);
       if ~isempty(local_energies)
        energy_per_dist(i,j)=min(local_energies);
       else
        energy_per_dist(i,j)=NaN;
       end
    end
end

%energy heatmap
figure
energy_heatmap=heatmap(freq_array,torque_array,energy_per_dist,"MissingDataColor",'r','Colormap',parula);
xlabel('Speed (Hz)')
ylabel('Torque (pN\cdot nm)')
heatmap_y_labels=strings(length(torque_array),1);
for i=1:length(heatmap_y_labels)
    if i==1
        heatmap_y_labels(i)=num2str(torque_array(i));
    elseif mod(i,10)==0
        heatmap_y_labels(i)=num2str(torque_array(i));
    else
        heatmap_y_labels(i)='';
    end
end
energy_heatmap.YDisplayLabels=heatmap_y_labels;
energy_heatmap.YDisplayData=flip(energy_heatmap.YDisplayData);

heatmap_x_labels=strings(length(freq_array),1);
for i=1:length(heatmap_x_labels)
    if i==1
        heatmap_x_labels(i)=num2str(freq_array(i));
    elseif mod(i,10)==0
        heatmap_x_labels(i)=num2str(freq_array(i));
    else
        heatmap_x_labels(i)='';
    end
end
energy_heatmap.XDisplayLabels=heatmap_x_labels;
exportgraphics(gca,[dir,'/energy_per_dist_heatmap.png'])

%contour plot
figure
contour(freq_array,torque_array,energy_per_dist)
xlabel('Speed (Hz)')
ylabel('Torque (pN\cdot nm)')
exportgraphics(gca,[dir,'/energy_per_dist_contour_plot.png'])