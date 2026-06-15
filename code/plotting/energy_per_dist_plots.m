clearvars
close all

full_sweep_data=table2array(readtable("../../data/fuller_sweep_table.txt"));
forces=full_sweep_data(:,4);
torques=full_sweep_data(:,5);
translational_speeds=full_sweep_data(:,6);

dir='../../reports/figures/energy_per_dist_plots';
if(~exist(dir,'dir'))
    mkdir(dir)
end
%% plot all load lines
opt_ineff_row=readmatrix("../../data/opt_ineff_row.txt");
opt_load_line_slope=opt_ineff_row(5);
[~,opt_ineff_index]=min(abs(torques-opt_load_line_slope));

freq_array=linspace(1,600,3);
load_lines=repmat(freq_array',1,length(torques))*diag(torques);

figure
hold on
for i=1:length(torques)
    if i==1
        plot(freq_array,load_lines(:,i),'Color','b','LineStyle','-')
    elseif i==opt_ineff_index
        plot(freq_array,load_lines(:,i),'Color','r','LineStyle','-','LineWidth',1,'ZData',ones(length(freq_array),1))
    else
        plot(freq_array,load_lines(:,i),'Color','b','LineStyle','-','HandleVisibility','off')
    end
end
legendlabels={'Load lines','Purcell inefficiency minimizing load line'};
legend(legendlabels,'Location','northwest')
xlabel('Speed (Hz)')
ylabel('Torque (pN\cdot nm)')
exportgraphics(gca,[dir,'/all_load_lines.png'])
hold off


%% Make energy heatmap in torque-speed space from closest point among "bundle"
freq_array=linspace(1,500,100);
torque_array=linspace(1,4000,100);

energy_per_dist=zeros(length(freq_array),length(torque_array));

for i=1:length(freq_array)
    for j=1:length(torque_array)
       local_torque_indices=find(abs(torques-torque_array(j)/freq_array(i))<0.2);
       if ~isempty(local_torque_indices)
        local_energies=torques(local_torque_indices)*2*pi*freq_array(i)./translational_speeds(local_torque_indices);
        energy_per_dist(i,j)=min(local_energies);
       else
        energy_per_dist(i,j)=NaN;
       end
    end
end

%energy heatmap
figure
energy_heatmap=heatmap(freq_array,torque_array,energy_per_dist',"MissingDataColor",'r','Colormap',parula);
%energy_per_dist is transposed since matlab puts matrix rows as heatmap
%rows, which makes the x-coordinate the column index rather than the row
%index which would fit the code better (row,column) -> (y,x)

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
contour(freq_array,torque_array,energy_per_dist')
xlabel('Speed (Hz)')
ylabel('Torque (pN\cdot nm)')
exportgraphics(gca,[dir,'/energy_per_dist_contour_plot.png'])

%% Make energy heatmap in torque-speed space only varying arc length, other two parameters fixed at purcell optimal values

%find minimum inefficiency location in full sweep data
purcell_ineff=2*pi*torques./(forces.*translational_speeds); %note these are all "per frequency" quantities, the 2pi accounts for the spare omega
[min_ineff,min_ineff_index]=min(purcell_ineff);

%find where wavelength and radius are at the minimum inefficiency values
arclen_indices=find(full_sweep_data(:,2)==full_sweep_data(min_ineff_index,2) & full_sweep_data(:,3)==full_sweep_data(min_ineff_index,3));
arclengths=full_sweep_data(arclen_indices,1);

%freq_array=linspace(1,500,100); %duplicated above

temp_tau_per_U=repmat(torques(arclen_indices)'./translational_speeds(arclen_indices)',length(freq_array),1);
energy_per_dist_arclen=2*pi*diag(freq_array)*temp_tau_per_U;

%{
% More understandable but slower version of the above two lines:
energy_per_dist_arclen=zeros(length(freq_array),length(arclen_indices));

for i=1:length(freq_array)
    for j=1:length(arclen_indices)
        energy_per_dist_arclen(i,j)=torques(arclen_indices(j))*2*pi*freq_array(i)./translational_speeds(arclen_indices(j));
    end
end
%}

freq_column=reshape(repmat(freq_array',1,length(arclen_indices))',[],1); %blocks of constant freq
torque_column=reshape(torques(arclen_indices)*freq_array,[],1); %torque changes within block of freq
energy_column=reshape(energy_per_dist_arclen',[],1); %transposed so ^ holds



energy_interpolator=scatteredInterpolant(freq_column,torque_column,energy_column,'natural','none');
[freq_grid,torque_grid]=meshgrid(freq_array,torque_array);
interpolated_energies=energy_interpolator(freq_grid,torque_grid);


figure
surf(freq_grid,torque_grid,interpolated_energies)
view(2) %top view
xlabel('Speed (Hz)')
ylabel('Torque (pN\cdot nm)')
exportgraphics(gca,[dir,'/energy_per_dist_surface_arclen_only.png'])

figure
contour(freq_grid,torque_grid,interpolated_energies);
xlabel('Speed (Hz)')
ylabel('Torque (pN\cdot nm)')
exportgraphics(gca,[dir,'/energy_per_dist_contour_arclen_only.png'])

%hold on
%scatter3(freq_column,torque_column,energy_column,'filled')
