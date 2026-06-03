clearvars
close all

full_sweep_data=table2array(readtable("../../data/full_sweep_table.txt"));
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
        plot(freq_array,load_lines(:,i),'Color','r','LineStyle','-','LineWidth',1)
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


%%
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
%why do I need to transpose energy_per_dist?? That makes the heatmap
%consistent with the range of available load line slopes, but shouldn't it
%be already???
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