close all
clearvars

%%
summarydir='../../data/nonlinearity_check_sweep/summary';
force_matrix=readmatrix([summarydir,'/force_table.txt']);
torque_matrix=readmatrix([summarydir,'/torque_table.txt']);
U_matrix=readmatrix([summarydir,'/U_table.txt']);
freq_array=readmatrix([summarydir,'/freq_array.txt']);
z_shift_array=readmatrix([summarydir,'/z_shift_array.txt']);
parameters=jsondecode(fileread([summarydir,'/geometric_parameters.json']));

freq_axis_label='Frequency (Hz)';
y_axis_labels_freq={'Drag force (nN)','Torque (pN{\cdot}nm)','Speed ({\mu}m/s)'};

plotdir='../../reports/figures/nonlinearity_check_plots';
if ~exist(plotdir,'dir')
    mkdir(plotdir)
end

slopes=zeros(length(z_shift_array),3); %force, torque, speed

%%

for i=1:length(z_shift_array)
    figure
    plot(freq_array,force_matrix(:,3*i))
    xlabel(freq_axis_label)
    ylabel(y_axis_labels_freq{1})
    linfit=polyfit(freq_array,force_matrix(:,3*i),1);
    slopes(i,1)=linfit(1);
    exportgraphics(gca,[plotdir,'/force_vs_freq_z_shift_',num2str(z_shift_array(i)),'.png'])

    figure
    plot(freq_array,torque_matrix(:,3*i))
    xlabel(freq_axis_label)
    ylabel(y_axis_labels_freq{2})
    linfit=polyfit(freq_array,torque_matrix(:,3*i),1);
    slopes(i,2)=linfit(1);
    exportgraphics(gca,[plotdir,'/torque_vs_freq_z_shift_',num2str(z_shift_array(i)),'.png'])

    figure
    plot(freq_array,U_matrix(:,3*i))
    xlabel(freq_axis_label);
    ylabel(y_axis_labels_freq{3})
    linfit=polyfit(freq_array,U_matrix(:,3*i),1);
    slopes(i,3)=linfit(1);
    exportgraphics(gca,[plotdir,'/U_vs_freq_z_shift_',num2str(z_shift_array(i)),'.png'])
end

%%
z_shift_axis_label='Distance to flagellum ({\mu}m)';
y_axis_labels_z={'Drag force per unit frequency (nN{\cdot}s)','Torque per unit frequency (pN{\cdot}nm{\cdot}s)','Speed per unit frequency ({\mu}m)'};
plot_titles={'drag_vs_z_shift','torque_vs_z_shift','U_vs_z_shift'};

for j=1:3
    figure
    plot(z_shift_array,slopes(:,j))
    xscale log
    xlabel(z_shift_axis_label)
    ylabel(y_axis_labels_z{j})
    exportgraphics(gca,[plotdir,'/',plot_titles{j},'.png'])
end

%%
figure
plot(z_shift_array,slopes(:,1)./slopes(:,3))
yline(6*pi*parameters.r*parameters.Viscosity,'LineStyle','--')
xscale log
xlabel(z_shift_axis_label)
legend('Simulation','Stokes drag')
ylabel('Drag coefficient (nN{\cdot}s/{\mu}m)')
exportgraphics(gca,[plotdir,'/drag_coeff_vs_z_shift.png'])