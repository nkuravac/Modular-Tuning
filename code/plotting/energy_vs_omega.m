%Construct load lines from the full sweep data, then compute the energy per
%distance at the intersection of each load line with a certain fixed-power
%curve and plot energy per distance vs. angular frequency.

%import data
close all
clearvars
format long

power=1000*200*2*pi; %pN nm /s, in terms of torque ~ 1000 pN nm and frequency ~ 200 Hz
full_sweep_data=table2array(readtable("../../data/full_sweep_table.txt"));
opt_ineff_row=readmatrix("../../data/opt_ineff_row.txt");
full_sweep_data=[full_sweep_data; opt_ineff_row];
full_sweep_data(:,4:end)=full_sweep_data(:,4:end)/(2*pi); %quantity/freq to quantity/angular freq
[full_sweep_data,torque_sort_indices]=sortrows(full_sweep_data,5,"ascend"); %sort by load line (torque vs. ang freq) slope
torque_slopes=full_sweep_data(:,5);
speed_slopes=full_sweep_data(:,6);
opt_ineff_index=find(torque_sort_indices==length(full_sweep_data(:,1)));

%%
%calculations
omega=sqrt(power./torque_slopes); %int between load line (linear with slope torque slopes) and fixed-power curve (hyperbola power/omega)
energy_per_dist=power./(full_sweep_data(:,6).*omega); %E per dist at that int
%energy_per_dist=torque_slopes.*omega./speed_slopes; %equivalent way to calculate E per dist at int


%%
%plot all of the load lines
omega_array=linspace(1,2000,3);
load_lines=repmat(omega_array',1,length(full_sweep_data(:,5)))*diag(full_sweep_data(:,5));

figure
hold on
theme light
for i=1:length(full_sweep_data(:,5))
    if i==1
        plot(omega_array,load_lines(:,i),'Color','b','LineStyle','-')
    elseif i==opt_ineff_index
        plot(omega_array,load_lines(:,i),'Color','r','LineStyle','-','LineWidth',1)
    else
        plot(omega_array,load_lines(:,i),'Color','b','LineStyle','-','HandleVisibility','off')
    end
end
legendlabels={'Load lines','Purcell inefficiency minimizing load line'};
legend(legendlabels,'Location','northwest')
xlabel('\Omega (2\pi/s)')
ylabel('\tau (pN\cdot nm)')
exportgraphics(gca,'../../reports/figures/all_load_lines.png')
hold off

%%
[~,min_diff_index]=min(diff(torque_slopes));

disp('---------------------')
disp(['The two nearest load line slopes are ',num2str(torque_slopes(min_diff_index),10),' and ',num2str(torque_slopes(min_diff_index+1),10)])
disp(['The energy per distance values at the intersections of these two load lines with a curve of constant power ',num2str(power), ' are ', num2str(energy_per_dist(min_diff_index)), ' and ', num2str(energy_per_dist(min_diff_index+1))])

min(torque_slopes)
torque_slopes(opt_ineff_index)
max(torque_slopes)


figure
hold on
for i=1:2
    plot(omega_array,torque_slopes(min_diff_index+i-1)*omega_array,'LineStyle','-','Color',[i-1,0,2-i])
end
hold off


%%
%plot E per dist against omega
figure
theme light
plot(omega, energy_per_dist,'Linestyle','none','Marker','o')
xlabel('\Omega (2\pi/s)')
ylabel('Energy per distance (pN)')
exportgraphics(gca,'../../reports/figures/energy_per_dist_vs_freq_questionable.png')