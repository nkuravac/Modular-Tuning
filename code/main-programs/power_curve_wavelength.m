close all
clearvars -except ang_rot
clc
format long

delete *.txt *.fig *.png

fsize = 20;
ang_rot = 0;
mu = 0.93; 
filament_radius = 0.012;
blob_size_on_flag = 2.139*filament_radius;
n_body=1000;
%arclen=10; %calculated "typical" arclen based on wavelength=1.8, axial length=8.3, R=0.2
%arclen_array=[20 18 16 14 12 10 8 6 4 2];

num_phase = 16;

rad_st = .44;
cell_a = rad_st;
pill_height=2.2;
ds_on_cell_body = sqrt(4*pi*cell_a^2/n_body);
opti_blob_size_on_cell_body = 0.375*ds_on_cell_body;

%%%build a cell body
        r=sqrt(cell_a*pill_height/2); %same surf area as pill
        [r_cb,con_pt] = FibonacciSphereCartesian(n_body,r); 
        new_body = r_cb;

opt_arclen_purcell_previous=9.648*r;
opt_wave_length_purcell_previous=3.687*r;
opt_hrad_purcell_previous=0.5695*r;

wave_length_array=[linspace(0.5*r,12*r,11) opt_wave_length_purcell_previous];
wave_length_array=sort(wave_length_array,"descend");

hrad_array=repmat(opt_hrad_purcell_previous,length(wave_length_array),1);
arclen_array=repmat(opt_arclen_purcell_previous,length(wave_length_array),1);

base_freq=154;

n_loadlines=length(arclen_array);
n_tslines=length(arclen_array);

energy_per_dist = zeros(length(arclen_array),1);
U_full_bacterium = energy_per_dist;
flag_torque = energy_per_dist;
E_per_d_per_freq=energy_per_dist;
purcell_ineff=energy_per_dist;
power_output=energy_per_dist;
ex_load_line_slopes=energy_per_dist;
freq_array=energy_per_dist;
ts_y_ints=energy_per_dist;
ts_x_ints=energy_per_dist;

bigdir = ['../../data/power_curve_wavelength',char(datetime('now','Format','MM-dd-yyyy_HH-mm'))];
mkdir(bigdir)

[~,~,~,base_torque_vec,~,~,~,~,~,~,~,~,~] = simulate_bacterium(bigdir,new_body,base_freq,opt_wave_length_purcell_previous,opt_hrad_purcell_previous,r,ds_on_cell_body,opti_blob_size_on_cell_body,opt_arclen_purcell_previous,blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);
base_torque=abs(base_torque_vec(3));


for repeat=1:2
for i=1:n_loadlines
for j=1:n_tslines

arclen=arclen_array(i);
wave_length=wave_length_array(i);
R=hrad_array(i);

if repeat==1
    freq=base_freq;
else
    freq=ts_y_ints(j)/(ex_load_line_slopes(i)+ts_y_ints(j)/ts_x_ints(j));
    freq_array(i,j)=freq;
end    
    %wave_length = wave_length_array(i);
    %R = helix_rad_array(i);

%run sim
[avg_sumF,avg_sumT,avg_sumFf,avg_sumTf,avg_sumFb,avg_sumTb,avg_U_model,avg_Omega_model,avg_Omega_cell_body_calc,avg_Omega_flag_calc,avg_Omega_net_calc,avg_U_net_calc,tmp] = simulate_bacterium(bigdir,new_body,freq,wave_length,R,r,ds_on_cell_body,opti_blob_size_on_cell_body,arclen,blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);


if repeat==1
    ex_load_line_slopes(i)=abs(avg_sumTf(3))/base_freq;
    ts_y_ints(i)=base_torque*base_freq/(sqrt(base_torque*base_freq/ex_load_line_slopes(i)))+ex_load_line_slopes(i)*sqrt(base_torque*base_freq/ex_load_line_slopes(i));
    ts_x_ints(i)=ts_y_ints(i)/ex_load_line_slopes(i);
    
    continue
end


%torque-speed relevant outputs
energy_per_dist(i,j)=abs((avg_sumTf(3)*(freq*2*pi))/avg_U_net_calc(3));
flag_torque(i,j)=avg_sumTf(3);
U_full_bacterium(i,j)=avg_U_net_calc(3);
E_per_d_per_freq(i,j)=energy_per_dist(i,j)/freq;
purcell_ineff(i,j)=abs(energy_per_dist(i,j)/avg_sumFf(3));
power_output(i,j)=abs(avg_sumTf(3)*freq*2*pi);

end
end
end

writematrix(energy_per_dist,[bigdir,'/energy_per_dist.txt'])
writematrix(flag_torque,[bigdir,'/flag_torque.txt'])
writematrix(U_full_bacterium,[bigdir,'/U_full_bacterium.txt'])
writematrix(E_per_d_per_freq,[bigdir,'/energy_per_dist_per_freq.txt'])
writematrix(purcell_ineff,[bigdir,'/purcell_ineff.txt'])
writematrix(power_output,[bigdir,'/power_output.txt'])

%% 

%figure('Theme','light')
hold on
brgrad=color_gradient([0 0 1],[1 0 0],n_loadlines);
for i=1:n_loadlines
   plot([0 freq_array(i,i)],[0 ex_load_line_slopes(i)*freq_array(i,i)],"Color",brgrad(i,:),"LineStyle","-")
end
for j=1:n_tslines
    plot([0 ts_x_ints(j)],[ts_y_ints(j) 0],"Color","#33ffaa")
end
fplot(@(x) base_torque*base_freq/x,[30 400],'Color','k')
xlabel("$\Omega_m/(2\pi)$ (Hz)","Interpreter","latex")
ylabel("$\tau$ (pN$\cdot$nm)","Interpreter","latex")
%saveas(gcf,[bigdir,'/ts_and_load_lines'])
%exportgraphics(gcf,[bigdir,'/ts_and_load_lines_wl.png'])
saveas(gcf,[bigdir,'/ts_and_load_lines_wl.png'])
hold off

[min_energy_per_dist_on_power_curve,min_energy_per_dist_on_power_curve_index]=min(diag(energy_per_dist));
e_per_d_wl=fopen([bigdir,'/min_energy_summary.txt'],"w");
fprintf(e_per_d_wl,'The minimum energy cost per distance on the fixed power curve is %12.4g\n. The corresponding load line has wavelength %12.4g micrometers or %12.4g *r.\n r=%12.4g',min_energy_per_dist_on_power_curve,wave_length_array(min_energy_per_dist_on_power_curve_index),wave_length_array(min_energy_per_dist_on_power_curve_index)/r,r);
fclose(e_per_d_wl);

[max_speed_on_power_curve,max_speed_on_power_curve_index]=min(-diag(U_full_bacterium));
max_speed_on_power_curve=abs(max_speed_on_power_curve);
U_max_wl=fopen([bigdir,'/max_speed_summary.txt'],"w");
fprintf(U_max_wl,'The maximum speed on the fixed power curve is %12.4g\n. The corresponding load line has wavelength %12.4g micrometers or %12.4g *r.\n r=%12.4g',max_speed_on_power_curve,wave_length_array(max_speed_on_power_curve_index),wave_length_array(max_speed_on_power_curve_index)/r,r);
fclose(U_max_wl);

[min_ineff_on_power_curve,min_ineff_on_power_curve_index]=min(diag(purcell_ineff));
min_ineff_on_power_curve=abs(min_ineff_on_power_curve);
min_ineff_wl=fopen([bigdir,'/ineff_summary.txt'],"w");
fprintf(min_ineff_wl,'The minimum Purcell inefficiency on the fixed power curve is %12.4g\n. The corresponding load line has wavelength %12.4g micrometers or %12.4g *r.\n The Purcell inefficiency at the intersection of the energy-per-distance-minimizing load line and the fixed power curve is %12.4g\n r=%12.4g',min_ineff_on_power_curve,wave_length_array(min_ineff_on_power_curve_index),wave_length_array(min_ineff_on_power_curve_index)/r,purcell_ineff(min_energy_per_dist_on_power_curve_index,min_energy_per_dist_on_power_curve_index),r);
fclose(min_ineff_wl);

%figure('Theme','light')
[min_energy_per_dist_col,min_energy_per_dist_col_index]=min(energy_per_dist);
plot(wave_length_array(min_energy_per_dist_col_index)/r,min_energy_per_dist_col,"LineStyle","none","Marker",".","Color","k")
xline(wave_length_array(min_ineff_on_power_curve_index)/r,"Color","g","Label","Optimal wavelength","LabelVerticalAlignment","middle")
xlabel("\lambda/r")
ylabel("Energy per distance (10^{-12} J/m)")
%saveas(gcf,[bigdir,'/min_energy_per_dist_vs_wl'])
%exportgraphics(gcf,[bigdir,'/min_energy_per_dist_vs_wl.png'])
saveas(gcf,[bigdir,'/min_energy_per_dist_vs_wl.png'])

plot(wave_length_array/r,diag(energy_per_dist))
xline(wave_length_array(min_ineff_on_power_curve_index)/r,"Color","g","Label","Purcell-inefficiency-minimizing wavelength","LabelVerticalAlignment","middle")
xlabel("\lambda/r")
ylabel("Energy per distance (10^{-12} J/m)")
saveas(gcf,[bigdir,'/energy_per_dist_vs_wl.png'])

%{
figure
hold on
for i=1:n_loadlines
    plot(y_ints./all_x_ints,energy_per_dist(i,:)/1e3)
end
xlabel("Absolute value of torque-speed curve slope (pN$\cdot$nm$\cdot$s)","Interpreter","latex")
ylabel("Energy per distance (10^{-12} J/m)")
legend("Higher torque geometry","Optimal geometry","Lower torque geometry")
saveas(gcf,[bigdir,'/e_per_dist_vs_ts_slope'])
saveas(gcf,[bigdir,'/e_per_dist_vs_ts_slope_fixed_x_int.png'])
%}

%{
figure
hold on
for j=1:n_tslines
    plot(arclen_array,energy_per_dist(:,j)/1e3,"Color",bggrad(j,:))
end
hold off
xlabel("Arclength (nm)")
ylabel("Energy per distance (10^{-12} J/m)")
saveas(gcf,[bigdir,'/e_per_dist_vs_arclen'])
exportgraphics(gcf,[bigdir,'/e_per_dist_vs_arclen.png'])


figure
hold on
for j=1:n_tslines
    plot(arclen_array,purcell_ineff(:,j),"Color",bggrad(j,:))
end
hold off
xlabel("Arclength (nm)")
ylabel("Purcell inefficiency")
saveas(gcf,[bigdir,'/purcell_ineff_vs_arclen'])
exportgraphics(gcf,[bigdir,'/purcell_ineff_vs_arclen.png'])

figure
hold on
for j=1:n_tslines
    plot(arclen_array,power_output(:,j),"Color",bggrad(j,:))

min_energy_per_dist = min(abs(energy_per_dist(:,j)));
min_energy_arclen_index=find(abs(energy_per_dist(:,j))==min_energy_per_dist);
opt_arclen_energy = arclen_array(min_energy_arclen_index)/r;

min_purcell_ineff = min(abs(purcell_ineff(:,j)));
min_ineff_arclen_index=find(abs(purcell_ineff(:,j))==min_purcell_ineff);
opt_arclen_purcell_ineff = arclen_array(min_ineff_arclen_index)/r;

opt_wl=fopen([bigdir,'/optimal_arclength_stall_torque_',num2str(stall_torque_vec(j)),'_max_speed_',num2str(all_x_ints(j)),'.txt'],"w");
fprintf(opt_wl,'The minimum energy cost per distance is %12.4g, and the corresponding optimal arclength is %12.4g micrometers.\n The minimum Purcell inefficiency is %12.4g, and the corresponding optimal arclength is %12.4g micrometers\n',min_energy_per_dist,opt_arclen_energy,min_purcell_ineff,opt_arclen_purcell_ineff);
fclose(opt_wl);
end
hold off
xlabel("Arclength (nm)")
ylabel("$\tau\Omega$ (pN$\cdot$ nm/s)",'Interpreter','latex')
saveas(gcf,[bigdir,'/power_output_vs_arclen'])
exportgraphics(gcf,[bigdir,'/power_output_vs_arclen.png'])
%}
