close all
clearvars -except ang_rot
clc
format long

delete *.txt *.fig *.png

%MAKE SURE TO SAVE MATRIX OF DATA!! NOT JUST PLOT!!

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
base_freq=154;

%%%build a cell body
        r=sqrt(cell_a*pill_height/2); %same surf area as pill
        [r_cb,con_pt] = FibonacciSphereCartesian(n_body,r); 
        new_body = r_cb;


opt_arclen_purcell_previous=9.648*r;
opt_wave_length_purcell_previous=3.687*r;
opt_hrad_purcell_previous=0.5695*r;

bigdir = ['./energy_optimize_',char(datetime('now','Format','MM-dd-yyyy_HH-mm'))];
mkdir(bigdir)

[~,~,~,base_torque_vec,~,~,~,~,~,~,~,~,~] = simulate_bacterium(bigdir,new_body,base_freq,opt_wave_length_purcell_previous,opt_hrad_purcell_previous,r,ds_on_cell_body,opti_blob_size_on_cell_body,opt_arclen_purcell_previous,blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);
base_torque=abs(base_torque_vec(3));     

fixed_power=base_torque*base_freq;
%% 

energy_per_dist=@(x) Energy_per_distance_fixed_power(base_freq,fixed_power,bigdir,new_body,x(1),x(2),r,ds_on_cell_body,opti_blob_size_on_cell_body,x(3),blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);
[x_energy_min,min_energy]=fmincon(energy_per_dist,[3.3938,0.485,10],[],[],[],[],[0.3*r,0.1*r,4],[12*r,0.8*r,20]);
opt_wl_energy=x_energy_min(1);
opt_hrad_energy=x_energy_min(2);
opt_arclen_energy=x_energy_min(3);

ineff_at_opt_energy=Purcell_inefficiency(bigdir,new_body,base_freq,opt_wl_energy,opt_hrad_energy,r,ds_on_cell_body,opti_blob_size_on_cell_body,opt_arclen_energy,blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);

opt_ineff=fopen([bigdir,'/optimal_parameters_energy.txt'],"w");
fprintf(opt_ineff,'The minimum energy per distance is %12.4g.\n The optimal wavelength is %12.4g micrometers or %12.4g*r.\n The optimal helical radius is %12.4g micrometers or %12.4g*r.\n The optimal arclength is %12.4g micrometers or %12.4g*r.\n The Purcell inefficiency at these parameters is %12.4g\n',min_energy,opt_wl_energy,opt_wl_energy/r,opt_hrad_energy,opt_hrad_energy/r,opt_arclen_energy,opt_arclen_energy/r,ineff_at_opt_energy);
fclose(opt_ineff);


%% 
function energy_per_dist = Energy_per_distance_fixed_power(base_freq,power,bigdir,new_body,wave_length,R,r,ds_on_cell_body,opti_blob_size_on_cell_body,arclen,blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase)
    [~,~,~,avg_sumTf,~,~,~,~,~,~,~,~,~]=simulate_bacterium(bigdir,new_body,base_freq,wave_length,R,r,ds_on_cell_body,opti_blob_size_on_cell_body,arclen,blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);
    freq=sqrt(abs(power*base_freq/avg_sumTf(3)));
    [~,~,~,avg_sumTf,~,~,~,~,~,~,~,avg_U_net_calc,~]=simulate_bacterium(bigdir,new_body,freq,wave_length,R,r,ds_on_cell_body,opti_blob_size_on_cell_body,arclen,blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);
    energy_per_dist=abs(avg_sumTf(3)*(freq*2*pi)/avg_U_net_calc(3));
end

function ineff = Purcell_inefficiency(bigdir,new_body,freq,wave_length,R,r,ds_on_cell_body,opti_blob_size_on_cell_body,arclen,blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase)
    [~,~,avg_sumFf,avg_sumTf,~,~,~,~,~,~,~,avg_U_net_calc,~]=simulate_bacterium(bigdir,new_body,freq,wave_length,R,r,ds_on_cell_body,opti_blob_size_on_cell_body,arclen,blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);
    ineff=abs(avg_sumTf(3)*(freq*2*pi)/(avg_sumFf(3)*avg_U_net_calc(3)));
end

