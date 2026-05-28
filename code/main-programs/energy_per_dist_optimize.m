close all
clearvars -except ang_rot
clc
format long

delete *.txt *.fig *.png

%% 


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

%optimal parameters for Purcell inefficiency from previous simulations
opt_arclen_purcell_previous=9.648*r;
opt_wl_purcell_previous=3.687*r;
opt_hrad_purcell_previous=0.5695*r;

bigdir = ['../../data/energy_optimize_',char(datetime('now','Format','MM-dd-yyyy_HH-mm'))];
mkdir(bigdir)

rawdir=[bigdir,'/raw'];
mkdir(rawdir)
summarydir=[bigdir,'/summary'];
mkdir(summarydir)

[~,~,~,base_torque_vec,~,~,~,~,~,~,~,~,~] = simulate_bacterium(rawdir,new_body,base_freq,opt_wl_purcell_previous,opt_hrad_purcell_previous,r,ds_on_cell_body,opti_blob_size_on_cell_body,opt_arclen_purcell_previous,blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);
base_torque=abs(base_torque_vec(3));     

fixed_power=base_torque*base_freq; %since I'm hard-coding base_freq anyway, this could probably just be hard-coded to be some nice number that's physically typical
%% 

energy_per_dist=@(x) Energy_per_distance_fixed_power(base_freq,fixed_power,rawdir,new_body,x(1),x(2),r,ds_on_cell_body,opti_blob_size_on_cell_body,x(3),blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);
[x_energy_min,min_energy]=fmincon(energy_per_dist,[3.3938,0.485,10],[],[],[],[],[0.3*r,0.1*r,4],[12*r,0.8*r,20]);
opt_wl_energy=x_energy_min(1);
opt_hrad_energy=x_energy_min(2);
opt_arclen_energy=x_energy_min(3);

ineff_at_opt_energy=Purcell_inefficiency(rawdir,new_body,base_freq,opt_wl_energy,opt_hrad_energy,r,ds_on_cell_body,opti_blob_size_on_cell_body,opt_arclen_energy,blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);
energy_per_dist_at_opt_Purcell=Energy_per_distance_fixed_power(base_freq,fixed_power,rawdir,new_body,opt_wl_purcell_previous,opt_hrad_purcell_previous,r,ds_on_cell_body,opti_blob_size_on_cell_body,opt_arclen_purcell_previous,blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);


opt_ineff=fopen([summarydir,'/optimal_parameters_energy.txt'],"w");
fprintf(opt_ineff,'The minimum energy per distance is %12.4g.\n The optimal wavelength is %12.4g micrometers or %12.4g*r.\n The optimal helical radius is %12.4g micrometers or %12.4g*r.\n The optimal arclength is %12.4g micrometers or %12.4g*r.\n The Purcell inefficiency at these parameters is %12.4g.\n The energy per distance at the Purcell-optimal parameters is %12.4g',min_energy,opt_wl_energy,opt_wl_energy/r,opt_hrad_energy,opt_hrad_energy/r,opt_arclen_energy,opt_arclen_energy/r,ineff_at_opt_energy,energy_per_dist_at_opt_Purcell);
fclose(opt_ineff);

%%
%{
wave_length_array=linspace(opt_wl_purcell_previous-3*r,opt_wl_purcell+3*r,21);
hrad_array=linspace(opt_hrad_purcell_previous-0.3*r,opt_hrad_purcell+0.3*r,21);
arclen_array=linspace(opt_arclen_purcell_previous-4*r,opt_arclen_purcell+4*r,5);

energy_array=zeros(length(wave_length_array),length(hrad_array),length(arclen_array));

for i_al=1:length(arclen_array)
    arclen=arclen_array(i_al);
    for i_wl=1:length(wave_length_array)
        for i_hr=1:length(hrad_array)

            energy_array(i_wl,i_hr,i_al)=Energy_per_distance_fixed_power(base_freq,power,rawdir,new_body,wave_length_array(i_wl),hrad_array(i_hr),r,ds_on_cell_body,opti_blob_size_on_cell_body,arclen_array(i_al),blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);

        end
    end
    writematrix(energy_array(:,:,i_al),[outputdir,'/energy_array_arclen_',num2str(arclen),'.txt'])
end
%}
%% 
function energy_per_dist = Energy_per_distance_fixed_power(base_freq,power,dir,new_body,wave_length,R,r,ds_on_cell_body,opti_blob_size_on_cell_body,arclen,blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase)
    [~,~,~,avg_sumTf,~,~,~,~,~,~,~,~,~]=simulate_bacterium(dir,new_body,base_freq,wave_length,R,r,ds_on_cell_body,opti_blob_size_on_cell_body,arclen,blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);
    freq=sqrt(abs(power*base_freq/avg_sumTf(3)));
    [~,~,~,avg_sumTf,~,~,~,~,~,~,~,avg_U_net_calc,~]=simulate_bacterium(dir,new_body,freq,wave_length,R,r,ds_on_cell_body,opti_blob_size_on_cell_body,arclen,blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);
    energy_per_dist=abs(avg_sumTf(3)*(freq*2*pi)/avg_U_net_calc(3));
end


