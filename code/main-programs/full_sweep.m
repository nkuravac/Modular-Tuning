close all
clearvars -except ang_rot
clc
format long

delete *.txt *.fig *.png

%% 

%values here come from the physical system (mu, filament_radius, etc.) or
%are based on previous experiments (blob sizes)
fsize = 20;
ang_rot = 0;
mu = 0.93; 
filament_radius = 0.012;
blob_size_on_flag = 2.139*filament_radius;
n_body=1000;
%arclen=10; %calculated "typical" arclen based on wavelength=1.8, axial length=8.3, R=0.2
num_phase = 16;
rad_st = .44;
cell_a = rad_st;
pill_height=2.2;
ds_on_cell_body = sqrt(4*pi*cell_a^2/n_body);
opti_blob_size_on_cell_body = 0.375*ds_on_cell_body;

%this value doesn't actually matter (the Stokes equations are linear, so all of
%the results are linear in frequency and we just divide this out) but we
%need to run the simulation at some value, and this is a typical one for E.
%coli.
base_freq=154;


%%%build a cell body
r=sqrt(cell_a*pill_height/2); %same surf area as pill from previous paper
[r_cb,con_pt] = FibonacciSphereCartesian(n_body,r); 
new_body = r_cb;

%optimal parameters for Purcell inefficiency from previous simulations
%{
opt_arclen_purcell_previous=9.648*r;
opt_wave_length_purcell_previous=3.687*r;
opt_hrad_purcell_previous=0.5695*r;
%}
%this ^ informs the following search ranges:
arclen_array=linspace(5*r,15*r,32);
wave_length_array=linspace(r,8*r,16);
hrad_array=linspace(0.3*r,0.8*r,16);

opt_arclen_purcell=9.648*r;
opt_wl_purcell=3.687*r;
opt_hrad_purcell=0.5695*r;

arclen_array=sort([arclen_array,opt_arclen_purcell]);
wave_length_array=sort([wave_length_array,opt_wl_purcell]);
hrad_array=sort([hrad_array,opt_hrad_purcell]);


bigdir = ['../../data/full_sweep_',char(datetime('now','Format','MM-dd-yyyy_HH-mm'))];
mkdir(bigdir)

rawdir=[bigdir,'/raw'];
mkdir(rawdir)
summarydir=[bigdir,'/summary'];
mkdir(summarydir)

full_sweep_data=zeros(length(arclen_array)*length(hrad_array)*length(wave_length_array),6);
%columns: arc length, wavelength, helix_radius, force (on body), torque (on body), translational speed (net)


%%

for i_al=1:length(arclen_array)
    for i_wl=1:length(wave_length_array)
        for i_hr=1:length(hrad_array)
            
            [~,~,~,~,avg_sumFb,avg_sumTb,~,~,~,~,~,avg_U_net_calc,~]=simulate_bacterium(rawdir,new_body,base_freq,wave_length_array(i_wl),hrad_array(i_hr),r,ds_on_cell_body,opti_blob_size_on_cell_body,arclen_array(i_al),blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);
            full_sweep_data((i_al-1)*length(wave_length_array)*length(hrad_array)+(i_wl-1)*length(hrad_array)+i_hr,:)=[arclen_array(i_al),wave_length_array(i_wl),hrad_array(i_hr),avg_sumFb(3),avg_sumTb(3),avg_U_net_calc(3)];
            %arranged in blocks for now: R varies within block of constant lambda, which varies within block of constant S. 
            %Purcell ineff needs F, T, and U, while energy per distance
            %just needs T and U. 
            
        end
    end
end
full_sweep_data(:,4:end)=full_sweep_data(:,4:end)/base_freq; %F, T and U are all linear in frequency (and go to zero at zero frequency) so we just want the slopes. 
sweep_column_names={'Arc_length_S','Wavelength_lambda','Helix_radius_R','Force_F_per_freq','Torque_Tau_per_freq','Speed_U_per_freq'};
full_sweep_table=array2table(full_sweep_data,"VariableNames",sweep_column_names);
writetable(full_sweep_table,[summarydir,'/full_sweep_table.txt']);

