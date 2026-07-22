close all
clearvars -except ang_rot
clc
format long

delete *.txt *.fig *.png

%% 
%This script simulates bacteria with a range of flagellum parameters near
%around the global maximum of Purcell efficiency

bigdir ='../../data/sweep_near_global_max_sphere'; %,char(datetime('now','Format','MM-dd-yyyy_HH-mm'))
rawdir=[bigdir,'/raw'];
summarydir=[bigdir,'/summary'];

if ~exist(bigdir,'dir')
    mkdir(bigdir)
end

if ~exist(rawdir,'dir')
    mkdir(rawdir)
end

if ~exist(summarydir,'dir')
    mkdir(summarydir)
end


%values here come from the physical system (mu, filament_radius, etc.) or
%are based on previous experiments (blob sizes)
fsize = 20;
ang_rot = 0;
mu = 0.93;
filament_radius = 0.005;
blob_size_on_flag = 2.139*filament_radius;
n_body=1000;
num_phase = 16;


%this value doesn't actually matter (the Stokes equations are linear, so all of
%the results are linear in frequency and we just divide this out) but we
%need to run the simulation at some value and this is in a typical range
%for bacteria.
base_freq=154;


%%%build a cell body
r=0.5; %nice number with similar fluid drag on body to typical Pseudomonas aeruginosa cell body
[r_cb,con_pt] = FibonacciSphereCartesian(n_body,r); 
new_body = r_cb;

ds_on_cell_body = sqrt(4*pi*r^2/n_body); %approximate distance between points
blob_to_ds_ratio=0.382; %calibrated by comparison to Stokes drag
opti_blob_size_on_cell_body = blob_to_ds_ratio*ds_on_cell_body;


arclen_array=linspace(2*r,30*r,40);
wavelength_array=linspace(r,15*r,30);
hrad_array=linspace(0.25*r,1.5*r,30);

parameter_struct=jsondecode(fileread('../../data/ineff_optimize_sphere/summary/config_parameters.json'));
opt_ineff_struct=jsondecode(fileread('../../data/ineff_optimize_sphere/summary/optimal_parameters_purcell_ineff.json'));

opt_arclen_purcell=opt_ineff_struct.arclength;
opt_wl_purcell=opt_ineff_struct.wavelength;
opt_hrad_purcell=opt_ineff_struct.R;

arclen_array=sort([arclen_array,opt_arclen_purcell]);
wavelength_array=sort([wavelength_array,opt_wl_purcell]);
hrad_array=sort([hrad_array,opt_hrad_purcell]);

arclen_column=[arclen_array';repmat(opt_arclen_purcell,length(wavelength_array)+length(hrad_array),1)];
wl_column=[repmat(opt_wl_purcell,length(arclen_array),1);wavelength_array';repmat(opt_wl_purcell,length(hrad_array),1)];
hrad_column=[repmat(opt_hrad_purcell,length(arclen_array)+length(wavelength_array),1);hrad_array'];

parameter_mat=[arclen_column,wl_column,hrad_column]; %varies one parameter at a time while keeping others fixed at optimal values
num_sims=length(parameter_mat(:,1));

full_sweep_data=zeros(num_sims,6);
%columns: arc length, wavelength, helix radius, force (on body), torque (on body), translational speed (net)
%arranged so that each row makes sense even out of context: the geometric parameters (first 3
%columns) lead to the physical outputs (second 3 columns) 


%%

for i=1:num_sims

        [~,~,~,~,avg_sumFb,avg_sumTb,~,~,~,~,~,avg_U_net_calc,~]=simulate_bacterium(rawdir,new_body,base_freq,wl_column(i),hrad_column(i),r,ds_on_cell_body,opti_blob_size_on_cell_body,arclen_column(i),blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);
        full_sweep_data(i,:)=[parameter_mat(i,:),avg_sumFb(3),avg_sumTb(3),avg_U_net_calc(3)];
   
end

full_sweep_data(:,4:end)=full_sweep_data(:,4:end)/base_freq; %F, T and U are all linear in frequency for a given geometry (and go to zero at zero frequency) so we just want the slopes. 
sweep_column_names={'Arc_length_S','Wavelength_lambda','Helix_radius_R','Force_F_per_freq','Torque_Tau_per_freq','Speed_U_per_freq'};
full_sweep_table=array2table(full_sweep_data,"VariableNames",sweep_column_names);
writetable(full_sweep_table,[summarydir,'/3_parameter_1d_sweeps.txt']);

row_counts=struct(...
    'num_arclens',length(arclen_array),...
    'num_wl',length(wavelength_array),...
    'num_hrad',length(hrad_array)...
    );

row_count_json=fopen([summarydir,'/row_counts.json'],'w');
fprintf(row_count_json,jsonencode(row_counts));
fclose(row_count_json);

config_output=struct(...
    'r',r,...
    'n_body',n_body,...
    'blob_to_ds_ratio',blob_to_ds_ratio,...
    'freq',base_freq,...
    'mu',mu,...
    'filament_radius',filament_radius...
    );

config_output_json=fopen([summarydir,'/config_parameters.json'],'w');
fprintf(config_output_json,jsonencode(config_output));
fclose(config_output_json);