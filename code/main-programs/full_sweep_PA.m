close all
clearvars -except ang_rot
clc
format long

delete *.txt *.fig *.png

%% 

%%%
%Edits:
%6/15: updated to include optimal parameters for PA in the sweep. 
%Then limited to just running arclength at optimal parameters for time
%purposes, as well as taking out the Wu parameters. Make sure to update
%optimal parameters if changed.
%
%6/12: Updated to simulate Pseudomonas aeruginosa (PA).
%Changed sphere to pill (simulation specs from OS+HN paper, biological from NCBI),
%filament radius to 5 nm (Wu torque-speed paper), and some file naming to be clear about the
%changes
%%%

bigdir ='../../data/arclen_sweep_PA_mid_Healy_body'; %,char(datetime('now','Format','MM-dd-yyyy_HH-mm'))
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
filament_radius = 0.005; %was 0.012 from E coli bundle
blob_size_on_flag = 2.139*filament_radius;
%n_body=1000;
num_phase = 16;


rad_st = 0.71/2; %midpoint of Healy measurement
cell_a = rad_st;
pill_height=1.91-0.71; %midpoint of Healy measurement %previous: 1.6+0.65=2.25, midpoint of 1.5 to 3.0 range. Was 2.2 for E. coli
%ds_on_cell_body = sqrt(4*pi*cell_a^2/n_body);
%opti_blob_size_on_cell_body = 0.375*ds_on_cell_body;

%this value doesn't actually matter (the Stokes equations are linear, so all of
%the results are linear in frequency and we just divide this out) but we
%need to run the simulation at some value, and this is a typical one for E.
%coli.
base_freq=154;

%{
%%%build a cell body
r=sqrt(cell_a*pill_height/2); %same surf area as pill from previous paper
[r_cb,con_pt] = FibonacciSphereCartesian(n_body,r); 
new_body = r_cb;
%}


ds_on_cell_body = 2*2*2*0.012;  
opti_blob_size_on_cell_body = ds_on_cell_body/6.4;
r=cell_a; %typical of PA
[r_cb,con_pt] = create_pill_shape_cylinder(r,pill_height,ds_on_cell_body,0,rawdir); 
new_body = r_cb;
n_body=length(new_body);


%optimal parameters for Purcell inefficiency from previous simulations for
%E coli
%{
opt_arclen_purcell_previous=9.648*r;
opt_wave_length_purcell_previous=3.687*r;
opt_hrad_purcell_previous=0.5695*r;
%}
%this ^ informs the following search ranges:
arclen_array=[linspace(r,30*r,60),linspace(5.4,6.3,20)];
%wave_length_array=[];%linspace(r,8*r,9);
%hrad_array=[];%linspace(0.3*r,1*r,9);

%parameter_struct=jsondecode(fileread('../../data/ineff_optimize_PA/config_parameters.json'));
%opt_ineff_struct=jsondecode(fileread('../../data/ineff_optimize_PA/optimal_parameters_PA_purcell_ineff.json'));
%opt_ineff_struct=jsondecode(fileread('../../data/ineff_optimize_sphere/summary/optimal_parameters_purcell_ineff.json'));

%prev_r=parameter_struct.r;

%opt_arclen_purcell=opt_ineff_struct.Arc_length*r/prev_r;
%opt_wl_purcell=opt_ineff_struct.Wavelength*r/prev_r;
%opt_hrad_purcell=opt_ineff_struct.Helix_radius*r/prev_r;

R_Wu=0.35;
L_Wu=6.281;
p_Wu=1.514;

arclen_array=sort([arclen_array,L_Wu]); %,L_Wu
wave_length_array=p_Wu;
hrad_array=R_Wu;
%wave_length_array=sort([wave_length_array,opt_wl_purcell]); %,p_Wu
%hrad_array=sort([hrad_array,opt_hrad_purcell]); %,R_Wu


%full_sweep_data=zeros(length(arclen_array)*length(hrad_array)*length(wave_length_array),6);
full_sweep_data=zeros(length(arclen_array)*length(wave_length_array),6);
%columns: arc length, wavelength, helix radius, force (on body), torque (on body), translational speed (net)


%

%%

for i_al=1:length(arclen_array)
    for i_wl=1:length(wave_length_array)
        for i_hr=1:length(hrad_array)
            
            [~,~,~,~,avg_sumFb,avg_sumTb,~,~,~,~,~,avg_U_net_calc,~]=simulate_bacterium(rawdir,new_body,base_freq,wave_length_array(i_wl),hrad_array(i_hr),r,ds_on_cell_body,opti_blob_size_on_cell_body,arclen_array(i_al),blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);
            full_sweep_data((i_al-1)*length(wave_length_array)*length(hrad_array)+(i_wl-1)*length(hrad_array)+i_hr,:)=[arclen_array(i_al),wave_length_array(i_wl),hrad_array(i_hr),avg_sumFb(3),avg_sumTb(3),avg_U_net_calc(3)];
            %arranged in blocks for now: R varies within block of constant lambda, which varies within block of constant S. 
            %Purcell ineff needs F, T, and U, while energy per distance
            %just needs T and U. 


            %[~,~,~,~,avg_sumFb,avg_sumTb,~,~,~,~,~,avg_U_net_calc,~]=simulate_bacterium(rawdir,new_body,base_freq,wave_length_array(i_wl),hrad_array(i_wl),r,ds_on_cell_body,opti_blob_size_on_cell_body,arclen_array(i_al),blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);
            %full_sweep_data((i_wl-1)*length(arclen_array)+i_al,:)=[arclen_array(i_al),wave_length_array(i_wl),hrad_array(i_wl),avg_sumFb(3),avg_sumTb(3),avg_U_net_calc(3)];
            %blocks - opt parameters all arc lengths, then Wu parameters all
            %arc lengths
            
        end
    end
end
full_sweep_data(:,4:end)=full_sweep_data(:,4:end)/base_freq; %F, T and U are all linear in frequency for a given geometry (and go to zero at zero frequency) so we just want the slopes. 
sweep_column_names={'Arc_length_S','Wavelength_lambda','Helix_radius_R','Force_F_per_freq','Torque_Tau_per_freq','Speed_U_per_freq'};
full_sweep_table=array2table(full_sweep_data,"VariableNames",sweep_column_names);
writetable(full_sweep_table,[summarydir,'/arclen_sweep_table_PA.txt']);

