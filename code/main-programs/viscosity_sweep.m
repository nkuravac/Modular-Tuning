close all
clearvars -except ang_rot
clc
format long

delete *.txt *.fig *.png

%% 

%%%


bigdir ='../../data/viscosity_sweep'; %,char(datetime('now','Format','MM-dd-yyyy_HH-mm'))
if ~exist(bigdir,'dir')
    mkdir(bigdir)
end

rawdir=[bigdir,'/raw'];
mkdir(rawdir)
summarydir=[bigdir,'/summary'];
mkdir(summarydir)


%values here come from the physical system (mu, filament_radius, etc.) or
%are based on previous experiments (blob sizes)
fsize = 20;
ang_rot = 0;
mu_array = logspace(-1,1,20);
mu_array=sort([mu_array,0.93]);
filament_radius = 0.005; %was 0.012 from E coli bundle
blob_size_on_flag = 2.139*filament_radius;
n_body=1000;
num_phase = 16;
%rad_st = .65/2;
%cell_a = rad_st;
%pill_height=1.6; %1.6+0.65=2.25, midpoint of 1.5 to 3.0 range. Was 2.2 for E. coli

%this value doesn't actually matter (the Stokes equations are linear, so all of
%the results are linear in frequency and we just divide this out) but we
%need to run the simulation at some value, and this is a typical one for E.
%coli.
base_freq=154;


%%%build a cell body
r=0.5;%0.61;%sqrt(cell_a*pill_height/2); %same surf area as pill from previous paper
[r_cb,con_pt] = FibonacciSphereCartesian(n_body,r); 
new_body = r_cb;

ds_on_cell_body = sqrt(4*pi*r^2/n_body);
opti_blob_size_on_cell_body = 0.382*ds_on_cell_body;

%{
ds_on_cell_body = 2*2*2*0.012;  
opti_blob_size_on_cell_body = ds_on_cell_body/6.4;
r=0.65/2; %typical of PA
[r_cb,con_pt] = create_pill_shape_cylinder(r,pill_height,ds_on_cell_body,0,rawdir); 
new_body = r_cb;
n_body=length(new_body);
%}


opt_ineff_struct=jsondecode(fileread('../../data/ineff_optimize_sphere/summary/optimal_parameters_purcell_ineff.json'));


opt_arclen_purcell=opt_ineff_struct.arclength;
opt_wl_purcell=opt_ineff_struct.wavelength;
opt_hrad_purcell=opt_ineff_struct.R;


full_sweep_data=zeros(length(mu_array),7);
%columns: arc length, wavelength, helix_radius, force (on body), torque (on
%body), translational speed (net), viscosity


%%

for i=1:length(mu_array)
            
            [~,~,~,~,avg_sumFb,avg_sumTb,~,~,~,~,~,avg_U_net_calc,~]=simulate_bacterium(rawdir,new_body,base_freq,opt_wl_purcell,opt_hrad_purcell,r,ds_on_cell_body,opti_blob_size_on_cell_body,opt_arclen_purcell,blob_size_on_flag,n_body,filament_radius,mu_array(i),ang_rot,fsize,num_phase);
            full_sweep_data(i,:)=[opt_arclen_purcell,opt_wl_purcell,opt_hrad_purcell,avg_sumFb(3),avg_sumTb(3),avg_U_net_calc(3),mu_array(i)];
            %repeating the first three parameters just to keep data
            %structure consistent.
            %Purcell ineff needs F, T, and U, while energy per distance
            %just needs T and U. 
            
end
full_sweep_data(:,4:6)=full_sweep_data(:,4:6)/base_freq; %F, T and U are all linear in frequency (and go to zero at zero frequency) so we just want the slopes. 
sweep_column_names={'Arc_length_S','Wavelength_lambda','Helix_radius_R','Force_F_per_freq','Torque_Tau_per_freq','Speed_U_per_freq','Viscosity_mu'};
full_sweep_table=array2table(full_sweep_data,"VariableNames",sweep_column_names);
writetable(full_sweep_table,[summarydir,'/viscosity_sweep_table.txt']);

