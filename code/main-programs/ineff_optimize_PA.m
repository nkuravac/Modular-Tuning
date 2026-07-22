close all
clearvars -except ang_rot
clc
format long

delete *.txt *.fig *.png

%%%
%Most recent edit 
%7/15: updated radius and pill height to match paper reference
%7/3: updated Purcell inefficiency calculation to reflect redefinition.
%6/15: Updated to simulate Pseudomonas aeruginosa (PA).
%Changed sphere to pill (simulation specs from OS+HN paper, biological from NCBI),
%filament radius to 5 nm (Wu torque-speed paper), and some file naming to be clear about the
%changes
%%%

bigdir = '../../data/ineff_optimize_PA_mid_Healy_body'; %char(datetime('now','Format','MM-dd-yyyy_HH-mm'))
mkdir(bigdir)

rawdir=[bigdir,'/raw'];
mkdir(rawdir)
summarydir=[bigdir,'/summary'];
mkdir(summarydir)

%values here come from the physical system (mu, filament_radius, etc.) or
%are based on previous experiments (blob sizes)
fsize = 20;
ang_rot = 0;
mu = 0.93;
filament_radius = 0.005; %was 0.012 from E coli bundle
blob_size_on_flag = 2.139*filament_radius;
num_phase = 16;
%rad_st = .44;
%cell_a = rad_st;
%pill_height=1.6; %1.6+0.65=2.25, midpoint of 1.5 to 3.0 range. Was 2.2 for E. coli


%this value doesn't actually matter (the Stokes equations are linear, so all of
%the results are linear in frequency and we just divide this out) but we
%need to run the simulation at some value, and this is a typical one for E.
%coli.
base_freq=154;


ds_on_cell_body = 2*2*2*0.012;  
opti_blob_size_on_cell_body = ds_on_cell_body/6.4;
r=0.71/2;
pill_height=1.91-2*r;
[r_cb,con_pt] = create_pill_shape_cylinder(r,pill_height,ds_on_cell_body,0,rawdir); 
new_body = r_cb;
n_body=length(new_body);

%%%find drag coefficient for body alone
u0=1;
U=repmat([0;0;u0],n_body,1);

M=MatrixStokeslets_3D(new_body,new_body,opti_blob_size_on_cell_body,mu);

F = M\U;

%calculate net force
Fnet=zeros(3,1);
for i=1:3
    Fnet(i)=sum(F(i:3:end));
end

drag_coeff=Fnet(3)/u0;

parameters=struct('mu',mu,'r',r,'pill_height',pill_height,'filament_radius',filament_radius,'base_freq',base_freq,'drag_coeff_calc',drag_coeff);
parameters_json=fopen([summarydir,'/config_parameters.json'],'w');
fprintf(parameters_json,jsonencode(parameters));
fclose(parameters_json);


%%

purcell_ineff=@(x) Purcell_inefficiency(drag_coeff,rawdir,new_body,base_freq,x(1),x(2),r,ds_on_cell_body,opti_blob_size_on_cell_body,x(3),blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);
[x_purcell_min,min_ineff]=fmincon(purcell_ineff,[1.5,0.35,6.1],[],[],[],[],[0.3*r,0.1*r,5*r],[12*r,1.2*r,40*r]);
opt_wl_purcell=x_purcell_min(1);
opt_hrad_purcell=x_purcell_min(2);
opt_arclen_purcell=x_purcell_min(3);

opt_ineff_struct.Minimum_Purcell_inefficiency=min_ineff;
opt_ineff_struct.Wavelength=opt_wl_purcell;
opt_ineff_struct.Helix_radius=opt_hrad_purcell;
opt_ineff_struct.Arc_length=opt_arclen_purcell;

opt_ineff_json=fopen([summarydir,'/optimal_parameters_PA_purcell_ineff.json'],'w');
fprintf(opt_ineff_json,jsonencode(opt_ineff_struct));
fclose(opt_ineff_json);