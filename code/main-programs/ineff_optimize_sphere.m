close all
clearvars -except ang_rot
clc
format long

delete *.txt *.fig *.png

%MAKE SURE TO SAVE MATRIX OF DATA!! NOT JUST PLOT!!

fsize = 20;
ang_rot = 0;
mu = 0.93;
filament_radius = 0.005; %0.012
blob_size_on_flag = 2.139*filament_radius;
n_body=1000;

num_phase = 16;

r=0.5; %0.61; %matches drag on big pill
ds_on_cell_body = sqrt(4*pi*r^2/n_body);
opti_blob_size_on_cell_body = 0.382*ds_on_cell_body;
freq=154;

%%%build a cell body
[r_cb,con_pt] = FibonacciSphereCartesian(n_body,r); 
new_body = r_cb;

bigdir = '../../data/ineff_optimize_sphere';
if ~exist(bigdir,'dir')
    mkdir(bigdir)
end

rawdir=[bigdir,'/raw'];
if ~exist(rawdir,'dir')
    mkdir(rawdir)
end
summarydir=[bigdir,'/summary'];
if ~exist(summarydir,'dir')
    mkdir(summarydir)
end

parameters=struct('r',r,'Sphere',true,'n_body',n_body,'Viscosity',mu,'Filament_radius',filament_radius,'Frequency',freq);
jsonencode(parameters)

parameters_json=fopen([summarydir,'/config_parameters.json'],'w');
fprintf(parameters_json,jsonencode(parameters));
fclose(parameters_json);

%%

purcell_ineff=@(x) Purcell_inefficiency(6*pi*mu*r,rawdir,new_body,freq,x(1),x(2),r,ds_on_cell_body,opti_blob_size_on_cell_body,x(3),blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);
[x_purcell_min,min_ineff]=fmincon(purcell_ineff,[3.3938,0.485,10],[],[],[],[],[0.3*r,0.1*r,4],[12*r,0.8*r,20]);
opt_wl_purcell=x_purcell_min(1);
opt_hrad_purcell=x_purcell_min(2);
opt_arclen_purcell=x_purcell_min(3);

opt_ineff_struct=struct(...
    'min_ineff',min_ineff,...
    'r',r,...
    'wavelength',x_purcell_min(1),...
    'R',x_purcell_min(2),...
    'arclength',x_purcell_min(3)...
    );

opt_ineff_json=fopen([summarydir,'/optimal_parameters_purcell_ineff.json'],'w');
fprintf(opt_ineff_json,jsonencode(opt_ineff_struct));
fclose(opt_ineff_json);