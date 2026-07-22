close all
clearvars -except ang_rot
clc
format long

delete *.txt *.fig *.png

%% 

bigdir ='../../data/nonlinearity_check_sweep'; %,char(datetime('now','Format','MM-dd-yyyy_HH-mm'))
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
mu = 0.93;
filament_radius = 0.005; %was 0.012 from E coli bundle
blob_size_on_flag = 2.139*filament_radius;
n_body=1000;
%arclen=10; %calculated "typical" arclen based on wavelength=1.8, axial length=8.3, R=0.2
num_phase = 16;
rad_st = .44;
cell_a = rad_st;
pill_height=1.6; %1.6+0.65=2.25, midpoint of 1.5 to 3.0 range. Was 2.2 for E. coli
ds_on_cell_body = sqrt(4*pi*cell_a^2/n_body);
opti_blob_size_on_cell_body = 0.375*ds_on_cell_body;

freq_array=linspace(10,1000,100);
z_shift_array=linspace(0,20,6);
z_shift_array=[z_shift_array,50,1000];

%%%build a cell body
r=sqrt(cell_a*pill_height/2); %same surf area as pill from previous paper
[r_cb,con_pt] = FibonacciSphereCartesian(n_body,r); 
new_body = r_cb;


%{
ds_on_cell_body = 2*2*2*0.012;  
opti_blob_size_on_cell_body = ds_on_cell_body/6.4;
r=0.65/2; %typical of PA
[r_cb,con_pt] = create_pill_shape_cylinder(r,pill_height,ds_on_cell_body,0,rawdir); 
new_body = r_cb;
n_body=length(new_body);
%}


opt_ineff_struct=jsondecode(fileread('../../data/optimal_parameters_purcell_ineff.json'));

prev_r=opt_ineff_struct.r;

opt_arclen_purcell=opt_ineff_struct.Arc_length*r/prev_r;
opt_wl_purcell=opt_ineff_struct.Wavelength*r/prev_r;
opt_hrad_purcell=opt_ineff_struct.Helix_radius*r/prev_r;

parameters=struct('r',r,'Arclength',opt_arclen_purcell,'Wavelength',opt_wl_purcell,'Helix_radius',opt_hrad_purcell,'Sphere',true,'Viscosity',mu);
jsonencode(parameters)

parameters_json=fopen([summarydir,'/geometric_parameters.json'],'w');
fprintf(parameters_json,jsonencode(parameters));
fclose(parameters_json);

forces=zeros(length(freq_array),3*length(z_shift_array));
torques=forces;
translational_speeds=forces;


%%

for i_z=1:length(z_shift_array)
    for i_freq=1:length(freq_array)
            
            [~,~,~,~,avg_sumFb,avg_sumTb,~,~,~,~,~,avg_U_net_calc,~]=simulate_bacterium_z_shift(z_shift_array(i_z),rawdir,new_body,freq_array(i_freq),opt_wl_purcell,opt_hrad_purcell,r,ds_on_cell_body,opti_blob_size_on_cell_body,opt_arclen_purcell,blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);
            column_slots=3*(i_z-1)+1:3*(i_z-1)+3;
            forces(i_freq,column_slots)=avg_sumFb;
            torques(i_freq,column_slots)=avg_sumTb;
            translational_speeds(i_freq,column_slots)=avg_U_net_calc;
            %frequency increases down each column. z-shift increases across
            %block of three columns - x,y,z components of each output. x,y
            %should be negligible, but worth checking.

            
    end
end

%%
writematrix(forces,[summarydir,'/force_table.txt'])
writematrix(torques,[summarydir,'/torque_table.txt'])
writematrix(translational_speeds,[summarydir,'/U_table.txt'])
writematrix(freq_array,[summarydir,'/freq_array.txt'])
writematrix(z_shift_array,[summarydir,'/z_shift_array.txt'])

data_explainer=fopen([summarydir,'/data_explainer.txt'],'w');
fprintf(data_explainer,'The force, torque, and U matrices are formatted such that frequency increases down each column while z-shift (distance from cell body end to flagellum start) increases across blocks of three columns.\n The blocks contain the x, y, and z components of the output. The x and y components should be negligible, but worth checking.');
fclose(data_explainer);