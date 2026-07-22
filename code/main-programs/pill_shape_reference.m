rawdir='../../data/pill_raw_outputs';
if ~exist(rawdir,'dir')
    mkdir(rawdir)
end

%cell body size from tian wu zhang yuan, ds and blob from nguyen
r = 0.8/2;
pill_height= 2.53-2*r;
ds_on_cell_body = 2*2*2*0.012;  
opti_blob_size_on_cell_body = ds_on_cell_body/6.4;
[r_cb,con_pt] = create_pill_shape_cylinder(r,pill_height,ds_on_cell_body,0,rawdir); 
new_body = r_cb;
n_body=length(new_body);

freq=154;
mu = 0.93;
u0 = 1;
u=[0; 0; u0];
U=repmat(u,n_body,1);

M = MatrixStokeslets_3D(new_body,new_body,opti_blob_size_on_cell_body,mu);

F = M\U;

%calculate net force
Fnet=zeros(3,1);
for i=1:3
    Fnet(i)=sum(F(i:3:end));
end

Fnet

r_eff=Fnet(3)/(6*pi*mu*u0) %effective radius
body_drag_coeff=Fnet(3)/u0

%%

fsize = 20;
ang_rot = 0;
mu = 0.93; 
filament_radius = 0.005;
blob_size_on_flag = 2.139*filament_radius;

num_phase = 16;

%from same paper as cell body measurements:
R=0.38;
wavelength=1.48;
arclen=6;

%from paper with torque-speed curve:
% R=0.350;
% wavelength=1.514;
% arclen=6.281;

[avg_sumF,avg_sumT,avg_sumFf,avg_sumTf,avg_sumFb,avg_sumTb,avg_U_model,avg_Omega_model,avg_Omega_cell_body_calc,avg_Omega_flag_calc,avg_Omega_net_calc,avg_U_net_calc,tmp]=simulate_bacterium(rawdir,new_body,freq,wavelength,R,r,ds_on_cell_body,opti_blob_size_on_cell_body,arclen,blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase)

purcell_efficiency=body_drag_coeff*avg_U_model(3)^2/(avg_sumTb(3)*(freq*2*pi))