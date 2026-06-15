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
freq=154;

%%%build a cell body
        r=sqrt(cell_a*pill_height/2); %same surf area as pill
        [r_cb,con_pt] = FibonacciSphereCartesian(n_body,r); 
        new_body = r_cb;


bigdir = ['../../data/ineff_optimize_',char(datetime('now','Format','MM-dd-yyyy_HH-mm'))];
mkdir(bigdir)

%% 

rawdir=[bigdir,'/raw'];
mkdir(rawdir)
summarydir=[bigdir,'/summary'];
mkdir(summarydir)

purcell_ineff=@(x) Purcell_inefficiency(rawdir,new_body,freq,x(1),x(2),r,ds_on_cell_body,opti_blob_size_on_cell_body,x(3),blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);
[x_purcell_min,min_ineff]=fmincon(purcell_ineff,[3.3938,0.485,10],[],[],[],[],[0.3*r,0.1*r,4],[12*r,0.8*r,20]);
opt_wl_purcell=x_purcell_min(1);
opt_hrad_purcell=x_purcell_min(2);
opt_arclen_purcell=x_purcell_min(3);

opt_ineff_struct.Minimum_Purcell_inefficiency=min_ineff;
opt_ineff_struct.r=r;
opt_ineff_struct.Wavelength=opt_wl_purcell;
opt_ineff_struct.Helix_radius=opt_hrad_purcell;
opt_ineff_struct.Arc_length=opt_arclen_purcell;

opt_ineff_json=fopen([summarydir,'/optimal_parameters_purcell_ineff.json'],'w');
fprintf(opt_ineff_json,jsonencode(opt_ineff_struct));
fclose(opt_ineff_json);

opt_ineff=fopen([summarydir,'/optimal_parameters_purcell_ineff.txt'],"w");
fprintf(opt_ineff,'The minimum Purcell inefficiency is %12.9g.\n The optimal wavelength is %12.9g micrometers or %12.9g*r.\n The optimal helical radius is %12.9g micrometers or %12.9g*r.\n The optimal arclength is %12.9g micrometers or %12.9g*r.\n',min_ineff,opt_wl_purcell,opt_wl_purcell/r,opt_hrad_purcell,opt_hrad_purcell/r,opt_arclen_purcell,opt_arclen_purcell/r);
fclose(opt_ineff);

%% 
%{
wave_length_array=linspace(opt_wl_purcell-3*r,opt_wl_purcell+3*r,21);
hrad_array=linspace(opt_hrad_purcell-0.3*r,opt_hrad_purcell+0.3*r,21);
arclen_array=linspace(opt_arclen_purcell-4*r,opt_arclen_purcell+4*r,5);

ineff_array=zeros(length(wave_length_array),length(hrad_array),length(arclen_array));

for i_al=1:length(arclen_array)
    arclen=arclen_array(i_al);
    for i_wl=1:length(wave_length_array)
        for i_hr=1:length(hrad_array)
    
            ineff_array(i_wl,i_hr,i_al)=purcell_ineff([wave_length_array(i_wl),hrad_array(i_hr),arclen_array(i_al)]);

        end
    end
writematrix(ineff_array(:,:,i_al),[outputdir,'/ineff_array_arclen_',num2str(arclen),'.txt'])
%}

%following code is for plotting heat maps of Purcell inefficiency vs.
%wavelength and helical radius at each arc length.
%{
surf(wave_length_array/r,hrad_array/r,ineff_array(:,:,i_al)');
view(0,90);
colorbar;
xlabel("\lambda/r")
ylabel("R/r")
saveas(gcf,[bigdir,'/purcell_ineff_heatmap_freq_',num2str(2*pi*freq),'_arclen_',num2str(arclen_array(i_al)),'.png'])

surf(wave_length_array/r,hrad_array/r,ineff_array(:,:,i_al)');
view(0,90);
colorbar;
clim([0,400]);
xlabel("\lambda/r")
ylabel("R/r")
saveas(gcf,[bigdir,'/modified_purcell_ineff_heatmap_freq_',num2str(2*pi*freq),'_arclen_',num2str(arclen_array(i_al)),'.png'])

ineff_cap=400;
truncated_ineff_array=ineff_array.*double(ineff_array<ineff_cap);
truncated_ineff_array(truncated_ineff_array==0)=1;

surf(wave_length_array/r,hrad_array/r,truncated_ineff_array(:,:,i_al)');
view(0,90);
colorbar;
xlabel("\lambda/r")
ylabel("R/r")
saveas(gcf,[bigdir,'/truncated_',num2str(ineff_cap),'_purcell_ineff_heatmap_freq_',num2str(2*pi*freq),'_arclen_',num2str(arclen_array(i_al)),'.png'])


end
%}