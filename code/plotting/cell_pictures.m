
bigdir ='../../reports/figures/sim_cell_pictures';
rawdir=[bigdir,'/helix_data'];

if ~exist(bigdir,'dir')
    mkdir(bigdir)
end

if ~exist(rawdir,'dir')
    mkdir(rawdir)
end


%values here come from the physical system (mu, filament_radius, etc.) or
%are based on previous experiments (blob sizes)
fsize = 20;
ang_rot = 0;
mu = 0.93;
filament_radius = 0.005; %was 0.012 from E coli bundle
blob_size_on_flag = 2.139*filament_radius;
n_body=1000;
num_phase = 16;


rad_st = 0.65/2; %0.44
cell_a = rad_st;
pill_height=1.6; %1.6+0.65=2.25, midpoint of 1.5 to 3.0 range. Was 2.2 for E. coli
ds_on_cell_body = sqrt(4*pi*cell_a^2/n_body);
opti_blob_size_on_cell_body = 0.375*ds_on_cell_body;

freq=154;


%%%build a cell body
r=sqrt(cell_a*pill_height/2); %same surf area as pill from previous paper
[r_cb,con_pt] = FibonacciSphereCartesian(n_body,r); 
new_body = r_cb;


parameter_struct=jsondecode(fileread('../../data/ineff_optimize_sphere/summary/config_parameters.json'));
opt_ineff_struct=jsondecode(fileread('../../data/ineff_optimize_sphere/summary/optimal_parameters_purcell_ineff.json'));

opt_arclen_purcell=opt_ineff_struct.arclength;
opt_wl_purcell=opt_ineff_struct.wavelength;
opt_hrad_purcell=opt_ineff_struct.R;

phase=1;

min_al=r;
max_al=20*r;

arclen_array=[min_al,(min_al+opt_arclen_purcell)/2,opt_arclen_purcell,(max_al+opt_arclen_purcell)/2,max_al];

geometric_parameters_struct=struct(...
    'r',r,...
    'Wavelength',opt_wl_purcell,...
    'Helix_radius',opt_hrad_purcell,...
    'Arc_length',arclen_array...
    );

geometric_parameters=fopen([bigdir,'/parameters.json'],'w');
fprintf(geometric_parameters,jsonencode(geometric_parameters_struct));
fclose(geometric_parameters);

%%
ll_slopes=zeros(length(arclen_array),1);

for i=1:length(arclen_array)

    [~,~,~,~,~,avg_sumTb,~,~,~,~,~,~,~]=simulate_bacterium(rawdir,new_body,freq,opt_wl_purcell,opt_hrad_purcell,r,ds_on_cell_body,opti_blob_size_on_cell_body,arclen_array(i),blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase);
    ll_slopes(i)=avg_sumTb(3)/freq;

end

%%
figure
hold on
freq_array=linspace(10,1350,100);
plot(freq_array,max_power/2/pi./freq_array,'Color','k');
for i=1:length(ll_slopes)
    if i==1
        handle_toggle='on';
    else
        handle_toggle='off';
    end
    freq_points=[0,sqrt(max_power/2/pi/ll_slopes(i))];
    plot(freq_points,ll_slopes(i)*freq_points,'HandleVisibility',handle_toggle,'Color','b')
end
hold off
xlabel('Frequency (Hz)')
ylabel('Torque (pN{\cdot}nm)')
legend('Fixed power curve','Load lines')
ylim([0,1400])
exportgraphics(gca,[bigdir,'/example_load_lines_arclen.png'])





%%
for i=1:length(arclen_array)

    [flagellum,UIB] = left_handed_arclength_helix(phase,freq,blob_size_on_flag,opt_hrad_purcell,opt_wl_purcell,arclen_array(i),rawdir);
    
    %%%pull the cell to the helix
    new_bodyt = move_body(new_body,flagellum);
    new_body = new_bodyt;

    %%%plot cell body + flagellum
    
    fig1=figure('Visible','off');
    view(3)
    axis equal
    axis off
    hold on
    
    plot3(new_body(:,1),new_body(:,2),new_body(:,3),'.r','MarkerFaceColor','r','markersize',2)
    plot3(flagellum(:,1),flagellum(:,2),flagellum(:,3),'-g','linewidth',2)
    saveas(fig1,[bigdir,'/full_cell_sphere_arclen_',num2str(arclen_array(i)),'.png'])
    saveas(fig1,[bigdir,'/full_cell_sphere_arclen_',num2str(arclen_array(i)),'.fig'])
end