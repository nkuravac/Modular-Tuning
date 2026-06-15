close all
clearvars -except ang_rot
clc
format shortG

delete *.txt *.fig *.png

%%%rotate the model using Rodrigues' rotation: left-hand rule, axis of
%%%rotation is the positive y axis, and the angle of rotation is ang_rot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fsize = 20;
ang_rot = 0;
mu = 0.93 %cP
ds_on_cell_body = 2*2*2*0.012;  
opti_blob_size_on_cell_body = ds_on_cell_body/6.4;
filament_radius = 0.012;
blob_size_on_flag = 2.139*filament_radius;

num_phase = 16;

Ls = 8.3; %8.3 +/- 2.0
R = 0.39; %radius of helix  0.39 +/- 0.05
freq = 226;  %154 +/- 30   %Table 3 in Darnton et al.
height_a = 2.5;  %2.5 +/- 0.6

%wave_length = 2.22 - 0.2; %2.22 +/- 0.2
%wave_length_array = [0.2000    0.5000    0.8000    1.1000    1.4000    1.7000    2.0200    2.22  2.3000 2.4200    2.6000 2.9000    3.2000    3.6000    4.0000    5.0000    7.0000    9.0000]; %microns
%wave_length_array = wave_length_array(4:12); %change this 1/27/2024
wave_length_array = 2.22

%d=[0.025 logspace(log10(.05), log10(2),20) logspace(log10(2.06), log10(6), 10) 8 10] + R;
%d=[d(12) d(14) d(17) d(22) d(end)]; %change this 1/27/2024
d = .8;  %wall at x = -d

rad_st = .44;
%cell_a = [rad_st-0.09/2 rad_st rad_st+0.09/2]; %0.88 */- 0.09
cell_a = rad_st;

for i_wl = 1:length(wave_length_array)
    wave_length = wave_length_array(i_wl);
    for i_cb = 1:length(cell_a)
        r = cell_a(i_cb);
        height = height_a;
        dir = ['./wlength_',num2str(wave_length),'_cell_rad_',num2str(r)];
        if(~exist(dir,'dir'))
            mkdir(dir)
        else
            rmdir(dir,'s')
            mkdir(dir)
        end   
        
        diary([dir,'/parameters.txt'])
        r,height,ds_on_cell_body,opti_blob_size_on_cell_body
        Ls,R,freq,wave_length,blob_size_on_flag
        diary off
        
        %%%build a cell body
        [r_cb,con_pt] = create_pill_shape_cylinder(r,height,ds_on_cell_body,0,dir); 
        new_body = r_cb;
        
        fig1=figure;
        axes('fontsize',fsize)
        view(3)
        axis equal
        hold all
        plot3(new_body(:,1),new_body(:,2),new_body(:,3),'ok','MarkerFaceColor','k','markersize',4)
        hold on
        plot3(con_pt(1),con_pt(2),con_pt(3),'og','MarkerFaceColor','g','markersize',6)

        % true_height = max(new_body(:,3)) - min(new_body(:,3))
        % height
        % true_radius = (max(new_body(:,2)) - min(new_body(:,2)))/2.0
        % r

        xlabel('x axis (\mum)');ylabel('y axis (\mum)');zlabel('z axis (\mum)')
        title(['cell height = ',num2str(height),' \mum'])
        saveas(fig1,[dir,'/cell_rad_',num2str(r),'_height_',num2str(height),'_con_pt.png'])
        %saveas(fig1,[dir,'/cell_rad_',num2str(r),'_height_',num2str(height),'_con_pt.fig'],'fig') 
        
        %%%Save for free space (force-free and torque-free)        
        avg_sumF=zeros(1,3);
        avg_sumT=zeros(1,3);
        avg_sumFf=zeros(1,3);
        avg_sumTf=zeros(1,3);
        avg_sumFb=zeros(1,3);
        avg_sumTb=zeros(1,3);
        avg_U_model=zeros(1,3);
        avg_Omega_model=zeros(1,3);
        avg_Omega_cell_body_calc=zeros(1,3);
        avg_Omega_flag_calc=zeros(1,3);
        avg_Omega_net_calc=zeros(1,3);
        avg_U_net_calc=zeros(1,3);

        %%%Save for near wall (force-free and torque-free)
        avg_sumF_wall=zeros(length(d),3);
        avg_sumT_wall=zeros(length(d),3);
        avg_sumFf_wall=zeros(length(d),3);
        avg_sumTf_wall=zeros(length(d),3);
        avg_sumFb_wall=zeros(length(d),3);
        avg_sumTb_wall=zeros(length(d),3);
        avg_U_model_wall=zeros(length(d),3);
        avg_Omega_model_wall=zeros(length(d),3);
        avg_Omega_cell_body_calc_wall=zeros(length(d),3);
        avg_Omega_flag_calc_wall=zeros(length(d),3);
        avg_Omega_net_calc_wall=zeros(length(d),3);
        avg_U_net_calc_wall=zeros(length(d),3);
        
        for phase = 1:num_phase
            clear flagellum
            close all
            disp('*****Phase*****')
            phase
            [flagellum,UIB] = left_handed_helix(phase,freq,blob_size_on_flag,R,wave_length,Ls,dir);
                        
            %%%pull the cell to the helix
            new_bodyt = move_body(new_body,flagellum);
            new_body = new_bodyt;
            con_pt = new_body(end,:); %last point on cell body is the same as the first point on flagellum
            first_point_on_flagellum = flagellum(1,:);
            num_cell = size(new_body,1);
            num_flag = size(flagellum,1);
            tmp = [num_cell;num_flag];

            %%%plot cell body + flagellum
            if phase == 1
                fig1=figure('Visible','off');
                axes('fontsize',fsize)
                view(3)
                axis equal
                hold all
    
                [X,Y,Z] = cylinder(r);
                Z = Z*height+r;
                surf(X,Y,Z,'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none')
    
                plot3(new_body(:,1),new_body(:,2),new_body(:,3),'.k','MarkerFaceColor','k','markersize',2)
                plot3(flagellum(:,1),flagellum(:,2),flagellum(:,3),'-g','linewidth',2)
            end

            %%%rotate cell body
            ind_top = find(new_body(:,3)> height-0.0001 & abs(new_body(:,1)) < 0.0001 & abs(new_body(:,2)) < 0.0001);
            if phase == 1
                plot3(new_body(ind_top,1),new_body(ind_top,2),new_body(ind_top,3),'ok','MarkerFaceColor','k','MarkerSize',2)
            end
            rot_axis = repmat([0;1;0],1,size(new_body,1));
            new_body_temp = (1 - cos(ang_rot))*repmat(dot(new_body',rot_axis),3,1).*rot_axis + cos(ang_rot)*new_body' - sin(ang_rot)*cross(rot_axis,new_body',1);
            new_body = new_body_temp' + new_body(ind_top,:)-new_body_temp(:,ind_top)';
            
            if (phase == 1)
                plot3(new_body(:,1),new_body(:,2),new_body(:,3),'.c','MarkerFaceColor','c','markersize',2)
            end

            %%%rotate flagellum
            rot_axis = repmat([0;1;0],1,size(flagellum,1));
            flag_temp = (1 - cos(ang_rot))*repmat(dot(flagellum',rot_axis),3,1).*rot_axis + cos(ang_rot)*flagellum' - sin(ang_rot)*cross(rot_axis,flagellum',1);
            flagellum = flag_temp' + new_body(end,:) - flagellum(1,:);
            
            if (phase == 1)
                plot3(flagellum(:,1),flagellum(:,2),flagellum(:,3),'.m','MarkerFaceColor','m','markersize',2)
            end

            %%%rotate UIB on flagellum
            rot_axis = repmat([0;1;0],1,size(UIB,1));
            UIB_temp = (1 - cos(ang_rot))*repmat(dot(UIB',rot_axis),3,1).*rot_axis + cos(ang_rot)*UIB' - sin(ang_rot)*cross(rot_axis,UIB',1);
            UIB = UIB_temp';

            if (phase == 1)
                quiver3(flagellum(:,1),flagellum(:,2),flagellum(:,3),UIB(:,1),UIB(:,2),UIB(:,3))
    
                xlabel('x axis (\mum)');ylabel('y axis (\mum)');zlabel('z axis (\mum)')
                title({['cell radius = ',num2str(r)],...
                    ['number of points on cell body = ',num2str(num_cell)],...
                    ['number of points on flagellum = ',num2str(num_flag)]})
                saveas(fig1,[dir,'/cell_rad_',num2str(r),'_phase_',num2str(phase),'.png'])
                saveas(fig1,[dir,'/cell_rad_',num2str(r),'_phase_',num2str(phase),'.fig'],'fig')
            end

            fid_FT_whole = fopen([dir,'/free_FT_whole_phase_',num2str(phase),'.txt'],'w');
            fid_FT_flag = fopen([dir,'/free_FT_flag_phase_',num2str(phase),'.txt'],'w');            
            fid_FT_cell = fopen([dir,'/free_FT_cell_phase_',num2str(phase),'.txt'],'w'); 
            fid_U_whole = fopen([dir,'/free_U_whole_phase_',num2str(phase),'.txt'],'w');
            fid_Omega_cell = fopen([dir,'/free_Omega_whole_phase_',num2str(phase),'.txt'],'w');

            %%%Free space
            flag_wall = 0;
            [sumF,sumT,sumFf,sumTf,sumFb,sumTb,U_model,Omega_model,FCMt,UIB_post] = force_free_and_torque_free_whole(dir,phase,new_body,flagellum,con_pt,opti_blob_size_on_cell_body,blob_size_on_flag,mu,Inf,UIB,flag_wall);
            fprintf(fid_FT_flag,'%12.4g %12.4g %12.4g %12.4g %12.4g %12.4g\n',sumFf,sumTf);
            fprintf(fid_FT_cell,'%12.4g %12.4g %12.4g %12.4g %12.4g %12.4g\n',sumFb,sumTb);
            fprintf(fid_FT_whole,'%12.4g %12.4g %12.4g %12.4g %12.4g %12.4g\n',sumF,sumT);
            fprintf(fid_U_whole,'%12.4g %12.4g %12.4g\n',U_model);
            fprintf(fid_Omega_cell,'%12.4g %12.4g %12.4g\n',Omega_model);
            avg_sumF = avg_sumF + sumF; 
            avg_sumT = avg_sumT + sumT;
            avg_sumFf = avg_sumFf + sumFf; 
            avg_sumTf = avg_sumTf + sumTf;  
            avg_sumFb = avg_sumFb + sumFb;
            avg_sumTb = avg_sumTb + sumTb;
            avg_U_model = avg_U_model + U_model;
            avg_Omega_model = avg_Omega_model + Omega_model;

            fclose(fid_FT_whole);
            fclose(fid_FT_flag); 
            fclose(fid_FT_cell);
            fclose(fid_U_whole);
            fclose(fid_Omega_cell);
            

            %%%calculate velocities at IB points (free space)
            IB_points = [new_body;flagellum(2:end,:)];
            sbody = 1; ebody = length(new_body); sflag = ebody + 1; eflag = length(IB_points);

            blob_size_array = [opti_blob_size_on_cell_body*ones(ebody,1);blob_size_on_flag*ones(eflag - sflag + 1,1)];
            [Omega_cell_body,Omega_flag,Omega_net,U_net,U_free] = velocities_at_IB_points(flag_wall,Inf,IB_points,blob_size_array,mu,sbody,ebody,sflag,eflag,FCMt,dir,phase);
            save([dir,'/free_U_calc_phase_',num2str(phase),'.txt'],'U_net','-ascii')
            save([dir,'/free_U_calc_IB_phase_',num2str(phase),'.txt'],'U_free','-ascii')
            save([dir,'/free_Omega_calc_phase_',num2str(phase),'.txt'],'Omega_cell_body','Omega_flag','Omega_net','-ascii')
            
            avg_Omega_cell_body_calc = avg_Omega_cell_body_calc + Omega_cell_body;
            avg_Omega_flag_calc = avg_Omega_flag_calc + Omega_flag;
            avg_Omega_net_calc = avg_Omega_net_calc + Omega_net;
            avg_U_net_calc = avg_U_net_calc + U_net;

            %with wall but with constraints: sum(F_y) = 0, sum(F_z) = 0, sum(T_x) = 0, sum(T_z) = 0, U_x = 0, Omega_f_x = 2*Omega_c_x = 2*2*pi*0.5
            %(Omega_model is Omega_flagellum)
            fid_FT_whole = fopen([dir,'/FT_whole_phase_',num2str(phase),'.txt'],'w');
            fid_FT_flag = fopen([dir,'/FT_flag_phase_',num2str(phase),'.txt'],'w');            
            fid_FT_cell = fopen([dir,'/FT_cell_phase_',num2str(phase),'.txt'],'w'); 
            fid_U_whole = fopen([dir,'/U_whole_phase_',num2str(phase),'.txt'],'w');
            fid_Omega_cell = fopen([dir,'/Omega_cell_phase_',num2str(phase),'.txt'],'w');
            fid_Omega_flag = fopen([dir,'/Omega_flag_phase_',num2str(phase),'.txt'],'w');
            fid_Omega_calc = fopen([dir,'/Omega_calc_phase_',num2str(phase),'.txt'],'w');
            fid_U_calc = fopen([dir,'/U_calc_phase_',num2str(phase),'.txt'],'w');
            
            %%%Wall at x = dist
            flag_wall = 1;
            for i = 1: length(d)
                dist = -d(i);%distance from the centerline (axis of symmetry) on a model to wall
                
                %%%force free & torque free (Omega_model is Omega_cell_body)
                [sumF,sumT,sumFf,sumTf,sumFb,sumTb,U_model,Omega_model,FCMt,UIB_post] = force_free_and_torque_free_whole(dir,phase,new_body,flagellum,con_pt,opti_blob_size_on_cell_body,blob_size_on_flag,mu,dist,UIB,flag_wall);
                fprintf(fid_FT_whole,'%12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g\n',dist,sumF,sumT);
                fprintf(fid_FT_flag,'%12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g\n',dist,sumFf,sumTf);
                fprintf(fid_FT_cell,'%12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g\n',dist,sumFb,sumTb);
                fprintf(fid_U_whole,'%12.4g %12.4g %12.4g %12.4g\n',dist,U_model);
                fprintf(fid_Omega_cell,'%12.4g %12.4g %12.4g %12.4g\n',dist,Omega_model);
                avg_sumF_wall(i,:) = avg_sumF_wall(i,:) + sumF; 
                avg_sumT_wall(i,:) = avg_sumT_wall(i,:) + sumT;
                avg_sumFf_wall(i,:) = avg_sumFf_wall(i,:) + sumFf; 
                avg_sumTf_wall(i,:) = avg_sumTf_wall(i,:) + sumTf;  
                avg_sumFb_wall(i,:) = avg_sumFb(i,:) + sumFb;
                avg_sumTb_wall(i,:) = avg_sumTb_wall(i,:) + sumTb;
                avg_U_model_wall(i,:) = avg_U_model_wall(i,:) + U_model;
                avg_Omega_model_wall(i,:) = avg_Omega_model_wall(i,:) + Omega_model;

                %%%calculate velocities at IB points (near wall)
                %save([dir,'/IB_points_phase_',num2str(phase),'_flag_wall_',num2str(flag_wall),'.mat'],'IB_points');
                [Omega_cell_body,Omega_flag,Omega_net,U_net_wall,U_wall] = velocities_at_IB_points(flag_wall,dist,IB_points,blob_size_array,mu,sbody,ebody,sflag,eflag,FCMt,dir,phase);
                fprintf(fid_Omega_calc,'%12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g\n',dist,Omega_cell_body,Omega_flag,Omega_net);
                fprintf(fid_U_calc,'%12.4g %12.4g %12.4g %12.4g \n',dist,U_net_wall);
                save([dir,'/U_calc_IB_phase_',num2str(phase),'.txt'],'U_wall','-ascii')%%compare with UIB_post at [dir,'/UIB_post_phase_',num2str(phase),'.txt']

                avg_Omega_cell_body_calc_wall(i,:) = avg_Omega_cell_body_calc_wall(i,:) + Omega_cell_body;
                avg_Omega_flag_calc_wall(i,:) = avg_Omega_flag_calc_wall(i,:) + Omega_flag;
                avg_Omega_net_calc_wall(i,:) = avg_Omega_net_calc_wall(i,:) + Omega_net;
                avg_U_net_calc_wall(i,:) = avg_U_net_calc_wall(i,:) + U_net_wall;
            end
            
            fclose(fid_FT_whole);
            fclose(fid_FT_flag); 
            fclose(fid_FT_cell);
            fclose(fid_U_whole);
            fclose(fid_Omega_cell);
            fclose(fid_Omega_flag);
            fclose(fid_Omega_calc);
            fclose(fid_U_calc);
        end
        avg_sumF = avg_sumF/num_phase; 
        avg_sumT = avg_sumT/num_phase;
        avg_sumFf = avg_sumFf/num_phase; 
        avg_sumTf = avg_sumTf/num_phase;  
        avg_sumFb = avg_sumFb/num_phase;
        avg_sumTb = avg_sumTb/num_phase;
        avg_U_model = avg_U_model/num_phase;
        avg_Omega_model = avg_Omega_model/num_phase;
        avg_Omega_cell_body_calc = avg_Omega_cell_body_calc/num_phase;
        avg_Omega_flag_calc = avg_Omega_flag_calc/num_phase;
        avg_Omega_net_calc = avg_Omega_net_calc/num_phase;
        avg_U_net_calc = avg_U_net_calc/num_phase;

        avg_sumF_wall = avg_sumF_wall/num_phase; 
        avg_sumT_wall = avg_sumT_wall/num_phase;
        avg_sumFf_wall = avg_sumFf_wall/num_phase; 
        avg_sumTf_wall = avg_sumTf_wall/num_phase;  
        avg_sumFb_wall = avg_sumFb_wall/num_phase;
        avg_sumTb_wall = avg_sumTb_wall/num_phase;
        avg_U_model_wall = avg_U_model_wall/num_phase;
        avg_U_net_calc_wall = avg_U_net_calc_wall/num_phase;
        avg_Omega_model_wall = avg_Omega_model_wall/num_phase;
        avg_Omega_cell_body_calc_wall = avg_Omega_cell_body_calc_wall/num_phase;
        avg_Omega_flag_calc_wall = avg_Omega_flag_calc_wall/num_phase;
        avg_Omega_net_calc_wall = avg_Omega_net_calc_wall/num_phase; 

        fid_FT_whole = fopen([dir,'/avg_FT_whole_force_free_torque_free.txt'],'w');
        fid_FT_flag = fopen([dir,'/avg_FT_flag_force_free_torque_free.txt'],'w');            
        fid_FT_cell = fopen([dir,'/avg_FT_cell_force_free_torque_free.txt'],'w'); 
        fid_U_whole = fopen([dir,'/avg_U_force_free_torque_free.txt'],'w');
        fid_Omega_whole = fopen([dir,'/avg_Omega_cell_body_force_free_torque_free.txt'],'w');
        fprintf(fid_FT_whole,'%12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g\n',Inf,avg_sumF,avg_sumT);
        fprintf(fid_FT_flag,'%12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g\n',Inf,avg_sumFf,avg_sumTf);
        fprintf(fid_FT_cell,'%12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g\n',Inf,avg_sumFb,avg_sumTb);
        fprintf(fid_U_whole,'%12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g\n',Inf,avg_U_model,avg_U_net_calc);
        fprintf(fid_Omega_whole,'%12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g\n',Inf,avg_Omega_model,avg_Omega_cell_body_calc,avg_Omega_flag_calc,avg_Omega_net_calc);
        for i = length(d):-1:1
            dist = -d(i);%distance from the centerline (axis of symmetry) on a model to wall

            fprintf(fid_FT_whole,'%12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g\n',dist,avg_sumF_wall(i,:),avg_sumT_wall(i,:));
            fprintf(fid_FT_flag,'%12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g\n',dist,avg_sumFf_wall(i,:),avg_sumTf_wall(i,:));
            fprintf(fid_FT_cell,'%12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g\n',dist,avg_sumFb_wall(i,:),avg_sumTb_wall(i,:));
            fprintf(fid_U_whole,'%12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g\n',dist,avg_U_model_wall(i,:),avg_U_net_calc_wall(i,:));
            fprintf(fid_Omega_whole,'%12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g\n',dist,avg_Omega_model_wall(i,:),avg_Omega_cell_body_calc_wall(i,:),avg_Omega_flag_calc_wall(i,:),avg_Omega_net_calc_wall(i,:));
        end
        fclose(fid_FT_whole);
        fclose(fid_FT_flag); 
        fclose(fid_FT_cell);
        fclose(fid_U_whole);
        fclose(fid_Omega_whole);  
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%
function [Omega_cell_body,Omega_flag,Omega_net,U_net,U] = velocities_at_IB_points(flag_wall,h,IB_points,d,mu,sbody,ebody,sflag,eflag,FCMt,dir,phase)

    %%%calculate velocities at IB points
    if flag_wall == 0
        M2t = MatrixStokeslet_with_blob_size_array(IB_points,IB_points,d,mu);
    else
        M2t = image_MatrixStokeslet_with_blob_size_array(IB_points,IB_points,d,mu,h);
    end
    %save([dir,'/M2t_phase_',num2str(phase),'_flag_wall_',num2str(flag_wall),'.mat'],'M2t')
    
    UCMt = M2t*FCMt;

    U = [UCMt(1:3:end),UCMt(2:3:end),UCMt(3:3:end)];
    %%%avg velocity of the model
    U_net = sum(U)/size(IB_points,1);

    r_vec = IB_points-IB_points(ebody,:); %con_pt is IB_points(ebody,:)
    %IB_points(ebody,:)is a referenced point
    %U = omega_x*cross(direction,IB_points-IB_points(ebody,:));
    direction =  [1,0,0]; %axis of rotation is the x axis
    Omega_vecs_tmp_x = angular_vel_in_direction(U,direction,r_vec,sflag);
    direction =  [0,1,0]; %axis of rotation is the y axis
    Omega_vecs_tmp_y = angular_vel_in_direction(U,direction,r_vec,sflag);
    direction =  [0,0,1]; %axis of rotation is the z axis
    Omega_vecs_tmp_z = angular_vel_in_direction(U,direction,r_vec,sflag);
    Omega_vecs = zeros(size(U));
    Omega_vecs(:,1) = Omega_vecs_tmp_x(:,1);
    Omega_vecs(:,2) = Omega_vecs_tmp_y(:,2);
    Omega_vecs(:,3) = Omega_vecs_tmp_z(:,3);

    Omega_cell_body = sum(Omega_vecs(sbody:ebody,:))/(ebody-sbody+1);
    Omega_flag = sum(Omega_vecs(sflag:eflag,:))/(eflag-sflag+1);
    Omega_net = sum(Omega_vecs)/size(IB_points,1);
end

function Omega_vecs_tmp = angular_vel_in_direction(U,direction,r_vec,sflag)

Omega_vecs_tmp = zeros(size(U));
for iv = 1:size(U,1)
    %iv
    if iv == sflag
        Omega_vecs_tmp(iv,:) = [0 0 0];
    else
        Omega_vec_dir = cross(direction,r_vec(iv,:));
        tmp = norm(Omega_vec_dir);
        if tmp < 1e-10
            Omega_vecs_tmp(iv,:) = [0 0 0];
        else
            Omega_vecs_tmp(iv,:) = dot(U(iv,:),Omega_vec_dir)/(tmp^2)*direction;
        end
    end
end

end
%%%%%%%%%%%%%%%%
function A = MatrixStokeslet_with_blob_size_array(x,y,d,mu)
%
% x = evaluation points
% y = location of the Stokeslets
% d2 = square of blob size
% 
%


nx = length(x(:,1));
ny = length(y(:,1));

A = zeros(3*nx,3*ny);

d2_array = d.^2;
for i = 1:ny
      d2 = d2_array(i);
      columnid = 3*(i-1);
      columnid1 = columnid+1;
      columnid2 = columnid+2;
      columnid3 = columnid+3;

      dx1 = x(:,1)-y(i,1);
      dx2 = x(:,2)-y(i,2);
      dx3 = x(:,3)-y(i,3);
      r2  = dx1.^2 + dx2.^2 + dx3.^2;
      R  = sqrt(r2+d2);
      H1 = 1/2./R + 1/2*d2./R.^3;
      H2 = 1/2./R.^3;

      A(1:3:end,columnid1) = H1 + dx1.*dx1.*H2;
      A(1:3:end,columnid2) = dx1.*dx2.*H2 ;
      A(1:3:end,columnid3) = dx1.*dx3.*H2 ;

      A(2:3:end,columnid1) = dx2.*dx1.*H2 ;
      A(2:3:end,columnid2) = H1 + dx2.*dx2.*H2 ;
      A(2:3:end,columnid3) = dx2.*dx3.*H2 ;

      A(3:3:end,columnid1) = dx3.*dx1.*H2 ;
      A(3:3:end,columnid2) = dx3.*dx2.*H2 ;
      A(3:3:end,columnid3) = H1 + dx3.*dx3.*H2 ;
end

A = A/(4*pi*mu);
end

%%%%%%%%%%%%%%%%%
function A = image_MatrixStokeslet_with_blob_size_array(x,y,d,mu,w)

%x,y,d,mu,w
%pause

%
% x = evaluation points
% y = location of the Stokeslets
% d = delta
% w = distances to the wall
%
h = y(:,1)-w;
%disp('image_MatrixStokeslet_with_blob_size_array')
%h
%pause
d2_array = d.^2;

nx = length(x(:,1));
ny = length(y(:,1));

A = zeros(3*nx,3*ny);

for i = 1 : ny
      d2 = d2_array(i);
      columnid = 3*(i-1);
      columnid1 = columnid+1;
      columnid2 = columnid+2;
      columnid3 = columnid+3;
      
      % Parameters set up
      dx1 = x(:,1)-y(i,1);            dx1_im = x(:,1) -(2*w-y(i,1)); 
      dx2 = x(:,2)-y(i,2);            dx2_im = dx2;
      dx3 = x(:,3)-y(i,3);            dx3_im = dx3;
      r2  = dx1.^2 + dx2.^2 + dx3.^2; r2_im = dx1_im.^2 + dx2_im.^2 + dx3_im.^2;
      R  = sqrt(r2+d2);               R_im = sqrt(r2_im +d2);
      
      H1 = 1/2./R + 1/2*d2./R.^3;     H1p = -1/2*sqrt(r2)./R.^3 - 3/2*d2*sqrt(r2)./R.^5;
      H2 = 1/2./R.^3;                 H2p = - 3/2*sqrt(r2)./R.^5;
      D1 = 1./R.^3-3*d2./R.^5; 
      D2 = -3./R.^5;
      
      H1_im = 1/2./R_im + 1/2*d2./R_im.^3;     H1p_im = -1/2*sqrt(r2_im)./R_im.^3 - 3/2*d2*sqrt(r2_im)./R_im.^5;
      H2_im = 1/2./R_im.^3;                 H2p_im = - 3/2*sqrt(r2_im)./R_im.^5;
      D1_im = 1./R_im.^3-3*d2./R_im.^5; 
      D2_im = -3./R_im.^5;
      

      A(1:3:end,columnid1) = H1 + dx1.*dx1.*H2 -...
                             H1_im - dx1_im.*dx1_im.*H2_im -...
                             h(i)^2*D1_im - h(i)^2*dx1_im.*dx1_im.*D2_im + ...
                             2*h(i)*dx1_im.*H2_im*2 + 2*h(i)*dx1_im.*H1p_im./sqrt(r2_im) + 2*h(i)*dx1_im.*dx1_im.*dx1_im.*H2p_im./sqrt(r2_im);
      A(1:3:end,columnid2) = dx1.*dx2.*H2 - dx2_im.*dx1_im.*H2_im +...
                             h(i)^2.*dx2_im.*dx1_im.*D2_im +... 
                             2*h(i)*dx2_im.*H1p_im./sqrt(r2_im)+2*h(i)*dx2_im.*H2_im - ...
                             2*h(i)*dx2_im.*H1p_im./sqrt(r2_im) - 2*h(i)*dx1_im.*dx1_im.*dx2_im.*H2p_im./sqrt(r2_im);
      A(1:3:end,columnid3) = dx1.*dx3.*H2 - dx1_im.*dx3_im.*H2_im + ...
                             h(i)^2.*dx3_im.*dx1_im.*D2_im +...
                             2*h(i)*dx3_im.*H1p_im./sqrt(r2_im)+2*h(i)*dx3_im.*H2_im - ...
                             2*h(i)*dx3_im.*H1p_im./sqrt(r2_im) - 2*h(i)*dx1_im.*dx1_im.*dx3_im.*H2p_im./sqrt(r2_im);


      A(2:3:end,columnid1) = dx1.*dx2.*H2 - dx2_im.*dx1_im.*H2_im -...
                             h(i)^2.*dx2_im.*dx1_im.*D2_im + ...
                             2*h(i)*dx2_im.*H2_im + 2*h(i)*dx1_im.*dx1_im.*dx2_im.*H2p_im./sqrt(r2_im);
      A(2:3:end,columnid2) = H1 + dx2.*dx2.*H2 -...
                             H1_im - dx2_im.*dx2_im.*H2_im +...
                             h(i)^2.*D1_im + h(i)^2.*dx2_im.*dx2_im.*D2_im - ...
                             2*h(i)*dx1_im.*H1p_im./sqrt(r2_im) - 2*h(i)*dx1_im.*H2_im - ...
                             2*h(i)*dx1_im.*H2_im - 2*h(i)*dx1_im.*dx2_im.*dx2_im.*H2p_im./sqrt(r2_im);
      A(2:3:end,columnid3) = dx2.*dx3.*H2 - dx2_im.*dx3_im.*H2_im + ...
                             h(i)^2.*dx3_im.*dx2_im.*D2_im -...
                             2*h(i)*dx1_im.*dx2_im.*dx3_im.*H2p_im./sqrt(r2_im);


      A(3:3:end,columnid1) = dx1.*dx3.*H2 - dx3_im.*dx1_im.*H2_im -...
                             h(i)^2.*dx3_im.*dx1_im.*D2_im + ...
                             2*h(i)*dx3_im.*H2_im + 2*h(i)*dx1_im.*dx1_im.*dx3_im.*H2p_im./sqrt(r2_im);
      A(3:3:end,columnid2) = dx2.*dx3.*H2 - dx2_im.*dx3_im.*H2_im + ...
                             h(i)^2.*dx3_im.*dx2_im.*D2_im -...
                             2*h(i)*dx1_im.*dx2_im.*dx3_im.*H2p_im./sqrt(r2_im);
      A(3:3:end,columnid3) = H1 + dx3.*dx3.*H2 -...
                             H1_im - dx3_im.*dx3_im.*H2_im +...
                             h(i)^2.*D1_im + h(i)^2.*dx3_im.*dx3_im.*D2_im - ...
                             2*h(i)*dx1_im.*H1p_im./sqrt(r2_im)-2*h(i)*dx1_im.*H2_im - ...
                             2*h(i)*dx1_im.*H2_im - 2*h(i)*dx1_im.*dx3_im.*dx3_im.*H2p_im./sqrt(r2_im);
end
A = A/(4*pi*mu);
end
