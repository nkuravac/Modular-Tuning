function [avg_sumF,avg_sumT,avg_sumFf,avg_sumTf,avg_sumFb,avg_sumTb,avg_U_model,avg_Omega_model,avg_Omega_cell_body_calc,avg_Omega_flag_calc,avg_Omega_net_calc,avg_U_net_calc,tmp] = simulate_bacterium(bigdir,new_body,freq,wave_length,R,r,ds_on_cell_body,opti_blob_size_on_cell_body,arclen,blob_size_on_flag,n_body,filament_radius,mu,ang_rot,fsize,num_phase)
    
    dir = [bigdir,'/freq_',num2str(freq),'_wl_',num2str(wave_length),'_hr_',num2str(R)];
        if(~exist(dir,'dir'))
            mkdir(dir)
        else
            rmdir(dir,'s')
            mkdir(dir)
        end   
        
        diary([dir,'/parameters.txt'])
        r,ds_on_cell_body,opti_blob_size_on_cell_body,
        freq,arclen,blob_size_on_flag,n_body,filament_radius,mu,wave_length,R
        diary off
  

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


        for phase = 1:num_phase
            clear flagellum
            close all
            disp('*****Phase*****')
            phase
            [flagellum,UIB] = left_handed_arclength_helix(phase,freq,blob_size_on_flag,R,wave_length,arclen,dir);
                        
            %%%pull the cell to the helix
            new_bodyt = move_body(new_body,flagellum);
            new_body = new_bodyt;
            con_pt = new_body(end,:); %last point on cell body is the same as the first point on flagellum
            %first_point_on_flagellum = flagellum(1,:);
            num_cell = size(new_body,1);
            num_flag = size(flagellum,1);
            tmp = [num_cell;num_flag];
            

            
            %%%plot cell body + flagellum
            %{
            if phase == 1
                fig1=figure('Visible','off'); %,'Theme','Light'
                axes('fontsize',fsize)
                view(3)
                axis equal
                %axis off
                hold all
    
                plot3(new_body(:,1),new_body(:,2),new_body(:,3),'.r','MarkerFaceColor','r','markersize',2)
                plot3(flagellum(:,1),flagellum(:,2),flagellum(:,3),'-g','linewidth',2)
                saveas(fig1,[dir,'/cell_rad_',num2str(r),'_wavelength_',num2str(wave_length),'_helix_rad_',num2str(R),'_full_cell.png'])
            end
            %}


            %%%rotate cell body
            ind_top = find(new_body(:,3)> 2*r-0.0001 & abs(new_body(:,1)) < 0.0001 & abs(new_body(:,2)) < 0.0001);
            %{
            if phase == 1
                plot3(new_body(ind_top,1),new_body(ind_top,2),new_body(ind_top,3),'ok','MarkerFaceColor','k','MarkerSize',2)
            end
            %}
            rot_axis = repmat([0;1;0],1,size(new_body,1));
            new_body_temp = (1 - cos(ang_rot))*repmat(dot(new_body',rot_axis),3,1).*rot_axis + cos(ang_rot)*new_body' - sin(ang_rot)*cross(rot_axis,new_body',1);
            new_body = new_body_temp' + new_body(ind_top,:)-new_body_temp(:,ind_top)';
            
            %{
            if (phase == 1)
                plot3(new_body(:,1),new_body(:,2),new_body(:,3),'.c','MarkerFaceColor','c','markersize',2)
            end
            %}

            %%%rotate flagellum
            rot_axis = repmat([0;1;0],1,size(flagellum,1));
            flag_temp = (1 - cos(ang_rot))*repmat(dot(flagellum',rot_axis),3,1).*rot_axis + cos(ang_rot)*flagellum' - sin(ang_rot)*cross(rot_axis,flagellum',1);
            flagellum = flag_temp' + new_body(end,:) - flagellum(1,:);
            
            %{
            if (phase == 1)
                plot3(flagellum(:,1),flagellum(:,2),flagellum(:,3),'.m','MarkerFaceColor','m','markersize',2)
            end
            %}

            %%%rotate UIB on flagellum
            rot_axis = repmat([0;1;0],1,size(UIB,1));
            UIB_temp = (1 - cos(ang_rot))*repmat(dot(UIB',rot_axis),3,1).*rot_axis + cos(ang_rot)*UIB' - sin(ang_rot)*cross(rot_axis,UIB',1);
            UIB = UIB_temp';
            
            %%%Free space
            flag_wall = 0;
            [sumF,sumT,sumFf,sumTf,sumFb,sumTb,U_model,Omega_model,FCMt,UIB_post] = force_free_and_torque_free_whole(dir,phase,new_body,flagellum,con_pt,opti_blob_size_on_cell_body,blob_size_on_flag,mu,Inf,UIB,flag_wall);
            
            avg_sumF = avg_sumF + sumF; 
            avg_sumT = avg_sumT + sumT;
            avg_sumFf = avg_sumFf + sumFf; 
            avg_sumTf = avg_sumTf + sumTf;  
            avg_sumFb = avg_sumFb + sumFb;
            avg_sumTb = avg_sumTb + sumTb;
            avg_U_model = avg_U_model + U_model;
            avg_Omega_model = avg_Omega_model + Omega_model;

            %%%calculate velocities at IB points (free space)
            IB_points = [new_body;flagellum(2:end,:)];
            sbody = 1; ebody = length(new_body); sflag = ebody + 1; eflag = length(IB_points);

            blob_size_array = [opti_blob_size_on_cell_body*ones(ebody,1);blob_size_on_flag*ones(eflag - sflag + 1,1)];
            
            [Omega_cell_body,Omega_flag,Omega_net,U_net,U_free] = velocities_at_IB_points(flag_wall,Inf,IB_points,blob_size_array,mu,sbody,ebody,sflag,eflag,FCMt,dir,phase);
            
            avg_Omega_cell_body_calc = avg_Omega_cell_body_calc + Omega_cell_body;
            avg_Omega_flag_calc = avg_Omega_flag_calc + Omega_flag;
            avg_Omega_net_calc = avg_Omega_net_calc + Omega_net;
            avg_U_net_calc = avg_U_net_calc + U_net;

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