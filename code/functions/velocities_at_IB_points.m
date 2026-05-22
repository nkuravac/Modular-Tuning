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