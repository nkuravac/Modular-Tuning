function A=NewMatrixStokeslet_new(new_body,flagellum,blob_size_on_cell_body,blob_size_on_flag,mu,con_pt,d,phase,flag_wall,dir)

    N=size(new_body,1)+size(flagellum,1); %already removed the first point on flagellum to avoid overlapped with the last point of the cell body
    I=eye(3);
    Z=zeros(3,3);
    if flag_wall == 0
        A = MatrixStokeslet_whole(new_body,flagellum,blob_size_on_cell_body,blob_size_on_flag,mu);
        %save([dir,'/M_A_phase_',num2str(phase),'_flag_wall_',num2str(flag_wall),'.mat'],'A');
        %%%%Without wall
    else
        %avg_x = con_pt(1,1);  %axis of symmetry is on this plane
        %w = avg_x - d;  %wall position
        A = image_MatrixStokeslet_whole(new_body,flagellum,blob_size_on_cell_body,blob_size_on_flag,mu,d);
        %save([dir,'/M_A_phase_',num2str(phase),'_flag_wall_',num2str(flag_wall),'.mat'],'A');
        %%%%With wall
    end
    
    
    A=[A,repmat(-I,N,1)];  %-I represents the minus sign in the equations containing translational velocites
    A=[A;repmat(-I,1,N),Z];
    temp1=zeros(3*N+3,3);


    model = [new_body;flagellum];
    for i=1:N
        Q=Rotation(model(i,:),con_pt);%- Omega x (x - x_c) = (x - x_c) x Omega
        temp1(3*(i-1)+1:3*i,1:3)=Q;
    end

    temp2=-temp1';%(x - x_c) x f
    A=[A,temp1;temp2,Z];

    % disp('In NewMatrixStokeslet_new')
    % con_pt
    % model(length(new_body),:) %con_pt
    % model(length(new_body)+1,:)
    % pause

end

function [Q]=Rotation(x,c) %
    Q=zeros(3,3); %minus signs in Q represents the minus sign in the equations containing rotational velocites
    dx = x(1) - c(1);
    dy = x(2) - c(2);
    dz = x(3) - c(3);
    Q(1,2)=-dz;
    Q(1,3)=dy;
    Q(2,1)=dz;
    Q(2,3)=-dx;
    Q(3,1)=-dy;
    Q(3,2)=dx;
end
