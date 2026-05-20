function [sumF,sumT,sumFf,sumTf,sumFb,sumTb,U_model,Omega_model,FCMt,UIB_post] = force_free_and_torque_free_whole(dir,phase,new_body,flagellum,con_pt,blob_size_on_cell_body,blob_size_on_flag,mu,d,UIB,flag_wall)

% disp('In force_free_and_torque_free_whole')
% con_pt
% pause

model = [new_body;flagellum(2:end,:)]; %new_body(end,:) = flagellum(1,:)
N = size(model,1);

%%%%%%%%%%%%%%%%%%%
%%%%%Given instantaneous velocity at IB Points, 
%%%%%impose extra conditions: sum of forces = 0, sum of torques = 0,
%%%%%then compute the resulting forces at the IB points
%%%%%%%%%%%%%%%%%%%
newU = zeros(N*3+6,1);
N1 = size(new_body,1);
Ulong = [zeros(N1,3);UIB(2:end,:)];
newU(1:N*3,1)=reshape(Ulong',N*3,1); 

%*******finding IB forces satisfying conservations of 
%*******linear/angular momenta
%*******i.e., set up the matrix system of regularized Stokeslets and two
%*******extra momentum conditions on forces and torques    

%%%free swimmer
%IB_points = [new_body;flagellum(2:end,:)];
%save([dir,'/IB_points_model_phase_',num2str(phase),'_flag_wall_',num2str(flag_wall),'.mat'],'IB_points');
M=NewMatrixStokeslet_new(new_body,flagellum(2:end,:),blob_size_on_cell_body,blob_size_on_flag,mu,con_pt,d,phase,flag_wall,dir);
% cond(M)
% pause


%*******given the instantaneous velocity  at time t, find the forces
%*******at the IB points, total translational velocity and total rotational velocity 
FIB = M\newU;
FCMt = FIB(1:end-6);

F1=FIB(1:3:end-6);%force in x
F2=FIB(2:3:end-6);%force in y
F3=FIB(3:3:end-6);%force in z

U_model = FIB(end-5: end-3)';
Omega_model = FIB(end-2: end)';

sumF = sum([F1 F2 F3]);
sumT = sum(cross(model-con_pt,[F1 F2 F3]));   
sumFf = sum([F1(N1+1:end) F2(N1+1:end) F3(N1+1:end)]);
sumTf = sum(cross(model(N1+1:end,:)-con_pt,[F1(N1+1:end) F2(N1+1:end) F3(N1+1:end)]));   
sumFb = sum([F1(1:N1) F2(1:N1) F3(1:N1)]);
sumTb = sum(cross(model(1:N1,:)-con_pt,[F1(1:N1) F2(1:N1) F3(1:N1)]));   

% sumT
% sumTf
% sumTb
% pause

%%%Check if Ulong + U_model + Omega_model x (xj - xr) = M*F
Omega_model_rep = repmat(Omega_model,N,1);
sum_Omega = cross(Omega_model_rep,model-con_pt); 

UIB_post = Ulong + repmat(U_model,N,1) + sum_Omega;

if (flag_wall == 0)
    disp('--- Free space (force free & torque free)')
    save([dir,'/free_UIB_post_phase_',num2str(phase),'_flag_wall_',num2str(flag_wall),'.txt'],'UIB_post','-ascii')
else
    disp('--- Near wall (force free & torque free)')
    save([dir,'/UIB_post_phase_',num2str(phase),'_flag_wall_',num2str(flag_wall),'.txt'],'UIB_post','-ascii')
end
fprintf('Net force: [%2.1g %2.1g %2.1g] \n', sumF)
fprintf('Net torque: [%2.1g %2.1g %2.1g] \n', sumT)

end
