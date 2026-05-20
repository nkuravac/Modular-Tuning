function [free_F,free_T,wall_F,wall_T,w] = force_and_torque_whole(new_body,flagellum,con_pt,blob_size_on_cell_body,blob_size_on_flag,mu,d,freq)

omega = freq*2*pi;

model = [new_body;flagellum(2:end,:)];

%%%free swimmer
N = size(model,1);
%U = repmat([0 0 1],N,1);
direction =  repmat([0,0,1],N,1);
U = omega*cross(direction,model);

U = reshape(U',N*3,1);
M = MatrixStokeslet_whole(new_body,flagellum(2:end,:),blob_size_on_cell_body,blob_size_on_flag,mu);
free_F = M\U;
free_F = [free_F(1:3:end) free_F(2:3:end) free_F(3:3:end)];
free_T = cross(model-con_pt,free_F);

%%%near wall
lowest_x = min(model(:,1));
avg_x = con_pt(1,1);  %axis of symmetry is on this plane
w = avg_x - d;  %wall position
M = image_MatrixStokeslet_whole(new_body,flagellum(2:end,:),blob_size_on_cell_body,blob_size_on_flag,mu,w);
%save('matrix_wall_whole.mat','M')
wall_F = M\U;
wall_F = [wall_F(1:3:end) wall_F(2:3:end) wall_F(3:3:end)];
wall_T = cross(model-con_pt,wall_F);

end