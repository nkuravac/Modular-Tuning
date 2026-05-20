function A = image_MatrixStokeslet_whole(new_body,flagellum,blob_size_on_cell_body,blob_size_on_flag,mu,w)

%new_body,flagellum,blob_size_on_cell_body,blob_size_on_flag,mu,w
%pause

model = [new_body;flagellum];
x = model; y = model;
ind_first_pt_flag = size(new_body,1) + 1;

nx = length(model(:,1));
ny = length(model(:,2));

A = zeros(3*nx,3*ny);

h = y(:,1)-w;
%disp('image_MatrixStokeslet_whole')
%h
%pause

%%%cell body
d = blob_size_on_cell_body;
d2 = d^2;


for i = 1 : ind_first_pt_flag-1

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

%%%flagellum
d = blob_size_on_flag;
d2 = d^2;

for i = ind_first_pt_flag : ny

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

% disp('In image_MatrixStokeslet_whole')
% disp('Cell body indices from 1 to ind_first_pt_flag - 1:')
% ind_first_pt_flag-1
% model(ind_first_pt_flag-1,:) %con_pt
% disp('Flagellum indices from ind_first_pt_flag to ')
% ny
% model(ind_first_pt_flag,:)
% pause

A = A/(4*pi*mu);

% A
% save('A.mat','A','model')


%pause
end