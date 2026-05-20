function A = MatrixStokeslet_whole(new_body,flagellum,blob_size_on_cell_body,blob_size_on_flag,mu)


model = [new_body;flagellum];
x = model; y = model;
ind_first_pt_flag = size(new_body,1) + 1;

nx = length(model(:,1));
ny = length(model(:,1));

A = zeros(3*nx,3*ny);

%%%cell body
d = blob_size_on_cell_body;
d2 = d^2;

for i = 1 : ind_first_pt_flag-1

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

%%%flagellum
d = blob_size_on_flag;
d2 = d^2;

for i = ind_first_pt_flag : ny

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

% disp('In MatrixStokeslet_whole')
% disp('Cell body indices from 1 to ind_first_pt_flag - 1:')
% ind_first_pt_flag-1
% model(ind_first_pt_flag-1,:)
% disp('Flagellum indices from ind_first_pt_flag to ')
% ny
% model(ind_first_pt_flag,:)
% pause

A = A/(4*pi*mu);
end