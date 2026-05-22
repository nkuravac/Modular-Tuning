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