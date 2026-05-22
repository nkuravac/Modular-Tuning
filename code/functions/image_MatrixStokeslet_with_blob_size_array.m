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