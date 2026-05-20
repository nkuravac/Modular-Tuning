function [P,con_pt] = FibonacciSphereCartesian(n,radius)
%n = # of points
fibphi=(1+sqrt(5))/2;

P=zeros(n,3);
indices=linspace(0,n-1,n);

phi=2*pi*indices/fibphi;
theta=acos(1-2*indices/(n-1));

P(:,1)=cos(phi).*sin(theta);
P(:,2)=sin(phi).*sin(theta);
P(:,3)=cos(theta);
P = radius*P;

con_pt=P(end,:);
end

%[appendix]{"version":"1.0"}
%---
