function [flagellum,UIB] = left_handed_helix(phase,freq,blob_size_on_flag,R,wave_length,Ls,dir)
    %R = 0.2;
    %wave_length = 10.9955*R;
    k = 2*pi/wave_length;
    num_phase = 16;
    ang = 2*pi/num_phase;
    %freq = 154;
    omega = -freq*2*pi;%%%ADD "-" to change the rotation direction of flagellum HN:2/28/20225
    
    %Ls = 8.3;
    
    %%%%%compute arclength
    f = @(t) sqrt(((1 - exp(-k^2*t.^2)).^2)*R^2*k^2.*cos(k*t+(phase-1)*ang).^2 +...
                    ((- exp(-k^2*t.^2)).^2)*R^2.*((-2*k^2*t).^2).*sin(k*t+(phase-1)*ang).^2 + ...
                ((1 - exp(-k^2*t.^2)).^2)*R^2*k^2.*sin(k*t+(phase-1)*ang).^2 +...
                    ((- exp(-k^2*t.^2)).^2)*R^2.*((-2*k^2*t).^2).*cos(k*t+(phase-1)*ang).^2 + ...
                    1);
    arclen = integral(f,0,Ls);
    
    numf = round(arclen/blob_size_on_flag)+1;
    
    s = linspace(0,Ls,numf)';
    pz = s;
    px = (1 - exp(-k^2*s.^2))*R.*sin(k*s+(phase-1)*ang); 
    py = (1 - exp(-k^2*s.^2))*R.*cos(k*s+(phase-1)*ang);
        
    pt = interparc(numf,px,py,pz); %Redistribute discrete points along flagellum so that they are equi-distant from each other   
%     POINT
%     pt(1,:)
    
    pt(1,:) = [0 0 0];
    
%     d1 = pt(1:end-1,:)-pt(2:end,:);d2 = (d1(:,1).^2+d1(:,2).^2+d1(:,3).^2).^0.5;
%     arclength = sum(d2)
%     pause
    
    transform = [1 0 0; 0 -1 0; 0 0 -1];
    flagellum = (transform*pt')';

    tmp = [R, wave_length, k, ang, phase, freq, blob_size_on_flag, numf, Ls, arclen, pt(1,:)];
    save([dir,'/helix_parameters_',num2str(phase),'.txt'],'tmp','-ascii')    
    
    %disp('Plot')
    
    %%%free swimmer
    N = size(flagellum,1);
    direction =  repmat([0,0,1],N,1);
    U = omega*cross(direction,flagellum);    
    
%     fig1=figure;
%     axes('Fontsize',24)
%     view(3)
%     hold all
%     plot3(flagellum(:,1),flagellum(:,2),flagellum(:,3),'-g','linewidth',4)
%     plot3(flagellum(1,1),flagellum(1,2),flagellum(1,3),'or','linewidth',4)
%     quiver3(flagellum(:,1),flagellum(:,2),flagellum(:,3),U(:,1),U(:,2),U(:,3),'b')
%     xlabel('x axis');ylabel('y axis');zlabel('z axis');
%     axis equal
% 
%     saveas(fig1,[dir,'/flag_phase_',num2str(phase),'.png'])
%     saveas(fig1,[dir,'/flag_phase_',num2str(phase),'.fig'],'fig') 

    
% %%%-----------------------Instantaneous Velocities----------------------%%%
%     UIB = zeros(N,3);
%     
%     utmp = [-w*(1 - exp(-k^2*s.^2))*R.*cos(k*s-w*t)...
%              w*(1 - exp(-k^2*s.^2))*R.*sin(k*s-w*t)...
%              zeros(size(s))];%Velocity
         
%    UIB = (transform*utmp')';
    
    UIB = U;
    save([dir,'/UIB_phase_',num2str(phase),'.txt'],'UIB','-ascii')

end
