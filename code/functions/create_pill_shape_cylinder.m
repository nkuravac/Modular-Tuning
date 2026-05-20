function [p_total,con_pt] = create_pill_shape_cylinder(radius,height,h,phase,dir)
    %%%load points.dat on the unit sphere (SCVT)
    file_id=fopen('points.dat','r');
    number_of_points=fscanf(file_id,'%d',[1 1]);
    A=fscanf(file_id, '%g %g %g', [3 inf]);
    A=radius*A';
    
    %{
    fig=figure
    axes('fontsize',16)
    axis equal
    view(3)
    hold all
    [x,y,z] = sphere(30);
    s1 = surf(radius*x,radius*y,radius*z);
    shading interp
    %}
    
    num_fixed = 28;
    %{
    plot3(A(:,1),A(:,2),A(:,3),'or','markerfacecolor','r')
    %first point: north pole; num_fixed point: south pole
    plot3(A(1:num_fixed,1),A(1:num_fixed,2),A(1:num_fixed,3),'ok','markerfacecolor','k')
    xlabel('x axis')
    ylabel('y axis')
    zlabel('z axis')
    saveas(fig,'points_on_scaled_unit_sphere.fig','fig')
    saveas(fig,'points_on_scaled_unit_sphere.png','png')
    %}

    p_equator = A(2:num_fixed-1,:);
    ind_n = find(A(:,3)>0);
    p_north = A(ind_n,:); length(p_north)
    p_south = [p_north(:,1) p_north(:,2) -p_north(:,3)];
    
    %{
    fig=figure
    axes('fontsize',16)
    axis equal
    view(3)
    hold all
    [x,y,z] = sphere(30);
    s1 = surf(radius*x,radius*y,radius*z);
    shading interp
    
    plot3(p_north(:,1),p_north(:,2),p_north(:,3),'or','markerfacecolor','r')
    plot3(p_south(:,1),p_south(:,2),p_south(:,3),'xr','markerfacecolor','r')
    plot3(p_equator(:,1),p_equator(:,2),p_equator(:,3),'ok','markerfacecolor','k')
    xlabel('x axis')
    ylabel('y axis')
    zlabel('z axis')
    saveas(fig,'points_on_flipped_sphere.fig','fig')
    saveas(fig,'points_on_flipped_sphere.png','png')
    %}

    %%%Build top/bottom dome for a cylinder
    p_bot = p_south; %bottom circle at z = 0
    p_top = [p_north(:,1) p_north(:,2) height+p_north(:,3)]; %top circle at z = height
    
    %{
    fig=figure
    axes('fontsize',16)
    axis equal
    view(3)
    hold all
    
    plot3(p_top(:,1),p_top(:,2),p_top(:,3),'or','markerfacecolor','r')
    plot3(p_bot(:,1),p_bot(:,2),p_bot(:,3),'xr','markerfacecolor','r')
    plot3(p_equator(:,1),p_equator(:,2),p_equator(:,3),'ok','markerfacecolor','k')
    plot3(p_equator(:,1),p_equator(:,2),p_equator(:,3)+height,'ok','markerfacecolor','k')
    xlabel('x axis')
    ylabel('y axis')
    zlabel('z axis')
    saveas(fig,'points_on_split_sphere.fig','fig')
    saveas(fig,'points_on_split_sphere.png','png')
    %}

    %%%Build a cylinder    
    num_side = ceil(height/h) + 1;
    tmp = linspace(0,height,num_side);
    p_side = [];
    for i = 1:length(tmp)
        p_side = [p_side;[p_equator(:,1),p_equator(:,2),p_equator(:,3)+tmp(i)]];
    end
    
    p_total = [p_top;p_side;p_bot];
    [val,ind_min] = min(p_total(:,3));
    con_pt = p_total(ind_min,:);
    p_total(ind_min,:) = [];
    p_total = [p_total;con_pt];
    %size(p_total)
    
    %{
    %%%plot cell body
    fsize = 24;
    fig1=figure;
    axes('fontsize',fsize)
    view(3)
    axis equal
    hold all

    [X,Y,Z] = cylinder(radius);
    Z = Z*height;
    surf(X,Y,Z,'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none')

    plot3(p_top(:,1),p_top(:,2),p_top(:,3),'ok','MarkerFaceColor','k','markersize',4)
    plot3(p_bot(:,1),p_bot(:,2),p_bot(:,3),'og','MarkerFaceColor','g','markersize',4)
    
%     plot3(p_top(ind_boundary,1),p_top(ind_boundary,2),p_top(ind_boundary,3),'sy','MarkerFaceColor','y','markersize',4)
%     plot3(p_bot(ind_boundary,1),p_bot(ind_boundary,2),p_bot(ind_boundary,3),'sy','MarkerFaceColor','y','markersize',4)
% 
    plot3(p_side(:,1),p_side(:,2),p_side(:,3),'ob','MarkerFaceColor','b','markersize',3)
    
    xlabel('x axis (\mum)');ylabel('y axis (\mum)');zlabel('z axis (\mum)')
    title(['cell radius = ',num2str(radius),' \mum'])
    saveas(fig1,['cell_rad_',num2str(radius),'_height_',num2str(height),'.png'])
    saveas(fig1,['cell_rad_',num2str(radius),'_height_',num2str(height),'.fig'],'fig') 
    %}
end