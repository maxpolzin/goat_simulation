function Rod_ploting(q_node,Duplication_M,L_co,C_coord)
%plot the rod
%q_node: the node position and oritention
Theta_plot=(-pi:pi/6:pi)';

for jj=1:size(Duplication_M,1)
    for ii=1:size(q_node,3)
        r_node_1=Duplication_M(jj,:)'.*q_node(1:3,:,ii);
        dr_dx_1=Duplication_M(jj,:)'.*q_node(4:6,:,ii);
        dr_dy_1=Duplication_M(jj,:)'.*q_node(7:9,:,ii);
        dr_dz_1=Duplication_M(jj,:)'.*q_node(10:12,:,ii);        
%         plot3(r_node_1(1,:),r_node_1(2,:),r_node_1(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
%         plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dy_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dy_1(2,:)], ...
%             [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dy_1(3,:)],'-','color',C_coord(1,:),'linewidth',1)
%         plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dz_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dz_1(2,:)], ...
%             [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dz_1(3,:)],'-','color',C_coord(2,:),'linewidth',1)

        surf1(:,:,1)=r_node_1(1,:)+L_co*( cos(Theta_plot)*dr_dy_1(1,:)+sin(Theta_plot)*dr_dz_1(1,:) );
        surf1(:,:,2)=r_node_1(2,:)+L_co*( cos(Theta_plot)*dr_dy_1(2,:)+sin(Theta_plot)*dr_dz_1(2,:) );
        surf1(:,:,3)=r_node_1(3,:)+L_co*( cos(Theta_plot)*dr_dy_1(3,:)+sin(Theta_plot)*dr_dz_1(3,:) );
        surf(surf1(:,:,1),surf1(:,:,2),surf1(:,:,3),10*abs(Theta_plot)+0*surf1(:,:,3),FaceColor="interp",EdgeColor=[0 0 0])
        shading interp

%         surf2(:,:,1)=r_node_1(1,:)-L_co*( cos(Theta_plot)*dr_dy_1(1,:)+sin(Theta_plot)*dr_dz_1(1,:) );
%         surf2(:,:,2)=r_node_1(2,:)-L_co*( cos(Theta_plot)*dr_dy_1(2,:)+sin(Theta_plot)*dr_dz_1(2,:) );
%         surf2(:,:,3)=r_node_1(3,:)-L_co*( cos(Theta_plot)*dr_dy_1(3,:)+sin(Theta_plot)*dr_dz_1(3,:) );
%         surf(surf2(:,:,1),surf2(:,:,2),surf2(:,:,3),0.7+0*surf2)
    end
end

end
