function Rod_ploting1(q_node,Duplication_M,L_co,Option)
%plot the rod
%q_node: the node position and oritention
Theta_plot=(-pi:pi/36:pi)';

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
        if Option>0
            surf1(:,:,1)=r_node_1(1,:)+L_co*( cos(Theta_plot)*dr_dy_1(1,:)+sin(Theta_plot)*dr_dz_1(1,:) );
            surf1(:,:,2)=r_node_1(2,:)+L_co*( cos(Theta_plot)*dr_dy_1(2,:)+sin(Theta_plot)*dr_dz_1(2,:) );
            surf1(:,:,3)=r_node_1(3,:)+L_co*( cos(Theta_plot)*dr_dy_1(3,:)+sin(Theta_plot)*dr_dz_1(3,:) );
            C_surf(:,:,1)=0.85+0*surf1(:,:,1)
            C_surf(:,:,2)=0.33+0*surf1(:,:,1)
            C_surf(:,:,3)=0.1+0*surf1(:,:,1)
            surf(surf1(:,:,1),surf1(:,:,2),surf1(:,:,3),C_surf,FaceColor="interp",EdgeColor=[0 0 0])

        end
        if Option<2
            plot3(q_node(1,:,ii),q_node(2,:,ii),q_node(3,:,ii),'--','color',[0.0 0.45 0.75])
            plot3(-q_node(1,:,ii),q_node(2,:,ii),q_node(3,:,ii),'--','color',[0.0 0.45 0.75])
            plot3(-q_node(1,:,ii),-q_node(2,:,ii),q_node(3,:,ii),'--','color',[0.0 0.45 0.75])
            plot3( q_node(1,:,ii),-q_node(2,:,ii),q_node(3,:,ii),'--','color',[0.0 0.45 0.75])
        end
        shading interp

        %         surf2(:,:,1)=r_node_1(1,:)-L_co*( cos(Theta_plot)*dr_dy_1(1,:)+sin(Theta_plot)*dr_dz_1(1,:) );
        %         surf2(:,:,2)=r_node_1(2,:)-L_co*( cos(Theta_plot)*dr_dy_1(2,:)+sin(Theta_plot)*dr_dz_1(2,:) );
        %         surf2(:,:,3)=r_node_1(3,:)-L_co*( cos(Theta_plot)*dr_dy_1(3,:)+sin(Theta_plot)*dr_dz_1(3,:) );
        %         surf(surf2(:,:,1),surf2(:,:,2),surf2(:,:,3),0.7+0*surf2)
    end
end

end
