% the script simulate the rod ring folding and twisting deformation 
% author: Qinghua Guan
%1： caculate the Bisymmetric configuration with different tendon lengths
%with a quarter of the ring
%2： analysis the stability of the Bisymmetric configuration by twist the
%ring circle as antisymmetric about the z-axis
%draw the map of the stability of Bisymmetric configurations
%% start the simulation
clear
clc
close all
% input mechinical information
R0_ring=50;%radius of the ring mm
L_sec=R0_ring*pi/2;
R_rod=5;%radius of the rod
E=20e6;%elastic modulus Pa
G=E/2*(1+0.4);%shearing modulus
Par_E=[E*pi*R_rod^2, E*pi*R_rod^4/64 0.2*G*pi*R_rod^4/32]';% the axis 
% bending and torsion stiffness

% initialization of the caculation
N_e=9;% the element number
L_e=L_sec/N_e;
N_node=5;%the sample node number to cacultate the elastic energy of each element
A=0.1e8;%the penalty parameter to constrain the frame vector Y Z as orthonormal
% it related to the elastic energy, when it too high the iteration would be too slow
% when it is too low, the frame vector Y Z won't be orthonormal anymore.

%% Evaluate the stability of every bisym config by applying twisting

% double the section to make it as anti-sym
% the origin config
load("Rod_BiSym_Folding.mat","Q_sec_M2")% load the twist or flat configuration as the original configuration
for mm=1:51
    Num_state(mm)=size(Q_sec_M2{mm},3);
end
for mm=45:51
%     figure(1);
%     clf
%     figure(2);
%     clf
%     figure(3);
%     clf
%     figure(4);
%     clf
    Q_sec_2=[];
    Fq_sec=[];
    U=[];
    for nn=1:Num_state(mm)
        mm
        nn
        Q_sec_1=Q_sec_M2{mm}(:,:,nn);
        N_e=2*(size(Q_sec_1,2)-1);% get the element number from the data
        Q_sec_2(1:3,:)=diag([-1,1,1])*Q_sec_1(1:3,:);
        Q_sec_2(4:6,:)=diag([-1,-1,-1])*diag([-1,1,1])*Q_sec_1(4:6,:);
        Q_sec_2(7:9,:)=diag([1,1,1])*diag([-1,1,1])*Q_sec_1(7:9,:);
        Q_sec_2(10:12,:)=diag([1,1,1])*diag([-1,1,1])*Q_sec_1(10:12,:);
        Q_sec=[Q_sec_1,flip(Q_sec_2(:,1:end-1),2)];

        L_sec_2=L_sec*2;
        L_e_2=L_sec_2/N_e;

        %plot the origin config
%         figure(1);
%         C_coord=[.0 .45 .74; .74 .45 .0];
%         L_co=5;
% 
%         r_node_1=Q_sec(1:3,:);
%         dr_dx_1=Q_sec(4:6,:);
%         dr_dy_1=Q_sec(7:9,:);
%         dr_dz_1=Q_sec(10:12,:);
% 
%         hold on
%         plot3(r_node_1(1,:),r_node_1(2,:),r_node_1(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
%         plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dy_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dy_1(2,:)], ...
%             [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dy_1(3,:)],'-','color',C_coord(1,:),'linewidth',1)
%         plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dz_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dz_1(2,:)], ...
%             [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dz_1(3,:)],'-','color',C_coord(2,:),'linewidth',1)
% 
%         r_node_2=[1 -1 1]'.*Q_sec(1:3,:);
%         dr_dx_2=[1 -1 1]'.*Q_sec(4:6,:);
%         dr_dy_2=[1 -1 1]'.*Q_sec(7:9,:);
%         dr_dz_2=[1 -1 1]'.*Q_sec(10:12,:);
% 
%         plot3(r_node_2(1,:),r_node_2(2,:),r_node_2(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
%         plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dy_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dy_2(2,:)], ...
%             [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dy_2(3,:)],'-','color',C_coord(1,:),'linewidth',1)
%         plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dz_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dz_2(2,:)], ...
%             [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dz_2(3,:)],'-','color',C_coord(2,:),'linewidth',1)
% 
%         box on
%         grid on
%         xlabel('X')
%         ylabel('Y')
%         zlabel('Z')
%         axis equal
%         axis([-70 70 -70 70 -70 70])
%         view([1 1 1])
%         hold off

        % Build the new boundary conditions for the stability analysis

        % start point x_disp in control and end point y_disp free
        %Using the de_demension matrix to define the constrain of the structure
        Proj_M_0=eye(12*(N_e+1),12*(N_e+1));%dedemesionalize the variables

        % start point
        D_fix=0*Q_sec;%fixed dimensions/varibales
        D_fix(1,1)=1;%control the X disp
        D_fix(2,1)=1;%fix the Y disp
        D_fix(3,1)=1;%fix the  Z disp
        % the local x y z axis in control
        D_fix(4:6,1)=1;%x axis in control
        D_fix(7:9,1)=1;%y axis  in control
        D_fix(10:12,1)=1;%z axis  in control

        % end point
        D_fix(1,end)=1;%fix the X disp
        D_fix(2,end)=1;%control the Y disp
        D_fix(3,end)=1;%control the Z disp
        % local x y z aixs in control
        D_fix([4 5 6],end)=1;%fix the rotation
        D_fix([7 8 9],end)=1;%fix the rotation
        D_fix([10 11 12],end)=1;%fix the rotation
        n_m=N_e/2+1;%the seq numb of middle point

        % the middle point position
        % based on anti sym condition
        % the distance of middle point to the z axis is constrained: x^2+y^2=C
        % introduce the R and theta to replace the x y disp of middle point
        r_m_0=Q_sec(2,n_m);
        theta_m_0=pi/2;
        % the present First-order expansion relationship between d(theta,R) and d(x,y)
        Con_M=[-sin(theta_m_0)/r_m_0 cos(theta_m_0)/r_m_0; cos(theta_m_0) sin(theta_m_0)];
        Proj_M_0((n_m-1)*12+[1:2],(n_m-1)*12+[1:2])=Con_M;
        %fix R: the distance to the z axis
        D_fix(2,n_m)=1;%fix the R
        T_angle=linspace(0,1*pi/180,11);
        r1_frame_0=reshape(Q_sec(4:12,1),3,3);
        r2_frame_0=reshape(Q_sec(4:12,end),3,3);
%         U=[];
%         E_elastic=[];
        %

        %
%         figure(2);
%         clf
%         box on
%         grid on
%         xlabel('X')
%         ylabel('Y')
%         zlabel('Z')
%         axis equal
%         axis([-70 70 -70 70 -70 70])
%         view([1 1 1])
%         hold off

        for ii=1:5
            mm
            nn
            ii
            if ii~=1
                Q_sec=Q_sec_M3{mm,nn}(:,:,ii-1);
            end
            %redefine the end points and frames
            alpha=T_angle(ii);
            RM_X=[1 0 0;
                0 cos(alpha) -sin(alpha);
                0 sin(alpha)  cos(alpha)];
            r1_frame=RM_X*r1_frame_0;
            Q_sec(4:12,1)=r1_frame(:);
            alpha=-T_angle(ii);
            RM_X=[1 0 0;
                0 cos(alpha) -sin(alpha);
                0 sin(alpha)  cos(alpha)];
            r2_frame=RM_X*r2_frame_0;
            Q_sec(4:12,end)=r2_frame(:);

            kk=0;
            flag1=1;
            flag3=0;
            L_step=0.01;
            while (kk==0||flag1>0.0001*L_step||kk<100)&&kk<400
                kk=kk+1;
                [dFq_dq_sec,Fq_sec(:,ii),q_node,U(ii),E_elastic{mm}(:,ii,nn)] = Jocob_rod_sec(Q_sec,N_e,L_e_2,N_node,Par_E,A);
                % update the boundary information
                if r_m_0<10e-4
                    r_m=0;
                    theta_m=theta_m_0;
                    D_fix(1,n_m)=1;
                else
                    D_fix(1,n_m)=0;
                    r_m=sqrt(sum(Q_sec(1:2,n_m).^2));
                    theta_m=atan2(Q_sec(2,n_m),Q_sec(1,n_m));
                    % the present First-order expansion relationship between d(theta,R) and d(x,y)
                    Con_M=[-sin(theta_m)/r_m cos(theta_m)/r_m; cos(theta_m) sin(theta_m)];
                    Proj_M_0((n_m-1)*12+[1:2],(n_m-1)*12+[1:2])=Con_M;
                end
                B_index=find(D_fix==0);% find the fixed demensions
                Proj_M=Proj_M_0(:,B_index);% extract the projection matrix of free paramters
                dFq_dq_sec_2=Proj_M'*dFq_dq_sec*Proj_M;
                dQ_sec_2=-pinv(dFq_dq_sec_2)*Proj_M'*Fq_sec(:,ii);
                dQ_sec=L_step*Proj_M*dQ_sec_2;
                dQ_sec=reshape(dQ_sec,12,[]);
                Q_sec=Q_sec+dQ_sec;
                % canonicalize the diplacement based on the varition of angle
                if r_m_0<10e-4
                    Q_sec(1:2,n_m)=[0,0]';
                else
                    d_theta=Con_M(1,:)*dQ_sec([1:2],n_m);
                    theta_m=theta_m+d_theta;
                    Q_sec(1:2,n_m)=r_m_0*[cos(theta_m);sin(theta_m)];
                end
                Q_sec_M3{mm,nn}(:,:,ii)=Q_sec;
                flag1=sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e+1);
                flag2=mean(abs(Proj_M'*Fq_sec(:,ii)));
%                 if kk==4000
%                     flag3=1
%                 elseif mod(kk,100)==0
%                     kk
%                     flag1
%                     for jj=1:N_e
% 
%                         figure(2);
%                         r_node_1=q_node(1:3,:,jj);
%                         dr_dx_1=q_node(4:6,:,jj);
%                         dr_dy_1=q_node(7:9,:,jj);
%                         dr_dz_1=q_node(10:12,:,jj);
% 
%                         hold on
%                         plot3(r_node_1(1,:),r_node_1(2,:),r_node_1(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
%                         plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dy_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dy_1(2,:)], ...
%                             [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dy_1(3,:)],'-','color',C_coord(1,:),'linewidth',1)
%                         plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dz_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dz_1(2,:)], ...
%                             [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dz_1(3,:)],'-','color',C_coord(2,:),'linewidth',1)
% 
%                         %                 r_node_2=[-1 -1 1]'.*q_node(1:3,:,jj);
%                         %                 dr_dx_2=[-1 -1 1]'.*q_node(4:6,:,jj);
%                         %                 dr_dy_2=[-1 -1 1]'.*q_node(7:9,:,jj);
%                         %                 dr_dz_2=[-1 -1 1]'.*q_node(10:12,:,jj);
%                         %
%                         %                 plot3(r_node_2(1,:),r_node_2(2,:),r_node_2(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
%                         %                 plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dy_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dy_2(2,:)], ...
%                         %                     [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dy_2(3,:)],'-','color',C_coord(1,:),'linewidth',1)
%                         %                 plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dz_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dz_2(2,:)], ...
%                         %                     [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dz_2(3,:)],'-','color',C_coord(2,:),'linewidth',1)
%                         %                 hold off
% 
%                     end
%                     figure(3);
%                     hold on
%                     plot(E_elastic{mm}(1,1:ii,nn),'-*')
%                     hold off
%                     pause(0.1)
%                 end
            end

%             for jj=1:N_e
% 
%                 figure(2);
%                 r_node_1=q_node(1:3,:,jj);
%                 dr_dx_1=q_node(4:6,:,jj);
%                 dr_dy_1=q_node(7:9,:,jj);
%                 dr_dz_1=q_node(10:12,:,jj);
% 
%                 hold on
%                 plot3(r_node_1(1,:),r_node_1(2,:),r_node_1(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
%                 plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dy_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dy_1(2,:)], ...
%                     [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dy_1(3,:)],'-','color',C_coord(1,:),'linewidth',1)
%                 plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dz_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dz_1(2,:)], ...
%                     [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dz_1(3,:)],'-','color',C_coord(2,:),'linewidth',1)
% 
%                 r_node_2=[-1 -1 1]'.*q_node(1:3,:,jj);
%                 dr_dx_2=[-1 -1 1]'.*q_node(4:6,:,jj);
%                 dr_dy_2=[-1 -1 1]'.*q_node(7:9,:,jj);
%                 dr_dz_2=[-1 -1 1]'.*q_node(10:12,:,jj);
% 
%                 plot3(r_node_2(1,:),r_node_2(2,:),r_node_2(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
%                 plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dy_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dy_2(2,:)], ...
%                     [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dy_2(3,:)],'-','color',C_coord(1,:),'linewidth',1)
%                 plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dz_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dz_2(2,:)], ...
%                     [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dz_2(3,:)],'-','color',C_coord(2,:),'linewidth',1)
%                 hold off
% 
%             end

%             figure(3);
%             hold on
%             plot(E_elastic{mm}(1,1:ii,nn),'-*')
%             hold off
%             pause(0.1)

        end

%         figure(4);
%         hold on
%         plot3(nn*ones(1,5),1:5,E_elastic{mm}(1,:,nn),'-*')
%         hold off
%         pause(0.1)
%         box on
%         grid on   
    end
%     save("Rod_BiSym_stability_1.mat",'E_elastic',"Q_sec_M3")
end


%%
load("Rod_BiSym_stability_1.mat",'E_elastic')
load("Rod_BiSym_Folding.mat","Q_sec_M2")% load the twist or flat configuration as the original configuration

hold on
for mm=1:51
    for nn=1:size(E_elastic{mm},3)
        Sign{mm}(nn,:)=E_elastic{mm}(1,2:end,nn)>E_elastic{mm}(1,1:end-1,nn);
    end
end
hold off
axis equal

figure(5)
hold on
for mm=1:51
   for nn=1:size(E_elastic{mm},3)
        if sum(Sign{mm}(nn,:))==4
            plot(Q_sec_M2{mm}(1,1,nn),Q_sec_M2{mm}(2,end,nn),'bo')
            plot(Q_sec_M2{mm}(2,end,nn),Q_sec_M2{mm}(1,1,nn),'bo')
        else 
            scatter(Q_sec_M2{mm}(1,1,nn),Q_sec_M2{mm}(2,end,nn),'r')
            scatter(Q_sec_M2{mm}(2,end,nn),Q_sec_M2{mm}(1,1,nn),'r')
        end
    end
end
hold off
axis equal
%%
load("Rod_BiSym_stability_1.mat",'E_elastic')
load("Rod_BiSym_Folding.mat","Q_sec_M2")% load the twist or flat configuration as the original configuration

Sign=[]
for mm=1:51
    for nn=1:size(E_elastic{mm},3)
        Sign{mm}(nn,:)=E_elastic{mm}(1,2:4,nn)>E_elastic{mm}(1,1:3,nn);
    end
end
hold off
axis equal

figure(6)
hold on
for mm=1:51
   for nn=1:size(E_elastic{mm},3)
        if sum(Sign{mm}(nn,:))==3
            plot(Q_sec_M2{mm}(1,1,nn),Q_sec_M2{mm}(2,end,nn),'bo')
            plot(Q_sec_M2{mm}(2,end,nn),Q_sec_M2{mm}(1,1,nn),'bo')
        else 
            scatter(Q_sec_M2{mm}(1,1,nn),Q_sec_M2{mm}(2,end,nn),'r')
            scatter(Q_sec_M2{mm}(2,end,nn),Q_sec_M2{mm}(1,1,nn),'r')
        end
    end
end
hold off
axis equal

%%
%%
load("Rod_BiSym_stability_1.mat",'E_elastic')
load("Rod_BiSym_Folding.mat","Q_sec_M2")% load the twist or flat configuration as the original configuration

Sign=[]
for mm=1:51
    for nn=1:size(E_elastic{mm},3)
        Sign{mm}(nn,:)=E_elastic{mm}(1,2,nn)>E_elastic{mm}(1,1,nn);
    end
end
hold off
axis equal

figure(8)
hold on
for mm=1:51
   for nn=1:size(E_elastic{mm},3)
        if sum(Sign{mm}(nn,:))==1
            plot(Q_sec_M2{mm}(1,1,nn),Q_sec_M2{mm}(2,end,nn),'bo')
            plot(Q_sec_M2{mm}(2,end,nn),Q_sec_M2{mm}(1,1,nn),'bo')
        else 
            scatter(Q_sec_M2{mm}(1,1,nn),Q_sec_M2{mm}(2,end,nn),'r')
            scatter(Q_sec_M2{mm}(2,end,nn),Q_sec_M2{mm}(1,1,nn),'r')
        end
    end
end
hold off
axis equal

%%
mm=43
figure(12)
for nn=1:68
    plot(E_elastic{mm}(1,:,nn))
    title(['mm=',num2str(mm),'nn=',num2str(nn)])
end

% Q_sec=Q_sec_M3{mm}(:,:,nn)



  %plot the origin config
%         figure(1);
%         C_coord=[.0 .45 .74; .74 .45 .0];
%         L_co=5;
% 
%         r_node_1=Q_sec(1:3,:);
%         dr_dx_1=Q_sec(4:6,:);
%         dr_dy_1=Q_sec(7:9,:);
%         dr_dz_1=Q_sec(10:12,:);
% 
%         hold on
%         plot3(r_node_1(1,:),r_node_1(2,:),r_node_1(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
%         plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dy_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dy_1(2,:)], ...
%             [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dy_1(3,:)],'-','color',C_coord(1,:),'linewidth',1)
%         plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dz_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dz_1(2,:)], ...
%             [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dz_1(3,:)],'-','color',C_coord(2,:),'linewidth',1)
% 
%         r_node_2=[1 -1 1]'.*Q_sec(1:3,:);
%         dr_dx_2=[1 -1 1]'.*Q_sec(4:6,:);
%         dr_dy_2=[1 -1 1]'.*Q_sec(7:9,:);
%         dr_dz_2=[1 -1 1]'.*Q_sec(10:12,:);
% 
%         plot3(r_node_2(1,:),r_node_2(2,:),r_node_2(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
%         plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dy_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dy_2(2,:)], ...
%             [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dy_2(3,:)],'-','color',C_coord(1,:),'linewidth',1)
%         plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dz_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dz_2(2,:)], ...
%             [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dz_2(3,:)],'-','color',C_coord(2,:),'linewidth',1)
% 
%         box on
%         grid on
%         xlabel('X')
%         ylabel('Y')
%         zlabel('Z')
%         axis equal
%         axis([-70 70 -70 70 -70 70])
%         view([1 1 1])
%         hold off








