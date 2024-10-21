%% build the whole rod ring to caculate the twisting folding
%% start the simulation
clear
clc
close all
%%
% input mechinical information
R0_ring=50;%radius of the ring mm
L_sec_0=R0_ring*pi/2;% set as the 1/4circle 
R_rod=5;%radius of the rod
E=20e6;%elastic modulus Pa
G=E/2*(1+0.4);%shearing modulus
Par_E=[E*pi*R_rod^2, E*pi*R_rod^4/64 0.2*G*pi*R_rod^4/32]';% the axis 
% bending and torsion stiffness

% initialization of the caculation
N_e_0=8;% the element number
L_e_0=L_sec_0/N_e_0;
N_node=5;%the sample node number to cacultate the elastic energy of each element
A=0.1e8;%the penalty parameter to constrain the frame vector Y Z as orthonormal
% it related to the elastic energy, when it too high the iteration would be too slow
% when it is too low, the frame vector Y Z won't be orthonormal anymore.
L_co=2;%the length of the coordinates
C_coord=[.0 .45 .74; .74 .45 .0];%the color of the coordinates
Duplication_M=[1 1 1;-1 -1 1];

%% build the whole rod ring to caculate the twisting folding

load("Rod_BiSym_Folding.mat","Q_sec_M3")% load the twist or flat configuration as the original configuration
Q_sec=[];
Q_sec_1=[];
Q_sec_2=[];
Fq_sec=[];
U=[];

mm=41
nn=5
ll=24

Q_sec_1=Q_sec_M3{mm}(:,:,ll,nn);
N_e_2=2*(size(Q_sec_1,2)-1);% get the element number from the data
Q_sec_2(1:3,:)=diag([-1,-1,1])*Q_sec_1(1:3,:);
Q_sec_2(4:6,:)=diag([-1,-1,1])*Q_sec_1(4:6,:);
Q_sec_2(7:9,:)=diag([-1,-1,1])*Q_sec_1(7:9,:);
Q_sec_2(10:12,:)=diag([-1,-1,1])*Q_sec_1(10:12,:);
Q_sec=[Q_sec_1,Q_sec_2(:,2:end)];

L_sec_2=L_sec_0*4;
L_e_2=L_sec_2/N_e_2;

n_1q=N_e_2/4+1;%the seq numb of 1/4 point
n_2q=2*N_e_2/4+1;%the seq numb of 2/4 point
n_3q=3*N_e_2/4+1;%the seq numb of 3/4 point

%plot the origin config
figure(8);
C_coord=[.0 .45 .74; .74 .45 .0];
L_co=5;

r_node_1=Q_sec(1:3,:);
dr_dx_1=Q_sec(4:6,:);
dr_dy_1=Q_sec(7:9,:);
dr_dz_1=Q_sec(10:12,:);

hold on
plot3(r_node_1(1,:),r_node_1(2,:),r_node_1(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dy_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dy_1(2,:)], ...
    [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dy_1(3,:)],'-','color',C_coord(1,:),'linewidth',1)
plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dz_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dz_1(2,:)], ...
    [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dz_1(3,:)],'-','color',C_coord(2,:),'linewidth',1)

plot3(r_node_1(1,[1,n_2q]),r_node_1(2,[1,n_2q]),r_node_1(3,[1,n_2q]),'-','color','b','linewidth',1)
plot3(r_node_1(1,[n_1q,n_3q]),r_node_1(2,[n_1q,n_3q]),r_node_1(3,[n_1q,n_3q]),'-','color','r','linewidth',1)
scatter3(r_node_1(1,1),r_node_1(2,1),r_node_1(3,1),'r')
scatter3(r_node_1(1,n_2q),r_node_1(2,n_2q),r_node_1(3,n_2q),'g')

box on
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
set(gca,'Zdir','reverse')
axis equal
axis([-70 70 -70 70 -70 70])
view([1 1 1])
hold off


%% Build the new boundary conditions for the stability analysis

A=0.05e8;


% start point x_disp in control and end point y_disp free
%Using the de_demension matrix to define the constrain of the structure
Proj_M_0=eye(12*(N_e_2+1),12*(N_e_2+1));%dedemesionalize the variables


% based on anti sym condition

% the start and end point 
% start point
D_fix=0*Q_sec;%fixed dimensions/varibales
D_fix(1,1)=1;%control the X disp
D_fix(2,1)=1;%fix the Y disp
D_fix(3,1)=1;%fix the  Z disp

%control the local rotation by constrain the local z axis

% D_fix(4:6,1)=1;%local z axis in control 

% the local x y z axis in control
% because the equalism relationship
%D_fix(1:3,1)=[1 1 1].*D_fix(4:6,end);%x axis in control 
%D_fix(4:6,1)=[1 1 1].*D_fix(4:6,end);%x axis in control 
%D_fix(7:9,1)=[1 1 1].*D_fix(4:6,end);%x axis in control 
%D_fix(10:12,1)=[1 1 1].*D_fix(4:6,end);%x axis in control 
% sus
Con_M_1=diag(repmat([1 1 1],1,4));
% Proj_M_0(1:12,end-11:end)=0
Proj_M_0(end-11:end,1:12)=Con_M_1;

% end point in constrain
D_fix(1,end)=1;%fix the X disp
D_fix(2,end)=1;%control the Y disp
D_fix(3,end)=1;%control the Z disp

% local x y z aixs in control
D_fix([4 5 6],end)=1;%fix the rotation
D_fix([7 8 9],end)=1;%fix the rotation
D_fix([10 11 12],end)=1;%fix the rotation
r1_frame_0=reshape(Q_sec(4:12,1),3,3);

%2/4 point constrain the local x in the plane A of X-[0 -1 0]
D_fix(1,n_2q)=1;%fix the X disp
D_fix(2,n_2q)=1;%control the Y disp
D_fix(3,n_2q)=1;%control the Z disp

% introduce the new local coordinate at n_2q
% r_frame_2q_0=reshape(Q_sec(4:12,n_2q),3,3);
% r_frame_2q_1(:,2)=[0 1 0]';
% V_x=r_frame_2q_0(:,1).*[1 0 1]';
% r_frame_2q_1(:,1)=V_x/norm(V_x);
% r_frame_2q_1(:,3)=cross(r_frame_2q_1(:,1),r_frame_2q_1(:,2));
% R_M_con3=r_frame_2q_1'
% Con_M_3=zeros(9,9);
% Con_M_3(1:3,1:3)=R_M_con3;
% Con_M_3(4:6,4:6)=R_M_con3;
% Con_M_3(7:9,7:9)=R_M_con3;
% D_fix(6,n_2q)=1;%constrain the new local x axis in the plane A


% the 1/4 and 3/4 point position
% based on constrain of the length
% the distance of 1/4 point to the 3/4 axis is constrained: x^2+y^2=C
% introduce the R and theta to replace the x y disp of 3/4 point
V_13q=Q_sec(1:3,n_3q)-Q_sec(1:3,n_1q);
r_xy_m_0=sqrt(sum(V_13q(1:2).^2));
r_m_0=sqrt(sum(V_13q.^2));
theta_m_0=atan2(V_13q(2),V_13q(1));
phi_m_0=atan2(V_13q(3),r_xy_m_0);
% the present First-order expansion relationship between d(xm,ym,zm,theta,phi,R) and d(x1,y1,z1,x3,y3,z3)
% x1=xm-cos(theta)*cos(phi)*R/2 
% y1=ym-sin(theta)*cos(phi)*R/2
% z1=zm-sin(phi)*R/2
% x3=xm+cos(theta)*cos(phi)*R/2 
% y3=ym+sin(theta)*cos(phi)*R/2
% z3=zm+sin(phi)*R/2
%dx3=dxm+(-sin(theta)*cos(phi)*R*dtheta-cos(theta)*sin(phi)*R*dphi+cos(theta)*cos(phi)*dR)
%dy3=dym+(cos(theta)*R*dtheta+sin(theta)*dR)
%dz3=dzm+(cos(phi)*R*dphi+sin(phi)dR)
Con_M_2_0=[1 0 0,   0,  0, 0;
           0 1 0,   0,  0, 0;
           0 0 1,   0,  0, 0;
           1 0 0,   -sin(theta_m_0)*cos(phi_m_0)*r_m_0, -cos(theta_m_0)*sin(phi_m_0)*r_m_0, cos(theta_m_0)*cos(phi_m_0);
           0 1 0,   cos(theta_m_0)*cos(phi_m_0)*r_m_0,  -sin(theta_m_0)*sin(phi_m_0)*r_m_0, sin(theta_m_0)*cos(phi_m_0);
           0 0 1,    0  cos(phi_m_0)*r_m_0, sin(phi_m_0);]
Con_M_2_1=[                                  1,                                    0,                   0,                                    0,                                    0,                  0;
                                             0,                                    1,                   0,                                    0,                                    0,                  0;
                                             0,                                    0,                   1,                                    0,                                    0,                  0;
           sin(theta_m_0)/(r_m_0*cos(phi_m_0)), -cos(theta_m_0)/(r_m_0*cos(phi_m_0)),                   0, -sin(theta_m_0)/(r_m_0*cos(phi_m_0)),  cos(theta_m_0)/(r_m_0*cos(phi_m_0)),                  0;
           (cos(theta_m_0)*sin(phi_m_0))/r_m_0,  (sin(phi_m_0)*sin(theta_m_0))/r_m_0, -cos(phi_m_0)/r_m_0, -(cos(theta_m_0)*sin(phi_m_0))/r_m_0, -(sin(phi_m_0)*sin(theta_m_0))/r_m_0, cos(phi_m_0)/r_m_0;
          -cos(phi_m_0)*cos(theta_m_0),         -cos(phi_m_0)*sin(theta_m_0),       -sin(phi_m_0),          cos(phi_m_0)*cos(theta_m_0),          cos(phi_m_0)*sin(theta_m_0),       sin(phi_m_0);];

Proj_M_0((n_1q-1)*12+[1:3],(n_1q-1)*12+[1:3])=Con_M_2_0(1:3,1:3);
Proj_M_0((n_1q-1)*12+[1:3],(n_3q-1)*12+[1:3])=Con_M_2_0(1:3,4:6);

Proj_M_0((n_3q-1)*12+[1:3],(n_1q-1)*12+[1:3])=Con_M_2_0(4:6,1:3);
Proj_M_0((n_3q-1)*12+[1:3],(n_3q-1)*12+[1:3])=Con_M_2_0(4:6,4:6);

%fix R: the distance to the z axis
D_fix(1,n_1q)=1;%fix the xm
D_fix(3,n_3q)=1;%fix the R

Disp_x_1q_0=Q_sec(1,n_1q);
%%
Disp_x_1q=linspace(0,5,11);
L_co=2;
% Q_sec(:,1)
% Q_sec(:,end);
% figure(11);
%     clf
figure(10);
clf
% T_angle=linspace(0,30*pi/180,11);
for ii=1
    ii    
    figure(9);
    clf
    if ii~=1
        Q_sec=Q_sec_M4{mm}(:,:,ii-1,ll,nn);
    end

%     %redefine the local z axis of start/end points
%     R_M_z=[cos(T_angle(ii))  0  -sin(T_angle(ii)) ;
%            0                 1               0;
%            sin(T_angle(ii))  0  cos(T_angle(ii)) ;];
%     r1_frame_1=r1_frame_0*R_M_z;
%     Q_sec(10:12,1)=r1_frame_1(:,3);
%     Q_sec(10:12,end)=r1_frame_1(:,3);  

    Q_sec(1,n_1q)=Disp_x_1q_0-Disp_x_1q(ii); 
    kk=0;
    flag1=1;
    L_step=0.02;
   
    while (kk==0||flag1>0.0001*L_step||kk<2/L_step)&&kk<30/L_step
        kk=kk+1;
        [dFq_dq_sec,Fq_sec(:,ii),q_node,U(ii),E_elastic_4{mm}(:,ii,ll,nn)] = Jocob_rod_sec(Q_sec,N_e_2,L_e_2,N_node,Par_E,A);
        % update the boundary information
        % update the boundary information
        if r_m_0<10e-2
            r_m=0;
            theta_m=theta_m_0;
            D_fix(1:3,n_3q)=1;
        else
            D_fix(1:2,n_3q)=0;
            V_13q=Q_sec(1:3,n_3q)-Q_sec(1:3,n_1q);
            r_xy_m=sqrt(sum(V_13q(1:2).^2));
            r_m=sqrt(sum(V_13q.^2));
            phi_m=atan2(V_13q(3),r_xy_m);
            if r_xy_m<10e-2
              D_fix(1:2,n_3q)=0;
              D_fix(3,n_3q)=1;
              Con_M_2_1=[ 0,   0,   0;
                          0,   0, 0;
                          0,  0,  1];
              Con_M_2_2=[ 1,   0,   0;
                          0,   1, 0;
                          0,  0,  1];
              Proj_M_0((n_3q-1)*12+[1:3],(n_1q-1)*12+[1:3])=Con_M_2_1;
              Proj_M_0((n_3q-1)*12+[1:3],(n_3q-1)*12+[1:3])=Con_M_2_2;
            else
                D_fix(1,n_3q)=0;
                D_fix(2,n_3q)=0;
                D_fix(3,n_3q)=1;
                theta_m=atan2(V_13q(2),V_13q(1));                
                % the present First-order expansion relationship between d(theta,phi,R) and d(x,y,z)
                Con_M_2_0=[ 1 0 0,   0,  0, 0;
                            0 1 0,   0,  0, 0;
                            0 0 1,   0,  0, 0;
                            1 0 0,   -sin(theta_m)*cos(phi_m)*r_m, -cos(theta_m)*sin(phi_m)*r_m, cos(theta_m)*cos(phi_m);
                            0 1 0,   cos(theta_m)*cos(phi_m)*r_m,  -sin(theta_m)*sin(phi_m)*r_m, sin(theta_m)*cos(phi_m);
                            0 0 1,    0,  cos(phi_m)*r_m, sin(phi_m);];
                Con_M_2_1=[                            1,                              0,               0,                              0,                              0,              0;
                                                       0,                              1,               0,                              0,                              0,              0;
                                                       0,                              0,               1,                              0,                              0,              0;
                           sin(theta_m)/(r_m*cos(phi_m)), -cos(theta_m)/(r_m*cos(phi_m)),               0, -sin(theta_m)/(r_m*cos(phi_m)),  cos(theta_m)/(r_m*cos(phi_m)),              0;
                           (cos(theta_m)*sin(phi_m))/r_m,  (sin(phi_m)*sin(theta_m))/r_m, -cos(phi_m)/r_m, -(cos(theta_m)*sin(phi_m))/r_m, -(sin(phi_m)*sin(theta_m))/r_m, cos(phi_m)/r_m;
                                -cos(phi_m)*cos(theta_m),       -cos(phi_m)*sin(theta_m),     -sin(phi_m),        cos(phi_m)*cos(theta_m),        cos(phi_m)*sin(theta_m),     sin(phi_m)];
                Proj_M_0((n_1q-1)*12+[1:3],(n_1q-1)*12+[1:3])=Con_M_2_0(1:3,1:3);
                Proj_M_0((n_1q-1)*12+[1:3],(n_3q-1)*12+[1:3])=Con_M_2_0(1:3,4:6);
                Proj_M_0((n_3q-1)*12+[1:3],(n_1q-1)*12+[1:3])=Con_M_2_0(4:6,1:3);
                Proj_M_0((n_3q-1)*12+[1:3],(n_3q-1)*12+[1:3])=Con_M_2_0(4:6,4:6);
            end
        end
        B_index=find(D_fix==0);% find the fixed demensions
        Proj_M=Proj_M_0(:,B_index);% extract the projection matrix of free paramters
        dFq_dq_sec_2=Proj_M'*dFq_dq_sec*Proj_M;
        dQ_sec_2=-inv(dFq_dq_sec_2)*Proj_M'*Fq_sec(:,ii);
        dQ_sec=L_step*Proj_M*dQ_sec_2;
        dQ_sec=reshape(dQ_sec,12,[]);        
        Q_sec=Q_sec+dQ_sec;
        V_13q=Q_sec(1:3,n_3q)-Q_sec(1:3,n_1q);
%         % canonicalize the diplacement based on the varition of angle
%         if r_m_0<10e-4
%             Q_sec(1:3,n_3q)=Q_sec(1:3,n_1q);
%         elseif r_xy_m<10e-4
%             V_13q=Q_sec(1:3,n_3q)-Q_sec(1:3,n_1q);
%             Q_sec(1:3,n_3q)=Q_sec(1:3,n_1q)+r_m_0*V_13q/norm(V_13q);
%         else
% %             d_theta=Con_M_2_2(1,1:3)*(dQ_sec([1:3],n_3q)-dQ_sec([1:3],n_1q));
% %             d_phi=Con_M_2_2(2,1:3)*(dQ_sec([1:3],n_3q)-dQ_sec([1:3],n_1q));
% %             theta_m=theta_m+d_theta;
% %             phi_m=phi_m+d_phi;
% %             Q_sec(1:3,n_3q)=Q_sec([1:3],n_1q)+r_m_0*[cos(theta_m)*cos(phi_m);sin(theta_m)*cos(phi_m);sin(phi_m)];
%         end
        Q_sec_M4{mm}(:,:,ii,ll,nn)=Q_sec;
        flag1=sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e_2+1);
        flag2=mean(abs(Proj_M'*Fq_sec(:,ii)));
        if kk==5000
            flag3=1
        elseif mod(kk,50)==0
            fprintf('mm=%u,nn=%u,ii=%u, kk=%u, flag1=%f, D13=%f \n',mm,nn,ii,kk,flag1,norm(V_13q)) 
            figure(9);
            hold on
            for jj=1:N_e_2
                r_node_1=q_node(1:3,:,jj);
                dr_dx_1=q_node(4:6,:,jj);
                dr_dy_1=q_node(7:9,:,jj);
                dr_dz_1=q_node(10:12,:,jj);                
                plot3(r_node_1(1,:),r_node_1(2,:),r_node_1(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
                plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dy_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dy_1(2,:)], ...
                    [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dy_1(3,:)],'-','color',C_coord(1,:),'linewidth',1)
                plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dz_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dz_1(2,:)], ...
                    [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dz_1(3,:)],'-','color',C_coord(2,:),'linewidth',1)         
                
            end
            
            plot3(Q_sec(1,[1,n_2q]),Q_sec(2,[1,n_2q]),Q_sec(3,[1,n_2q]),'-','color','b','linewidth',1)            
            plot3(Q_sec(1,[n_1q,n_3q]),Q_sec(2,[n_1q,n_3q]),Q_sec(3,[n_1q,n_3q]),'-','color','r','linewidth',1)
            scatter3(Q_sec(1,1),Q_sec(2,1),Q_sec(3,1),'r')
            scatter3(Q_sec(1,n_2q),Q_sec(2,n_2q),Q_sec(3,n_2q),'g')
            hold off
            set(gca,'Zdir','reverse')
            axis equal
            box on
            grid on
            view([1 1 1])

            figure(10);
            hold on
            plot(E_elastic_4{mm}(1,1:ii,ll,nn),'-*')
            hold off
            xlabel('x')
            ylabel('y')
            zlabel('z')
            pause(0.01)
        end

    end


    figure(11);
    hold on
    Rod_ploting(q_node,Duplication_M,L_co,C_coord)
    plot3(Q_sec(1,[1,n_2q]),Q_sec(2,[1,n_2q]),Q_sec(3,[1,n_2q]),'-','color','b','linewidth',1)
    plot3(Q_sec(1,[n_1q,n_3q]),Q_sec(2,[n_1q,n_3q]),Q_sec(3,[n_1q,n_3q]),'-','color','r','linewidth',1)
     scatter3(Q_sec(1,1),Q_sec(2,1),Q_sec(3,1),'r')
            scatter3(Q_sec(1,n_2q),Q_sec(2,n_2q),Q_sec(3,n_2q),'g')
    box on
    grid on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    axis equal
    axis([-70 70 -70 70 -50 50])
    set(gca,'Zdir','reverse')
    view([1 1 1])
    hold off

end


ii=1
jj=1
figure(12);
hold on
plot3(ll*ones(jj-ii+1),ii:jj,E_elastic_4{mm}(1,ii:jj,ll,nn),'-*')

hold off
pause(0.1)
view([1 1 1])
box on
grid on

 save('Rod_BiSym_Folding.mat','E_elastic_4','Q_sec_M4','-append')