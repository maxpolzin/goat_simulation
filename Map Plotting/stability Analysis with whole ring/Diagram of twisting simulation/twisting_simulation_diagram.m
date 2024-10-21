% twisting simulation diagram
%% build the whole rod ring to caculate the twisting folding
%try to fold the ring without too much control of the boundary
%% start the simulation
clear
clc
close all
%
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
L_co=1.5;%the length of the coordinates
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
ll=24/1


Q_sec_1=Q_sec_M3{mm}(:,:,ll,nn);
N_e_2=2*(size(Q_sec_1,2)-1);% get the element number from the data
Q_sec_2(1:3,:)=diag([-1,-1,1])*Q_sec_1(1:3,:);
Q_sec_2(4:6,:)=diag([-1,-1,1])*Q_sec_1(4:6,:);
Q_sec_2(7:9,:)=diag([-1,-1,1])*Q_sec_1(7:9,:);
Q_sec_2(10:12,:)=diag([-1,-1,1])*Q_sec_1(10:12,:);
Q_sec=[Q_sec_1,Q_sec_2(:,2:end)];

Duplication_M=[1,1,1]

L_sec_2=L_sec_0*4;
L_e_2=L_sec_2/N_e_2;

n_1q=N_e_2/4+1;%the seq numb of 1/4 point
n_2q=2*N_e_2/4+1;%the seq numb of 2/4 point
n_3q=3*N_e_2/4+1;%the seq numb of 3/4 point

r_node_1=Q_sec(1:3,:);
dr_dx_1=Q_sec(4:6,:);
dr_dy_1=Q_sec(7:9,:);
dr_dz_1=Q_sec(10:12,:);

Duplication_M=[1 1 1;-1,-1,1]


%plot the origin config
figure(7);
hold on
for ii=1:size(Q_sec,2)-1
    q_node = Curve_interp(Q_sec(:,ii),Q_sec(:,ii+1),L_e_2,5);
    Rod_ploting3(q_node,Duplication_M,2*L_co,2)
end

plot3(r_node_1(1,[1,n_2q]),r_node_1(2,[1,n_2q]),r_node_1(3,[1,n_2q]),'-','color','g','linewidth',1)
plot3(r_node_1(1,[n_1q,n_3q]),r_node_1(2,[n_1q,n_3q]),r_node_1(3,[n_1q,n_3q]),'-','color','r','linewidth',1)
scatter3(r_node_1(1,1),r_node_1(2,1),r_node_1(3,1),'r')
scatter3(r_node_1(1,n_2q),r_node_1(2,n_2q),r_node_1(3,n_2q),'g')
l = light;
l.Color = [1 1 1];
l.Position = [1 0 -1];
box off
grid off

xlabel('X')
ylabel('Y')
zlabel('Z')
set(gca,'Zdir','reverse')
set(gca,'Ydir','reverse')
axis equal
% axis off
axis([-70 70 -70 70 -20 30])
view([1 1 0.5])
hold off

%%



figure(8);
C_coord=[.0 .45 .74; .74 .45 .0];
L_co=1.5;



Duplication_M=[1 1 1]

hold on
for ii=1:size(Q_sec,2)-1
    q_node = Curve_interp(Q_sec(:,ii),Q_sec(:,ii+1),L_e_2,5)
    Rod_ploting3(q_node,Duplication_M,L_co,2)
end





plot3(r_node_1(1,[1,n_2q]),r_node_1(2,[1,n_2q]),r_node_1(3,[1,n_2q]),'-','color','b','linewidth',1)
plot3(r_node_1(1,[n_1q,n_3q]),r_node_1(2,[n_1q,n_3q]),r_node_1(3,[n_1q,n_3q]),'-','color','r','linewidth',1)
scatter3(r_node_1(1,1),r_node_1(2,1),r_node_1(3,1),'r')
scatter3(r_node_1(1,n_2q),r_node_1(2,n_2q),r_node_1(3,n_2q),'g')
l = light;
l.Color = [1 1 1];
l.Position = [1 0 -1];
box off
grid off

xlabel('X')
ylabel('Y')
zlabel('Z')
set(gca,'Zdir','reverse')
axis([-70 70 -70 70 -70 70])
axis equal
axis off
view([1 1 1])
hold off

%%
% Build the new boundary conditions for the stability analysis

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
% Con_M_2_0=[1 0 0,   sin(theta_m_0)*cos(phi_m_0)*r_m_0/2, cos(theta_m_0)*sin(phi_m_0)*r_m_0/2, -cos(theta_m_0)*cos(phi_m_0)/2;
%            0 1 0,   sin(theta_m_0)*cos(phi_m_0)*r_m_0/2, cos(theta_m_0)*sin(phi_m_0)*r_m_0/2, -cos(theta_m_0)*cos(phi_m_0)/2;
%            0 0 1,   0  -cos(phi_m_0)*r_m_0/2, -sin(phi_m_0)/2;
%            1 0 0,   -sin(theta_m_0)*cos(phi_m_0)*r_m_0/2, -cos(theta_m_0)*sin(phi_m_0)*r_m_0/2, cos(theta_m_0)*cos(phi_m_0)/2;
%            0 1 0,   cos(theta_m_0)*cos(phi_m_0)*r_m_0/2,  -sin(theta_m_0)*sin(phi_m_0)*r_m_0/2, sin(theta_m_0)*cos(phi_m_0)/2;
%            0 0 1,    0  cos(phi_m_0)*r_m_0/2, sin(phi_m_0)/2;];
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

Proj_M_0((n_1q-1)*12+(1:3),(n_1q-1)*12+(1:3))=Con_M_2_0(1:3,1:3);
Proj_M_0((n_1q-1)*12+[1:3],(n_3q-1)*12+[1:3])=Con_M_2_0(1:3,4:6);

Proj_M_0((n_3q-1)*12+[1:3],(n_1q-1)*12+[1:3])=Con_M_2_0(4:6,1:3);
Proj_M_0((n_3q-1)*12+[1:3],(n_3q-1)*12+[1:3])=Con_M_2_0(4:6,4:6);

%fix R: the distance to the z axis
D_fix(1,n_1q)=0;%fix the xm
D_fix(3,n_3q)=1;%fix the R

Disp_y_0q_0=Q_sec(2,1);
% Disp_x_1q_0=Q_sec(1,n_1q);
%
Disp_y_0q=linspace(0,60,121);
L_co=2;
% Q_sec(:,1)
% Q_sec(:,end);
% figure(11);
%     clf
load('Rod_Free_Folding.mat')
figure(10);
clf
% T_angle=linspace(0,30*pi/180,11);
pic_num_Ydisp=1

%%
for ii=[1:5:107 107]
%      pic_num_Ydisp=ii
    ii;    
    figure(9);
    clf
    if ii~=0
        Q_sec=Q_sec_M4{mm}(:,:,ii,ll,nn);
    end 

    Q_sec(2,1)=Disp_y_0q_0+Disp_y_0q(ii);
    Q_sec(2,end)=Disp_y_0q_0+Disp_y_0q(ii);
    Q_sec(2,n_2q)=-Disp_y_0q_0-Disp_y_0q(ii); 
    L_step=0.02;
    if ismember(ii,[1:21])
        Q_sec(1,n_1q)=Q_sec(1,n_1q)+0.5;
    elseif ismember(ii,[22:41])
        Q_sec(1,n_1q)=Q_sec(1,n_1q)+0.1;
        Q_sec(2,n_1q)=Q_sec(2,n_1q)+0.1;
    elseif ismember(ii,[42:61])
        Q_sec(1,n_1q)=Q_sec(1,n_1q)+0.1;
        Q_sec(2,n_1q)=Q_sec(2,n_1q)+0.1;
    elseif ismember(ii,[62:81])
        Q_sec(1,n_1q)=Q_sec(1,n_1q)+0.1;
        Q_sec(2,n_1q)=Q_sec(2,n_1q)+0.1;
    elseif ismember(ii,[82:93])
        Q_sec(1,n_1q)=Q_sec(1,n_1q)+0.05;
        Q_sec(2,n_1q)=Q_sec(2,n_1q)+0.05;
    elseif ismember(ii,[94:101])
        Q_sec(1,n_1q)=Q_sec(1,n_1q)+0;
        Q_sec(2,n_1q)=Q_sec(2,n_1q)+0;
        L_step=0.01
    elseif ismember(ii,[102:121])
        Q_sec(1,n_1q)=Q_sec(1,n_1q)+0;
        Q_sec(2,n_1q)=Q_sec(2,n_1q)+0;
        L_step=0.01
    end
    kk=0;
    flag1=1;
    flag2=1;
    
   
    while (kk==0||flag1>0.0001*L_step||kk<2/L_step)&&flag2&&kk<1000
        kk=kk+1;
        [dFq_dq_sec,Fq_sec(:,ii),q_node,U_flag(kk),E_elastic_4{mm}(:,ii,ll,nn)] = Jocob_rod_sec(Q_sec,N_e_2,L_e_2,N_node,Par_E,A);
        % update the boundary information
        % update the boundary information
        if r_m_0<10e-2
            r_m=0;
            theta_m=theta_m_0;
            
            D_fix(1:3,n_3q)=1;
            %use the x1=x2 y1=y2 z1=z2 reduced the demension directly
            Con_M_2_0=[ 1 0 0,   0,  0, 0;
                        0 1 0,   0,  0, 0;
                        0 0 1,   0,  0, 0;
                        1 0 0,   0,  0, 0;
                        0 1 0,   0,  0, 0;
                        0 0 1,   0,  0, 0;];
            Proj_M_0((n_1q-1)*12+(1:3),(n_1q-1)*12+(1:3))=Con_M_2_0(1:3,1:3);
            Proj_M_0((n_1q-1)*12+[1:3],(n_3q-1)*12+[1:3])=Con_M_2_0(1:3,4:6);

            Proj_M_0((n_3q-1)*12+[1:3],(n_1q-1)*12+[1:3])=Con_M_2_0(4:6,1:3);
            Proj_M_0((n_3q-1)*12+[1:3],(n_3q-1)*12+[1:3])=Con_M_2_0(4:6,4:6);
        else
            D_fix(1:2,n_3q)=0;
            V_13q=Q_sec(1:3,n_3q)-Q_sec(1:3,n_1q);
            r_xy_m=sqrt(sum(V_13q(1:2).^2));
            r_m=sqrt(sum(V_13q.^2));
            phi_m=atan2(V_13q(3),r_xy_m);
            if r_xy_m<10e-2
                % using the (xm ym zm, x3 y3 z3) repersent x1,x2..., 
                % xm=(x2+x1)/2   ym=(y2+y1)/2  zm=(z2+z1)/2
                % x3=(x2-x1)/2   y3=(y2-y1)/2  z3=(z2-z1)/2
                D_fix(1:2,n_3q)=0;
                D_fix(3,n_3q)=1;
                Con_M_2_0=[1 0 0,   -1,  0, 0;
                           0 1 0,   0,  -1, 0;
                           0 0 1,   0,  0, -1;
                           1 0 0,   1, 0, 0;
                           0 1 0,   0, 1, 0;
                           0 0 1,   0  0, 1;];                
                Proj_M_0((n_1q-1)*12+(1:3),(n_1q-1)*12+(1:3))=Con_M_2_0(1:3,1:3);
                Proj_M_0((n_1q-1)*12+[1:3],(n_3q-1)*12+[1:3])=Con_M_2_0(1:3,4:6);

                Proj_M_0((n_3q-1)*12+[1:3],(n_1q-1)*12+[1:3])=Con_M_2_0(4:6,1:3);
                Proj_M_0((n_3q-1)*12+[1:3],(n_3q-1)*12+[1:3])=Con_M_2_0(4:6,4:6);
            else
                D_fix(1,n_3q)=0;
                D_fix(2,n_3q)=0;
                D_fix(3,n_3q)=1;
                theta_m=atan2(V_13q(2),V_13q(1));
%                 the present First-order expansion relationship between d(theta,phi,R) and d(x,y,z)
%                 Con_M_2_0=[ 1 0 0,   sin(theta_m)*cos(phi_m)*r_m/2, cos(theta_m)*sin(phi_m)*r_m/2, -cos(theta_m)*cos(phi_m)/2;
%                             0 1 0,   sin(theta_m)*cos(phi_m)*r_m/2, cos(theta_m)*sin(phi_m)*r_m/2, -cos(theta_m)*cos(phi_m)/2;
%                             0 0 1,   0  -cos(phi_m)*r_m/2, -sin(phi_m)/2;
%                             1 0 0,   -sin(theta_m)*cos(phi_m)*r_m/2, -cos(theta_m)*sin(phi_m)*r_m/2, cos(theta_m)*cos(phi_m)/2;
%                             0 1 0,   cos(theta_m)*cos(phi_m)*r_m/2,  -sin(theta_m)*sin(phi_m)*r_m/2, sin(theta_m)*cos(phi_m)/2;
%                             0 0 1,    0  cos(phi_m)*r_m/2, sin(phi_m)/2;];
%                 Con_M_2_1=inv(Con_M_2_0);
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
        %canonicalize the diplacement based on the varition of angle
        V_13q=Q_sec(1:3,n_3q)-Q_sec(1:3,n_1q);
        if r_m_0<10e-2            
            Q_sec(1:3,n_1q)=(Q_sec(1:3,n_1q)+Q_sec(1:3,n_3q))/2;
            Q_sec(1:3,n_3q)=Q_sec(1:3,n_1q);
        elseif r_xy_m<10e-2 
            Q_sec(1:3,n_1q)=(Q_sec(1:3,n_1q)+Q_sec(1:3,n_3q))/2-0.5*r_m_0*V_13q/norm(V_13q);
            Q_sec(1:3,n_3q)=(Q_sec(1:3,n_1q)+Q_sec(1:3,n_3q))/2+0.5*r_m_0*V_13q/norm(V_13q);
        else
%             d_theta=Con_M_2_1(4,1:6)*[dQ_sec([1:3],n_1q);dQ_sec([1:3],n_3q)];
%             d_phi  =Con_M_2_1(5,1:6)*[dQ_sec([1:3],n_1q);dQ_sec([1:3],n_3q)];
%             theta_m=theta_m+d_theta;
%             phi_m=phi_m+d_phi;
%             Q_sec(1:3,n_1q)=(Q_sec(1:3,n_1q)+Q_sec(1:3,n_3q))/2-0.5*r_m_0*[cos(theta_m)*cos(phi_m);sin(theta_m)*cos(phi_m);sin(phi_m)];
%             Q_sec(1:3,n_3q)=(Q_sec(1:3,n_1q)+Q_sec(1:3,n_3q))/2+0.5*r_m_0*[cos(theta_m)*cos(phi_m);sin(theta_m)*cos(phi_m);sin(phi_m)];

%             Q_sec(1:3,n_1q)=(Q_sec(1:3,n_1q)+Q_sec(1:3,n_3q))/2-0.5*r_m_0*V_13q/norm(V_13q);
%             Q_sec(1:3,n_3q)=(Q_sec(1:3,n_1q)+Q_sec(1:3,n_3q))/2+0.5*r_m_0*V_13q/norm(V_13q);
              Q_sec(1:3,n_3q)=Q_sec(1:3,n_1q)+r_m_0*V_13q/norm(V_13q);
        end
        Q_sec_M4{mm}(:,:,ii,ll,nn)=Q_sec;
        q_node4{mm}(:,:,:,ii,ll,nn)=q_node;
        flag1=sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e_2+1);
        if kk<=1/L_step
            flag2=1;
        elseif (max(U_flag((kk-1/L_step):kk))-min(U_flag((kk-1/L_step):kk)))<0.1e7%energy variation
            flag2=0;
        end
        X_1(kk)=reshape(Q_sec(1,n_1q),1,[]);
        Y_2(kk)=reshape(Q_sec(2,end),1,[]);
        Z_1(kk)=reshape(Q_sec(3,n_1q),1,[]);
        if kk==1||mod(kk,50)==0
            fprintf('mm=%u,nn=%u,ii=%u, kk=%u, flag1=%f, D13=%f \n',mm,nn,ii,kk,flag1,norm(V_13q))            
            figure(9);
            hold on
            Rod_ploting(q_node,Duplication_M,L_co,C_coord)            
            plot3(Q_sec(1,[1,n_2q]),Q_sec(2,[1,n_2q]),Q_sec(3,[1,n_2q]),'-','color','b','linewidth',1)            
            plot3(Q_sec(1,[n_1q,n_3q]),Q_sec(2,[n_1q,n_3q]),Q_sec(3,[n_1q,n_3q]),'-','color','r','linewidth',1)
            scatter3(Q_sec(1,1),Q_sec(2,1),Q_sec(3,1),'r')
            scatter3(Q_sec(1,n_2q),Q_sec(2,n_2q),Q_sec(3,n_2q),'g')
            hold off
            title(['X_1=',num2str(Q_sec(1,n_1q),'%.0f'), ...
            ', Y_2=', num2str(-Q_sec(2,1),'%.0f'), ...
            ', Z_1=', num2str(Q_sec(3,n_1q),'%.0f')],FontSize=18)
            set(gca,'Zdir','reverse')
            axis equal
            box on
            grid on
            view([1 1 1])
        end
        if mod(kk,10)==1
            figure(5);
            if kk==1
                clf
            end
            subplot(3,1,1)
            plot(1:kk,U_flag(1:kk))
            ylabel('U_{flag}')
            subplot(3,1,2)
            if kk>1/L_step
            plot((kk-1/L_step):kk,U_flag((kk-1/L_step):kk)-U_flag((kk-1/L_step)))
            end
            box on
            grid on
            ylabel('U_{flag}')
            set(gcf, 'Position', [0 50 300 750]);
            drawnow
            subplot(3,1,3)
            hold on
            plot(kk,E_elastic_4{mm}(1,ii,ll,nn),'-*')
            hold off            
            ylabel('E')             
            pause(0.01)
        end

    end


    figure(11);
    clf
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
%     title(['X_1=',num2str(Q_sec(1,n_1q),'%.0f'), ...
%             ', Y_2=', num2str(-Q_sec(2,1),'%.0f'), ...
%             ', Z_1=', num2str(Q_sec(3,n_1q),'%.0f')],FontSize=18)
    title(['T_1=',num2str(norm(V_13q),'%.0f'), ...
            ', T_2=', num2str(-2*Q_sec(2,1),'%.0f')],FontSize=18)
    axis equal
    axis([-70 70 -90 90 -50 70])
    set(gca,'Zdir','reverse')
    view([1 1 1])
    hold off
    set(gcf, 'Position', [305 50 400 360]);
    F = getframe(gcf);  % 获取当前绘图窗口的图片
    Im = frame2im(F);   % 返回与电影帧相关的图片数据
    [Amap, map] = rgb2ind(Im, 256); % 将RGB图片转化为索引图片
    if pic_num_Ydisp == 1
        imwrite(Amap, map,['Free_Folding_Ydis_4','_mm',num2str(mm) '.gif'],'gif', 'Loopcount', Inf, 'DelayTime', 0.1);
        % 将第一张图片写入‘sin.gif’文件中，并且将GIF播放次数设置成无穷，即保存的GIF图会一直动下去
    else
        imwrite(Amap, map,['Free_Folding_Ydis_4','_mm',num2str(mm) '.gif'],'gif','WriteMode','append','DelayTime',0.1);
        % 依次将其他的图片写入到GIF文件当中
    end
    pic_num_Ydisp =pic_num_Ydisp +1
end


% ii=1
% jj=1
% figure(12);
% hold on
% plot3(ll*ones(jj-ii+1),ii:jj,E_elastic_4{mm}(1,ii:jj,ll,nn),'-*')

hold off
pause(0.1)
view([1 1 1])
box on
grid on
save('Rod_Free_Folding_2.mat','E_elastic_4','Q_sec_M4','q_node4')
 save('Rod_Free_Folding_2.mat','E_elastic_4','Q_sec_M4','-append')
 %%
 load('Rod_Free_Folding_2.mat')
 L_co=3
 ii=1
 q_node=q_node4{mm}(:,:,:,ii,ll,nn)
 Q_sec=Q_sec_M4{mm}(:,:,ii,ll,nn)
 figure(11);
%  clf
 hold on
 Rod_ploting3(q_node,Duplication_M,L_co,2)
 plot3(Q_sec(1,[1,n_2q]),Q_sec(2,[1,n_2q]),Q_sec(3,[1,n_2q]),'-','color','g','linewidth',2)
 plot3(Q_sec(1,[n_1q,n_3q]),Q_sec(2,[n_1q,n_3q]),Q_sec(3,[n_1q,n_3q]),'-','color','r','linewidth',2)
 scatter3(Q_sec(1,1),Q_sec(2,1),Q_sec(3,1),'r')
 scatter3(Q_sec(1,n_2q),Q_sec(2,n_2q),Q_sec(3,n_2q),'g')
 box on
 grid on
 xlabel('X')
 ylabel('Y')
 zlabel('Z')
 %     title(['X_1=',num2str(Q_sec(1,n_1q),'%.0f'), ...
 %             ', Y_2=', num2str(-Q_sec(2,1),'%.0f'), ...
 %             ', Z_1=', num2str(Q_sec(3,n_1q),'%.0f')],FontSize=18)
%  title(['T_1=',num2str(norm(V_13q),'%.0f'), ...
%      ', T_2=', num2str(-2*Q_sec(2,1),'%.0f')],FontSize=18)
 axis equal
 axis([-30 50 -90 90 -30 70])
 axis off
 set(gca,'Zdir','reverse')
 set(gca,'Ydir','reverse')
 view([1 -1 1])
 hold off
%  set(gcf, 'Position', [305 50 400 360]);
% l = light;
% l.Color = [1 1 1];
% l.Position = [1 0 -1];





