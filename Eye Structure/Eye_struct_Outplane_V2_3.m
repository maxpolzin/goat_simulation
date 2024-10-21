% the script simulate the eye structure
% author: Qinghua Guan 
% Anyone modify or regenerate based on this code should attached this Authorization statement

%1： caculate the Bisymmetric configuration with different tendon lengths
%with a quarter of the ring
%2： analysis the stability of the Bisymmetric configuration by twist the
%ring circle as antisymmetric about the z-axis
%draw the map of the stability of Bisymmetric configurations
%% start the simulation
clear
clc
close all
%
% input mechinical information
R0_ring=500*1.00;%radius of the ring mm
%the cross angle
Angle_cross=10*pi/180;
% the middle distance 
L_mid=R0_ring*0.3;


L_sec_0=1*R0_ring*pi/4*1.00;% set as the 1/4circle
R_rod=5;%radius of the rod
L_co=2*R_rod;%the length of the coordinates
E=2e4;%elastic modulus 10kPa
G=E/2*(1+0.4);%shearing modulus
Par_E=[E*pi*R_rod^2, E*pi*R_rod^4/4 0.2*G*pi*R_rod^4/2]';% the axis 
% Par_E=1.0e+07 *[ 1.5708; 0.6136; 0.1718]*1e2
% bending and torsion stiffness

% initialization of the caculation
N_e_0=5;% the element number
L_e_0=L_sec_0/N_e_0;
N_node=4;%the sample node number to cacultate the elastic energy of each element
A=1e3;%the penalty parameter to constrain the frame vector Y Z as orthonormal
% it related to the elastic energy, when it too high the iteration would be too slow
% when it is too low, the frame vector Y Z won't be orthonormal anymore.


% caculate the inplane deformation with a 8th of the structure

%define the origin configuration of the two rod of the eye structure


%define the start point and frame of two rod
r1_frame=[0 1 0; -1 0 0; 0 0 1]';
% Tune the start frame
alpha1=0*pi/180;%X axis
beta1=0*pi/180;%Y axis
theta1=0*pi/180;%Z axis
RM_A=[1 0 0;
      0 cos(alpha1) -sin(alpha1);
      0 sin(alpha1)  cos(alpha1)];

RM_B=[cos(beta1)  0  sin(beta1);
      0           1  0;
      -sin(beta1) 0  cos(beta1)];

RM_C=[cos(theta1) -sin(theta1)  0;
      sin(theta1)  cos(theta1)  0
      0            0            1];
r1_frame=RM_C*RM_B*RM_A*r1_frame;% rotaion seq X Y Z 

%generate the start frame for two rods
RM_rd1=[cos(Angle_cross), 0, -sin(Angle_cross);
    0,                1, 0;
    sin(Angle_cross), 0, cos(Angle_cross)];
RM_rd2=[cos(-Angle_cross), 0, -sin(-Angle_cross);
    0,                1, 0;
    sin(-Angle_cross), 0, cos(-Angle_cross)];
r1_frame_rd1=r1_frame*RM_rd1;
r1_frame_rd2=r1_frame*RM_rd2;

r1_P=[R0_ring 0 0]';
q1_rd1=[r1_P;r1_frame_rd1(:)];
q1_rd2=[r1_P;r1_frame_rd2(:)];


%define the end points and frames
r2_frame=[0 1 0; -1 0 0; 0 0 1]';
% Tune the start frame
alpha2=0*pi/180;%X axis
beta2=0*pi/180;%Y axis
theta2=45*pi/180;%Z axis

RM_A=[1 0 0;
      0 cos(alpha2) -sin(alpha2);
      0 sin(alpha2)  cos(alpha2)];

RM_B=[cos(beta2)  0  sin(beta2);
      0           1  0;
      -sin(beta2) 0  cos(beta2)];

RM_C=[cos(theta2) -sin(theta2)  0;
      sin(theta2)  cos(theta2)  0
      0            0            1];
r2_frame=RM_C*RM_A*RM_B*r2_frame;% rotaion seq X Y Z 

%generate the start frame for two rods
r2_frame_rd1=r2_frame;
r2_frame_rd2=r2_frame;
r2_p_rd1=[R0_ring*sqrt(2)/2 R0_ring*sqrt(2)/2, L_mid/2]';
r2_p_rd2=[R0_ring*sqrt(2)/2 R0_ring*sqrt(2)/2,-L_mid/2]';

q2_rd1=[r2_p_rd1;r2_frame_rd1(:)];
q2_rd2=[r2_p_rd2;r2_frame_rd2(:)];

%generate the original solution by interplotion
Q_sec_rd1=Curve_interp(q1_rd1,q2_rd1,L_sec_0,N_e_0);
Q_sec_rd2=Curve_interp(q1_rd2,q2_rd2,L_sec_0,N_e_0);

C_coord=[.0 .45 .74; .74 .45 .0];%the color of the coordinates
Duplication_M=[1 1 1;-1 1 1;-1 -1 1;1 -1 1;];


hold on
for nn=1:(size(Q_sec_rd1,2)-1)
    q_node=Curve_interp(Q_sec_rd1(:,nn),Q_sec_rd1(:,nn+1),L_e_0,N_node);
    Rod_ploting1(q_node,Duplication_M,L_co,C_coord);
    q_node=Curve_interp(Q_sec_rd2(:,nn),Q_sec_rd2(:,nn+1),L_e_0,N_node);
    Rod_ploting1(q_node,Duplication_M,L_co,C_coord);
end
hold off
box on
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
% axis([-600 600 -600 600 -600 600])
view([1 1 1])
% caculate the origin equilibrium state
% Build the boundary conditions for the 8th of structure
%Using the de_demension matrix to define the constrain of the structure
Proj_M_rd1=eye(12*(N_e_0+1),12*(N_e_0+1));%dedemesionalize the variables
% start point
D_fix_rd1=0*Q_sec_rd1;%fixed dimensions/varibales
D_fix_rd1(2,1)=1;%fix the Y disp
D_fix_rd1(3,1)=1;%fix the  Z disp
% fix the local frame
D_fix_rd1(4:12,1)=1;%fix the rotation

%end point
D_fix_rd1(1,end)=0;%X in control
D_fix_rd1(2,end)=1;%Y is constrained with X
D_fix_rd1(3,end)=1;%fix the Z disp
Proj_M_rd1(12*N_e_0+2,12*N_e_0+1)=1;

% Fix local frame
D_fix_rd1([4 5 6],end)=1;%fix the rotation
D_fix_rd1([7 8 9],end)=1;%fix the rotation
D_fix_rd1([10 11 12],end)=1;%fix the rotation

B_index_rd1=find(D_fix_rd1==0);% find the fixed demensions
Proj_M_rd1=Proj_M_rd1(:,B_index_rd1);% extract the projection matrix of free paramters


N_x=0;
% L_co=2;%the length of the coordinates
C_coord=[.0 .45 .74; .74 .45 .0];%the color of the coordinates

Duplication_M=[1 1 1];%-1 1 1;-1 -1 1;1 -1 1;];


Damp_1=0.05;%reduce the iteration speed
mm=1;
kk=0;
flag1=1;
flag3=0;
while (kk==0||flag1>0.000001*R0_ring||flag2>0.10)&&kk<2000
    kk=kk+1;
    [dFq_dq_sec,Fq_sec_M(:,mm),q_node,U_M(mm)] =Jocob_rod_sec(Q_sec_rd1,N_e_0,L_e_0,N_node,Par_E,A);
    dFq_dq_sec_2=Proj_M_rd1'*dFq_dq_sec*Proj_M_rd1;
    dQ_sec_1=-pinv(dFq_dq_sec_2)*Proj_M_rd1'*Fq_sec_M(:,mm);
    dQ_sec=Proj_M_rd1*dQ_sec_1;
    dQ_sec=reshape(dQ_sec,12,[]);
    Q_sec_rd1=Q_sec_rd1+Damp_1*dQ_sec;
    Q_sec_M(:,:,mm)=Q_sec_rd1;
    flag1=sqrt(sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e_0+1))
    flag2=0;%mean(abs(Proj_M_rd1'*Fq_sec_M(:,mm)));%
    % AA=reshape(Fq_sec_M(:,mm),12,[]);
    if kk==2000
        flag3=1;
    end

end





%double the structure to quater
Duplication_RM=[0 1 0;
                1 0 0
                0 0 1];
Q_sec_rd1_1=Q_sec_rd1;
Q_sec_rd1_2=reshape(Duplication_RM*reshape(Q_sec_rd1,3,[]),12,[]);
Q_sec_rd1_2(4:6,:)=-Q_sec_rd1_2(4:6,:);
Q_sec_rd1_2=flip(Q_sec_rd1_2,2);
Q_sec_rd1=[Q_sec_rd1_1(:,1:end),Q_sec_rd1_2(:,2:end)];
Fq_sec_M=[]
[dFq_dq_sec,Fq_sec_M(:,1),q_node,U_M(1)] =Jocob_rod_sec(Q_sec_rd1_1,N_e_0,L_e_0,N_node,Par_E,A);
[dFq_dq_sec,Fq_sec_M(:,2),q_node,U_M(2)] =Jocob_rod_sec(Q_sec_rd1_2,N_e_0,L_e_0,N_node,Par_E,A);
reshape(Fq_sec_M,12,[],2)

Duplication_RM=[1 0 0;
                0 1 0
                0 0 -1]; 
Q_sec_rd2=Q_sec_rd1;
Q_sec_rd2=reshape(Duplication_RM*reshape(Q_sec_rd2,3,[]),12,[]);
Q_sec_rd2(10:12,:)=-Q_sec_rd2(10:12,:);
Fq_sec_M=[]
[dFq_dq_sec,Fq_sec_M(:,1),q_node,U_M(1)] =Jocob_rod_sec(Q_sec_rd1,N_e_0,L_e_0,N_node,Par_E,A);
[dFq_dq_sec,Fq_sec_M(:,2),q_node,U_M(2)] =Jocob_rod_sec(Q_sec_rd2,N_e_0,L_e_0,N_node,Par_E,A);
Fq_sec_M=reshape(Fq_sec_M,12,[],2)
Fq_sec_M(:,N_e_0+1,1)+Fq_sec_M(:,N_e_0+1,2)


Duplication_M=[1 1 1;-1 1 1;-1 -1 1;1 -1 1;];
figure(2);
hold on
Rod_ploting1(Q_sec_rd1,Duplication_M,L_co,C_coord)
Rod_ploting1(Q_sec_rd2,Duplication_M,L_co,C_coord)
Rod_ploting2(Q_sec_rd1,Duplication_M,L_co*3,C_coord)
Rod_ploting2(Q_sec_rd2,Duplication_M,L_co*3,C_coord)
hold off
box on
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
% axis([-600 600 -600 600 -600 600])
view([1 1 1])

N_e_1=N_e_0*2;
L_sec_1=L_sec_0*2;
N_e_1=N_e_0*2;% the element number
L_e_1=L_sec_1/N_e_1;

R_wheel=120;% The radius of the wheel
N_spoke=12;% The spoke num of the wheel
Q_sec=[Q_sec_rd1,Q_sec_rd2]
% save("Eye_Origin.mat",'Q_sec','Q_sec_rd1','Q_sec_rd2')
%% caculate the inplane deformation of bisymmetric configuration
load("Eye_Origin.mat")
load('Eye_Inplane_V2_3.mat')
A=5e3;
% Q_sec=Q_sec_M1(:,:,1);
% Q_sec_rd1=Q_sec(:,1:(N_e_1+1));
% Q_sec_rd2=Q_sec(:,(N_e_1+2):(2*(N_e_1+1)));
% Build the boundary conditions
%Using the de_demension matrix to define the constrain of the structure
Proj_M_0=eye(2*12*(N_e_1+1),2*12*(N_e_1+1));%dedemesionalize the variables

%rod 1
D_fix_rd1=0*Q_sec_rd1;%fixed dimensions/varibales signed as default 0
% start point
% constrain the start point of rd1
%  constarin X Y Z disp
D_fix_rd1(1,1)=1;%control  X disp
D_fix_rd1(2,1)=1;%fix the Y disp
D_fix_rd1(3,1)=1;%fix the  Z disp
%  constarin the local frame
D_fix_rd1(4:12,1)=1;%fix the rotation

%end point
% constrain the end point of rd1
%  constarin X Y Z disp
D_fix_rd1(1,end)=1;%fix the X disp
D_fix_rd1(2,end)=0;%free Y disp
D_fix_rd1(3,end)=1;%fix the Z disp
%  constarin the local frame
D_fix_rd1(4:12,end)=1;%fix the rotation

%mid point
% couple middle point with the middle-middle point of rd1 and rd2
%  couple  X Y Z disp 
%define the projection relationship between mid point of rd1 and mid-mid
%piont: [X_1 Y_1 Z_1]'=[X Y Z]'+RM_local*[0 0 L_mid]'=[X Y Z]'+RM_z*L_mid/2
Con_M_1=eye(3);
Con_M_2=L_mid/2*eye(3);
Proj_M_0(12*(N_e_0)+(1:3),12*(N_e_0)+(1:3))=Con_M_1;
Proj_M_0(12*(N_e_0)+(1:3),12*(N_e_0)+(10:12))=Con_M_2;
%  couple local frame
Proj_M_0(12*N_e_0+(4:12),12*N_e_0+(4:12))=eye(9);

%%rod 2
D_fix_rd2=0*Q_sec_rd1;%fixed dimensions/varibales

% start point
% couple rd2 start point with the start point rd1
%  couple X Y Z disp
Con_M_1=eye(3);
Proj_M_0(12*(N_e_1+1)+(1:3),(1:3))=Con_M_1;
%  couple local frame 
Proj_M_0(12*(N_e_1+1)+(4:12),(4:12))=eye(9);
% constrain the rd2 start point with the start point of rd1
% contrain  X Y Z disp 
D_fix_rd2(1,1)=1;%constrain Y disp
D_fix_rd2(2,1)=1;%constrain Y disp
D_fix_rd2(3,1)=1;%constrain  Z disp
%  constrain the local frame
D_fix_rd2(4:12,1)=1;%fix the rotation

%end point
% couple rd2 end point with the end point of rd1
%  couple X Y Z disp 
Con_M_1=eye(3);
Proj_M_0(12*(2*N_e_1+1)+(1:3),12*N_e_1+(1:3))=Con_M_1;
%  couple local frame 
Proj_M_0(12*(2*N_e_1+1)+(4:12),12*N_e_1+(4:12))=eye(9);
% constrain the rd2 end point the end point of rd1
%  constrain  X Y Z disp
D_fix_rd2(1,end)=1;%constrain X disp
D_fix_rd2(2,end)=1;%constrain Y disp
D_fix_rd2(3,end)=1;%constrain Z disp
%  constrain the local frame
D_fix_rd2(4:12,end)=1;%fix the rotation

%middle point
%couple rd2 middle point with the middle-middle point rd1 and rd2
% couple X Y Z disp 
% [X_2 Y_2 Z_2]'=[X Y Z]'-RM_local*[0 0 L_mid]'=[X Y Z]'-RM_z*L_mid
Con_M_1=eye(3);
Con_M_2=-L_mid/2*eye(3);
Proj_M_0(12*(N_e_1+1+N_e_0)+(1:3),12*(N_e_0)+(1:3))=Con_M_1;
Proj_M_0(12*(N_e_1+1+N_e_0)+(1:3),12*(N_e_0)+(10:12))=Con_M_2;
% couple local frame 
Proj_M_0(12*(N_e_1+1+N_e_0)+(4:12),12*(N_e_0)+(4:12))=eye(9);
% constrain the rd2 middle point the middle-middle point rd1 and rd2
%  constrain  X Y Z disp
D_fix_rd2(1,N_e_0+1)=1;%constrain X disp
D_fix_rd2(2,N_e_0+1)=1;%constrain Y disp
D_fix_rd2(3,N_e_0+1)=1;%constrain Z disp
%contrain local frame
D_fix_rd2(4:12,N_e_0+1)=1;%constrain local frame



D_fix=[D_fix_rd1,D_fix_rd2];
B_index=D_fix==0;% find the fixed demensions
Proj_M=Proj_M_0(:,B_index);% extract the projection matrix of free paramters


Q_sec=[Q_sec_rd1,Q_sec_rd2];
X0=Q_sec_rd1(1,1);
disp_x=0:5:X0;% sampling displacement x
N_x=length(disp_x)-1
% L_co=20;%the length of the coordinates
C_coord=[.0 .45 .74; .74 .45 .0];%the color of the coordinates

Duplication_M=[1 1 1;-1 1 1;-1 -1 1;1 -1 1;];
Damp_1=0.5;%reduce the iteration speed
pic_num=1
P_wheel=[];

%%
for mm=1%:N_x+1
    mm
    if mm~=1
        Q_sec=Q_sec_M1(:,:,mm-1);
        Q_sec_rd1=Q_sec(:,1:(N_e_1+1));
        Q_sec_rd2=Q_sec(:,(N_e_1+2):(2*(N_e_1+1)));
    end  
    Q_sec(1,1)=X0-disp_x(mm);
    Q_sec(1,N_e_1+2)=Q_sec(1,1);
    
    %define the end points and frames
    kk=0;
    flag1=1;
    flag3=0;
    while (kk<=50||flag1>0.000001*R0_ring||flag2>0.10)&&kk<2000
        kk=kk+1;
        [dFq_dq_sec1_rd1,Fq_sec_M1_rd1(:,mm),q_node1_rd1(:,:,:,mm),U_Cal_sec1_rd1(mm),U_v_sec1_rd1(:,mm),U_E_ele1_rd1(mm,:),Dens_E_node1_rd1(:,:,mm)] =Jocob_rod_sec(Q_sec_rd1,N_e_1,L_e_1,N_node,Par_E,A);
        [dFq_dq_sec1_rd2,Fq_sec_M1_rd2(:,mm),q_node1_rd2(:,:,:,mm),U_Cal_sec1_rd2(mm),U_v_sec1_rd2(:,mm),U_E_ele1_rd2(mm,:),Dens_E_node1_rd2(:,:,mm)] =Jocob_rod_sec(Q_sec_rd2,N_e_1,L_e_1,N_node,Par_E,A);
        dFq_dq_sec=[dFq_dq_sec1_rd1,0*dFq_dq_sec1_rd1;0*dFq_dq_sec1_rd1,dFq_dq_sec1_rd2];
        Fq_sec_M1(:,mm)=[Fq_sec_M1_rd1(:,mm);Fq_sec_M1_rd2(:,mm)];
        dFq_dq_sec_2=Proj_M'*dFq_dq_sec*Proj_M;
        Fq_sec_M1_2=Proj_M_0'*Fq_sec_M1(:,mm);
        F_tendon1(:,mm)=Fq_sec_M1_2(1:2,1);
        %         Fq_sec_M1_2(12*(N_e_0+1)+(1:12)-12)
        dQ_sec_1=-pinv(dFq_dq_sec_2)*Proj_M'*Fq_sec_M1(:,mm);
        dQ_sec=Proj_M*dQ_sec_1;
        dQ_sec=reshape(dQ_sec,12,[]);        
        Q_sec=Q_sec+Damp_1*dQ_sec;

        Q_sec_rd1=Q_sec(:,1:(N_e_1+1));
        Q_sec_rd2=Q_sec(:,(N_e_1+2):(2*(N_e_1+1)));
        P_wheel(:,1)=Q_sec_rd1(1:3,N_e_0+1);
        P_wheel(:,2)=Q_sec_rd2(1:3,N_e_0+1);
        P_wheel(:,3)=mean(P_wheel(1:3,1:2),2);
        V_wheel=P_wheel(:,1)-P_wheel(:,2);
        P_wheel(:,1:2)=P_wheel(:,3)+L_mid/2*V_wheel/norm(V_wheel)*[1,-1];
        D_mid=sqrt(sum((P_wheel(:,1)-P_wheel(:,2)).^2,1));
        Q_sec_rd1(1:3,N_e_0+1)=P_wheel(:,1);
        Q_sec_rd2(1:3,N_e_0+1)=P_wheel(:,2);
        Q_sec=[Q_sec_rd1,Q_sec_rd2];
        Q_sec_M1(:,:,mm)=Q_sec;

        flag1=sqrt(sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e_0+1));
        fprintf('kk=%u, flag1=%f, D_mid=%.2f\n',kk,flag1, D_mid)
        flag2=0;%mean(abs(Proj_M_rd1'*Fq_sec_M1(:,mm)));
        if kk==2000
            flag3=1;
        end
%         if mod(kk,20)==1
%             figure(2);
%             clf
%             hold on
%             Rod_ploting1(q_node_rd1,Duplication_M,L_co,C_coord)
%             Rod_ploting1(q_node_rd2,Duplication_M,L_co,C_coord)
%             hold off
%             box on
%             grid on
%             xlabel('X')
%             ylabel('Y')
%             zlabel('Z')
%             axis equal
%             view([1 1 1])
%         end

    end

end

mm=51
figure(3);
subplot(2,2,1)
disp_xy_B=[reshape(Q_sec_M1(1,1,1:mm),1,[]);reshape(Q_sec_M1(2,N_e_1+1,1:mm),1,[])];
plot(disp_xy_B(1,:),disp_xy_B(2,:))
hold on
plot(disp_xy_B(2,:),disp_xy_B(1,:))
axis equal
hold off
box on
grid on
xlabel('X')
ylabel('Y')
subplot(2,2,2)
plot(disp_x(1:mm),U_v_sec1_rd1(1:2,1:mm)+U_v_sec1_rd2(1:2,1:mm))
subplot(2,2,3)
plot(disp_x(1:mm),-F_tendon1(1,1:mm))
subplot(2,2,4)
plot(disp_x(1:mm),-F_tendon1(2,1:mm))

% save('Eye_inplane_V2_3.mat',"Q_sec_M1",'q_node1_rd1','F_tendon1', ...
%     'U_v_sec1_rd1','U_E_ele1_rd1','Dens_E_node1_rd1' , ...
%     'q_node1_rd2','U_v_sec1_rd2','U_E_ele1_rd2', ...
%     'Dens_E_node1_rd2')
%%
% pic_num=1
% for mm=1:51
%     Q_sec=Q_sec_M1(:,:,mm);
%     Q_sec_rd1=Q_sec(:,1:(N_e_1+1));
%     Q_sec_rd2=Q_sec(:,(N_e_1+2):(2*(N_e_1+1)));
%     figure(2);
%     clf
%     hold on
%     Rod_ploting1(q_node1_rd1(:,:,:,mm),Duplication_M,L_co,C_coord,U_E_ele1_rd1(mm,:),Dens_E_node1_rd1(:,:,mm))
%     Rod_ploting1(q_node1_rd2(:,:,:,mm),Duplication_M,L_co,C_coord,U_E_ele1_rd2(mm,:),Dens_E_node1_rd2(:,:,mm))
%     %%plotting wheels
%     P_wheel=[];
%     P_wheel(:,1)=Q_sec_rd1(1:3,N_e_0+1);
%     P_wheel(:,2)=Q_sec_rd2(1:3,N_e_0+1);
%     P_wheel(:,3)=mean(P_wheel(1:3,:),2);
%     P_wheel(:,4)=Q_sec_rd1(1:3,1);
%     P_wheel(:,5)=Q_sec_rd2(1:3,1);
%     P_wheel(:,6)=mean(P_wheel(:,4:5),2);
%     P_wheel(:,7)=Q_sec_rd1(1:3,end);
%     P_wheel(:,8)=Q_sec_rd2(1:3,end);
%     P_wheel(:,9)=mean(P_wheel(:,7:8),2);
%     RM_wheel=reshape(Q_sec_rd1(4:12,N_e_0+1),3,3);
%     Wheel_Connect_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M)
%     plot3(P_wheel(1,6)*[1,-1],P_wheel(2,6)*[1,1],P_wheel(3,6)*[1,1],'r-','linewidth',2)
%     plot3(P_wheel(1,9)*[1,1],P_wheel(2,9)*[1,-1],P_wheel(3,9)*[1,1],'g-','linewidth',2)
%     c=colorbar
%     if mm==1
%         caxis([1.95e1 2.2e1])
%     end
%     c.Label.String = 'Energy Density (mJ/mm)';
%     hold off
%     box on
%     grid on
%     xlabel('X')
%     ylabel('Y')
%     zlabel('Z')
%     axis equal
%     axis([-600 600 -700 700 -100 100])
%     view([1 1 1])
%     title(['X_1=',num2str((X0-disp_x(mm)),'%.0f'),', Y_2=',num2str(Q_sec_M1(2,N_e_1+1,mm),'%.0f')])
%     drawnow;% 立即刷新当前绘图窗口，这是matlab绘图中动态展示的关键
%     F = getframe(gcf);  % 获取当前绘图窗口的图片
%     Im = frame2im(F);   % 返回与电影帧相关的图片数据
%     [Amap, map] = rgb2ind(Im, 256); % 将RGB图片转化为索引图片
%     if pic_num == 1
%         imwrite(Amap, map, 'Eye_inplane_V2_3_2.gif', 'gif', 'Loopcount', Inf, 'DelayTime', 0.1);
%         % 将第一张图片写入‘sin.gif’文件中，并且将GIF播放次数设置成无穷，即保存的GIF图会一直动下去
%     else
%         imwrite(Amap, map,'Eye_inplane_V2_3_3.gif','gif','WriteMode','append','DelayTime',0.1);
%         % 依次将其他的图片写入到GIF文件当中
%     end
%     pic_num = pic_num + 1;
%     pause(0.1)
% end
load('Eye_Outplane_V2_3.mat')
load('Eye_Bisym_outplane_V2_3.mat')
v = VideoWriter('Eye_inplane_2.mp4', 'MPEG-4');
v.FrameRate = 30;  % Set the frame rate (frames per second)
v.Quality = 90;    % Set the video quality (0 to 100)
open(v);
% writeVideo(v, F);
% close(v)
pic_num=1

for mm=1:51%N_x+1
    Q_sec=Q_sec_M1(:,:,mm);
    Q_sec_rd1=Q_sec(:,1:(N_e_1+1));
    Q_sec_rd2=Q_sec(:,(N_e_1+2):(2*(N_e_1+1)));
    figure(2);
    clf
    hold on
    Rod_ploting1(q_node1_rd1(:,:,:,mm),Duplication_M,L_co,C_coord,U_E_ele1_rd1(mm,:),Dens_E_node1_rd1(:,:,mm))
    Rod_ploting1(q_node1_rd2(:,:,:,mm),Duplication_M,L_co,C_coord,U_E_ele1_rd2(mm,:),Dens_E_node1_rd2(:,:,mm))
    %%plotting wheels
    P_wheel=[];
    P_wheel(:,1)=Q_sec_rd1(1:3,N_e_0+1);
    P_wheel(:,2)=Q_sec_rd2(1:3,N_e_0+1);
    P_wheel(:,3)=mean(P_wheel(1:3,:),2);
    P_wheel(:,4)=Q_sec_rd1(1:3,1);
    P_wheel(:,5)=Q_sec_rd2(1:3,1);
    P_wheel(:,6)=mean(P_wheel(:,4:5),2);
    P_wheel(:,7)=Q_sec_rd1(1:3,end);
    P_wheel(:,8)=Q_sec_rd2(1:3,end);
    P_wheel(:,9)=mean(P_wheel(:,7:8),2);
    RM_wheel=reshape(Q_sec_rd1(4:12,N_e_0+1),3,3);
    Wheel_Connect_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M)
    plot3(P_wheel(1,6)*[1,-1],P_wheel(2,6)*[1,1],P_wheel(3,6)*[1,1],'r-','linewidth',2)
    plot3(P_wheel(1,9)*[1,1],P_wheel(2,9)*[1,-1],P_wheel(3,9)*[1,1],'g-','linewidth',2)
    c=colorbar   
    if mm==1
        caxis([1.95e1 3.7e1])
    end
    c.Label.String = 'Energy Density (mJ/mm)';
    hold off
    box on
    grid on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    axis equal
    axis([-500 500 -700 700 -700 200])
    view([1 1 1])
    title(['D_1=',num2str(2*(X0-disp_x(mm)+9.5),'%.0f'),', D_2=',num2str(2*Q_sec_M1(2,N_e_1+1,mm)+19.5,'%.0f')])
    set(gca,'FontSize',16,'fontname','Times')
    drawnow;% 立即刷新当前绘图窗口，这是matlab绘图中动态展示的关键
    F = getframe(gcf);  % 获取当前绘图窗口的图片
    Im = frame2im(F);   % 返回与电影帧相关的图片数据
    [Amap, map] = rgb2ind(Im, 256); % 将RGB图片转化为索引图片
    if pic_num == 1
        imwrite(Amap, map, 'Eye_inplane_2.gif', 'gif', 'Loopcount', Inf, 'DelayTime', 0.1);
        % 将第一张图片写入‘sin.gif’文件中，并且将GIF播放次数设置成无穷，即保存的GIF图会一直动下去
    else
        imwrite(Amap, map,'Eye_inplane_2.gif','gif','WriteMode','append','DelayTime',0.1);
        % 依次将其他的图片写入到GIF文件当中
    end
    pic_num = pic_num + 1;
    writeVideo(v, F);

    pause(0.1)
end
close(v)




%% caculate outplane deformation when contract tendon2
load('Eye_Inplane_V2_3.mat',"Q_sec_M1")
% load('Eye_Biaxis_Outplane_V2.mat',"Q_sec_M2")
A=1e3;
Damp_1=0.1;%reduce the iteration speed
% Build the boundary conditions for the boundary
% start point x_disp fixed and end point y_disp in control
%Using the de_demension matrix to define the constrain of the structure
Proj_M_0=eye(2*12*(N_e_1+1),2*12*(N_e_1+1));%dedemesionalize the variables

%rod 1
D_fix_rd1=0*Q_sec_rd1;
% start point
% couple start point with the end-middle point of rd1 and rd2
%  couple X Y Z disp 
Con_M_1=eye(3);
Proj_M_0(1:3,1:3)=Con_M_1;
%  couple local frame 
Proj_M_0(4:12,4:12)=eye(9);
% constrain the start-middle point of rd1 and rd2
%  constrain  X Y Z disp
D_fix_rd1(1,1)=1;%control  the X disp
D_fix_rd1(2,1)=1;%fix the Y disp
D_fix_rd1(3,1)=0;%free the  Z disp
%  constrain/fix the local frame
D_fix_rd1(4:6,1)=1;%fix the X axis
D_fix_rd1(7:12,1)=0;%free the X ratation of Y Z axis
D_fix_rd1([8 11],1)=1;%free the X ratation of Y Z axis

%end point
% couple end point with the end-middle point of rd1 and rd2
%  couple X Y Z disp
Con_M_1=eye(3);
Proj_M_0(12*N_e_1+(1:3),12*N_e_1+(1:3))=Con_M_1;
%  couple local frame
Proj_M_0(12*N_e_1+(4:12),12*N_e_1+(4:12))=eye(9);
% constrain the end-middle point of rd1 and rd2
%  constrain  X Y Z disp 
D_fix_rd1(1,end)=1;%fix the X disp
D_fix_rd1(2,end)=1;%control the Y disp
D_fix_rd1(3,end)=1;%fix the Z disp
%  constrain/Fix local frame
D_fix_rd1(4:6,end)=1;%fix the X axis
D_fix_rd1(7:12,end)=0;%free the X ratation of Y Z axis
D_fix_rd1([7 10],end)=1;%free the X ratation of Y Z axis

%mid point
% couple rd1 middle point with the middle-middle point rd1 and rd2
%define the projection relationship between piont
%  couple X Y Z disp
% [X_1 Y_1 Z_1]'=[X Y Z]'+RM_local*[0 0 L_mid]'=[X Y Z]'+RM_z*L_mid
Con_M_1=eye(3);
Con_M_2=L_mid/2*eye(3);
Proj_M_0(12*N_e_0+(1:3),12*N_e_0+(1:3))=Con_M_1;
Proj_M_0(12*N_e_0+(1:3),12*N_e_0+(10:12))=Con_M_2;
%  couple local frame  
Proj_M_0(12*N_e_0+(4:12),12*N_e_0+(4:12))=eye(9);
% non-constrain of middle-middle point

%%rod 2
D_fix_rd2=0*Q_sec_rd1;%fixed dimensions/varibales

%start point
% couple rd2 start point with the start-middle point rd1 and rd2
%  couple X Y Z disp 
Con_M_1=eye(3);
Proj_M_0(12*(N_e_1+1)+(1:3),(1:3))=Con_M_1;
%  couple local frame 
Proj_M_0(12*(N_e_1+1)+(4:12),(4:12))=eye(9);
% constrain the rd2 start point with the start-middle point rd1 and rd2 
%  constrain  X Y Z disp
D_fix_rd2(1,1)=1;%control  X disp
D_fix_rd2(2,1)=1;%fix the Y disp
D_fix_rd2(3,1)=1;%fix the  Z disp
%  constrain the local frame
D_fix_rd2(4:12,1)=1;%fix the rotation

%end point
% couple rd2 end point with the end-middle point rd1 and rd2
%  couple X Y Z disp 
Con_M_1=eye(3);
Proj_M_0(12*(2*N_e_1+1)+(1:3),12*N_e_1+(1:3))=Con_M_1;
%  couple local frame 
Proj_M_0(12*(2*N_e_1+1)+(4:12),12*N_e_1+(4:12))=eye(9);
% constrain the rd2 end point the end-middle point rd1 and rd2
%  constrain  X Y Z disp
D_fix_rd2(1,end)=1;%control  X disp
D_fix_rd2(2,end)=1;%fix the Y disp
D_fix_rd2(3,end)=1;%fix the  Z disp
%  constrain the local frame
D_fix_rd2(4:12,end)=1;%fix the rotation

%middle point
%couple rd2 middle point with the middle-middle point rd1 and rd2
% couple X Y Z disp 
Con_M_1=eye(3);
Con_M_2=-L_mid/2*eye(3);
Proj_M_0(12*(N_e_1+1+N_e_0)+(1:3),12*N_e_0+(1:3))=Con_M_1;
Proj_M_0(12*(N_e_1+1+N_e_0)+(1:3),12*N_e_0+(10:12))=Con_M_2;
% couple local frame 
Proj_M_0(12*(N_e_1+1+N_e_0)+(4:12),12*+N_e_0+(4:12))=eye(9);
% constrain the rd2 middle point the middle-middle point rd1 and rd2
%  constrain  X Y Z disp
D_fix_rd2(1,N_e_0+1)=1;%control  X disp
D_fix_rd2(2,N_e_0+1)=1;%fix the Y disp
D_fix_rd2(3,N_e_0+1)=1;%fix the  Z disp
%  constrain the local frame
D_fix_rd2(4:12,N_e_0+1)=1;%fix the rotation

D_fix=[D_fix_rd1,D_fix_rd2];

B_index=D_fix==0;% find the fixed demensions
Proj_M=Proj_M_0(:,B_index);% extract the projection matrix of free paramters


% since the symmetry of the structure, just a quater of the structure need to be caculate
Fq_sec=[];
U=[];
P_wheel=[]
U_M1_rd1=[]
U_M1_rd2=[]
% for mm=50 
mm=50 
    figure(3);
    clf
    Q_sec=Q_sec_M1(:,:,mm);
    Q_sec_rd1=Q_sec(:,1:(N_e_1+1));
    Q_sec_rd2=Q_sec(:,(N_e_1+2):(2*(N_e_1+1)));
    % caculate the range of disp y based on the boundary 
    disp_x_0=Q_sec(1,1);
    disp_y_0=Q_sec(2,end);
    N_y=fix((Q_sec(2,end)-Q_sec(1,1))/1);
    disp_y=flip([disp_x_0:5:disp_y_0, disp_y_0]);% sampling displacement y
    %%
    for nn=1:length(disp_y)
        mm
        nn
        Damp_1=0.0001
        Damp_2=0.01
        if ismember(nn,1:37)
            Damp_1=0.005
            Damp_2=0.5
        elseif ismember(nn,38:42)
            Damp_1=0.002
            Damp_2=0.2
        elseif ismember(nn,43:45)
            Damp_1=0.001
            Damp_2=0.1
        elseif ismember(nn,46:51)
            Damp_1=0.001
            Damp_2=0.05
        elseif ismember(nn,52:62)
            Damp_1=0.001
            Damp_2=0.02
        elseif ismember(nn,63:106)
            Damp_1=0.001
            Damp_2=0.01
%         elseif ismember(nn,39:91)
%             Damp_1=0.00001
%             Damp_2=0.005
%         elseif ismember(nn,92:100)
%             Damp_1=0.0001
%             Damp_2=0.002
%         elseif ismember(nn,19:39)
%             Damp_1=0.0001
%             Damp_2=0.1
%         elseif ismember(nn,40:44)
%             Damp_1=0.0001
%             Damp_2=0.05
%         elseif ismember(nn,45:length(disp_y))
%             Damp_1=0.0001
%             Damp_2=0.025
%         elseif ismember(nn,8:10)
%             Damp=0.1
%         elseif ismember(nn,10:length(disp_y))
%             Damp=0.01
        end
        if nn~=1
            Q_sec=Q_sec_M2{mm}(:,:,nn-1);
        end
%         Q_sec(1,1)=disp_y(nn);
%         Q_sec(1,N_e_1+2)=disp_y(nn);        
        Q_sec(2,N_e_1+1)=disp_y(nn);        
        Q_sec(2,2*(N_e_1+1))=disp_y(nn);
        if ismember(nn,1:20)
            Q_sec(3,1)=Q_sec(3,1)-4;% add the z axis distorion
            Q_sec(3,N_e_1+2)=Q_sec(3,N_e_1+2)-4;% add the z axis distorion
        else
            Q_sec(3,1)=Q_sec(3,1);
            Q_sec(3,N_e_1+2)=Q_sec(3,N_e_1+2);
        end
        Q_sec_rd1=Q_sec(:,1:(N_e_1+1));
        Q_sec_rd2=Q_sec(:,(N_e_1+2):(2*(N_e_1+1)));

        %define the end points and frames
        kk=0;
        flag1=1;
        flag2=1;
        Damp=Damp_1
        while (kk<=100||flag1>0.000001*R0_ring)&&flag2&&(kk<5/Damp||kk<6000)
            kk=kk+1;
            Damp=0.995^(kk)*Damp_1+(1-0.995^(kk))*Damp_2;
            [dFq_dq_sec2_rd1,Fq_sec_M2_rd1(:,mm),q_node2_rd1{mm}(:,:,:,nn),U_Cal_sec2_rd1{mm}(nn),U_v_sec2_rd1{mm}(:,nn),U_E_ele2_rd1{mm}(nn,:),Dens_E_node2_rd1{mm}(:,:,nn)] =Jocob_rod_sec(Q_sec_rd1,N_e_1,L_e_1,N_node,Par_E,A);
            [dFq_dq_sec2_rd2,Fq_sec_M2_rd2(:,mm),q_node2_rd2{mm}(:,:,:,nn),U_Cal_sec2_rd2{mm}(nn),U_v_sec2_rd2{mm}(:,nn),U_E_ele2_rd2{mm}(nn,:),Dens_E_node2_rd2{mm}(:,:,nn)] =Jocob_rod_sec(Q_sec_rd2,N_e_1,L_e_1,N_node,Par_E,A);
            U_Cal_sec2_flag(kk)=U_Cal_sec2_rd1{mm}(nn)+U_Cal_sec2_rd2{mm}(nn);
            dFq_dq_sec=[dFq_dq_sec2_rd1,0*dFq_dq_sec2_rd1;0*dFq_dq_sec2_rd1,dFq_dq_sec2_rd2];
            Fq_sec_M2(:,mm)=[Fq_sec_M2_rd1(:,mm);Fq_sec_M2_rd2(:,mm)];
            dFq_dq_sec_2=Proj_M'*dFq_dq_sec*Proj_M;
            Fq_sec_M1_2=Proj_M_0'*Fq_sec_M2(:,mm);
            F_tendon2{mm}(1:2,nn)=Fq_sec_M1_2(1:2);
            dQ_sec_1=-pinv(dFq_dq_sec_2)*Proj_M'*Fq_sec_M2(:,mm);
            dQ_sec=Proj_M*dQ_sec_1;
            dQ_sec=reshape(dQ_sec,12,[]);
            Q_sec=Q_sec+Damp*dQ_sec;            
            Q_sec_rd1=Q_sec(:,1:(N_e_1+1));
            Q_sec_rd2=Q_sec(:,(N_e_1+2):(2*(N_e_1+1)));
            P_wheel(:,1)=Q_sec_rd1(1:3,N_e_0+1);
            P_wheel(:,2)=Q_sec_rd2(1:3,N_e_0+1);
            P_wheel(:,3)=mean(P_wheel(1:3,1:2),2);
            V_wheel=P_wheel(:,1)-P_wheel(:,2);
            P_wheel(:,1:2)=P_wheel(:,3)+L_mid/2*V_wheel/norm(V_wheel)*[1,-1];
            D_mid=sqrt(sum((P_wheel(:,1)-P_wheel(:,2)).^2,1));
            Q_sec_rd1(1:3,N_e_0+1)=P_wheel(:,1);
            Q_sec_rd2(1:3,N_e_0+1)=P_wheel(:,2);
            Q_sec=[Q_sec_rd1,Q_sec_rd2];
            Q_sec_M2{mm}(:,:,nn)=Q_sec;   

            flag1=sqrt(sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e_0+1));
            flag2=1;
            if kk<300
                flag2=1;
            elseif (max(U_Cal_sec2_flag(kk-100:kk))-min(U_Cal_sec2_flag(kk-100:kk)))<2e1;%energy variation
                flag2=0;
            end
            if kk==2000
                flag3=1;
            end
            if mod(kk,20)==1
                fprintf('kk=%u, flag1=%f, D_mid=%.2f\n',kk,flag1, D_mid)
                figure(4);
                subplot(2,2,1)
                plot(U_Cal_sec2_flag(1:kk))
                subplot(2,2,2)
                if kk>100
                    plot(U_Cal_sec2_flag(kk-100:kk)-U_Cal_sec2_flag(kk-100))
                end
            end

            if mod(kk,100)==2
                figure(3);
                clf
                hold on
                Rod_ploting1(q_node2_rd1{mm}(:,:,:,nn),Duplication_M,L_co,C_coord,U_E_ele2_rd1{mm}(nn,:),Dens_E_node2_rd1{mm}(:,:,nn))
                Rod_ploting1(q_node2_rd2{mm}(:,:,:,nn),Duplication_M,L_co,C_coord,U_E_ele2_rd2{mm}(nn,:),Dens_E_node2_rd2{mm}(:,:,nn))
                %%plotting wheels
                %%plotting wheels
                P_wheel=[];
                P_wheel(:,1)=Q_sec_rd1(1:3,N_e_0+1);
                P_wheel(:,2)=Q_sec_rd2(1:3,N_e_0+1);
                P_wheel(:,3)=mean(P_wheel(1:3,:),2);
                P_wheel(:,4)=Q_sec_rd1(1:3,1);
                P_wheel(:,5)=Q_sec_rd2(1:3,1);
                P_wheel(:,6)=mean(P_wheel(:,4:5),2);
                P_wheel(:,7)=Q_sec_rd1(1:3,end);
                P_wheel(:,8)=Q_sec_rd2(1:3,end);
                P_wheel(:,9)=mean(P_wheel(:,7:8),2);
                RM_wheel=reshape(Q_sec_rd1(4:12,N_e_0+1),3,3);
                Wheel_Connect_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M)
                plot3(P_wheel(1,6)*[1,-1],P_wheel(2,6)*[1,1],P_wheel(3,6)*[1,1],'r-','linewidth',2)
                plot3(P_wheel(1,9)*[1,1],P_wheel(2,9)*[1,-1],P_wheel(3,9)*[1,1],'g-','linewidth',2)
                c=colorbar;
%                 if mm==1
%                     caxis([1.95e1 2.2e1]);
%                 end
                c.Label.String = 'Energy Density (mJ/mm)';
                hold off
                box on
                grid on
                xlabel('X')
                ylabel('Y')
                zlabel('Z')                
                axis([-600 600 -700 700 -700 700])
                axis equal
                view([1 1 1])
                title(['X_1=',num2str(Q_sec(1,1),'%.0f'), ...
                    ', X_2=', num2str(Q_sec(2,N_e_1+1),'%.0f'), ...
                    ', Z_1=', num2str(Q_sec(3,1),'%.0f')])
                figure(4);
                subplot(2,2,1)
                plot(U_Cal_sec2_flag(1:kk))
                subplot(2,2,2)
                if kk>100
                    plot(U_Cal_sec2_flag(kk-100:kk)-U_Cal_sec2_flag(kk-100))
                end
                subplot(2,2,3)
                X_1=reshape(sum(Q_sec_M2{mm}(1,[1,N_e_1+2],1:nn))/2,1,[]);
                Y_2=reshape(sum(Q_sec_M2{mm}(2,[N_e_1+1,2*N_e_1+2],1:nn))/2,1,[]);
                Z_2=reshape(sum(Q_sec_M2{mm}(3,[1,N_e_1+2],1:nn))/2,1,[]);
                yyaxis left
                hline1=plot(Y_2(1)-Y_2(1:nn),-F_tendon2{mm}(1,1:nn),'-o');
                ylim([0 inf])
                yyaxis right
                hline2=plot(Y_2(1)-Y_2(1:nn),F_tendon2{mm}(2,1:nn),'-o');
                ylim([0 inf])
                pause(0.01)
            end
        end
%         save('Eye_outplane_V2_3.mat',"Q_sec_M2",'q_node2_rd1','F_tendon2', ...
%             'U_v_sec2_rd1','U_E_ele2_rd1','Dens_E_node2_rd1' , ...
%             'q_node2_rd2','U_v_sec2_rd2','U_E_ele2_rd2', ...
%             'Dens_E_node2_rd2')
        X_1=reshape(Q_sec_M2{mm}(1,1,1:nn),1,[]);
        Y_2=reshape(Q_sec_M2{mm}(2,N_e_1+1,1:nn),1,[]);
        Z_2=reshape(Q_sec_M2{mm}(3,1,1:nn),1,[]);

        figure(4);
        subplot(2,2,4)
        plot3(X_1,Y_2,Z_2,'-o')
        box on
        grid on
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        view([1 0 0])
    end

% end

mm=50
X_1=reshape(sum(Q_sec_M2{mm}(1,[1,N_e_1+2],:))/2,1,[]);
Y_2=reshape(sum(Q_sec_M2{mm}(2,[N_e_1+1,2*N_e_1+2],:))/2,1,[]);
Z_2=reshape(sum(Q_sec_M2{mm}(3,[1,N_e_1+2],:))/2,1,[]);
U_M=U_Cal_sec2_rd1{mm}+U_Cal_sec2_rd2{mm};
figure(5);
plotyy(Y_2(1)-Y_2,U_M,Y_2(1)-Y_2,Z_2)
box on
grid on
xlabel('X')
ylabel('Y')
%%
load('Eye_Outplane_V2_3.mat')
load('Eye_Bisym_outplane_V2_3.mat')
v = VideoWriter('Eye_Outplane_2_2.mp4', 'MPEG-4');
v.FrameRate = 30;  % Set the frame rate (frames per second)
v.Quality = 90;    % Set the video quality (0 to 100)
open(v);
% writeVideo(v, F);
% close(v)
pic_num=1 
for nn=[1:length(disp_y) length(disp_y)]
      Q_sec=Q_sec_M2{mm}(:,:,nn);
            Q_sec_rd1=Q_sec(:,1:(N_e_1+1));
            Q_sec_rd2=Q_sec(:,(N_e_1+2):(2*(N_e_1+1)));
            figure(3);
            clf
            hold on
            Rod_ploting1(q_node2_rd1{mm}(:,:,:,nn),Duplication_M,L_co,C_coord,U_E_ele2_rd1{mm}(nn,:),Dens_E_node2_rd1{mm}(:,:,nn))
            Rod_ploting1(q_node2_rd2{mm}(:,:,:,nn),Duplication_M,L_co,C_coord,U_E_ele2_rd2{mm}(nn,:),Dens_E_node2_rd2{mm}(:,:,nn))
            
            %%plotting wheels
            P_wheel=[];
            P_wheel(:,1)=Q_sec_rd1(1:3,N_e_0+1);
            P_wheel(:,2)=Q_sec_rd2(1:3,N_e_0+1);
            P_wheel(:,3)=mean(P_wheel(1:3,:),2);
            P_wheel(:,4)=Q_sec_rd1(1:3,1);
            P_wheel(:,5)=Q_sec_rd2(1:3,1);
            P_wheel(:,6)=mean(P_wheel(:,4:5),2);
            P_wheel(:,7)=Q_sec_rd1(1:3,end);
            P_wheel(:,8)=Q_sec_rd2(1:3,end);
            P_wheel(:,9)=mean(P_wheel(:,7:8),2);
            RM_wheel=reshape(Q_sec_rd1(4:12,N_e_0+1),3,3);
            Wheel_Connect_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M)
            plot3(P_wheel(1,6)*[1,-1],P_wheel(2,6)*[1,1],P_wheel(3,6)*[1,1],'r-','linewidth',2)
            plot3(P_wheel(1,9)*[1,1],P_wheel(2,9)*[1,-1],P_wheel(3,9)*[1,1],'g-','linewidth',2)
            c=colorbar;
%             if mm==1
%                 caxis([1.95e1 2.2e1]);
%             end
            c.Label.String = 'Energy Density (mJ/mm)';
            hold off
            box on
            grid on
            h=gca;
            %                 set(h, 'ZDir', 'reverse');
            axis equal
            axis([-500 500 -700 700 -700 200])
            view([1 1 1])
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            title(['D_1=',num2str(2*P_wheel(1,6)+9.5,'%.0f'), ...
                ', D_2=', num2str(2*P_wheel(2,9)+9.5,'%.0f'), ...
                ', Z_1=', num2str(P_wheel(3,6),'%.0f')])
            set(gca,'FontSize',16,'fontname','Times')
            drawnow;% 立即刷新当前绘图窗口，这是matlab绘图中动态展示的关键
            F = getframe(gcf);  % 获取当前绘图窗口的图片
            Im = frame2im(F);   % 返回与电影帧相关的图片数据
            [Amap, map] = rgb2ind(Im, 256); % 将RGB图片转化为索引图片
            if pic_num == 1
                imwrite(Amap, map, 'Eye_Outplane_2_2.gif', 'gif', 'Loopcount', Inf, 'DelayTime', 0.1);
                % 将第一张图片写入‘sin.gif’文件中，并且将GIF播放次数设置成无穷，即保存的GIF图会一直动下去
            else
                imwrite(Amap, map,'Eye_Outplane_2_2.gif','gif','WriteMode','append','DelayTime',0.1);
                % 依次将其他的图片写入到GIF文件当中
            end
            pic_num = pic_num + 1;
            writeVideo(v, F);
            pause(0.1)

end
close(v)
%% caculate bisymmetric configurations in the range inside the boundary
load('Eye_Outplane_V2_3.mat')
% load('Eye_Biaxis_Outplane_V2.mat',"Q_sec_M2")
A=0.7e3;
Damp_1=0.1;%reduce the iteration speed
% Build the boundary conditions for the boundary
% start point x_disp fixed and end point y_disp in control
%Using the de_demension matrix to define the constrain of the structure
Proj_M_0=eye(2*12*(N_e_1+1),2*12*(N_e_1+1));%dedemesionalize the variables

%rod 1
D_fix_rd1=0*Q_sec_rd1;
% start point
% couple start point with the end-middle point of rd1 and rd2
%  couple X Y Z disp 
Con_M_1=eye(3);
Proj_M_0(1:3,1:3)=Con_M_1;
%  couple local frame 
Proj_M_0(4:12,4:12)=eye(9);
% constrain the start-middle point of rd1 and rd2
%  constrain  X Y Z disp
D_fix_rd1(1,1)=1;%control  the X disp
D_fix_rd1(2,1)=1;%fix the Y disp
D_fix_rd1(3,1)=0;%free the  Z disp
%  constrain/fix the local frame
D_fix_rd1(4:6,1)=1;%fix the X axis
D_fix_rd1(7:12,1)=0;%free the X ratation of Y Z axis
D_fix_rd1([8 11],1)=1;%free the X ratation of Y Z axis

%end point
% couple end point with the end-middle point of rd1 and rd2
%  couple X Y Z disp
Con_M_1=eye(3);
Proj_M_0(12*N_e_1+(1:3),12*N_e_1+(1:3))=Con_M_1;
%  couple local frame
Proj_M_0(12*N_e_1+(4:12),12*N_e_1+(4:12))=eye(9);
% constrain the end-middle point of rd1 and rd2
%  constrain  X Y Z disp 
D_fix_rd1(1,end)=1;%fix the X disp
D_fix_rd1(2,end)=1;%control the Y disp
D_fix_rd1(3,end)=1;%fix the Z disp
%  constrain/Fix local frame
D_fix_rd1(4:6,end)=1;%fix the X axis
D_fix_rd1(7:12,end)=0;%free the X ratation of Y Z axis
D_fix_rd1([7 10],end)=1;%free the X ratation of Y Z axis

%mid point
% couple rd1 middle point with the middle-middle point rd1 and rd2
%define the projection relationship between piont
%  couple X Y Z disp
% [X_1 Y_1 Z_1]'=[X Y Z]'+RM_local*[0 0 L_mid]'=[X Y Z]'+RM_z*L_mid
Con_M_1=eye(3);
Con_M_2=L_mid/2*eye(3);
Proj_M_0(12*N_e_0+(1:3),12*N_e_0+(1:3))=Con_M_1;
Proj_M_0(12*N_e_0+(1:3),12*N_e_0+(10:12))=Con_M_2;
%  couple local frame  
Proj_M_0(12*N_e_0+(4:12),12*N_e_0+(4:12))=eye(9);
% non-constrain of middle-middle point

%%rod 2
D_fix_rd2=0*Q_sec_rd1;%fixed dimensions/varibales

%start point
% couple rd2 start point with the start-middle point rd1 and rd2
%  couple X Y Z disp 
Con_M_1=eye(3);
Proj_M_0(12*(N_e_1+1)+(1:3),(1:3))=Con_M_1;
%  couple local frame 
Proj_M_0(12*(N_e_1+1)+(4:12),(4:12))=eye(9);
% constrain the rd2 start point with the start-middle point rd1 and rd2 
%  constrain  X Y Z disp
D_fix_rd2(1,1)=1;%control  X disp
D_fix_rd2(2,1)=1;%fix the Y disp
D_fix_rd2(3,1)=1;%fix the  Z disp
%  constrain the local frame
D_fix_rd2(4:12,1)=1;%fix the rotation

%end point
% couple rd2 end point with the end-middle point rd1 and rd2
%  couple X Y Z disp 
Con_M_1=eye(3);
Proj_M_0(12*(2*N_e_1+1)+(1:3),12*N_e_1+(1:3))=Con_M_1;
%  couple local frame 
Proj_M_0(12*(2*N_e_1+1)+(4:12),12*N_e_1+(4:12))=eye(9);
% constrain the rd2 end point the end-middle point rd1 and rd2
%  constrain  X Y Z disp
D_fix_rd2(1,end)=1;%control  X disp
D_fix_rd2(2,end)=1;%fix the Y disp
D_fix_rd2(3,end)=1;%fix the  Z disp
%  constrain the local frame
D_fix_rd2(4:12,end)=1;%fix the rotation

%middle point
%couple rd2 middle point with the middle-middle point rd1 and rd2
% couple X Y Z disp 
Con_M_1=eye(3);
Con_M_2=-L_mid/2*eye(3);
Proj_M_0(12*(N_e_1+1+N_e_0)+(1:3),12*N_e_0+(1:3))=Con_M_1;
Proj_M_0(12*(N_e_1+1+N_e_0)+(1:3),12*N_e_0+(10:12))=Con_M_2;
% couple local frame 
Proj_M_0(12*(N_e_1+1+N_e_0)+(4:12),12*+N_e_0+(4:12))=eye(9);
% constrain the rd2 middle point the middle-middle point rd1 and rd2
%  constrain  X Y Z disp
D_fix_rd2(1,N_e_0+1)=1;%control  X disp
D_fix_rd2(2,N_e_0+1)=1;%fix the Y disp
D_fix_rd2(3,N_e_0+1)=1;%fix the  Z disp
%  constrain the local frame
D_fix_rd2(4:12,N_e_0+1)=1;%fix the rotation

D_fix=[D_fix_rd1,D_fix_rd2];

B_index=D_fix==0;% find the fixed demensions
Proj_M=Proj_M_0(:,B_index);% extract the projection matrix of free paramters


% since the symmetry of the structure, just a quater of the structure need to be caculate

P_wheel=[]

for mm=50 
    figure(3);
    clf
    Q_sec=Q_sec_M2{50}(:,:,end);
    Q_sec_rd1=Q_sec(:,1:(N_e_1+1));
    Q_sec_rd2=Q_sec(:,(N_e_1+2):(2*(N_e_1+1)));
    % caculate the range of disp y based on the boundary 
%     disp_x_0=Q_sec(1,1);
    disp_y_0=Q_sec(2,end);
    N_y=fix((Q_sec(2,end)-Q_sec(1,1))/1);
    disp_y=flip([0:5:disp_y_0, disp_y_0]);% sampling displacement y
    for nn=51:length(disp_y)
        mm
        nn
        Damp_1=0.0001
        Damp_2=0.01
        if ismember(nn,1:37)
            Damp_1=0.002
            Damp_2=0.2
        elseif ismember(nn,38:40)
            Damp_1=0.001
            Damp_2=0.1
        elseif ismember(nn,41:49)
            Damp_1=0.001
            Damp_2=0.05
        elseif ismember(nn,50:60)
            Damp_1=0.001
            Damp_2=0.02
%         elseif ismember(nn,46:51)
%             Damp_1=0.001
%             Damp_2=0.05
%         elseif ismember(nn,52:62)
%             Damp_1=0.001
%             Damp_2=0.02
%         elseif ismember(nn,63:106)
%             Damp_1=0.001
%             Damp_2=0.01
        end
        if nn~=1
            Q_sec=Q_sec_M3{mm}(:,:,nn-1);
        end
        Q_sec(1,1)=disp_y(nn);
        Q_sec(1,N_e_1+2)=disp_y(nn);        
        Q_sec(2,N_e_1+1)=disp_y(nn);        
        Q_sec(2,2*(N_e_1+1))=disp_y(nn);

        Q_sec_rd1=Q_sec(:,1:(N_e_1+1));
        Q_sec_rd2=Q_sec(:,(N_e_1+2):(2*(N_e_1+1)));

        %define the end points and frames
        kk=0;
        flag1=1;
        flag2=1;
        Damp=Damp_1
        while (kk<=100||flag1>0.000001*R0_ring)&&flag2&&(kk<5/Damp||kk<6000)
            kk=kk+1;
            Damp=0.995^(kk)*Damp_1+(1-0.995^(kk))*Damp_2;
            [dFq_dq_sec3_rd1,Fq_sec_M3_rd1(:,mm),q_node3_rd1{mm}(:,:,:,nn),U_Cal_sec3_rd1{mm}(nn),U_v_sec3_rd1{mm}(:,nn),U_E_ele3_rd1{mm}(nn,:),Dens_E_node3_rd1{mm}(:,:,nn)] =Jocob_rod_sec(Q_sec_rd1,N_e_1,L_e_1,N_node,Par_E,A);
            [dFq_dq_sec3_rd2,Fq_sec_M3_rd2(:,mm),q_node3_rd2{mm}(:,:,:,nn),U_Cal_sec3_rd2{mm}(nn),U_v_sec3_rd2{mm}(:,nn),U_E_ele3_rd2{mm}(nn,:),Dens_E_node3_rd2{mm}(:,:,nn)] =Jocob_rod_sec(Q_sec_rd2,N_e_1,L_e_1,N_node,Par_E,A);
            U_Cal_sec2_flag(kk)=U_Cal_sec3_rd1{mm}(nn)+U_Cal_sec3_rd2{mm}(nn);
            dFq_dq_sec=[dFq_dq_sec3_rd1,0*dFq_dq_sec3_rd1;0*dFq_dq_sec3_rd1,dFq_dq_sec3_rd2];
            Fq_sec_M2(:,mm)=[Fq_sec_M3_rd1(:,mm);Fq_sec_M3_rd2(:,mm)];
            dFq_dq_sec_2=Proj_M'*dFq_dq_sec*Proj_M;
            Fq_sec_M1_2=Proj_M_0'*Fq_sec_M2(:,mm);
            F_tendon3{mm}(1:2,nn)=Fq_sec_M1_2(1:2);
            dQ_sec_1=-pinv(dFq_dq_sec_2)*Proj_M'*Fq_sec_M2(:,mm);
            dQ_sec=Proj_M*dQ_sec_1;
            dQ_sec=reshape(dQ_sec,12,[]);
            Q_sec=Q_sec+Damp*dQ_sec;            
            Q_sec_rd1=Q_sec(:,1:(N_e_1+1));
            Q_sec_rd2=Q_sec(:,(N_e_1+2):(2*(N_e_1+1)));
            P_wheel(:,1)=Q_sec_rd1(1:3,N_e_0+1);
            P_wheel(:,2)=Q_sec_rd2(1:3,N_e_0+1);
            P_wheel(:,3)=mean(P_wheel(1:3,1:2),2);
            V_wheel=P_wheel(:,1)-P_wheel(:,2);
            P_wheel(:,1:2)=P_wheel(:,3)+L_mid/2*V_wheel/norm(V_wheel)*[1,-1];
            D_mid=sqrt(sum((P_wheel(:,1)-P_wheel(:,2)).^2,1));
            Q_sec_rd1(1:3,N_e_0+1)=P_wheel(:,1);
            Q_sec_rd2(1:3,N_e_0+1)=P_wheel(:,2);
            Q_sec=[Q_sec_rd1,Q_sec_rd2];
            Q_sec_M3{mm}(:,:,nn)=Q_sec;   

            flag1=sqrt(sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e_0+1));
            flag2=1;
            if kk<=10/Damp_2
                flag2=1;
            elseif (max(U_Cal_sec2_flag((kk-10/Damp_2):kk))-min(U_Cal_sec2_flag((kk-10/Damp_2):kk)))<2e1;%energy variation
                flag2=0;
            end
            if kk==2000
                flag3=1;
            end
            if mod(kk,20)==1
                fprintf('kk=%u, flag1=%f, D_mid=%.2f\n',kk,flag1, D_mid)
                figure(4);
                subplot(2,2,1)
                plot(U_Cal_sec2_flag(1:kk))
                subplot(2,2,2)
                if kk>10/Damp_2
                    plot(U_Cal_sec2_flag((kk-10/Damp_2):kk)-U_Cal_sec2_flag(kk-10/Damp_2))
                end
            end

            if mod(kk,100)==2
                figure(3);
                clf
                hold on
                Rod_ploting1(q_node3_rd1{mm}(:,:,:,nn),Duplication_M,L_co,C_coord,U_E_ele3_rd1{mm}(nn,:),Dens_E_node3_rd1{mm}(:,:,nn))
                Rod_ploting1(q_node3_rd2{mm}(:,:,:,nn),Duplication_M,L_co,C_coord,U_E_ele3_rd2{mm}(nn,:),Dens_E_node3_rd2{mm}(:,:,nn))
                %%plotting wheels
                %%plotting wheels
                P_wheel=[];
                P_wheel(:,1)=Q_sec_rd1(1:3,N_e_0+1);
                P_wheel(:,2)=Q_sec_rd2(1:3,N_e_0+1);
                P_wheel(:,3)=mean(P_wheel(1:3,:),2);
                P_wheel(:,4)=Q_sec_rd1(1:3,1);
                P_wheel(:,5)=Q_sec_rd2(1:3,1);
                P_wheel(:,6)=mean(P_wheel(:,4:5),2);
                P_wheel(:,7)=Q_sec_rd1(1:3,end);
                P_wheel(:,8)=Q_sec_rd2(1:3,end);
                P_wheel(:,9)=mean(P_wheel(:,7:8),2);
                RM_wheel=reshape(Q_sec_rd1(4:12,N_e_0+1),3,3);
                Wheel_Connect_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M)
                plot3(P_wheel(1,6)*[1,-1],P_wheel(2,6)*[1,1],P_wheel(3,6)*[1,1],'r-','linewidth',2)
                plot3(P_wheel(1,9)*[1,1],P_wheel(2,9)*[1,-1],P_wheel(3,9)*[1,1],'g-','linewidth',2)
                c=colorbar;
%                 if mm==1
%                     caxis([1.95e1 2.2e1]);
%                 end
                c.Label.String = 'Energy Density (mJ/mm)';
                hold off
                box on
                grid on
                xlabel('X')
                ylabel('Y')
                zlabel('Z')                
                axis([-600 600 -700 700 -700 700])
                axis equal
                view([1 1 1])
                title(['X_1=',num2str(Q_sec(1,1),'%.0f'), ...
                    ', X_2=', num2str(Q_sec(2,N_e_1+1),'%.0f'), ...
                    ', Z_1=', num2str(Q_sec(3,1),'%.0f')])
                figure(4);
                subplot(2,2,1)
                plot(U_Cal_sec2_flag(1:kk))
                subplot(2,2,2)
                if kk>10/Damp_2
                    plot(U_Cal_sec2_flag(kk-10/Damp_2:kk)-U_Cal_sec2_flag(kk-10/Damp_2))
                end
                subplot(2,2,3)
                X_1=reshape(sum(Q_sec_M3{mm}(1,[1,N_e_1+2],1:nn))/2,1,[]);
                Y_2=reshape(sum(Q_sec_M3{mm}(2,[N_e_1+1,2*N_e_1+2],1:nn))/2,1,[]);
                Z_2=reshape(sum(Q_sec_M3{mm}(3,[1,N_e_1+2],1:nn))/2,1,[]);                
                hline1=plot(Y_2(1)-Y_2(1:nn),-F_tendon3{mm}(1,1:nn),'-o');
                hold on
                hline2=plot(Y_2(1)-Y_2(1:nn),F_tendon3{mm}(2,1:nn),'-o');
                hold off
                ylim([150 inf])
                pause(0.01)
            end
        end
        save('Eye_Bisym_outplane_V2_3.mat',"Q_sec_M3",'q_node3_rd1','F_tendon3', ...
            'U_v_sec3_rd1','U_E_ele3_rd1','Dens_E_node3_rd1' , ...
            'q_node3_rd2','U_v_sec3_rd2','U_E_ele3_rd2', ...
            'Dens_E_node3_rd2')
        X_1=reshape(Q_sec_M3{mm}(1,1,1:nn),1,[]);
        Y_2=reshape(Q_sec_M3{mm}(2,N_e_1+1,1:nn),1,[]);
        Z_2=reshape(Q_sec_M3{mm}(3,1,1:nn),1,[]);

        figure(4);
        subplot(2,2,4)
        plot3(X_1,Y_2,Z_2,'-o')
        box on
        grid on
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        view([1 0 0])
    end

end

mm=50
X_1=reshape(sum(Q_sec_M3{mm}(1,[1,N_e_1+2],:))/2,1,[]);
Y_2=reshape(sum(Q_sec_M3{mm}(2,[N_e_1+1,2*N_e_1+2],:))/2,1,[]);
Z_2=reshape(sum(Q_sec_M3{mm}(3,[1,N_e_1+2],:))/2,1,[]);
U_M=U_Cal_sec3_rd1{mm}+U_Cal_sec3_rd2{mm};
figure(5);
plotyy(Y_2(1)-Y_2,U_M,Y_2(1)-Y_2,Z_2)
box on
grid on
xlabel('X')
ylabel('Y')

%%

load('Eye_Bisym_outplane_V2_3.mat')
v = VideoWriter('Eye_Bisym_Outplane_2.mp4', 'MPEG-4');
v.FrameRate = 30;  % Set the frame rate (frames per second)
v.Quality = 90;    % Set the video quality (0 to 100)
open(v);
pic_num=1
pic_num=1 
for nn=1:length(Q_sec_M3{mm}(1,1,:))
            Q_sec=Q_sec_M3{mm}(:,:,nn);
            Q_sec_rd1=Q_sec(:,1:(N_e_1+1));
            Q_sec_rd2=Q_sec(:,(N_e_1+2):(2*(N_e_1+1)));
            figure(3);
            clf
            hold on
            Rod_ploting1(q_node3_rd1{mm}(:,:,:,nn),Duplication_M,L_co,C_coord,U_E_ele3_rd1{mm}(nn,:),Dens_E_node3_rd1{mm}(:,:,nn))
            Rod_ploting1(q_node3_rd2{mm}(:,:,:,nn),Duplication_M,L_co,C_coord,U_E_ele3_rd2{mm}(nn,:),Dens_E_node3_rd2{mm}(:,:,nn))
            %%plotting wheels
            %%plotting wheels
            P_wheel=[];
            P_wheel(:,1)=Q_sec_rd1(1:3,N_e_0+1);
            P_wheel(:,2)=Q_sec_rd2(1:3,N_e_0+1);
            P_wheel(:,3)=mean(P_wheel(1:3,:),2);
            P_wheel(:,4)=Q_sec_rd1(1:3,1);
            P_wheel(:,5)=Q_sec_rd2(1:3,1);
            P_wheel(:,6)=mean(P_wheel(:,4:5),2);
            P_wheel(:,7)=Q_sec_rd1(1:3,end);
            P_wheel(:,8)=Q_sec_rd2(1:3,end);
            P_wheel(:,9)=mean(P_wheel(:,7:8),2);
            RM_wheel=reshape(Q_sec_rd1(4:12,N_e_0+1),3,3);
            Wheel_Connect_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M)
            plot3(P_wheel(1,6)*[1,-1],P_wheel(2,6)*[1,1],P_wheel(3,6)*[1,1],'r-','linewidth',2)
            plot3(P_wheel(1,9)*[1,1],P_wheel(2,9)*[1,-1],P_wheel(3,9)*[1,1],'g-','linewidth',2)
            c=colorbar;
            %                 if mm==1
            %                     caxis([1.95e1 2.2e1]);
            %                 end
            c.Label.String = 'Energy Density (mJ/mm)';
            hold off
            box on
            grid on
            h=gca;
            %                 set(h, 'ZDir', 'reverse');
            axis equal
            axis([-500 500 -700 700 -700 200])
            view([1 1 1])
            xlabel('X')
            ylabel('Y')
            zlabel('Z')

            title(['X_1=',num2str(2*P_wheel(1,6)-0.2,'%.0f'), ...
                ', Y_2=', num2str(2*P_wheel(2,9)-0.2,'%.0f'), ...
                ', Z_1=', num2str(P_wheel(3,6),'%.0f')])
            set(gca,'FontSize',16,'fontname','Times')
            drawnow;% 立即刷新当前绘图窗口，这是matlab绘图中动态展示的关键
            F = getframe(gcf);  % 获取当前绘图窗口的图片
            Im = frame2im(F);   % 返回与电影帧相关的图片数据
            [Amap, map] = rgb2ind(Im, 256); % 将RGB图片转化为索引图片
            if pic_num == 1
                imwrite(Amap, map, 'Eye_Bisym_Outplane_2.gif', 'gif', 'Loopcount', Inf, 'DelayTime', 0.1);
                % 将第一张图片写入‘sin.gif’文件中，并且将GIF播放次数设置成无穷，即保存的GIF图会一直动下去
            else
                imwrite(Amap, map,'Eye_Bisym_Outplane_2.gif','gif','WriteMode','append','DelayTime',0.1);
                % 依次将其他的图片写入到GIF文件当中
            end
            writeVideo(v, F);
            pic_num = pic_num + 1;
            pause(0.1)

end
close(v);

%%
clear
clc
close all
load('Eye_Inplane_V2_3.mat')
load('Eye_Outplane_V2_3.mat')
load('Eye_Bisym_outplane_V2_3.mat')



% input mechinical information
R0_ring=500;%radius of the ring mm
%the cross angle
Angle_cross=0*pi/180;
% the middle distance 
L_mid=R0_ring*0.3;


L_sec_0=1*R0_ring*pi/2;% set as the 1/4circle
R_rod=5;%radius of the rod
L_co=2*R_rod;%the length of the coordinates

% initialization of the caculation
N_e_0=10;% the element number
L_e_0=L_sec_0/N_e_0;
N_node=4;%the sample node number to cacultate the elastic energy of each element

N_e_1=N_e_0;
L_sec_1=L_sec_0;
N_e_1=N_e_0;% the element number
L_e_1=L_sec_1/N_e_1;
N_e_0=N_e_1/2
L_sec_0=L_sec_1/2;
Duplication_M=[1 1 1;-1 1 1;-1 -1 1;1 -1 1;];
C_coord=[.0 .45 .74; .74 .45 .0];%the color of the coordinates
R_wheel=120;% The radius of the wheel
N_spoke=12;% The spoke num of the wheel





mm=1
Q_sec=Q_sec_M1(:,:,mm)
Q_sec_rd1=Q_sec(:,1:(N_e_1+1));
Q_sec_rd2=Q_sec(:,(N_e_1+2):(2*(N_e_1+1)));
X0=Q_sec_rd1(1,1);
disp_x=[0:5:X0,X0];% sampling displacement x
figure(3);
clf
hold on
Rod_ploting1(q_node1_rd1(:,:,:,mm),Duplication_M,L_co,C_coord,U_E_ele1_rd1(mm,:),Dens_E_node1_rd1(:,:,mm));
Rod_ploting1(q_node1_rd2(:,:,:,mm),Duplication_M,L_co,C_coord,U_E_ele1_rd2(mm,:),Dens_E_node1_rd2(:,:,mm));
%%plotting wheels
P_wheel(:,1)=Q_sec_rd1(1:3,N_e_0+1);
P_wheel(:,2)=Q_sec_rd2(1:3,N_e_0+1);
P_wheel(:,3)=mean(P_wheel(:,1:2),2);
P_wheel(:,4)=Q_sec_rd1(1:3,1);
P_wheel(:,5)=Q_sec_rd2(1:3,1);
P_wheel(:,6)=mean(P_wheel(:,4:5),2);
P_wheel(:,7)=Q_sec_rd1(1:3,end);
P_wheel(:,8)=Q_sec_rd2(1:3,end);
P_wheel(:,9)=mean(P_wheel(:,7:8),2);
RM_wheel=reshape(Q_sec(4:12,N_e_0+1),3,3);
Wheel_Connect_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M,1)
plot3(P_wheel(1,6)*[1,-1],P_wheel(2,6)*[1,1],P_wheel(3,6)*[1,1],'r-','linewidth',2)
plot3(P_wheel(1,9)*[1,1],P_wheel(2,9)*[1,-1],P_wheel(3,9)*[1,1],'g-','linewidth',2)
c=colorbar
caxis([0 350])
c.Label.String = 'Energy Density (mJ/mm)';
colorbar('off')
hold off
box off
grid off
xlabel('X')
ylabel('Y')
zlabel('Z')
axis([-600 600 -700 700 -150 150])
axis off
axis equal
view([1 1 1])
title(['\itL_{T1}\rm=',num2str(2*(X0-disp_x(mm))+19,'%.0f'),', \itL_{T2}\rm=',num2str(2*Q_sec_M1(2,N_e_1+1,mm)+19,'%.0f')])
set(gca,'FontSize',16,'fontname','Time','linewidth',1)
set(gcf, 'Position', [100 100 380 280]*0.7);
set(gca, 'Position', [0.0  0.0  1 0.8])   %[x y width height]

%
mm=50
Q_sec=Q_sec_M1(:,:,mm)
Q_sec_rd1=Q_sec(:,1:(N_e_1+1));
Q_sec_rd2=Q_sec(:,(N_e_1+2):(2*(N_e_1+1)));
figure(4);
clf
hold on
Rod_ploting1(q_node1_rd1(:,:,:,mm),Duplication_M,L_co,C_coord,U_E_ele1_rd1(mm,:),Dens_E_node1_rd1(:,:,mm));
Rod_ploting1(q_node1_rd2(:,:,:,mm),Duplication_M,L_co,C_coord,U_E_ele1_rd2(mm,:),Dens_E_node1_rd2(:,:,mm));
%%plotting wheels
P_wheel(:,1)=Q_sec_rd1(1:3,N_e_0+1);
P_wheel(:,2)=Q_sec_rd2(1:3,N_e_0+1);
P_wheel(:,3)=mean(P_wheel(:,1:2),2);
P_wheel(:,4)=Q_sec_rd1(1:3,1);
P_wheel(:,5)=Q_sec_rd2(1:3,1);
P_wheel(:,6)=mean(P_wheel(:,4:5),2);
P_wheel(:,7)=Q_sec_rd1(1:3,end);
P_wheel(:,8)=Q_sec_rd2(1:3,end);
P_wheel(:,9)=mean(P_wheel(:,7:8),2);
RM_wheel=reshape(Q_sec(4:12,N_e_0+1),3,3);
Wheel_Connect_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M,1)
plot3(P_wheel(1,6)*[1,-1],P_wheel(2,6)*[1,1],P_wheel(3,6)*[1,1],'r-','linewidth',2)
plot3(P_wheel(1,9)*[1,1],P_wheel(2,9)*[1,-1],P_wheel(3,9)*[1,1],'g-','linewidth',2)
c=colorbar
caxis([0 350])
c.Label.String = 'Energy Density (mJ/mm)';
colorbar('off')
hold off
box off
grid off
xlabel('X')
ylabel('Y')
zlabel('Z')
axis([-600 600 -700 700 -150 150])
axis off
axis equal
view([1 1 1])
title(['\itL_{T1}\rm=',num2str(2*(X0-disp_x(mm))+9,'%.0f'),', \itL_{T2}\rm=',num2str(2*Q_sec_M1(2,N_e_1+1,mm),'%.0f')])
set(gca,'FontSize',16,'fontname','Time','linewidth',1)
set(gcf, 'Position', [100 100 380 280]*0.7);
set(gca, 'Position', [0.0  0.0  1 0.8])   %[x y width height]


%
mm=50
nn=79
Q_sec=Q_sec_M2{mm}(:,:,nn);
Q_sec_rd1=Q_sec(:,1:(N_e_1+1));
Q_sec_rd2=Q_sec(:,(N_e_1+2):(2*(N_e_1+1)));
figure(5);
clf
hold on
Rod_ploting1(q_node2_rd1{mm}(:,:,:,nn),Duplication_M,L_co,C_coord,U_E_ele2_rd1{mm}(nn,:),Dens_E_node2_rd1{mm}(:,:,nn))
Rod_ploting1(q_node2_rd2{mm}(:,:,:,nn),Duplication_M,L_co,C_coord,U_E_ele2_rd2{mm}(nn,:),Dens_E_node2_rd2{mm}(:,:,nn))
%%plotting wheels
P_wheel(:,1)=Q_sec_rd1(1:3,N_e_0+1);
P_wheel(:,2)=Q_sec_rd2(1:3,N_e_0+1);
P_wheel(:,3)=mean(P_wheel(:,1:2),2);
P_wheel(:,4)=Q_sec_rd1(1:3,1);
P_wheel(:,5)=Q_sec_rd2(1:3,1);
P_wheel(:,6)=mean(P_wheel(:,4:5),2);
P_wheel(:,7)=Q_sec_rd1(1:3,end);
P_wheel(:,8)=Q_sec_rd2(1:3,end);
P_wheel(:,9)=mean(P_wheel(:,7:8),2);
RM_wheel=reshape(Q_sec(4:12,N_e_0+1),3,3);
Wheel_Connect_ploting(P_wheel,RM_wheel,0.8*R_wheel,N_spoke,Duplication_M,1)
plot3(P_wheel(1,6)*[1,-1],P_wheel(2,6)*[1,1],P_wheel(3,6)*[1,1],'r-','linewidth',2)
plot3(P_wheel(1,9)*[1,1],P_wheel(2,9)*[1,-1],P_wheel(3,9)*[1,1],'g-','linewidth',2)
c=colorbar
caxis([0 350])
c.Label.String = 'Energy Density (mJ/mm)';
colorbar('off')
%%plotting wheels
P_wheel=Q_sec(1:3,N_e_0+1);
RM_wheel=reshape(Q_sec(4:12,N_e_0+1),3,3);
% Wheel_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M)
hold off
box off
grid off
h=gca;
%                 set(h, 'ZDir', 'reverse');

% axis([-350 350 -600 600 -400 200])
axis equal
axis off
view([1 1 1])
title(['\itL_{T1}\rm=',num2str(2*(X0-disp_x(mm)+4.5),'%.0f'),', \itL_{T2}\rm=',num2str(2*Q_sec(2,N_e_1+1)+9,'%.0f')])
xlabel('X')
ylabel('Y')
zlabel('Z')
set(gca,'FontSize',16,'fontname','Time','linewidth',1)
set(gcf, 'Position', [100 100 380 280]*0.7);
set(gca, 'Position', [0.0  0.0  1 0.8])   %[x y width height]

%
%
mm=50
nn=51
Q_sec=Q_sec_M3{mm}(:,:,nn);
Q_sec_rd1=Q_sec(:,1:(N_e_1+1));
Q_sec_rd2=Q_sec(:,(N_e_1+2):(2*(N_e_1+1)));
figure(6);
clf
hold on
Rod_ploting1(q_node3_rd1{mm}(:,:,:,nn),Duplication_M,L_co,C_coord,U_E_ele3_rd1{mm}(nn,:),Dens_E_node3_rd1{mm}(:,:,nn))
Rod_ploting1(q_node3_rd2{mm}(:,:,:,nn),Duplication_M,L_co,C_coord,U_E_ele3_rd2{mm}(nn,:),Dens_E_node3_rd2{mm}(:,:,nn))
%%plotting wheels
P_wheel(:,1)=Q_sec_rd1(1:3,N_e_0+1);
P_wheel(:,2)=Q_sec_rd2(1:3,N_e_0+1);
P_wheel(:,3)=mean(P_wheel(:,1:2),2);
P_wheel(:,4)=Q_sec_rd1(1:3,1);
P_wheel(:,5)=Q_sec_rd2(1:3,1);
P_wheel(:,6)=mean(P_wheel(:,4:5),2);
P_wheel(:,7)=Q_sec_rd1(1:3,end);
P_wheel(:,8)=Q_sec_rd2(1:3,end);
P_wheel(:,9)=mean(P_wheel(:,7:8),2);
RM_wheel=reshape(Q_sec(4:12,N_e_0+1),3,3);
Wheel_Connect_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M,0)
plot3(P_wheel(1,6)*[1,-1],P_wheel(2,6)*[1,1],P_wheel(3,6)*[1,1],'r-','linewidth',2)
plot3(P_wheel(1,9)*[1,1],P_wheel(2,9)*[1,-1],P_wheel(3,9)*[1,1],'g-','linewidth',2)
c=colorbar
caxis([0 350])
c.Label.String = 'Nominal Energy Density (mJ/mm)';
% colorbar('off')
%%plotting wheels
P_wheel=Q_sec(1:3,N_e_0+1);
RM_wheel=reshape(Q_sec(4:12,N_e_0+1),3,3);
% Wheel_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M)
hold off
box off
grid off
h=gca;
axis([-250 250 -250 250 -500 10])
axis equal
axis off
view([1 1 1])
% title(['\itL_{T1}\rm=',num2str(2*(X0-disp_x(mm)+4.5),'%.0f'),', \itL_{T2}\rm=',num2str(2*Q_sec(2,N_e_1+1),'%.0f')])
xlabel('X')
ylabel('Y')
zlabel('Z')
set(gca,'FontSize',16,'fontname','Time','linewidth',1)
set(gcf, 'Position', [100 100 380 280]*0.7);
% set(gcf, 'Position', [100 100 380 280]*2);
set(gca, 'Position', [0.0  0.1  0.65 0.9])   %[x y width height]

% 

F_tendon=F_tendon1(:,1:50)
F_tendon=[F_tendon,F_tendon2{50},F_tendon3{50}]/2
figure(7)
plot(-F_tendon(1,:),'LineWidth',2);
hold on
plot(F_tendon(2,:),'LineWidth',2);
hold off
% set (gca,'XDir','reverse');
xlabel('Tendon disp')
ylabel('Tension force (N)')
% ylabel(haxes(1), 'F_{T1}','linewidth',2);
% ylabel(haxes(2), 'F_{T2}','linewidth',2);

set(gca,'FontSize',16,'fontname','Times','linewidth',1)
legend('Tendon1','Tendon2','FontSize',14,'fontname','times new roman')
set(gcf, 'Position', [100 100 380 280].*[1 1 1 0.8]);
set(gca, 'Position', [0.2  0.3  0.75  0.6])   %[x y width height]

%%
v = VideoWriter('Eye_tendon_force.mp4', 'MPEG-4');
v.FrameRate = 30;  % Set the frame rate (frames per second)
v.Quality = 90;    % Set the video quality (0 to 100)
open(v);
pic_num=1
for ii=[1:130 130:length(F_tendon(1,:)) length(F_tendon(1,:))]
    figure(7)
    plot(-F_tendon(1,:),'r','LineWidth',2);
    hold on
    plot(F_tendon(2,:),'g','LineWidth',2);
    plot(ii,-F_tendon(1,ii),'ko','LineWidth',2)
    plot(ii,F_tendon(2,ii),'ko','LineWidth',2)
    hold off
    % set (gca,'XDir','reverse');
    xlabel('Simulation step k')
    ylabel('Tension force (N)')
    % ylabel(haxes(1), 'F_{T1}','linewidth',2);
    % ylabel(haxes(2), 'F_{T2}','linewidth',2);
    set(gca,'FontSize',16,'fontname','Times','linewidth',1)
%     legend('Tendon1','Tendon2','FontSize',14,'fontname','times new roman')
    set(gcf, 'Position', [100 100 380 280].*[1 1 1 0.8]);
    set(gca, 'Position', [0.2  0.3  0.75  0.6])   %[x y width height]
    drawnow;% 立即刷新当前绘图窗口，这是matlab绘图中动态展示的关键
    F = getframe(gcf);  % 获取当前绘图窗口的图片
    Im = frame2im(F);   % 返回与电影帧相关的图片数据
    [Amap, map] = rgb2ind(Im, 256); % 将RGB图片转化为索引图片
    if pic_num == 1
        imwrite(Amap, map, 'Eye_tendon_force.gif', 'gif', 'Loopcount', Inf, 'DelayTime', 0.1);
        % 将第一张图片写入‘sin.gif’文件中，并且将GIF播放次数设置成无穷，即保存的GIF图会一直动下去
    else
        imwrite(Amap, map,'Eye_tendon_force.gif','gif','WriteMode','append','DelayTime',0.1);
        % 依次将其他的图片写入到GIF文件当中
    end
    pic_num = pic_num + 1;
    writeVideo(v, F);
    pause(0.1)
end
close(v);
