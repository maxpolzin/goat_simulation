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


%% caculate the Bisymmetric configuration with different tendon 1 length with a quarter of the ring

%define the start point and frame
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

r1_P=[50 0 0]';
q1=[r1_P;r1_frame(:)];

%define the end points and frames
r2_frame=[-1 0 0; 0 -1 0;  0 0 1]';

% Tune the start frame
alpha2=0*pi/180;%X axis
beta2=0*pi/180;%Y axis
theta2=0*pi/180;%Z axis

r2_frame=[-1 0 0; 0 -1 0;  0 0 1]';
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

r2_p=[0 50 0]';
q2=[r2_p;r2_frame(:)];

%generate the original solution by interplotion
Q_sec=Curve_interp(q1,q2,L_sec_0,N_e_0);

% 1 caculate the boundary of bisymmetric configuration
% Build the boundary conditions for the boundary
% start point x_disp in control and end point y_disp free
%Using the de_demension matrix to define the constrain of the structure
Proj_M=eye(12*(N_e_0+1),12*(N_e_0+1));%dedemesionalize the variables
% start point
D_fix=0*Q_sec;%fixed dimensions/varibales
D_fix(1,1)=1;%control the X disp
D_fix(2,1)=1;%fix the Y disp
D_fix(3,1)=1;%fix the  Z disp
% fix the local x axis
D_fix(4:6,1)=1;%fix the rotation

%end point
D_fix(1,end)=1;%fix the X disp
D_fix(3,end)=1;%fix the Z disp

% Fix local x aixs
D_fix([4 5 6],end)=1;%fix the rotation
D_fix([7 8 9],end)=1;%fix the rotation
D_fix([10 11 12],end)=1;%fix the rotation

B_index=find(D_fix==0);% find the fixed demensions
Proj_M=Proj_M(:,B_index);% extract the projection matrix of free paramters


N_x=50;
disp_x=linspace(0,R0_ring,N_x+1);% sampling displacement x
L_co=2;%the length of the coordinates
C_coord=[.0 .45 .74; .74 .45 .0];%the color of the coordinates


for mm=1:N_x+1
    mm
    if mm~=1
        Q_sec=Q_sec_M1(:,:,mm-1);
    end    
    Q_sec(1,1)=R0_ring-disp_x(mm);
    %define the end points and frames
    kk=0;
    flag1=1;
    flag3=0;
    while (kk==0||flag1>0.00001||flag2>0.10)&&kk<2000
        kk=kk+1;
        [dFq_dq_sec,Fq_sec_M1(:,mm),q_node,U_M1(mm)] =Jocob_rod_sec(Q_sec,N_e_0,L_e_0,N_node,Par_E,A);
        dFq_dq_sec_2=Proj_M'*dFq_dq_sec*Proj_M;
        dQ_sec_2=-pinv(dFq_dq_sec_2)*Proj_M'*Fq_sec_M1(:,mm);
        dQ_sec=Proj_M*dQ_sec_2;
        dQ_sec=reshape(dQ_sec,12,[]);
        Q_sec=Q_sec+dQ_sec;
        Q_sec_M1(:,:,mm)=Q_sec;
        flag1=sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e_0+1);
        flag2=mean(abs(Proj_M'*Fq_sec_M1(:,mm)));
        if kk==2000
            flag3=1;
        end
    end
    if mod(mm,5)==1
        figure(1);
        for ii=1:N_e_0
            r_node_1=q_node(1:3,:,ii);
            dr_dx_1=q_node(4:6,:,ii);
            dr_dy_1=q_node(7:9,:,ii);
            dr_dz_1=q_node(10:12,:,ii);            
            hold on
            plot3(r_node_1(1,:),r_node_1(2,:),r_node_1(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
            plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dy_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dy_1(2,:)], ...
                [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dy_1(3,:)],'-','color',C_coord(1,:),'linewidth',1)
            plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dz_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dz_1(2,:)], ...
                [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dz_1(3,:)],'-','color',C_coord(2,:),'linewidth',1)

            r_node_2=[-1 1 1]'.*r_node_1;
            dr_dx_2=[-1 1 1]'.*dr_dx_1;
            dr_dy_2=[-1 1 1]'.*dr_dy_1;
            dr_dz_2=[-1 1 1]'.*dr_dz_1;
            plot3(r_node_2(1,:),r_node_2(2,:),r_node_2(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
            plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dy_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dy_2(2,:)], ...
                [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dy_2(3,:)],'-','color',C_coord(1,:),'linewidth',1)
            plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dz_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dz_2(2,:)], ...
                [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dz_2(3,:)],'-','color',C_coord(2,:),'linewidth',1)

            r_node_3=[-1 -1 1]'.*r_node_1;
            dr_dx_3=[-1 -1 1]'.*dr_dx_1;
            dr_dy_3=[-1 -1 1]'.*dr_dy_1;
            dr_dz_3=[-1 -1 1]'.*dr_dz_1;
            plot3(r_node_3(1,:),r_node_3(2,:),r_node_3(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
            plot3([r_node_3(1,:);r_node_3(1,:)+L_co*dr_dy_3(1,:)],[r_node_3(2,:);r_node_3(2,:)+L_co*dr_dy_3(2,:)], ...
                [r_node_3(3,:);r_node_3(3,:)+L_co*dr_dy_3(3,:)],'-','color',C_coord(1,:),'linewidth',1)
            plot3([r_node_3(1,:);r_node_3(1,:)+L_co*dr_dz_3(1,:)],[r_node_3(2,:);r_node_3(2,:)+L_co*dr_dz_3(2,:)], ...
                [r_node_3(3,:);r_node_3(3,:)+L_co*dr_dz_3(3,:)],'-','color',C_coord(2,:),'linewidth',1)

            r_node_4=[1 -1 1]'.*r_node_1;
            dr_dx_4=[1 -1 1]'.*dr_dx_1;
            dr_dy_4=[1 -1 1]'.*dr_dy_1;
            dr_dz_4=[1 -1 1]'.*dr_dz_1;
            plot3(r_node_4(1,:),r_node_4(2,:),r_node_4(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
            plot3([r_node_4(1,:);r_node_4(1,:)+L_co*dr_dy_4(1,:)],[r_node_4(2,:);r_node_4(2,:)+L_co*dr_dy_4(2,:)], ...
                [r_node_4(3,:);r_node_4(3,:)+L_co*dr_dy_4(3,:)],'-','color',C_coord(1,:),'linewidth',1)
            plot3([r_node_4(1,:);r_node_4(1,:)+L_co*dr_dz_4(1,:)],[r_node_4(2,:);r_node_4(2,:)+L_co*dr_dz_4(2,:)], ...
                [r_node_4(3,:);r_node_4(3,:)+L_co*dr_dz_4(3,:)],'-','color',C_coord(2,:),'linewidth',1)
        end
        box on
        grid on
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        axis equal
        view([1 1 1])
        pause(0.1)
    end
end

figure(2);
disp_xy_B=[reshape(Q_sec_M1(1,1,:),1,[]);reshape(Q_sec_M1(2,end,:),1,[])];
plot(disp_xy_B(1,:),disp_xy_B(2,:))
hold on
plot(disp_xy_B(2,:),disp_xy_B(1,:))
hold off
box on
grid on
xlabel('X')
ylabel('Y')
axis equal
axis([0 70 0 70])

save('Rod_BiSym_Folding.mat',"Q_sec_M1",'-append')
%% caculate bisymmetric configurations in the range inside the boundary
load('Rod_BiSym_Folding.mat',"Q_sec_M1")
A=0.1e8;
%%32-51:0.1e8
%31-29:0.05e8
%28-25:0.02e8
%24-14:0.05e8
%13-11:0.1e8
%10-9:0.05e8
%8-1:0.02e8
% Build the boundary conditions for the boundary
% start point x_disp in control and end point y_disp free
%Using the de_demension matrix to define the constrain of the structure
Proj_M=eye(12*(N_e_0+1),12*(N_e_0+1));%dedemesionalize the variables
% start point
D_fix=0*Q_sec;%fixed dimensions/varibales
D_fix(1,1)=1;%control the X disp
D_fix(2,1)=1;%fix the Y disp
% D_fix(3,1)=1;%fix the  Z disp

% fix the local x y z axis
D_fix(4:6,1)=1;%fix the rotation
D_fix(7:9,1)=1;%fix the rotation
D_fix(10:11,1)=1;%fix the rotation


%end point
D_fix(1,end)=1;%fix the X disp
D_fix(2,end)=1;%control the Y disp
D_fix(3,end)=1;%control the Z disp

% Fix local x aixs
D_fix([4 5 6],end)=1;%fix the rotation
% D_fix([7 8 9],end)=1;%fix the rotation
% D_fix([10 11 12],end)=1;%fix the rotation

B_index=find(D_fix==0);% find the fixed demensions
Proj_M=Proj_M(:,B_index);% extract the projection matrix of free paramters


% since the symmetry of the structure, just a half area need to be caculate
Fq_sec=[]
U=[]
for mm=51%N_x+1  
    figure(3);
    clf
    Q_sec=Q_sec_M1(:,:,mm);
    % caculate the range of disp y based on the boundary 
    disp_y_0=Q_sec(2,end);
    N_y=fix((Q_sec(2,end)-Q_sec(1,1))/1)
    disp_y=flip([Q_sec(1,1):1:disp_y_0, disp_y_0]);% sampling displacement y
    for nn=1:length(disp_y)
        mm
        nn
        if nn~=1
            Q_sec=Q_sec_M2{mm}(:,:,nn-1);
        end
        Q_sec(2,end)=disp_y(nn);
        if ismember(nn,[1:20])
            Q_sec(3,1)=Q_sec(3,1)+1;
        else
            Q_sec(3,1)=Q_sec(3,1);
        end
        %define the end points and frames
        kk=0;
        flag1=1;
        while (kk==0||flag1>0.00001||flag2>0.10)&&kk<5000
            kk=kk+1;
            [dFq_dq_sec,Fq_sec(:,mm),q_node,U(mm)] = Jocob_rod_sec(Q_sec,N_e_0,L_e_0,N_node,Par_E,A);
            dFq_dq_sec_2=Proj_M'*dFq_dq_sec*Proj_M;
            dQ_sec_2=-pinv(dFq_dq_sec_2)*Proj_M'*Fq_sec(:,mm);
            dQ_sec=Proj_M*dQ_sec_2;
            dQ_sec=reshape(dQ_sec,12,[]);
            Q_sec=Q_sec+dQ_sec;
            Q_sec_M2{mm}(:,:,nn)=Q_sec;
            flag1=sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e_0+1);
            flag2=mean(abs(Proj_M'*Fq_sec(:,mm)));
            if kk==5000
                flag3=1;
            end
        end
        if mod(nn,5)==1
            figure(3);
            for ii=1:N_e_0
                r_node_1=q_node(1:3,:,ii);
                dr_dx_1=q_node(4:6,:,ii);
                dr_dy_1=q_node(7:9,:,ii);
                dr_dz_1=q_node(10:12,:,ii);
                
                hold on
                plot3(r_node_1(1,:),r_node_1(2,:),r_node_1(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
                plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dy_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dy_1(2,:)], ...
                    [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dy_1(3,:)],'-','color',C_coord(1,:),'linewidth',1)
                plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dz_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dz_1(2,:)], ...
                    [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dz_1(3,:)],'-','color',C_coord(2,:),'linewidth',1)

                r_node_2=[-1 1 1]'.*r_node_1;
                dr_dx_2=[-1 1 1]'.*dr_dx_1;
                dr_dy_2=[-1 1 1]'.*dr_dy_1;
                dr_dz_2=[-1 1 1]'.*dr_dz_1;
                plot3(r_node_2(1,:),r_node_2(2,:),r_node_2(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
                plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dy_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dy_2(2,:)], ...
                    [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dy_2(3,:)],'-','color',C_coord(1,:),'linewidth',1)
                plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dz_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dz_2(2,:)], ...
                    [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dz_2(3,:)],'-','color',C_coord(2,:),'linewidth',1)

                r_node_3=[-1 -1 1]'.*r_node_1;
                dr_dx_3=[-1 -1 1]'.*dr_dx_1;
                dr_dy_3=[-1 -1 1]'.*dr_dy_1;
                dr_dz_3=[-1 -1 1]'.*dr_dz_1;
                plot3(r_node_3(1,:),r_node_3(2,:),r_node_3(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
                plot3([r_node_3(1,:);r_node_3(1,:)+L_co*dr_dy_3(1,:)],[r_node_3(2,:);r_node_3(2,:)+L_co*dr_dy_3(2,:)], ...
                    [r_node_3(3,:);r_node_3(3,:)+L_co*dr_dy_3(3,:)],'-','color',C_coord(1,:),'linewidth',1)
                plot3([r_node_3(1,:);r_node_3(1,:)+L_co*dr_dz_3(1,:)],[r_node_3(2,:);r_node_3(2,:)+L_co*dr_dz_3(2,:)], ...
                    [r_node_3(3,:);r_node_3(3,:)+L_co*dr_dz_3(3,:)],'-','color',C_coord(2,:),'linewidth',1)

                r_node_4=[1 -1 1]'.*r_node_1;
                dr_dx_4=[1 -1 1]'.*dr_dx_1;
                dr_dy_4=[1 -1 1]'.*dr_dy_1;
                dr_dz_4=[1 -1 1]'.*dr_dz_1;
                plot3(r_node_4(1,:),r_node_4(2,:),r_node_4(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
                plot3([r_node_4(1,:);r_node_4(1,:)+L_co*dr_dy_4(1,:)],[r_node_4(2,:);r_node_4(2,:)+L_co*dr_dy_4(2,:)], ...
                    [r_node_4(3,:);r_node_4(3,:)+L_co*dr_dy_4(3,:)],'-','color',C_coord(1,:),'linewidth',1)
                plot3([r_node_4(1,:);r_node_4(1,:)+L_co*dr_dz_4(1,:)],[r_node_4(2,:);r_node_4(2,:)+L_co*dr_dz_4(2,:)], ...
                    [r_node_4(3,:);r_node_4(3,:)+L_co*dr_dz_4(3,:)],'-','color',C_coord(2,:),'linewidth',1)

            end
            box on
            grid on
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            axis equal
            view([1 1 1])
            pause(0.1)
        end

    end
    X_1=reshape(Q_sec_M2{mm}(1,1,:),1,[]);
    Y_2=reshape(Q_sec_M2{mm}(2,end,:),1,[]);
    Z_2=reshape(Q_sec_M2{mm}(3,1,:),1,[]);

    figure(2);
    hold on
    plot3(X_1,Y_2,Z_2,'-o')
    box on
    grid on
    axis equal
    hold off
end
save('Rod_BiSym_Folding.mat',"Q_sec_M2",'-append')

%%

for mm=1:51
    figure(4);
    clf

nn=size(Q_sec_M2{mm}(:,:,:),3)
C_coord=[.0 .45 .74; .74 .45 .0];
L_co=5;
Q_sec=Q_sec_M2{mm}(:,:,nn);
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

r_node_2=[-1 1 1]'.*Q_sec(1:3,:);
dr_dx_2=[-1 1 1]'.*Q_sec(4:6,:);
dr_dy_2=[-1 1 1]'.*Q_sec(7:9,:);
dr_dz_2=[-1 1 1]'.*Q_sec(10:12,:);

plot3(r_node_2(1,:),r_node_2(2,:),r_node_2(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dy_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dy_2(2,:)], ...
    [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dy_2(3,:)],'-','color',C_coord(1,:),'linewidth',1)
plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dz_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dz_2(2,:)], ...
    [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dz_2(3,:)],'-','color',C_coord(2,:),'linewidth',1)

r_node_3=[-1 -1 1]'.*r_node_1;
dr_dx_3=[-1 -1 1]'.*dr_dx_1;
dr_dy_3=[-1 -1 1]'.*dr_dy_1;
dr_dz_3=[-1 -1 1]'.*dr_dz_1;
plot3(r_node_3(1,:),r_node_3(2,:),r_node_3(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
plot3([r_node_3(1,:);r_node_3(1,:)+L_co*dr_dy_3(1,:)],[r_node_3(2,:);r_node_3(2,:)+L_co*dr_dy_3(2,:)], ...
    [r_node_3(3,:);r_node_3(3,:)+L_co*dr_dy_3(3,:)],'-','color',C_coord(1,:),'linewidth',1)
plot3([r_node_3(1,:);r_node_3(1,:)+L_co*dr_dz_3(1,:)],[r_node_3(2,:);r_node_3(2,:)+L_co*dr_dz_3(2,:)], ...
    [r_node_3(3,:);r_node_3(3,:)+L_co*dr_dz_3(3,:)],'-','color',C_coord(2,:),'linewidth',1)

r_node_4=[1 -1 1]'.*r_node_1;
dr_dx_4=[1 -1 1]'.*dr_dx_1;
dr_dy_4=[1 -1 1]'.*dr_dy_1;
dr_dz_4=[1 -1 1]'.*dr_dz_1;
plot3(r_node_4(1,:),r_node_4(2,:),r_node_4(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
plot3([r_node_4(1,:);r_node_4(1,:)+L_co*dr_dy_4(1,:)],[r_node_4(2,:);r_node_4(2,:)+L_co*dr_dy_4(2,:)], ...
    [r_node_4(3,:);r_node_4(3,:)+L_co*dr_dy_4(3,:)],'-','color',C_coord(1,:),'linewidth',1)
plot3([r_node_4(1,:);r_node_4(1,:)+L_co*dr_dz_4(1,:)],[r_node_4(2,:);r_node_4(2,:)+L_co*dr_dz_4(2,:)], ...
    [r_node_4(3,:);r_node_4(3,:)+L_co*dr_dz_4(3,:)],'-','color',C_coord(2,:),'linewidth',1)

box on
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
axis([-70 70 -70 70 -70 70])
view([1 1 1])
hold off
end
%% Evaluate the stability of every bisym config by applying twisting

% double the section to make it as anti-sym
% the origin config
load("Rod_BiSym_Folding.mat","Q_sec_M2")% load the twist or flat configuration as the original configuration
Q_sec=[];
Q_sec_1=[];
Q_sec_2=[];
Fq_sec=[];
U=[];

mm=46
nn=30

Q_sec_1=Q_sec_M2{mm}(:,:,nn);
N_e_1=2*(size(Q_sec_1,2)-1);% get the element number from the data
Q_sec_2(1:3,:)=diag([1,-1,1])*Q_sec_1(1:3,:);
Q_sec_2(4:6,:)=diag([-1,-1,-1])*diag([1,-1,1])*Q_sec_1(4:6,:);
Q_sec_2(7:9,:)=diag([1,1,1])*diag([1,-1,1])*Q_sec_1(7:9,:);
Q_sec_2(10:12,:)=diag([1,1,1])*diag([1,-1,1])*Q_sec_1(10:12,:);
Q_sec=[flip(Q_sec_2(:,2:end),2),Q_sec_1];

L_sec_1=L_sec_0*2;
L_e_1=L_sec_1/N_e_1;

%plot the origin config
figure(4);
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

r_node_2=[-1 1 1]'.*Q_sec(1:3,:);
dr_dx_2=[-1 1 1]'.*Q_sec(4:6,:);
dr_dy_2=[-1 1 1]'.*Q_sec(7:9,:);
dr_dz_2=[-1 1 1]'.*Q_sec(10:12,:);

plot3(r_node_2(1,:),r_node_2(2,:),r_node_2(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dy_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dy_2(2,:)], ...
    [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dy_2(3,:)],'-','color',C_coord(1,:),'linewidth',1)
plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dz_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dz_2(2,:)], ...
    [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dz_2(3,:)],'-','color',C_coord(2,:),'linewidth',1)

box on
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
axis([-70 70 -70 70 -50 50])
view([1 1 1])
hold off

% Build the new boundary conditions for the stability analysis

A=0.05e8;


% start point x_disp in control and end point y_disp free
%Using the de_demension matrix to define the constrain of the structure
Proj_M_0=eye(12*(N_e_1+1),12*(N_e_1+1));%dedemesionalize the variables


% based on anti sym condition

% the start and end point 
% start point
D_fix=0*Q_sec;%fixed dimensions/varibales
D_fix(1,1)=1;%control the X disp
D_fix(2,1)=1;%fix the Y disp
D_fix(3,1)=1;%fix the  Z disp

%control the local rotation by constrain the local z axis
D_fix(10:12,1)=1;%local z axis in control 
r1_frame_0=reshape(Q_sec(4:12,1),3,3);

% the local x y z axis in control
% because the equalism relationship
%D_fix(1:3,1)=[-1 -1 1].*D_fix(4:6,end);%x axis in control 
%D_fix(4:6,1)=[-1 -1 1].*D_fix(4:6,end);%x axis in control 
%D_fix(7:9,1)=[-1 -1 1].*D_fix(4:6,end);%x axis in control 
%D_fix(10:12,1)=[-1 -1 1].*D_fix(4:6,end);%x axis in control 
% sus
Con_M_1=diag(repmat([-1 -1 1],1,4));
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

n_m=N_e_1/2+1;%the seq numb of middle point

% the middle point position
% based on anti sym condition
% the distance of middle point to the z axis is constrained: x^2+y^2=C
% introduce the R and theta to replace the x y disp of middle point
r_m_0=Q_sec(1,n_m);
theta_m_0=0;
% the present First-order expansion relationship between d(theta,R) and d(x,y)
Con_M_2_0=[-sin(theta_m_0)*r_m_0 cos(theta_m_0); cos(theta_m_0)*r_m_0 sin(theta_m_0)];
Con_M_2_2=[-sin(theta_m_0)/r_m_0 cos(theta_m_0)/r_m_0; cos(theta_m_0) sin(theta_m_0)];
Proj_M_0((n_m-1)*12+[1:2],(n_m-1)*12+[1:2])=Con_M_2_0;
%fix R: the distance to the z axis
D_fix(2,n_m)=1;%fix the R
D_fix(3,n_m)=1;%fix the Z



T_angle=linspace(0,90*pi/180,46);

% Q_sec(:,1)
% Q_sec(:,end);
figure(4);
    clf
for ii=1:46
    ii    
    figure(3);
    clf
    if ii~=1
        Q_sec=Q_sec_M3{mm}(:,:,ii-1,nn);
    end
    %redefine the local z axis of start/end points
    R_M_z=[cos(T_angle(ii))  0  -sin(T_angle(ii)) ;
           0                 1               0;
           sin(T_angle(ii))  0  cos(T_angle(ii)) ;];
    r1_frame_1=r1_frame_0*R_M_z;
    Q_sec(10:12,1)=r1_frame_1(:,3);
    Q_sec(10:12,end)=[-1 -1 1]'.*r1_frame_1(:,3);


    kk=0;
    flag1=1;
    L_step=0.05;
   
    while (kk==0||flag1>0.0001*L_step||kk<2/L_step)&&kk<20/L_step
        kk=kk+1;
        [dFq_dq_sec,Fq_sec(:,ii),q_node,U(ii),E_elastic{mm}(:,ii,nn)] = Jocob_rod_sec(Q_sec,N_e_1,L_e_1,N_node,Par_E,A);
        % update the boundary information
        if r_m_0<10e-4
            r_m=0;
            theta_m=theta_m_0;
            D_fix(1:2,n_m)=1;
        else
            D_fix(1,n_m)=0;
            r_m=sqrt(sum(Q_sec(1:2,n_m).^2));
            theta_m=atan2(Q_sec(2,n_m),Q_sec(1,n_m));
            % the present First-order expansion relationship between d(theta,R) and d(x,y)
            Con_M_2_0=[-sin(theta_m)*r_m cos(theta_m); cos(theta_m)*r_m sin(theta_m)];
            Con_M_2_2=[-sin(theta_m)/r_m cos(theta_m)/r_m; cos(theta_m) sin(theta_m)];
            Proj_M_0((n_m-1)*12+[1:2],(n_m-1)*12+[1:2])=Con_M_2_0;
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
            d_theta=Con_M_2_2(1,:)*dQ_sec([1:2],n_m);
            theta_m=theta_m+d_theta;
            Q_sec(1:2,n_m)=r_m_0*[cos(theta_m);sin(theta_m)];
        end
        Q_sec_M3{mm}(:,:,ii,nn)=Q_sec;
        flag1=sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e_1+1);
        flag2=mean(abs(Proj_M'*Fq_sec(:,ii)));
        if kk==5000
            flag3=1
        elseif mod(kk,100)==0
            kk
            flag1
            figure(3);
            for jj=1:N_e_1
                r_node_1=q_node(1:3,:,jj);
                dr_dx_1=q_node(4:6,:,jj);
                dr_dy_1=q_node(7:9,:,jj);
                dr_dz_1=q_node(10:12,:,jj);

                hold on
                plot3(r_node_1(1,:),r_node_1(2,:),r_node_1(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
                plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dy_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dy_1(2,:)], ...
                    [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dy_1(3,:)],'-','color',C_coord(1,:),'linewidth',1)
                plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dz_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dz_1(2,:)], ...
                    [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dz_1(3,:)],'-','color',C_coord(2,:),'linewidth',1)

                r_node_2=[-1 -1 1]'.*q_node(1:3,:,jj);
                dr_dx_2=[-1 -1 1]'.*q_node(4:6,:,jj);
                dr_dy_2=[-1 -1 1]'.*q_node(7:9,:,jj);
                dr_dz_2=[-1 -1 1]'.*q_node(10:12,:,jj);

                plot3(r_node_2(1,:),r_node_2(2,:),r_node_2(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
                plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dy_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dy_2(2,:)], ...
                    [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dy_2(3,:)],'-','color',C_coord(1,:),'linewidth',1)
                plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dz_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dz_2(2,:)], ...
                    [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dz_2(3,:)],'-','color',C_coord(2,:),'linewidth',1)
                hold off

            end
            axis equal
            box on
            grid on
            view([1 1 1])
            figure(5);
            hold on
            plot(E_elastic{mm}(1,1:ii,nn),'-*')
            hold off
            pause(0.01)
        end

    end


    figure(4);
    hold on
    for jj=1:N_e_1

        r_node_1=q_node(1:3,:,jj);
        dr_dx_1=q_node(4:6,:,jj);
        dr_dy_1=q_node(7:9,:,jj);
        dr_dz_1=q_node(10:12,:,jj);

        plot3(r_node_1(1,:),r_node_1(2,:),r_node_1(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
        plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dy_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dy_1(2,:)], ...
            [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dy_1(3,:)],'-','color',C_coord(1,:),'linewidth',1)
        plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dz_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dz_1(2,:)], ...
            [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dz_1(3,:)],'-','color',C_coord(2,:),'linewidth',1)

        r_node_2=[-1 -1 1]'.*q_node(1:3,:,jj);
        dr_dx_2=[-1 -1 1]'.*q_node(4:6,:,jj);
        dr_dy_2=[-1 -1 1]'.*q_node(7:9,:,jj);
        dr_dz_2=[-1 -1 1]'.*q_node(10:12,:,jj);

        plot3(r_node_2(1,:),r_node_2(2,:),r_node_2(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
        plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dy_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dy_2(2,:)], ...
            [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dy_2(3,:)],'-','color',C_coord(1,:),'linewidth',1)
        plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dz_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dz_2(2,:)], ...
            [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dz_2(3,:)],'-','color',C_coord(2,:),'linewidth',1)

    end
    plot3(Q_sec(1,[1,N_e_1+1]),Q_sec(2,[1,N_e_1+1]),Q_sec(3,[1,N_e_1+1]),'-','color','b','linewidth',1)
    plot3([1,-1]*Q_sec(1,n_m),[1,-1]*Q_sec(2,n_m),[1,1]*Q_sec(3,n_m),'-','color','r','linewidth',1)
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
jj=15
figure(6);
hold on
plot3(nn*ones(jj-ii+1),ii:jj,E_elastic{mm}(1,ii:jj,nn),'-*')

hold off
pause(0.1)
view([1 1 1])
box on
grid on

 save('Rod_BiSym_Folding.mat','E_elastic','Q_sec_M3','-append')

%%

%plot the origin config
figure(7)
ll=13
Q_sec=Q_sec_M3{mm}(:,:,ll,nn)
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

r_node_2=[-1 -1 1]'.*Q_sec(1:3,:);
dr_dx_2=[-1 -1 1]'.*Q_sec(4:6,:);
dr_dy_2=[-1 -1 1]'.*Q_sec(7:9,:);
dr_dz_2=[-1 -1 1]'.*Q_sec(10:12,:);

plot3(r_node_2(1,:),r_node_2(2,:),r_node_2(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dy_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dy_2(2,:)], ...
    [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dy_2(3,:)],'-','color',C_coord(1,:),'linewidth',1)
plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dz_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dz_2(2,:)], ...
    [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dz_2(3,:)],'-','color',C_coord(2,:),'linewidth',1)

plot3(r_node_2(1,[1,N_e_1+1]),r_node_2(2,[1,N_e_1+1]),r_node_2(3,[1,N_e_1+1]),'-','color','b','linewidth',1)
plot3([1,-1]*r_node_2(1,n_m),[1,-1]*r_node_2(2,n_m),[1,1]*r_node_2(3,n_m),'-','color','r','linewidth',1)
set(gca,'Zdir','reverse')
box on
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
axis([-70 70 -70 70 -50 50])
view([1 1 1])
hold off



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
            kk
            flag1
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

