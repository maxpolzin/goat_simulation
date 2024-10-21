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
A=0.01e8;%the penalty parameter to constrain the frame vector Y Z as orthonormal
% it related to the elastic energy, when it too high the iteration would be too slow
% when it is too low, the frame vector Y Z won't be orthonormal anymore.
%% caculate the Bisymmetric configuration with different tendon lengths with a quarter of the ring

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
Q_sec=Curve_interp(q1,q2,L_sec,N_e);

% 1 caculate the boundary of bisymmetric configuration
% Build the boundary conditions for the boundary
% start point x_disp in control and end point y_disp free
%Using the de_demension matrix to define the constrain of the structure
Proj_M=eye(12*(N_e+1),12*(N_e+1));%dedemesionalize the variables
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
pic_num=1
for mm=1:N_x+1
    mm
    if mm~=1
        Q_sec=Q_sec_M1(:,:,mm-1)
    end    
    Q_sec(1,1)=R0_ring-disp_x(mm);
    %define the end points and frames
    kk=0;
    flag1=1;
    flag3=0;
    while (kk==0||flag1>0.00001||flag2>0.10)&&kk<2000
        kk=kk+1;
        [dFq_dq_sec,Fq_sec_M1(:,mm),q_node,U_M1(mm)] =Jocob_rod_sec(Q_sec,N_e,L_e,N_node,Par_E,A);
        dFq_dq_sec_2=Proj_M'*dFq_dq_sec*Proj_M;
        dQ_sec_2=-pinv(dFq_dq_sec_2)*Proj_M'*Fq_sec_M1(:,mm);
        dQ_sec=Proj_M*dQ_sec_2;
        dQ_sec=reshape(dQ_sec,12,[]);
        Q_sec=Q_sec+dQ_sec;
        Q_sec_M1(:,:,mm)=Q_sec;
        flag1=sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e+1);
        flag2=mean(abs(Proj_M'*Fq_sec_M1(:,mm)));
        if kk==2000
            flag3=1;
        end
    end
    if mod(mm,5)==1
        for ii=1:N_e
            r_node_1=q_node(1:3,:,ii);
            dr_dx_1=q_node(4:6,:,ii);
            dr_dy_1=q_node(7:9,:,ii);
            dr_dz_1=q_node(10:12,:,ii);

            figure(1)
            plot3(r_node_1(1,:),r_node_1(2,:),r_node_1(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
            hold on
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
        drawnow;% 立即刷新当前绘图窗口，这是matlab绘图中动态展示的关键
        F = getframe(gcf);  % 获取当前绘图窗口的图片
        Im = frame2im(F);   % 返回与电影帧相关的图片数据
        [Amap, map] = rgb2ind(Im, 256); % 将RGB图片转化为索引图片
        if pic_num == 1
            imwrite(Amap, map, 'Inplane.gif', 'gif', 'Loopcount', Inf, 'DelayTime', 0.1);
            % 将第一张图片写入‘sin.gif’文件中，并且将GIF播放次数设置成无穷，即保存的GIF图会一直动下去
        else
            imwrite(Amap, map,'Inplane.gif','gif','WriteMode','append','DelayTime',0.1);
            % 依次将其他的图片写入到GIF文件当中
        end
        pic_num = pic_num + 1;
        pause(0.1)

    end
end

figure(10)
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

%% caculate bisymmetric configurations in the range inside the boundary

% Build the boundary conditions for the boundary
% start point x_disp in control and end point y_disp free
%Using the de_demension matrix to define the constrain of the structure
Proj_M=eye(12*(N_e+1),12*(N_e+1));%dedemesionalize the variables
% start point
D_fix=0*Q_sec;%fixed dimensions/varibales
D_fix(1,1)=1;%control the X disp
D_fix(2,1)=1;%fix the Y disp
D_fix(3,1)=1;%fix the  Z disp
% fix the local x axis
D_fix(4:6,1)=1;%fix the rotation

%end point
D_fix(1,end)=1;%fix the X disp
D_fix(2,end)=1;%control the Y disp

% Fix local x aixs
D_fix([4 5 6],end)=1;%fix the rotation
D_fix([7 8 9],end)=1;%fix the rotation
D_fix([10 11 12],end)=1;%fix the rotation

B_index=find(D_fix==0);% find the fixed demensions
Proj_M=Proj_M(:,B_index);% extract the projection matrix of free paramters


% since the symmetry of the structure, just a half area need to be caculate
pic_num=1
for mm=30%N_x+1  
    figure(2)
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
            Q_sec(3,end)=Q_sec(3,end)+2;
        else
            Q_sec(3,end)=Q_sec(3,end);
        end
        %define the end points and frames
        kk=0;
        flag1=1;
        while (kk==0||flag1>0.00001||flag2>0.10)&&kk<5000
            kk=kk+1;
            [dFq_dq_sec,Fq_sec(:,mm),q_node,U(mm)] = Jocob_rod_sec(Q_sec,N_e,L_e,N_node,Par_E,A);
            dFq_dq_sec_2=Proj_M'*dFq_dq_sec*Proj_M;
            dQ_sec_2=-pinv(dFq_dq_sec_2)*Proj_M'*Fq_sec(:,mm);
            dQ_sec=Proj_M*dQ_sec_2;
            dQ_sec=reshape(dQ_sec,12,[]);
            Q_sec=Q_sec+dQ_sec;
            Q_sec_M2{mm}(:,:,nn)=Q_sec;
            flag1=sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e+1);
            flag2=mean(abs(Proj_M'*Fq_sec(:,mm)));
            if kk==5000
                flag3=1;
            end
        end
        if mod(nn,2)==1
            for ii=1:N_e
                r_node_1=q_node(1:3,:,ii);
                dr_dx_1=q_node(4:6,:,ii);
                dr_dy_1=q_node(7:9,:,ii);
                dr_dz_1=q_node(10:12,:,ii);

                figure(2)
                plot3(r_node_1(1,:),r_node_1(2,:),r_node_1(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
                hold on                
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
            drawnow;% 立即刷新当前绘图窗口，这是matlab绘图中动态展示的关键
            F = getframe(gcf);  % 获取当前绘图窗口的图片
            Im = frame2im(F);   % 返回与电影帧相关的图片数据
            [Amap, map] = rgb2ind(Im, 256); % 将RGB图片转化为索引图片
            if pic_num == 1
                imwrite(Amap, map, 'Outplane.gif', 'gif', 'Loopcount', Inf, 'DelayTime', 0.1);
                % 将第一张图片写入‘sin.gif’文件中，并且将GIF播放次数设置成无穷，即保存的GIF图会一直动下去
            else
                imwrite(Amap, map,'Outplane.gif','gif','WriteMode','append','DelayTime',0.1);
                % 依次将其他的图片写入到GIF文件当中
            end
            pic_num = pic_num + 1;
            pause(0.1)
            hold off
        end

    end
    X_1=reshape(Q_sec_M2{mm}(1,1,:),1,[])
    Y_2=reshape(Q_sec_M2{mm}(2,end,:),1,[])
    Z_2=reshape(Q_sec_M2{mm}(3,end,:),1,[])
    figure(20)
    hold on
    plot3(X_1,Y_2,Z_2,'-o')
    box on
    grid on
    axis equal
    hold off
end
% save('Rod_BiSym_Folding.mat',"Q_sec_M1","Q_sec_M2")

%% Evaluate the stability of every bisym config by applying twisting

% double the section to make it as anti-sym
% the origin config
load("Rod_BiSym_Folding.mat","Q_sec_M2")% load the twist or flat configuration as the original configuration

mm=51
nn=1
Q_sec_1=Q_sec_M2{mm}(:,:,nn);
N_e=2*(size(Q_sec_1,2)-1);% get the element number from the data
Q_sec_2(1:3,:)=diag([-1,1,1])*Q_sec_1(1:3,:);
Q_sec_2(4:6,:)=diag([-1,-1,-1])*diag([-1,1,1])*Q_sec_1(4:6,:);
Q_sec_2(7:9,:)=diag([1,1,1])*diag([-1,1,1])*Q_sec_1(7:9,:);
Q_sec_2(10:12,:)=diag([1,1,1])*diag([-1,1,1])*Q_sec_1(10:12,:);
Q_sec=[Q_sec_1,flip(Q_sec_2(:,1:end-1),2)];

L_sec=L_sec*2;
L_e=L_sec/N_e;

%plot the origin config
figure(1)
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

r_node_2=[1 -1 1]'.*Q_sec(1:3,:);
dr_dx_2=[1 -1 1]'.*Q_sec(4:6,:);
dr_dy_2=[1 -1 1]'.*Q_sec(7:9,:);
dr_dz_2=[1 -1 1]'.*Q_sec(10:12,:);

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
r_m=Q_sec(2,n_m);
theta_m=pi/2;
% the present First-order expansion relationship between d(theta,R) and d(x,y)
Con_M=[-sin(theta_m)/r_m cos(theta_m)/r_m; cos(theta_m) sin(theta_m)];
Proj_M_0((n_m-1)*12+[1:2],(n_m-1)*12+[1:2])=Con_M;
%fix R: the distance to the z axis
D_fix(2,n_m)=1;%fix the R

T_angle=linspace(0,10*pi/180,11);
r1_frame_0=reshape(Q_sec(4:12,1),3,3);
r2_frame_0=reshape(Q_sec(4:12,end),3,3);

figure(2)
clf
pic_num=1
dFq_dq_sec=[]
Fq_sec=[]
q_node=[]
U=[]
for ii=1:10
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
    % update the boundary information
    r_m=Q_sec(2,n_m);
    theta_m=pi/2;
    % the present First-order expansion relationship between d(theta,R) and d(x,y)
    Con_M=[-sin(theta_m)/r_m cos(theta_m)/r_m; cos(theta_m) sin(theta_m)];
    Proj_M_0((n_m-1)*12+[1:2],(n_m-1)*12+[1:2])=Con_M;
    B_index=find(D_fix==0);% find the fixed demensions
    Proj_M=Proj_M_0(:,B_index);% extract the projection matrix of free paramters


    kk=0;
    flag1=1;
    while (kk==0||flag1>0.00001||flag2>0.10)&&kk<5000
        kk=kk+1;
        [dFq_dq_sec,Fq_sec(:,ii),q_node,U(ii)] = Jocob_rod_sec(Q_sec,N_e,L_e,N_node,Par_E,A);
        dFq_dq_sec_2=Proj_M'*dFq_dq_sec*Proj_M;
        dQ_sec_2=-pinv(dFq_dq_sec_2)*Proj_M'*Fq_sec(:,ii);
        dQ_sec=Proj_M*dQ_sec_2;
        dQ_sec=reshape(dQ_sec,12,[]);
        Q_sec=Q_sec+dQ_sec;
        Q_sec_M3{mm,nn}(:,:,ii)=Q_sec;
        flag1=sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e+1);
        flag2=mean(abs(Proj_M'*Fq_sec(:,ii)));
        if kk==5000
            flag3=1;
        end

    end

    for ii=1:N_e
        r_node_1=q_node(1:3,:,ii);
        dr_dx_1=q_node(4:6,:,ii);
        dr_dy_1=q_node(7:9,:,ii);
        dr_dz_1=q_node(10:12,:,ii);

        figure(2)
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

        r_node_2=[1 -1 1]'.*Q_sec(1:3,:);
        dr_dx_2=[1 -1 1]'.*Q_sec(4:6,:);
        dr_dy_2=[1 -1 1]'.*Q_sec(7:9,:);
        dr_dz_2=[1 -1 1]'.*Q_sec(10:12,:);

        plot3(r_node_2(1,:),r_node_2(2,:),r_node_2(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
        plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dy_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dy_2(2,:)], ...
            [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dy_2(3,:)],'-','color',C_coord(1,:),'linewidth',1)
        plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dz_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dz_2(2,:)], ...
            [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dz_2(3,:)],'-','color',C_coord(2,:),'linewidth',1)

    end
    box on
    grid on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    axis equal
    axis([-70 70 -70 70 -50 50])
    view([1 1 1])
    drawnow;% 立即刷新当前绘图窗口，这是matlab绘图中动态展示的关键
    F = getframe(gcf);  % 获取当前绘图窗口的图片
    Im = frame2im(F);   % 返回与电影帧相关的图片数据
    [Amap, map] = rgb2ind(Im, 256); % 将RGB图片转化为索引图片
    if pic_num == 1
        imwrite(Amap, map, 'Twist.gif', 'gif', 'Loopcount', Inf, 'DelayTime', 0.1);
        % 将第一张图片写入‘sin.gif’文件中，并且将GIF播放次数设置成无穷，即保存的GIF图会一直动下去
    else
        imwrite(Amap, map,'Twist.gif','gif','WriteMode','append','DelayTime',0.1);
        % 依次将其他的图片写入到GIF文件当中
    end
    pic_num = pic_num + 1;
    hold off

end


figure(3)

plot(U)








%% define the caculation parameters of one section
%the element number for caculation
% the length of the element
% the initial solution of the section


% load set its original state
load("Rod_twist_2.mat","Q_sec_MM","U_M")% load the twist or flat configuration as the original configuration
N_disp=48%the disp number of the state 1:51 related to 0:50
[a b]=min(U_M(:,N_disp))
Q_sec_1=Q_sec_MM{N_disp}(:,:,b)
N_e=2*(size(Q_sec_1,2)-1);% get the element number from the data
Q_sec_2(1:3,:)=diag([-1,1,-1])*Q_sec_1(1:3,:)
Q_sec_2(4:6,:)=diag([-1,-1,-1])*diag([-1,1,-1])*Q_sec_1(4:6,:)
Q_sec_2(7:9,:)=diag([1,1,1])*diag([-1,1,-1])*Q_sec_1(7:9,:)
Q_sec_2(10:12,:)=diag([-1,-1,-1])*diag([-1,1,-1])*Q_sec_1(10:12,:)
Q_sec=[Q_sec_1,flip(Q_sec_2(:,1:end-1),2)]

L_sec=100*pi/2;
L_sec_1=L_sec;
L_sec_2=L_sec;
L_e=L_sec_1/N_e;

R_rod=5;
E=20e6
G=E/2*(1+0.4)
[Par_E]=[E*pi*R_rod^2, E*pi*R_rod^4/64 0.2*G*pi*R_rod^4/32]';
A=0.04e8;%can be to high or too low, depending on overall elatic energy level. the constrain parameter for 

% Q_sec=interp1(0:N_e,Q_sec',0:0.5:N_e)'
% 
% N_e=2*N_e
% L_e=L_e/2

% Q_sec=Q_sec_M(:,:,66); 
figure(1)
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


r_node_2=[1 -1 -1]'.*Q_sec(1:3,:);
dr_dx_2=[-1 1 1]'.*Q_sec(4:6,:);
dr_dy_2=[1 -1 -1]'.*Q_sec(7:9,:);
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





%plot its original configuration
N_node=3;

%%
%set the boundary conditional

%Using the de_demension matrix to define the constrain of the structure
Proj_M=eye(12*(N_e+1),12*(N_e+1));%dedemesionalize the variables

D_fix=0*Q_sec;%fixed dimensions/varibales
%BC1 Antisymmetric twisted folding

% r1 start point

% position:
% control: X
% Fixed: Y Z
D_fix(1,1)=1;%control the X disp
D_fix(2,1)=1;%fix the Y disp
D_fix(3,1)=1;%fix the  Z disp

% Rotation
% fix the local y axis
D_fix(7:9,1)=1;%fix the rotation

% r2
% position:
% Free: X
% Fixed: Y Z
D_fix(1,end)=1;%control the Y disp
D_fix(2,end)=1;%fix the X disp
D_fix(3,end)=1;%fix the  Z disp

% Rotation
% Fix: local y axis 
D_fix(7:9,end)=1;%fix the rotation

%r3
% position: 
N_m=N_e/2+1
% control Y disp: 
% fixed Z disp:
D_fix(2,N_m)=1;% Y disp in control
D_fix(3,N_m)=1;% Z disp fixed
B_index=find(D_fix==0);
%
Proj_M=Proj_M(:,B_index);

% displacement x

% Angle_y=linspace(0,90,19)*pi/180;
%define the start points and frames
P_m_0=Q_sec(1:3,N_m);
disp_x=linspace(0,50,101);
U=0*disp_x;
for mm=1:101
    mm
    if mm~=1
       Q_sec=Q_sec_M(:,:,mm-1);    
    end
    Q_sec(2,N_m)=P_m_0(2)-disp_x(mm);
    if ismember(mm,[1:10])
        Q_sec(1,N_m)=Q_sec(1,N_m)+2;
%         Q_sec(3,N_m)=Q_sec(3,N_m)+1;
    else
        Q_sec(1,N_m)=Q_sec(1,N_m);
    end
    %define the end points and frames
    kk=0;
    flag1=1;
    flag3=0;
    while (kk==0||flag1>0.00001||flag2>0.10)&&kk<4000
        kk=kk+1;
        [dFq_dq_sec,Fq_sec_M1(:,mm),q_node,U(mm)] = Jocob_rod_sec(Q_sec,N_e,L_e,N_node,Par_E,A);
        dFq_dq_sec_2=Proj_M'*dFq_dq_sec*Proj_M;
        dQ_sec_2=-pinv(dFq_dq_sec_2)*Proj_M'*Fq_sec_M1(:,mm);
        dQ_sec=Proj_M*dQ_sec_2;
        dQ_sec=reshape(dQ_sec,12,[]);
        Q_sec=Q_sec+0.5*dQ_sec;
        flag1=sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e+1);
        flag2=mean(abs(Proj_M'*Fq_sec_M1(:,mm)));
        if kk==4000
            flag3=1
        end
    end
    
    Angle_T_1(mm)=atan2d(Q_sec(9,1),Q_sec(7,1));
    Angle_T_2(mm)=atan2d(Q_sec(9,end),-Q_sec(8,end));
    Q_sec_M(:,:,mm)=Q_sec;
    
    % cla
    if 1==1 %mod(mm,1)==0
        figure(1)
        C_coord=[.0 .45 .74; .74 .45 .0];
        for ii=1:N_e
            r_node_1=q_node(1:3,:,ii);
            dr_dx_1=q_node(4:6,:,ii);
            dr_dy_1=q_node(7:9,:,ii);
            dr_dz_1=q_node(10:12,:,ii);
            L_co=2;
            figure(1)
            hold on
            plot3(r_node_1(1,:),r_node_1(2,:),r_node_1(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
            plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dy_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dy_1(2,:)], ...
                [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dy_1(3,:)],'-','color',C_coord(1,:),'linewidth',1)
            plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dz_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dz_1(2,:)], ...
                [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dz_1(3,:)],'-','color',C_coord(2,:),'linewidth',1)
            r_node_2=[1 -1 -1]'.*r_node_1;
            dr_dx_2=[-1 1 1]'.*dr_dx_1;
            dr_dy_2=[1 -1 -1]'.*dr_dy_1;
            dr_dz_2=[-1 1 1]'.*dr_dz_1;
            plot3(r_node_2(1,:),r_node_2(2,:),r_node_2(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
            plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dy_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dy_2(2,:)], ...
                [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dy_2(3,:)],'-','color',C_coord(1,:),'linewidth',1)
            plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dz_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dz_2(2,:)], ...
                [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dz_2(3,:)],'-','color',C_coord(2,:),'linewidth',1)

        end
        box on
        grid on
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        axis equal
        view([1 1 1])
        axis([-70 70 -100 100 -50 50])
        title(num2str(mm))
        pause(0.1)
    end

end


% save("Rod_Folding_Twisting_2.mat","Q_sec_M",'U_M')
% reshape(Q_sec_M([10 12],1,:),2,[])
% for ii=1:51
%     Angle_T_1(ii)=atan2d(Q_sec_M(10,1,ii),Q_sec_M(12,1,ii));
%     Angle_T_2(ii)=atan2d(-Q_sec_M(10,end,ii),Q_sec_M(12,end,ii));
% end
% Angle_T=Angle_T_1-Angle_T_2
% figure (2)
% hold on
% % plot3(disp_x(1:7)*0+Angle_T(1:7),disp_x(1:7),U(1:7),'-o')
% plot3(Angle_T,disp_x,U,'-o')
% grid on
% box on
% xlabel('TAngle')
% ylabel('disp_x')
% zlabel('U')
% 
% figure (3)
% hold on
% % plot3(disp_x(1:7)*0+Angle_T(1:7),disp_x(1:7),-Fq_sec(1,1:7),'-o')
% plot3(Angle_T,disp_x,-Fq_sec(1,:),'-o')
% grid on
% box on
% xlabel('TAngle')
% ylabel('disp_x')
% zlabel('F_x')



%%








