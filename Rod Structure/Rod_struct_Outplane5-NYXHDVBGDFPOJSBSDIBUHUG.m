% the script simulate the eye structure
% author: Qinghua Guan 
% Anyone modify or regenerate based on this code should attached this Authorization statement

%1： caculate the Bisymmetric configuration with different tendon lengths
%with a quarter of the ring
%2： analysis the stability of the Bisymmetric configuration by twist the
%ring circle as antisymmetric about the z-axis
%draw the map of the stability of Bisymmetric configurations
%V_5: 1st inplane deformation, 2nd out plane to Symmetrical, 3
%% start the simulation
clear
clc
close all
%
% input mechinical information
R0_ring=500;%radius of the ring mm
%the cross angle
Angle_cross=0*pi/180;
% the middle distance 
L_mid=0*R0_ring*0.3;


L_sec_0=1*R0_ring*pi/2;% set as the 1/4circle
R_rod=5;%radius of the rod
L_co=2*R_rod;%the length of the coordinates
E=2e4;%elastic modulus 10kPa
G=E/2*(1+0.4);%shearing modulus
Par_E=[E*pi*R_rod^2, E*pi*R_rod^4/4 0.2*G*pi*R_rod^4/2]';% the axis 
% Par_E=1.0e+07 *[ 1.5708; 0.6136; 0.1718]*1e2
% bending and torsion stiffness

% initialization of the caculation
N_e_0=10;% the element number
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

r1_p=[R0_ring, 0 ,0]';
q1=[r1_p;r1_frame(:)];

%define the end points and frames
r2_frame=[0 1 0; -1 0 0; 0 0 1]';
% Tune the start frame
alpha2=0*pi/180;%X axis
beta2=0*pi/180;%Y axis
theta2=90*pi/180;%Z axis

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


r2_p=[0 R0_ring,0]';

q2=[r2_p;r2_frame(:)];

%generate the original solution by interplotion
Q_sec=Curve_interp(q1,q2,L_sec_0,N_e_0);


C_coord=[.0 .45 .74; .74 .45 .0];%the color of the coordinates
Duplication_M=[1 1 1;-1 1 1;-1 -1 1;1 -1 1;];


hold on
for nn=1:(size(Q_sec,2)-1)
    q_node2=Curve_interp(Q_sec(:,nn),Q_sec(:,nn+1),L_e_0,N_node);
    Rod_ploting1(q_node2,Duplication_M,L_co,C_coord);
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
Proj_M=eye(12*(N_e_0+1),12*(N_e_0+1));%dedemesionalize the variables
% start point
D_fix=0*Q_sec;%fixed dimensions/varibales
D_fix(1,1)=0;%free the X disp
D_fix(2,1)=1;%fix the Y disp
D_fix(3,1)=1;%fix the  Z disp
% fix the local frame
D_fix(4:12,1)=1;%fix the rotation

%end point
D_fix(1,end)=1;%fix the X disp
D_fix(2,end)=0;%free the Y disp
D_fix(3,end)=1;%fix the Z disp
% Proj_M(12*N_e_0+2,12*N_e_0+1)=1;

% Fix local frame
D_fix([4 5 6],end)=1;%fix the rotation
D_fix([7 8 9],end)=1;%fix the rotation
D_fix([10 11 12],end)=1;%fix the rotation

B_index=find(D_fix==0);% find the fixed demensions
Proj_M=Proj_M(:,B_index);% extract the projection matrix of free paramters


N_x=0;
% L_co=2;%the length of the coordinates
C_coord=[.0 .45 .74; .74 .45 .0];%the color of the coordinates


Damp=0.5;%reduce the iteration speed
mm=1;
kk=0;
flag1=1;
flag3=0;
while (kk==0||flag1>0.000001*R0_ring||flag2>0.10)&&kk<2000
    kk=kk+1;
    [dFq_dq_sec2,Fq_sec_M(:,mm),q_node2,U_Cal_sec(mm)] =Jocob_rod_sec(Q_sec,N_e_0,L_e_0,N_node,Par_E,A);
    dFq_dq_sec_2=Proj_M'*dFq_dq_sec2*Proj_M;
    dQ_sec_1=-pinv(dFq_dq_sec_2)*Proj_M'*Fq_sec_M(:,mm);
    dQ_sec=Proj_M*dQ_sec_1;
    dQ_sec=reshape(dQ_sec,12,[]);
    Q_sec=Q_sec+Damp*dQ_sec;
    Q_sec_M(:,:,mm)=Q_sec;
    flag1=sqrt(sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e_0+1))
    flag2=0;%mean(abs(Proj_M'*Fq_sec_M(:,mm)));%
    % AA=reshape(Fq_sec_M(:,mm),12,[]);
    if kk==2000
        flag3=1;
    end
end

Duplication_M=[1 1 1;-1 1 1;-1 -1 1;1 -1 1;];
figure(2);
hold on
Rod_ploting2(Q_sec,Duplication_M,L_co*5,C_coord)

for nn=1:(size(Q_sec,2)-1)
    q_node2=Curve_interp(Q_sec(:,nn),Q_sec(:,nn+1),L_e_0,N_node); 
    Rod_ploting1(q_node2,Duplication_M,L_co,C_coord);
end


N_e_1=N_e_0;
L_sec_1=L_sec_0;
N_e_1=N_e_0;% the element number
L_e_1=L_sec_1/N_e_1;
N_e_0=N_e_1/2
L_sec_0=L_sec_1/2;


P_wheel=Q_sec(1:3,N_e_0+1);
RM_wheel=reshape(Q_sec(4:12,N_e_0+1),3,3);
R_wheel=120;% The radius of the wheel
N_spoke=12;% The spoke num of the wheel
Wheel_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M)

hold off
box on
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
% axis([-600 600 -600 600 -600 600])
view([1 1 1])

%% caculate the inplane deformation
load('Rod_Inplane5.mat',"Q_sec_M1")
Q_sec=Q_sec_M1(:,:,1);
% Build the boundary conditions for the boundary
%Using the de_demension matrix to define the constrain of the structure
Proj_M=eye(12*(N_e_1+1),12*(N_e_1+1));%dedemesionalize the variables

%rod 1
D_fix=0*Q_sec;
% start point

% constrain the start point
%  constrain  X Y Z disp
D_fix(1,1)=1;%control  X disp
D_fix(2,1)=1;%fix the Y disp
D_fix(3,1)=1;%fix the  Z disp
%  constrain/fix the local frame
D_fix(4:12,1)=1;%fix the rotation

%end point

% constrain the end point
%  constrain  X Y Z disp 
D_fix(1,end)=1;%fix the X disp
D_fix(2,end)=0;%free Y disp
D_fix(3,end)=1;%fix the Z disp
%  constrain/Fix local frame
D_fix(4:12,end)=1;%fix the rotation

B_index=find(D_fix==0);% find the fixed demensions
Proj_M=Proj_M(:,B_index);% extract the projection matrix of free paramters


X0=Q_sec(1,1);
disp_x=[0:5:X0,X0];% sampling displacement x
N_x=length(disp_x)-1
% L_co=20;%the length of the coordinates
C_coord=[.0 .45 .74; .74 .45 .0];%the color of the coordinates

Duplication_M=[1 1 1;-1 1 1;-1 -1 1;1 -1 1;];
Damp=0.5;%reduce the iteration speed
pic_num=1
for mm=1:N_x+1
    mm
    if mm~=1
        Q_sec=Q_sec_M1(:,:,mm-1);        
    end  
    Q_sec(1,1)=X0-disp_x(mm);    
    
    %define the end points and frames
    kk=0;
    flag1=1;
    flag3=0;
    while (kk<=100||flag1>0.000001*R0_ring||flag2>0.10)&&kk<2000
        kk=kk+1;
        [dFq_dq_sec1,Fq_sec_M1(:,mm),q_node1(:,:,:,mm),U_Cal_sec1(mm),U_v_sec1(:,mm),U_E_ele1(mm,:),U_E_node1(:,:,mm)] =Jocob_rod_sec(Q_sec,N_e_1,L_e_1,N_node,Par_E,A);
        F_tendon_1(mm)=Fq_sec_M1(1:2,mm);
        dFq_dq_sec_2=Proj_M'*dFq_dq_sec1*Proj_M;
        Fq_sec_M1_2=Proj_M'*Fq_sec_M1(:,mm);
        dQ_sec_1=-pinv(dFq_dq_sec_2)*Proj_M'*Fq_sec_M1(:,mm);
        dQ_sec=Proj_M*dQ_sec_1;
        dQ_sec=reshape(dQ_sec,12,[]);        
        Q_sec=Q_sec+Damp*dQ_sec;
        Q_sec_M1(:,:,mm)=Q_sec;        
        flag1=sqrt(sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e_0+1));
%         fprintf('kk=%u, flag1=%f\n',kk,flag1)
        flag2=0;%mean(abs(Proj_M'*Fq_sec_M1(:,mm)));
        if kk==2000
            flag3=1;
        end
        %         if mod(mm,5)==1||mod(kk,20)==1

        %         end
    end

    figure(2);
    clf
    hold on
    Rod_ploting1(q_node1(:,:,:,mm),Duplication_M,L_co,C_coord,U_E_ele1(mm,:),U_E_node1(:,:,mm));    
    %%plotting wheels
    P_wheel=Q_sec(1:3,N_e_0+1);
    RM_wheel=reshape(Q_sec(4:12,N_e_0+1),3,3);
    Wheel_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M)
    plot3(Q_sec(1,1)*[1,-1],Q_sec(2,1)*[1,1],Q_sec(3,1)*[1,1],'r-')
    plot3(Q_sec(1,end)*[1,1],Q_sec(2,end)*[1,-1],Q_sec(3,end)*[1,1],'g-')
    colorbar     
    hold off
    box on
    grid on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    axis equal
    axis([-600 600 -700 700 -150 150])
    view([1 1 1])
    title(['X_1=',num2str((X0-disp_x(mm)),'%.0f'),', Y_2=',num2str(Q_sec_M1(2,N_e_1+1,mm),'%.0f')])
    drawnow;% 立即刷新当前绘图窗口，这是matlab绘图中动态展示的关键    
    pause(0.1)
end

for mm=1:N_x+1
    Q_sec=Q_sec_M1(:,:,mm)
    figure(2);
    clf
    hold on
    Rod_ploting1(q_node1(:,:,:,mm),Duplication_M,L_co,C_coord,U_E_ele1(mm,:),U_E_node1(:,:,mm));    
    %%plotting wheels
    P_wheel=Q_sec(1:3,N_e_0+1);
    RM_wheel=reshape(Q_sec(4:12,N_e_0+1),3,3);
    Wheel_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M)
    plot3(Q_sec(1,1)*[1,-1],Q_sec(2,1)*[1,1],Q_sec(3,1)*[1,1],'r-')
    plot3(Q_sec(1,end)*[1,1],Q_sec(2,end)*[1,-1],Q_sec(3,end)*[1,1],'g-')
    colorbar   
    if mm==1
        caxis([1.92e8 2.0e8])
    end
    hold off
    box on
    grid on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    axis equal
    axis([-600 600 -700 700 -150 150])
    view([1 1 1])
    title(['X_1=',num2str((X0-disp_x(mm)),'%.0f'),', Y_2=',num2str(Q_sec_M1(2,N_e_1+1,mm),'%.0f')])
    drawnow;% 立即刷新当前绘图窗口，这是matlab绘图中动态展示的关键
    F = getframe(gcf);  % 获取当前绘图窗口的图片
    Im = frame2im(F);   % 返回与电影帧相关的图片数据
    [Amap, map] = rgb2ind(Im, 256); % 将RGB图片转化为索引图片
    if pic_num == 1
        imwrite(Amap, map, 'Rod_inplane5.gif', 'gif', 'Loopcount', Inf, 'DelayTime', 0.1);
        % 将第一张图片写入‘sin.gif’文件中，并且将GIF播放次数设置成无穷，即保存的GIF图会一直动下去
    else
        imwrite(Amap, map,'Rod_inplane5.gif','gif','WriteMode','append','DelayTime',0.1);
        % 依次将其他的图片写入到GIF文件当中
    end
    pic_num = pic_num + 1;
    pause(0.1)
end


figure(3);
disp_xy_B=[reshape(Q_sec_M1(1,1,:),1,[]);reshape(Q_sec_M1(2,N_e_1+1,:),1,[])];
subplot(2,2,1)
plot(disp_xy_B(1,:),disp_xy_B(2,:))
hold on
plot(disp_xy_B(2,:),disp_xy_B(1,:))
axis equal
hold off
box on
grid on
xlabel('X')
ylabel('Y')

F_tendon_1=Fq_sec_M1(1:2,:);
subplot(2,2,3)
plot(disp_x,-F_tendon_1(1,:))
subplot(2,2,4)
plot(disp_x,-F_tendon_1(2,:))
save('Rod_Inplane5.mat',"Q_sec_M1",'q_node1','F_tendon_1','U_v_sec1','U_E_ele1','U_E_node1')

%% caculate bisymmetric configurations in the range inside the boundary
load('Rod_Inplane5.mat',"Q_sec_M1")
% load('Rod_BiSym_Outplane5.mat')
A=1e9;
Damp=0.1;%reduce the iteration speed
% Build the boundary conditions for the boundary
% start point x_disp fixed and end point y_disp in control
%Using the de_demension matrix to define the constrain of the structure
Proj_M=eye(12*(N_e_1+1),12*(N_e_1+1));%dedemesionalize the variables

%rod 1
D_fix=0*Q_sec;
% start point

% constrain the start-middle point of rd1 and rd2
%  constrain  X Y Z disp
D_fix(1,1)=1;%fix  the X disp
D_fix(2,1)=1;%fix the Y disp
D_fix(3,1)=0;%free the  Z disp
%  constrain/fix the local frame
D_fix(4:6,1)=1;%fix the X axis
D_fix(7:12,1)=0;%free the X ratation of Y Z axis
D_fix([8 11],1)=1;%free the X ratation of Y Z axis

%end point

% constrain the end-middle point of rd1 and rd2
%  constrain  X Y Z disp 
D_fix(1,end)=1;%fix the X disp
D_fix(2,end)=1;%control the Y disp
D_fix(3,end)=1;%fix the Z disp
%  constrain/Fix local frame
D_fix(4:12,end)=1;%fix the X axis
% D_fix(7:12,end)=0;%free the X ratation of Y Z axis
% D_fix([7 10],end)=1;%free the X ratation of Y Z axis



B_index=D_fix==0;% find the fixed demensions
Proj_M=Proj_M(:,B_index);% extract the projection matrix of free paramters


% since the symmetry of the structure, just a quater of the structure need to be caculate

% U=[];
P_wheel=[]


for mm=51%N_x+1  
    figure(4);
    clf
    Q_sec=Q_sec_M1(:,:,mm);    
    % caculate the range of disp y based on the boundary 
    disp_x_0=Q_sec(1,1);
    disp_y_0=Q_sec(2,end);    
    disp_y=flip([disp_x_0:5:disp_y_0, disp_y_0]);% sampling displacement y

    for nn=1:80%length(disp_y)
        mm
        nn
        Damp_1=0.01
        Damp_2=0.1
        if ismember(nn,1:3)
            Damp_1=0.01
            Damp_2=0.05
        elseif ismember(nn,4:94)
            Damp_1=0.005
            Damp_2=0.05
        elseif ismember(nn,95:102)
            Damp_1=0.001
            Damp_2=0.02
%         elseif ismember(nn,37:83)
%             Damp_1=0.0001
%             Damp_2=0.01
% %         elseif ismember(nn,40:63)
% %             Damp_1=0.0001
% %             Damp_2=0.05
%         elseif ismember(nn,84:length(disp_y))
%             Damp_1=0.0001
%             Damp_2=0.025
        end
        if nn~=1
            Q_sec=Q_sec_M2{mm}(:,:,nn-1);
        end
%         Q_sec(1,1)=disp_y(nn);        
        Q_sec(2,N_e_1+1)=disp_y(nn);     
        
        if ismember(nn,1:30)
            Q_sec(3,1)=Q_sec(3,1)+0.2*(15-0.5*nn);% add the z axis distorion
        end
        

        %define the end points and frames
        kk=0;
        flag1=1;
        flag2=1;
        Damp=Damp_1
        while (kk<=100||flag1>0.000001*R0_ring)&&flag2&&(kk<5/Damp||kk<3000)
            kk=kk+1;
            Damp=0.995^(kk)*Damp_1+(1-0.995^(kk))*Damp_2;
            [dFq_dq_sec2,Fq_sec_M2(:,nn),q_node2{mm}(:,:,:,nn),U_Cal_sec2{mm}(nn),U_v_sec2{mm}(:,nn),U_E_ele2{mm}(nn,:),U_E_node2{mm}(:,:,nn)] =Jocob_rod_sec(Q_sec,N_e_1,L_e_1,N_node,Par_E,A);
            U_M2_flag(kk)=U_Cal_sec2{mm}(nn);            
            F_tendon_2{mm}(:,nn)=Fq_sec_M2(1:2,nn);
            dFq_dq_sec_2=Proj_M'*dFq_dq_sec2*Proj_M;
            Fq_sec_M1_2=Proj_M'*Fq_sec_M2(:,nn);
            dQ_sec_1=-pinv(dFq_dq_sec_2)*Proj_M'*Fq_sec_M2(:,nn);
            dQ_sec=Proj_M*dQ_sec_1;
            dQ_sec=reshape(dQ_sec,12,[]);
            %the distance between two mid points of rod1 and rod2
            Q_sec=Q_sec+Damp*dQ_sec;    
            P_wheel(:,1)=Q_sec(1:3,N_e_0+1);            
            Q_sec_M2{mm}(:,:,nn)=Q_sec;
            flag1=sqrt(sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e_0+1));
            flag2=1;
            if kk<=10/Damp_2
                flag2=1;
            elseif (max(U_M2_flag((kk-10/Damp_2):kk))-min(U_M2_flag((kk-10/Damp_2):kk)))<1e7;%energy variation
                flag2=0;
            end

            if kk==2000
                flag3=1;
            end
            if mod(kk,20)==1
                fprintf('kk=%u, flag1=%f \n',kk,flag1)
            end
            if mod(kk,100)==2
                figure(3);
                clf
                hold on
                Rod_ploting1(q_node2{mm}(:,:,:,nn),Duplication_M,L_co,C_coord,U_E_ele2{mm}(nn,:),U_E_node2{mm}(:,:,nn))
                plot3(Q_sec(1,1)*[1,-1],Q_sec(2,1)*[1,1],Q_sec(3,1)*[1,1],'r-')
                plot3(Q_sec(1,end)*[1,1],Q_sec(2,end)*[1,-1],Q_sec(3,end)*[1,1],'g-')
                %%plotting wheels
                P_wheel=Q_sec(1:3,N_e_0+1);
                RM_wheel=reshape(Q_sec(4:12,N_e_0+1),3,3);
                Wheel_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M)
                colorbar
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
                    ', X_2=', num2str(Q_sec(2,end),'%.0f'), ...
                    ', Z_1=', num2str(Q_sec(3,1),'%.0f')])
                figure(4);
                subplot(2,2,1)
                plot(U_M2_flag(1:kk))
                subplot(2,2,2)
                if kk>10/Damp_2
                    plot(U_M2_flag(kk-10/Damp_2:kk)-U_M2_flag(kk-10/Damp_2))
                end
                pause(0.01)
            end
        end
        save('Rod_Outplane5.mat',"Q_sec_M2",'q_node2',"U_Cal_sec2","F_tendon_2",'U_E_ele2','U_E_node2')
        X_1=reshape(Q_sec_M2{mm}(1,1,1:nn),1,[]);
        Y_2=reshape(Q_sec_M2{mm}(2,end,1:nn),1,[]);
        Z_1=reshape(Q_sec_M2{mm}(3,1,1:nn),1,[]);
        
        figure(4);
%         hold on
        subplot(2,1,2)
        plot3(X_1,Y_2,Z_1,'-o')
        box on
        grid on
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        view([1 0 0])
        
%         hold off
    end

end

mm=51
X_1=reshape(Q_sec_M2{mm}(1,1,1:nn),1,[]);
Y_2=reshape(Q_sec_M2{mm}(2,end,1:nn),1,[]);
Z_1=reshape(Q_sec_M2{mm}(3,1,1:nn),1,[]);
U_Cal_sec=U_Cal_sec2{mm};
figure(5);
plotyy(Y_2(1)-Y_2,U_Cal_sec,Y_2(1)-Y_2,Z_1)
box on
grid on
xlabel('X')
ylabel('Y')
figure(6)
plot(Y_2(1)-Y_2,F_tendon_2{mm}(1,:))

% save('Eye_BiSym_Outplane1.mat',"Q_sec_M2",'-append')
%%
load('Rod_Outplane5.mat')
pic_num=1
for nn=1:length(disp_y)
    Q_sec=Q_sec_M2{mm}(:,:,nn);

    figure(3);
    clf
    hold on
    Rod_ploting1(q_node2{mm}(:,:,:,nn),Duplication_M,L_co,C_coord,U_E_ele2{mm}(nn,:),U_E_node2{mm}(:,:,nn))
    plot3(Q_sec(1,1)*[1,-1],Q_sec(2,1)*[1,1],Q_sec(3,1)*[1,1],'r-')
    plot3(Q_sec(1,end)*[1,1],Q_sec(2,end)*[1,-1],Q_sec(3,end)*[1,1],'g-')
    %%plotting wheels
    P_wheel=Q_sec(1:3,N_e_0+1);
    RM_wheel=reshape(Q_sec(4:12,N_e_0+1),3,3);
    Wheel_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M)
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

    title(['X_1=',num2str(Q_sec(1,1),'%.0f'), ...
        ', X_2=', num2str(Q_sec(2,N_e_1+1),'%.0f'), ...
        ', Z_1=', num2str(Q_sec(3,1),'%.0f')])
    drawnow;% 立即刷新当前绘图窗口，这是matlab绘图中动态展示的关键
    F = getframe(gcf);  % 获取当前绘图窗口的图片
    Im = frame2im(F);   % 返回与电影帧相关的图片数据
    [Amap, map] = rgb2ind(Im, 256); % 将RGB图片转化为索引图片
    if pic_num == 1
        imwrite(Amap, map, 'Rod_Outplane5.gif', 'gif', 'Loopcount', Inf, 'DelayTime', 0.1);
        % 将第一张图片写入‘sin.gif’文件中，并且将GIF播放次数设置成无穷，即保存的GIF图会一直动下去
    else
        imwrite(Amap, map,'Rod_Outplane5.gif','gif','WriteMode','append','DelayTime',0.1);
        % 依次将其他的图片写入到GIF文件当中
    end
    pic_num = pic_num + 1;
    pause(0.1)

end

%%

for mm=1:51
    figure(4);
    clf

nn=size(Q_sec_M2{mm}(:,:,:),3);
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

mm=41
nn=5

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



% the local x y z axis in control
% because the equalism relationship
%D_fix(1:3,1)=[-1 -1 1].*D_fix(4:6,end);%x axis in control 
%D_fix(4:6,1)=[-1 -1 1].*D_fix(4:6,end);%x axis in control 
%D_fix(7:9,1)=[-1 -1 1].*D_fix(4:6,end);%x axis in control 
%D_fix(10:12,1)=[-1 -1 1].*D_fix(4:6,end);%x axis in control 
% sus
Con_M_1=diag(repmat([-1 -1 1],1,4))
% Proj_M_0(1:12,end-11:end)=0
Proj_M_0(end-11:end,1:12)=Con_M_1

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
Con_M_2_2=[-sin(theta_m_0)/r_m_0 cos(theta_m_0)/r_m_0; cos(theta_m_0) sin(theta_m_0)];
Proj_M_0((n_m-1)*12+[1:2],(n_m-1)*12+[1:2])=Con_M_2_2;
%fix R: the distance to the z axis
D_fix(2,n_m)=1;%fix the R
D_fix(3,n_m)=1;%fix the Z
%the angle between local x axis and the plane X-Y in control
D_fix(6,n_m)=1;
%local y axis in the plane X-Y
D_fix(9,n_m)=1;

T_angle=linspace(0,30*pi/180,31);

% Q_sec(:,1)
% Q_sec(:,end);
figure(4);
    clf
for ii=1:31
    ii    
    figure(3);
    clf
    if ii~=1
        Q_sec=Q_sec_M3{mm}(:,:,ii-1,nn);
    end
    %redefine the middle points and frames    
    Q_sec(6,n_m)=sin(T_angle(ii));  

    kk=0;
    flag1=1;
    L_step=0.1;
   
    while (kk==0||flag1>0.0001*L_step||kk<2/L_step)&&kk<800
        kk=kk+1;
        [dFq_dq_sec2,Fq_sec(:,ii),q_node2,U(ii),E_elastic{mm}(:,ii,nn)] = Jocob_rod_sec(Q_sec,N_e_1,L_e_1,N_node,Par_E,A);
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
            Con_M_2_2=[-sin(theta_m)/r_m cos(theta_m)/r_m; cos(theta_m) sin(theta_m)];
            Proj_M_0((n_m-1)*12+[1:2],(n_m-1)*12+[1:2])=Con_M_2_2;
        end
        B_index=find(D_fix==0);% find the fixed demensions
        Proj_M=Proj_M_0(:,B_index);% extract the projection matrix of free paramters
        dFq_dq_sec_2=Proj_M'*dFq_dq_sec2*Proj_M;
        dQ_sec_1=-pinv(dFq_dq_sec_2)*Proj_M'*Fq_sec(:,ii);
        dQ_sec=L_step*Proj_M*dQ_sec_1;
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
                r_node_1=q_node2(1:3,:,jj);
                dr_dx_1=q_node2(4:6,:,jj);
                dr_dy_1=q_node2(7:9,:,jj);
                dr_dz_1=q_node2(10:12,:,jj);

                hold on
                plot3(r_node_1(1,:),r_node_1(2,:),r_node_1(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
                plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dy_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dy_1(2,:)], ...
                    [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dy_1(3,:)],'-','color',C_coord(1,:),'linewidth',1)
                plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dz_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dz_1(2,:)], ...
                    [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dz_1(3,:)],'-','color',C_coord(2,:),'linewidth',1)

                r_node_2=[-1 -1 1]'.*q_node2(1:3,:,jj);
                dr_dx_2=[-1 -1 1]'.*q_node2(4:6,:,jj);
                dr_dy_2=[-1 -1 1]'.*q_node2(7:9,:,jj);
                dr_dz_2=[-1 -1 1]'.*q_node2(10:12,:,jj);

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

        r_node_1=q_node2(1:3,:,jj);
        dr_dx_1=q_node2(4:6,:,jj);
        dr_dy_1=q_node2(7:9,:,jj);
        dr_dz_1=q_node2(10:12,:,jj);

        plot3(r_node_1(1,:),r_node_1(2,:),r_node_1(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
        plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dy_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dy_1(2,:)], ...
            [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dy_1(3,:)],'-','color',C_coord(1,:),'linewidth',1)
        plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dz_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dz_1(2,:)], ...
            [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dz_1(3,:)],'-','color',C_coord(2,:),'linewidth',1)

        r_node_2=[-1 -1 1]'.*q_node2(1:3,:,jj);
        dr_dx_2=[-1 -1 1]'.*q_node2(4:6,:,jj);
        dr_dy_2=[-1 -1 1]'.*q_node2(7:9,:,jj);
        dr_dz_2=[-1 -1 1]'.*q_node2(10:12,:,jj);

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
jj=31
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

Q_sec=Q_sec_M3{mm}(:,:,11,nn)
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

L_sec_2=L_sec_1*2;
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
axis([-70 70 -70 70 -50 50])
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



%2/4 point
D_fix(1,n_2q)=1;%fix the X disp
D_fix(2,n_2q)=1;%control the Y disp
D_fix(3,n_2q)=1;%control the Z disp


% the 1/4 and 3/4 point position
% based on constrain of the length
% the distance of 1/4 point to the 3/4 axis is constrained: x^2+y^2=C
% introduce the R and theta to replace the x y disp of 3/4 point
V_12q=Q_sec(1:3,n_3q)-Q_sec(1:3,n_1q);
r_m_0=sqrt(sum(V_12q.^2));
theta_m_0=atan2(V_12q(2),V_12q(1));
% the present First-order expansion relationship between d(theta,R,x1,y1) and d(x3,y3)
% x3=x1+cos(theta)*R y3=y1+sin(theta)*R
% z3=z1
%dx3=dx1-sin(theta)*R*dtheta+cos(theta)*dR
%dy3=dy1+cos(theta)*R*dtheta+sin(theta)*dR
%dz3=dz1
Con_M_2_1=eye(2);
Con_M_2_2=[-sin(theta_m_0)/r_m_0 cos(theta_m_0)/r_m_0; cos(theta_m_0) sin(theta_m_0)];
Con_M_2_3=1;
Proj_M_0((n_3q-1)*12+[1:2],(n_1q-1)*12+[1:2])=Con_M_2_1;
Proj_M_0((n_3q-1)*12+[1:2],(n_3q-1)*12+[1:2])=Con_M_2_2;
Proj_M_0((n_3q-1)*12+3,(n_1q-1)*12+3)=Con_M_2_3;
%fix R: the distance to the z axis
D_fix(2,n_3q)=1;%fix the R
D_fix(3,n_3q)=1;%fix the Z control by z1

D_fix(1,n_1q)=1;%x disp in control 

%the angle between local x axis and the plane X-Y in control
D_fix(6,n_1q)=1;
D_fix(9,n_1q)=1;
%local y axis in the plane X-Y
D_fix(6,n_3q)=1;
D_fix(9,n_3q)=1;

Disp_x_1q_0=Q_sec(1,n_1q);
%%

Disp_y_1q=linspace(0,5,11);

% Q_sec(:,1)
% Q_sec(:,end);
% figure(11);
%     clf
    figure(10);
    clf
for ii=1
    ii    
    figure(9);
    clf
    if ii~=1
        Q_sec=Q_sec_M4{mm}(:,:,ii-1,ll,nn);
    end
    %redefine the middle points and frames    
    Q_sec(1,n_1q)=Disp_x_1q_0+Disp_y_1q(ii);  

    kk=0;
    flag1=1;
    L_step=0.01;
   
    while (kk==0||flag1>0.0001*L_step||kk<2/L_step)&&kk<800
        kk=kk+1;
        [dFq_dq_sec2,Fq_sec(:,ii),q_node2,U(ii),E_elastic_4{mm}(:,ii,ll,nn)] = Jocob_rod_sec(Q_sec,N_e_2,L_e_1,N_node,Par_E,A);
        % update the boundary information
        if r_m_0<10e-4
            r_m=0;
            theta_m=theta_m_0;
            D_fix(1:2,n_3q)=1;
        else
            D_fix(1,n_3q)=0;
            V_12q=Q_sec(1:3,n_3q)-Q_sec(1:3,n_1q);
            r_m=sqrt(sum(V_12q.^2));
            theta_m=atan2(V_12q(2),V_12q(1));
            % the present First-order expansion relationship between d(theta,R) and d(x,y)
            Con_M_2_2=[-sin(theta_m)/r_m cos(theta_m)/r_m; cos(theta_m) sin(theta_m)];
            Proj_M_0((n_3q-1)*12+[1:2],(n_3q-1)*12+[1:2])=Con_M_2_2;
        end
        B_index=find(D_fix==0);% find the fixed demensions
        Proj_M=Proj_M_0(:,B_index);% extract the projection matrix of free paramters
        dFq_dq_sec_2=Proj_M'*dFq_dq_sec2*Proj_M;
        dQ_sec_1=-pinv(dFq_dq_sec_2)*Proj_M'*Fq_sec(:,ii);
        dQ_sec=L_step*Proj_M*dQ_sec_1;
        dQ_sec=reshape(dQ_sec,12,[]);
        Q_sec=Q_sec+dQ_sec;
        % canonicalize the diplacement based on the varition of angle
        if r_m_0<10e-4
            Q_sec(1:2,n_3q)=[theta_m_0,0]';
        else
            d_theta=Con_M_2_2(1,:)*(dQ_sec([1:2],n_3q)-dQ_sec([1:2],n_1q));
            theta_m=theta_m+d_theta;
            Q_sec(1:2,n_3q)=Q_sec([1:2],n_1q)+r_m_0*[cos(theta_m);sin(theta_m)];
        end
        Q_sec_M4{mm}(:,:,ii,ll,nn)=Q_sec;
        flag1=sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e_2+1);
        flag2=mean(abs(Proj_M'*Fq_sec(:,ii)));
        if kk==5000
            flag3=1
        elseif mod(kk,20)==0
            kk
            flag1
            figure(9);
            hold on
            for jj=1:N_e_2
                r_node_1=q_node2(1:3,:,jj);
                dr_dx_1=q_node2(4:6,:,jj);
                dr_dy_1=q_node2(7:9,:,jj);
                dr_dz_1=q_node2(10:12,:,jj);

                
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

        r_node_1=q_node2(1:3,:,jj);
        dr_dx_1=q_node2(4:6,:,jj);
        dr_dy_1=q_node2(7:9,:,jj);
        dr_dz_1=q_node2(10:12,:,jj);

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

