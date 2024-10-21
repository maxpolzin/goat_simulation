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
    q_node3=Curve_interp(Q_sec(:,nn),Q_sec(:,nn+1),L_e_0,N_node);
    Rod_ploting1(q_node3,Duplication_M,L_co,C_coord);
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
    [dFq_dq_sec3,Fq_sec_M(:,mm),q_node3,U_Cal_sec(mm)] =Jocob_rod_sec(Q_sec,N_e_0,L_e_0,N_node,Par_E,A);
    dFq_dq_sec_2=Proj_M'*dFq_dq_sec3*Proj_M;
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
    q_node3=Curve_interp(Q_sec(:,nn),Q_sec(:,nn+1),L_e_0,N_node); 
    Rod_ploting1(q_node3,Duplication_M,L_co,C_coord);
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

%% caculate the inplane deformation when contract tendon1
load('Rod_Inplane5.mat')
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


%%
for mm=1:N_x+1
    mm
    if mm~=1
        Q_sec=Q_sec_M1(:,:,mm-1);   
        Q_sec=Q_sec_M1(:,:,mm); %when the result has been caculated
    end  
    Q_sec(1,1)=X0-disp_x(mm);    
    
    %define the end points and frames
    kk=0;
    flag1=1;
    flag3=0;
    while (kk<=100||flag1>0.000001*R0_ring||flag2>0.10)&&kk<2000
        kk=kk+1;
        [dFq_dq_sec1,Fq_sec_M1(:,mm),q_node1(:,:,:,mm),U_Cal_sec1(mm),U_v_sec1(:,mm),U_E_ele1(mm,:),Dens_E_node1(:,:,mm)] =Jocob_rod_sec(Q_sec,N_e_1,L_e_1,N_node,Par_E,A);
        F_tendon1(1:2,mm)=Fq_sec_M1(1:2,mm);
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
    Rod_ploting1(q_node1(:,:,:,mm),Duplication_M,L_co,C_coord,U_E_ele1(mm,:),Dens_E_node1(:,:,mm));    
    %%plotting wheels
    P_wheel=Q_sec(1:3,N_e_0+1);
    RM_wheel=reshape(Q_sec(4:12,N_e_0+1),3,3);
    Wheel_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M)
    plot3(Q_sec(1,1)*[1,-1],Q_sec(2,1)*[1,1],Q_sec(3,1)*[1,1],'r-')
    plot3(Q_sec(1,end)*[1,1],Q_sec(2,end)*[1,-1],Q_sec(3,end)*[1,1],'g-')
    c=colorbar     
    c.Label.String = 'Energy Density (mJ/mm)';
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

save('Rod_Inplane5.mat',"Q_sec_M1",'q_node1','F_tendon1','U_v_sec1','U_E_ele1','Dens_E_node1')



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

F_tendon1=Fq_sec_M1(1:2,:);
subplot(2,2,3)
plot(disp_x,-F_tendon1(1,:))
subplot(2,2,4)
plot(disp_x,-F_tendon1(2,:))

figure(3)
[haxes,hline1,hline2]=plotyy(disp_x(1:51),-F_tendon1(1,1:51),disp_x(1:51),F_tendon1(2,1:51))

set(hline1,'LineWidth',2);
set(hline2,'LineWidth',2);
% set (gca,'XDir','reverse');
xlabel('X disp')
ylabel(haxes(1), 'F_{T1}','linewidth',2);
ylabel(haxes(2), 'F_{T2}','linewidth',2);
set(haxes,'FontSize',16,'fontname','Time','linewidth',1)
set(gcf, 'Position', [100 100 540 400]);

save('Rod_Inplane5.mat',"Q_sec_M1",'q_node1','F_tendon1','U_v_sec1','U_E_ele1','Dens_E_node1')

%%
v = VideoWriter('Rod_inplane5_2.mp4', 'MPEG-4');
v.FrameRate = 30;  % Set the frame rate (frames per second)
v.Quality = 90;    % Set the video quality (0 to 100)
open(v);
% writeVideo(v, F);
% close(v)
pic_num=1
for mm=1:51%N_x+1
    Q_sec=Q_sec_M1(:,:,mm);
    figure(2);
    clf
    hold on
    Rod_ploting1(q_node1(:,:,:,mm),Duplication_M,L_co,C_coord,U_E_ele1(mm,:),Dens_E_node1(:,:,mm));    
    %%plotting wheels
    P_wheel=Q_sec(1:3,N_e_0+1);
    RM_wheel=reshape(Q_sec(4:12,N_e_0+1),3,3);
    Wheel_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M)
    plot3(Q_sec(1,1)*[1,-1],Q_sec(2,1)*[1,1],Q_sec(3,1)*[1,1],'r-','linewidth',2)
    plot3(Q_sec(1,end)*[1,1],Q_sec(2,end)*[1,-1],Q_sec(3,end)*[1,1],'g-','linewidth',2)
    c=colorbar;   
    if mm==1
        caxis([1.95e1 2.2e1])
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
    title(['D_1=',num2str(2*(X0-disp_x(mm)),'%.0f'),', D_2=',num2str(2*Q_sec_M1(2,N_e_1+1,mm),'%.0f')])
    set(gca,'FontSize',16,'fontname','Times')
    drawnow;% 立即刷新当前绘图窗口，这是matlab绘图中动态展示的关键
    F = getframe(gcf);  % 获取当前绘图窗口的图片
    Im = frame2im(F);   % 返回与电影帧相关的图片数据
    [Amap, map] = rgb2ind(Im, 256); % 将RGB图片转化为索引图片
    if pic_num == 1
        imwrite(Amap, map, 'Rod_inplane5_2.gif', 'gif', 'Loopcount', Inf, 'DelayTime', 0.1);
        % 将第一张图片写入‘sin.gif’文件中，并且将GIF播放次数设置成无穷，即保存的GIF图会一直动下去
    else
        imwrite(Amap, map,'Rod_inplane5_2.gif','gif','WriteMode','append','DelayTime',0.1);
        % 依次将其他的图片写入到GIF文件当中
    end
    pic_num = pic_num + 1
    writeVideo(v, F);
    pause(0.1)
end
close(v)



%% caculate outplane configurations when contract tendon2
load('Rod_Inplane5.mat')
load('Rod_Outplane5.mat')
A=1e3;
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

mm=51%N_x+1  
    
    Q_sec=Q_sec_M1(:,:,mm);    
    % caculate the range of disp y based on the boundary 
    disp_x_0=Q_sec(1,1);
    disp_y_0=Q_sec(2,end);    
    disp_y=flip([disp_x_0:5:disp_y_0, disp_y_0]);% sampling displacement y
%%
figure(4);
    clf
    for nn=1:80%length(disp_y)
        mm
        nn
        Damp_1=0.01
        Damp_2=0.1
        if ismember(nn,1:3)
            Damp_1=0.005
            Damp_2=0.05
        elseif ismember(nn,4:80)
            Damp_1=0.004
            Damp_2=0.04
       
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
            Q_sec(3,1)=Q_sec(3,1)-0.2*(15-0.5*nn);% add the z axis distorion
        end
        

        %define the end points and frames
        kk=0;
        flag1=1;
        flag2=1;
        Damp=Damp_1
        while (kk<=100||flag1>0.000001*R0_ring)&&flag2&&(kk<5/Damp||kk<3000)
            kk=kk+1;
            Damp=0.995^(kk)*Damp_1+(1-0.995^(kk))*Damp_2;
            [dFq_dq_sec2,Fq_sec_M2(:,nn),q_node2{mm}(:,:,:,nn),U_Cal_sec2{mm}(nn),U_v_sec2{mm}(:,nn),U_E_ele2{mm}(nn,:),Dens_E_node2{mm}(:,:,nn)] =Jocob_rod_sec(Q_sec,N_e_1,L_e_1,N_node,Par_E,A);
            U_M2_flag(kk)=U_Cal_sec2{mm}(nn);            
            F_tendon2{mm}(:,nn)=Fq_sec_M2(1:2,nn);
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
            elseif (max(U_M2_flag((kk-10/Damp_2):kk))-min(U_M2_flag((kk-10/Damp_2):kk)))<1e1;%energy variation
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
                Rod_ploting1(q_node2{mm}(:,:,:,nn),Duplication_M,L_co,C_coord,U_E_ele2{mm}(nn,:),Dens_E_node2{mm}(:,:,nn))
                plot3(Q_sec(1,1)*[1,-1],Q_sec(2,1)*[1,1],Q_sec(3,1)*[1,1],'r-')
                plot3(Q_sec(1,end)*[1,1],Q_sec(2,end)*[1,-1],Q_sec(3,end)*[1,1],'g-')
                %%plotting wheels
                P_wheel=Q_sec(1:3,N_e_0+1);
                RM_wheel=reshape(Q_sec(4:12,N_e_0+1),3,3);
                Wheel_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M)
                c=colorbar
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
                    ', X_2=', num2str(Q_sec(2,end),'%.0f'), ...
                    ', Z_1=', num2str(Q_sec(3,1),'%.0f')])
                figure(4);
                subplot(2,2,1)
                plot(U_M2_flag(1:kk))
                subplot(2,2,2)
                if kk>10/Damp_2
                    plot(U_M2_flag(kk-10/Damp_2:kk)-U_M2_flag(kk-10/Damp_2))
                end
                subplot(2,2,3)
                X_1=reshape(Q_sec_M2{mm}(1,1,1:nn),1,[]);
                Y_2=reshape(Q_sec_M2{mm}(2,end,1:nn),1,[]);
                Z_1=reshape(Q_sec_M2{mm}(3,1,1:nn),1,[]);
                yyaxis left
                hline1=plot(Y_2(1)-Y_2(1:nn),-F_tendon2{mm}(1,1:nn));
                ylim([80 inf])
                yyaxis right
                hline2=plot(Y_2(1)-Y_2(1:nn),F_tendon2{mm}(2,1:nn));
                ylim([0 inf])
                pause(0.01)
            end
        end
        save('Rod_Outplane5.mat',"Q_sec_M2",'q_node2',"U_Cal_sec2","F_tendon2",'U_E_ele2','Dens_E_node2')
        
        
        figure(4);
%         hold on
        subplot(2,2,4)
        plot3(X_1,Y_2,Z_1,'-o')
        box on
        grid on
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        view([1 0 0])
        
%         hold off
    end



mm=51
nn=80
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
% subplot(1,2,1)
[haxes,hline1,hline2]=plotyy(Y_2(1)-Y_2,-F_tendon2{mm}(1,:),Y_2(1)-Y_2,F_tendon2{mm}(2,:));

set(hline1,'LineWidth',2);
set(hline2,'LineWidth',2);
% set (gca,'XDir','reverse');
xlabel('Y disp')
ylabel(haxes(1), 'F_{T1}','linewidth',2);
ylabel(haxes(2), 'F_{T2}','linewidth',2);
set(haxes,'FontSize',16,'fontname','Time','linewidth',1)
set(gcf, 'Position', [100 100 540 400]);
set(gca, 'Position', [0.17  0.18  0.65  0.74])   %[x y width height]


%%
load('Rod_Inplane5.mat')
load('Rod_Outplane5.mat')
v = VideoWriter('Rod_Outplane5_2.mp4', 'MPEG-4');
v.FrameRate = 30;  % Set the frame rate (frames per second)
v.Quality = 90;    % Set the video quality (0 to 100)
open(v);
% writeVideo(v, F);
% close(v)
mm=51
pic_num=1
for nn=1:length(disp_y)
    Q_sec=Q_sec_M2{mm}(:,:,nn);
    figure(3);
    clf
    hold on
    Rod_ploting1(q_node2{mm}(:,:,:,nn),Duplication_M,L_co,C_coord,U_E_ele2{mm}(nn,:),Dens_E_node2{mm}(:,:,nn))
    plot3(Q_sec(1,1)*[1,-1],Q_sec(2,1)*[1,1],Q_sec(3,1)*[1,1],'r-',linewidth=2)
    plot3(Q_sec(1,end)*[1,1],Q_sec(2,end)*[1,-1],Q_sec(3,end)*[1,1],'g-',linewidth=2)
    c=colorbar
    c.Label.String = 'Energy Density (mJ/mm)';
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
    title(['D_1=',num2str(2*Q_sec(1,1),'%.0f'), ...
        ', D_2=', num2str(2*Q_sec(2,N_e_1+1),'%.0f'), ...
        ', Z_1=', num2str(Q_sec(3,1),'%.0f')])
    set(gca,'FontSize',16,'fontname','Times')
    drawnow;% 立即刷新当前绘图窗口，这是matlab绘图中动态展示的关键
    F = getframe(gcf);  % 获取当前绘图窗口的图片
    Im = frame2im(F);   % 返回与电影帧相关的图片数据
    [Amap, map] = rgb2ind(Im, 256); % 将RGB图片转化为索引图片
    if pic_num == 1
        imwrite(Amap, map, 'Rod_Outplane5_2.gif', 'gif', 'Loopcount', Inf, 'DelayTime', 0.1);
        % 将第一张图片写入‘sin.gif’文件中，并且将GIF播放次数设置成无穷，即保存的GIF图会一直动下去
    else
        imwrite(Amap, map,'Rod_Outplane5_2.gif','gif','WriteMode','append','DelayTime',0.1);
        % 依次将其他的图片写入到GIF文件当中
    end
    pic_num = pic_num + 1;
    writeVideo(v, F); 
    pause(0.1)
end
close(v)
% set(gcf, 'Position', [100  100  600  500])   %[x y width height]

%% caculate bisymmetric configurations in the range inside the boundary
load('Rod_Inplane5.mat')
load('Rod_BiSym_Outplane5.mat')
A=1e3;
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
P_wheel=[]

for mm=51%N_x+1  
    figure(4);
    clf
    Q_sec=Q_sec_M2{mm}(:,:,end);    
    % caculate the range of disp y based on the boundary 
    disp_x_0=Q_sec(1,1);
    disp_y_0=Q_sec(2,end);    
    disp_y=flip([0:5:disp_y_0, disp_y_0]);% sampling displacement y

    for nn=48:51%length(disp_y)
        mm
        nn
        Damp_1=0.01
        Damp_2=0.1
        if ismember(nn,1:43)
            Damp_1=0.005
            Damp_2=0.05
        elseif ismember(nn,44:46)
            Damp_1=0.001
            Damp_2=0.01
       
        elseif ismember(nn,47:83)
            Damp_1=0.0005
            Damp_2=0.005
%         elseif ismember(nn,48:83)
%             Damp_1=0.0002
%             Damp_2=0.002
% %         elseif ismember(nn,40:63)
% %             Damp_1=0.0001
% %             Damp_2=0.05
%         elseif ismember(nn,84:length(disp_y))
%             Damp_1=0.0001
%             Damp_2=0.025
        end
        if nn~=1
            Q_sec=Q_sec_M3{mm}(:,:,nn-1);
        end
        Q_sec(1,1)=disp_y(nn);        
        Q_sec(2,N_e_1+1)=disp_y(nn);     
        
        if ismember(nn,1:0)
            Q_sec(3,1)=Q_sec(3,1)-0.2*(15-0.5*nn);% add the z axis distorion
        end        

        %define the end points and frames
        kk=0;
        flag1=1;
        flag2=1;
        Damp=Damp_1
        while (kk<=100||flag1>0.000001*R0_ring)&&flag2&&(kk<5/Damp||kk<4000)
            kk=kk+1;
            Damp=0.995^(kk)*Damp_1+(1-0.995^(kk))*Damp_2;
            [dFq_dq_sec3,Fq_sec_M3(:,nn),q_node3{mm}(:,:,:,nn),U_Cal_sec3{mm}(nn),U_v_sec3{mm}(:,nn),U_E_ele3{mm}(nn,:),Dens_E_node3{mm}(:,:,nn)] =Jocob_rod_sec(Q_sec,N_e_1,L_e_1,N_node,Par_E,A);
            U_Cal_flag(kk)=U_Cal_sec3{mm}(nn);            
            F_tendon3{mm}(:,nn)=Fq_sec_M3(1:2,nn);
            dFq_dq_sec_2=Proj_M'*dFq_dq_sec3*Proj_M;
            Fq_sec_M1_2=Proj_M'*Fq_sec_M3(:,nn);
            dQ_sec_1=-pinv(dFq_dq_sec_2)*Proj_M'*Fq_sec_M3(:,nn);
            dQ_sec=Proj_M*dQ_sec_1;
            dQ_sec=reshape(dQ_sec,12,[]);
            %the distance between two mid points of rod1 and rod2
            Q_sec=Q_sec+Damp*dQ_sec;    
            P_wheel(:,1)=Q_sec(1:3,N_e_0+1);            
            Q_sec_M3{mm}(:,:,nn)=Q_sec;
            flag1=sqrt(sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e_0+1));
            flag2=1;
            Samp_L=10/Damp_2;
            if kk<=1000
                flag2=1;
                Samp_L=10/Damp_2;
            elseif kk>1000&kk<=1500
                flag2=1;
                Samp_L=7/Damp_2;
            elseif kk>1500&kk<=2000
                flag2=1;
                Samp_L=5/Damp_2;
            elseif kk>2000&kk<=2500
                flag2=1;
                Samp_L=3/Damp_2;
            elseif kk>2500&kk<=3000;%energy variation
                Samp_L=2/Damp_2;
            elseif kk>3000;%energy variation
                Samp_L=1/Damp_2;
            end

            if kk<=Samp_L
            elseif (max(U_Cal_flag((kk-Samp_L):kk))-min(U_Cal_flag((kk-Samp_L):kk)))<1e1;%energy variation
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
                Rod_ploting1(q_node3{mm}(:,:,:,nn),Duplication_M,L_co,C_coord,U_E_ele3{mm}(nn,:),Dens_E_node3{mm}(:,:,nn))
                plot3(Q_sec(1,1)*[1,-1],Q_sec(2,1)*[1,1],Q_sec(3,1)*[1,1],'r-')
                plot3(Q_sec(1,end)*[1,1],Q_sec(2,end)*[1,-1],Q_sec(3,end)*[1,1],'g-')
                %%plotting wheels
                P_wheel=Q_sec(1:3,N_e_0+1);
                RM_wheel=reshape(Q_sec(4:12,N_e_0+1),3,3);
                Wheel_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M)
                c=colorbar
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
                    ', X_2=', num2str(Q_sec(2,end),'%.0f'), ...
                    ', Z_1=', num2str(Q_sec(3,1),'%.0f')])
                figure(4);
                subplot(2,2,1)
                plot(U_Cal_flag(1:kk))
                subplot(2,2,2)
                if kk>Samp_L
                    plot(U_Cal_flag(kk-Samp_L:kk)-U_Cal_flag(kk-Samp_L))
                end
                X_1=reshape(Q_sec_M3{mm}(1,1,1:nn),1,[]);
                Y_2=reshape(Q_sec_M3{mm}(2,end,1:nn),1,[]);
                Z_1=reshape(Q_sec_M3{mm}(3,1,1:nn),1,[]);
                subplot(2,2,4)
                plotyy(Y_2(1)-Y_2(1:nn),-F_tendon3{mm}(1,1:nn),Y_2(1)-Y_2(1:nn),F_tendon3{mm}(2,1:nn))
                box on
                grid on
                xlabel('X')
                ylabel('Y')
                pause(0.01)
            end
        end
        save('Rod_BiSym_Outplane5.mat',"Q_sec_M3",'q_node3',"U_Cal_sec3","F_tendon3",'U_E_ele3','Dens_E_node3')
        
        figure(4);
%         hold on
        subplot(2,2,3)
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
% % nn=80
X_1=reshape(Q_sec_M3{mm}(1,1,1:nn),1,[]);
Y_2=reshape(Q_sec_M3{mm}(2,end,1:nn),1,[]);
Z_1=reshape(Q_sec_M3{mm}(3,1,1:nn),1,[]);
U_Cal_sec=U_Cal_sec3{mm};
figure(5);
plotyy(Y_2(1)-Y_2,U_Cal_sec,Y_2(1)-Y_2,Z_1)
box on
grid on
xlabel('X')
ylabel('Y')

figure(6)
% subplot(1,2,1)
[haxes,hline1,hline2]=plotyy(Y_2(1)-Y_2,-F_tendon3{mm}(1,:),Y_2(1)-Y_2,F_tendon3{mm}(2,:));

set(hline1,'LineWidth',2);
set(hline2,'LineWidth',2);
% set (gca,'XDir','reverse');
xlabel('Y disp')
ylabel(haxes(1), 'F_{T1}','linewidth',2);
ylabel(haxes(2), 'F_{T2}','linewidth',2);
set(haxes,'FontSize',16,'fontname','Time','linewidth',1)
set(gcf, 'Position', [100 100 540 400]);
set(gca, 'Position', [0.17  0.18  0.65  0.74])   %[x y width height]
%%
load('Rod_Inplane5.mat')
load('Rod_Outplane5.mat')
load('Rod_BiSym_Outplane5.mat')
v = VideoWriter('Rod_BiSym_Outplane5_2.mp4', 'MPEG-4');
v.FrameRate = 30;  % Set the frame rate (frames per second)
v.Quality = 90;    % Set the video quality (0 to 100)
open(v);
% writeVideo(v, F);
% close(v)
mm=51
pic_num=1
for nn=1:51
    Q_sec=Q_sec_M3{mm}(:,:,nn);
    figure(3);
    clf
    hold on
    Rod_ploting1(q_node3{mm}(:,:,:,nn),Duplication_M,L_co,C_coord,U_E_ele3{mm}(nn,:),Dens_E_node3{mm}(:,:,nn))
    plot3(Q_sec(1,1)*[1,-1],Q_sec(2,1)*[1,1],Q_sec(3,1)*[1,1],'r-',LineWidth=2)
    plot3(Q_sec(1,end)*[1,1],Q_sec(2,end)*[1,-1],Q_sec(3,end)*[1,1],'g-',LineWidth=2)
    c=colorbar
    c.Label.String = 'Energy Density (mJ/mm)';
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

    title(['D_1=',num2str(2*Q_sec(1,1),'%.0f'), ...
        ', D_2=', num2str(2*Q_sec(2,N_e_1+1),'%.0f'), ...
        ', Z_1=', num2str(Q_sec(3,1),'%.0f')])
     set(gca,'FontSize',16,'fontname','Times')
    drawnow;% 立即刷新当前绘图窗口，这是matlab绘图中动态展示的关键
    F = getframe(gcf);  % 获取当前绘图窗口的图片
    Im = frame2im(F);   % 返回与电影帧相关的图片数据
    [Amap, map] = rgb2ind(Im, 256); % 将RGB图片转化为索引图片
    if pic_num == 1
        imwrite(Amap, map, 'Rod_BiSym_Outplane5_2.gif', 'gif', 'Loopcount', Inf, 'DelayTime', 0.1);
        % 将第一张图片写入‘sin.gif’文件中，并且将GIF播放次数设置成无穷，即保存的GIF图会一直动下去
    else
        imwrite(Amap, map,'Rod_BiSym_Outplane5_2.gif','gif','WriteMode','append','DelayTime',0.1);
        % 依次将其他的图片写入到GIF文件当中
    end
    pic_num = pic_num + 1;
    writeVideo(v, F);
    pause(0.1)
end
close(v)
%%
load('Rod_Inplane5.mat')
load('Rod_Outplane5.mat')
load('Rod_Bisym_Outplane5.mat')

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
X0=Q_sec(1,1);
disp_x=[0:5:X0,X0];% sampling displacement x
figure(1);
clf
hold on
Rod_ploting1(q_node1(:,:,:,mm),Duplication_M,L_co,C_coord,U_E_ele1(mm,:),Dens_E_node1(:,:,mm));
%%plotting wheels
P_wheel=Q_sec(1:3,N_e_0+1);
RM_wheel=reshape(Q_sec(4:12,N_e_0+1),3,3);
% Wheel_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M)
plot3(Q_sec(1,1)*[1,-1],Q_sec(2,1)*[1,1],Q_sec(3,1)*[1,1],'r-',LineWidth=2)
plot3(Q_sec(1,end)*[1,1],Q_sec(2,end)*[1,-1],Q_sec(3,end)*[1,1],'g-',LineWidth=2)
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
axis([-600 600 -700 700 -10 10])
axis off
axis equal
view([1 1 1])
title(['\itL_{T1}\rm=',num2str(2*(X0-disp_x(mm)),'%.0f'),', \itL_{T2}\rm=',num2str(2*Q_sec_M1(2,N_e_1+1,mm),'%.0f')])
set(gca,'FontSize',16,'fontname','Time','linewidth',1)
set(gcf, 'Position', [100 100 380 280]*0.7);
set(gca, 'Position', [0.0  0.0  1 0.8])   %[x y width height]



mm=51
Q_sec=Q_sec_M1(:,:,mm)
X0=Q_sec(1,1);
disp_x=[0:5:X0,X0];% sampling displacement x
figure(2);
clf
hold on
Rod_ploting1(q_node1(:,:,:,mm),Duplication_M,L_co,C_coord,U_E_ele1(mm,:),Dens_E_node1(:,:,mm));
%%plotting wheels
P_wheel=Q_sec(1:3,N_e_0+1);
RM_wheel=reshape(Q_sec(4:12,N_e_0+1),3,3);
% Wheel_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M)
plot3(Q_sec(1,1)*[1,-1],Q_sec(2,1)*[1,1],Q_sec(3,1)*[1,1],'r-',LineWidth=2)
plot3(Q_sec(1,end)*[1,1],Q_sec(2,end)*[1,-1],Q_sec(3,end)*[1,1],'g-',LineWidth=2)
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
axis([-600 600 -700 700 -10 10])
axis off
axis equal
view([1 1 1])
title(['\itL_{T1}\rm=',num2str(2*(X0-disp_x(mm)),'%.0f'),', \itL_{T2}\rm=',num2str(2*Q_sec_M1(2,N_e_1+1,mm),'%.0f')])
set(gca,'FontSize',16,'fontname','Time','linewidth',1)
set(gcf, 'Position', [100 100 380 280]*0.7);
set(gca, 'Position', [0.0  0.0  1 0.8])   %[x y width height]

%
mm=51
nn=80
Q_sec=Q_sec_M2{mm}(:,:,nn);
figure(3);
clf
hold on
Rod_ploting1(q_node2{mm}(:,:,:,nn),Duplication_M,L_co,C_coord,U_E_ele2{mm}(nn,:),Dens_E_node2{mm}(:,:,nn))
plot3(Q_sec(1,1)*[1,-1],Q_sec(2,1)*[1,1],Q_sec(3,1)*[1,1],'r-',LineWidth=2)
plot3(Q_sec(1,end)*[1,1],Q_sec(2,end)*[1,-1],Q_sec(3,end)*[1,1],'g-',LineWidth=2)
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

axis([-300 300 -300 300 -500 50])
axis equal
axis off
view([1 1 1])
title(['\itL_{T1}\rm=',num2str(2*(X0-disp_x(mm)),'%.0f'),', \itL_{T2}\rm=',num2str(2*Q_sec(2,N_e_1+1),'%.0f')])
xlabel('X')
ylabel('Y')
zlabel('Z')
set(gca,'FontSize',16,'fontname','Time','linewidth',1)
set(gcf, 'Position', [100 100 380 280]*0.7);
set(gca, 'Position', [0.0  0.0  1 0.8])   %[x y width height]




%
mm=51
nn=51
Q_sec=Q_sec_M3{mm}(:,:,nn);
figure(4);
clf
hold on
Rod_ploting1(q_node3{mm}(:,:,:,nn),Duplication_M,L_co,C_coord,U_E_ele3{mm}(nn,:),Dens_E_node3{mm}(:,:,nn))
plot3(Q_sec(1,1)*[1,-1],Q_sec(2,1)*[1,1],Q_sec(3,1)*[1,1],'r-',LineWidth=2)
plot3(Q_sec(1,end)*[1,1],Q_sec(2,end)*[1,-1],Q_sec(3,end)*[1,1],'g-',LineWidth=2)
c=colorbar
caxis([0 350])
% c.Label.String = 'Energy Density (mJ/mm)';

% % To hide the colorbar
% colorbar('visible', 'off')

%%plotting wheels
P_wheel=Q_sec(1:3,N_e_0+1);
RM_wheel=reshape(Q_sec(4:12,N_e_0+1),3,3);
theta=150*pi/180
RM_wheel=RM_wheel*[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]
Wheel_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M)
hold off
box off
grid off
axis([-300 300 -300 300 -500 10])
axis equal
axis off
view([1 1 1])
% title(['\itL_{T1}\rm=',num2str(2*Q_sec(1,1),'%.0f'),', \itL_{T2}\rm=',num2str(2*Q_sec(2,N_e_1+1),'%.0f')])
xlabel('X')
ylabel('Y')
zlabel('Z')
set(gca,'FontSize',16,'fontname','Time','linewidth',1)
set(gcf, 'Position', [100 100 380 280]*0.7);
set(gca, 'Position', [0.0  0.0  1 0.8])   %[x y width height]
set(gca, 'Position', [0.0  0.1  0.65 0.9])   %[x y width height]

F_tendon=F_tendon1(:,1:51)
F_tendon=[F_tendon,F_tendon2{51}]
F_tendon=[F_tendon,F_tendon3{51}]

figure(6)
plot(-F_tendon(1,:),'r','LineWidth',2);
hold on
plot(F_tendon(2,:),'g','LineWidth',2);
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

clear L_tendon_1 L_tendon_2
L_tendon_1=reshape(Q_sec_M1(1,1,1:51),1,[])*2/1000*4/pi
L_tendon_2=reshape(Q_sec_M1(2,end,1:51),1,[])*2/1000*4/pi
L_tendon_1=[L_tendon_1,reshape(Q_sec_M2{51}(1,1,1:80),1,[])*2/1000*4/pi]
L_tendon_2=[L_tendon_2,reshape(Q_sec_M2{51}(2,end,1:80),1,[])*2/1000*4/pi]
L_tendon_1=[L_tendon_1(1:end-1),reshape(Q_sec_M3{51}(1,1,1:51),1,[])*2/1000*4/pi]
L_tendon_2=[L_tendon_2(1:end-1),reshape(Q_sec_M3{51}(2,end,1:51),1,[])*2/1000*4/pi]

figure(7)
plot(0:length(L_tendon_1)-1,L_tendon_1,'linewidth',2)
hold on
plot(0:length(L_tendon_1)-1,L_tendon_2,'linewidth',2)
hold off
ylim([0,1.8])
xlabel('Simulation step k')
ylabel('Tension Length (m)')
set(gca,'FontSize',16,'fontname','Times')
legend('Tendon1','Tendon2','FontSize',14,'fontname','times new roman')
set(gcf, 'Position', [100 100 380 320].*[1 1 1 0.8]);
set(gca, 'Position', [0.2  0.3  0.75  0.6])   %[x y width height]
%%
v = VideoWriter('Rod_tendon_force.mp4', 'MPEG-4');
v.FrameRate = 30;  % Set the frame rate (frames per second)
v.Quality = 90;    % Set the video quality (0 to 100)
open(v);
% writeVideo(v, F);
% close(v)
pic_num=1
for ii=1:length(F_tendon(1,:))
    figure(6)
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
        imwrite(Amap, map, 'Rod_tendon_force.gif', 'gif', 'Loopcount', Inf, 'DelayTime', 0.1);
        % 将第一张图片写入‘sin.gif’文件中，并且将GIF播放次数设置成无穷，即保存的GIF图会一直动下去
    else
        imwrite(Amap, map,'Rod_tendon_force.gif','gif','WriteMode','append','DelayTime',0.1);
        % 依次将其他的图片写入到GIF文件当中
    end
    pic_num = pic_num + 1;
    writeVideo(v, F);
% close(v)
    pause(0.1)
end
% writeVideo(v, F);
close(v)


























