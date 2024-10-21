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
    q_node0=Curve_interp(Q_sec(:,nn),Q_sec(:,nn+1),L_e_0,N_node);
    Rod_ploting1(q_node0,Duplication_M,L_co,C_coord);
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
    [dFq_dq_sec0,Fq_sec_M0(:,mm),q_node0,U_Cal_sec0(mm)] =Jocob_rod_sec(Q_sec,N_e_0,L_e_0,N_node,Par_E,A);
    dFq_dq_sec_2=Proj_M'*dFq_dq_sec0*Proj_M;
    dQ_sec_1=-pinv(dFq_dq_sec_2)*Proj_M'*Fq_sec_M0(:,mm);
    dQ_sec=Proj_M*dQ_sec_1;
    dQ_sec=reshape(dQ_sec,12,[]);
    Q_sec=Q_sec+Damp*dQ_sec;
    Q_sec_M0(:,:,mm)=Q_sec;
    flag1=sqrt(sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e_0+1))
    flag2=0;%mean(abs(Proj_M'*Fq_sec_M(:,mm)));%
    % AA=reshape(Fq_sec_M(:,mm),12,[]);
    if kk==2000
        flag3=1;
    end
end

%plot the circle
Duplication_M=[1 1 1;-1 1 1;-1 -1 1;1 -1 1;];
figure(2);
hold on
Rod_ploting2(Q_sec,Duplication_M,L_co*5,C_coord)

for nn=1:(size(Q_sec,2)-1)
    q_node0=Curve_interp(Q_sec(:,nn),Q_sec(:,nn+1),L_e_0,N_node); 
    Rod_ploting1(q_node0,Duplication_M,L_co,C_coord);
end
hold off
box on
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
view([1 1 1])

% Redefine the element number
N_e_1=N_e_0;
L_sec_1=L_sec_0;
N_e_1=N_e_0;% the element number
L_e_1=L_sec_1/N_e_1;
N_e_0=N_e_1/2
L_sec_0=L_sec_1/2;



save('Rod_Origin',"Q_sec_M0")
%% caculate the inplane deformation when contract tendon1
load('Rod_Origin.mat')
load('Rod_Inplane.mat')
Q_sec=Q_sec_M0(:,:,1);
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

for mm=51%:N_x+1
    Q_sec=Q_sec_M1(:,:,mm)
    figure(2);
    clf
    hold on
    Rod_ploting1(q_node1(:,:,:,mm),Duplication_M,L_co,C_coord,U_E_ele1(mm,:),Dens_E_node1(:,:,mm));
    plot3(Q_sec(1,1)*[1,-1],Q_sec(2,1)*[1,1],Q_sec(3,1)*[1,1],'r-')
    plot3(Q_sec(1,end)*[1,1],Q_sec(2,end)*[1,-1],Q_sec(3,end)*[1,1],'g-')
    c=colorbar   
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

save('Rod_Inplane.mat',"Q_sec_M1",'q_node1','F_tendon1','U_Cal_sec1','U_v_sec1','U_E_ele1','Dens_E_node1')

%% caculate outplane configurations when contract tendon2
load('Rod_Inplane.mat')
load('Rod_Outplane.mat')
A=1e3;
Damp=0.1;%reduce the iteration speed
% Build the boundary conditions for the boundary
% start point x_disp fixed and end point y_disp in control
%Using the de_demension matrix to define the constrain of the structure
Proj_M=eye(12*(N_e_1+1),12*(N_e_1+1));%dedemesionalize the variables

%rod 1
D_fix=0*Q_sec;
% start point

% constrain the start point
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


for mm=1:N_x+1  
    figure(3);
    clf
    Q_sec=Q_sec_M1(:,:,mm);    
    % caculate the range of disp y based on the boundary 
    disp_x_0=Q_sec(1,1);
    disp_y_0=Q_sec(2,end);    
    disp_y=flip([disp_x_0:5:disp_y_0, disp_y_0]);% sampling displacement y
    Fq_sec_M2=[];%initializing
    for nn=133:length(disp_y) 
        Damp_1=0.001; 
        Damp_2=0.2; 
        if ismember(mm,[43 45 46 47 49 57])
            Damp_1=0.001
            Damp_2=0.4
        elseif ismember(mm,[58:60])
            Damp_1=0.001
            Damp_2=0.4
        elseif ismember(mm,61:63)
            Damp_1=0.001
            Damp_2=0.4
        elseif ismember(mm,64:67)
            Damp_1=0.001
            Damp_2=0.4
        elseif ismember(mm,68:70)
            Damp_1=0.001
            Damp_2=0.4
            A=0.7e3;
        elseif ismember(mm,71:74)
            Damp_1=0.001
            Damp_2=0.4
            A=0.7e3;
        elseif ismember(mm,75:78)
            Damp_1=0.001
            Damp_2=0.4
            A=0.7e3;
        elseif ismember(mm,79:99)
            Damp_1=0.001
            Damp_2=0.4
            A=0.5e3;
        elseif ismember(mm,100)
            Damp_1=0.001
            Damp_2=0.2
            A=0.5e3;
        elseif ismember(mm,101)
            Damp_1=0.001
            Damp_2=0.2
            A=0.5e3;
                      
% %         elseif ismember(nn,37:83)
% %             Damp_1=0.0001
% %             Damp_2=0.01
% % %         elseif ismember(nn,40:63)
% % %             Damp_1=0.0001
% % %             Damp_2=0.05
% %         elseif ismember(nn,84:length(disp_y))
% %             Damp_1=0.0001
% %             Damp_2=0.025
        end
        if nn~=1
            Q_sec=Q_sec_M2{mm}(:,:,nn-1);
        end
%         Q_sec(1,1)=disp_y(nn);        
        Q_sec(2,N_e_1+1)=disp_y(nn); 
        
        if ismember(mm,1:11)&&ismember(nn,1:10)
            Q_sec(3,1)=Q_sec(3,1)-0.5;% add the z axis distorion
        elseif ismember(mm,12:18)&&ismember(nn,1:10)
            Q_sec(3,1)=Q_sec(3,1)+1;% add the z axis distorion
        elseif ismember(mm,19:23)&&ismember(nn,1:10)
            Q_sec(3,1)=Q_sec(3,1)+2;% add the z axis distorion
        elseif ismember(mm,[28 29 34])&&ismember(nn,1:10)
            Q_sec(3,1)=Q_sec(3,1)+5;% add the z axis distorion
        elseif ismember(mm,24:40)&&ismember(nn,1:10)
            Q_sec(3,1)=Q_sec(3,1)+4;% add the z axis distorion        
        elseif ismember(mm,[45])&&ismember(nn,1:15)
            Q_sec(3,1)=Q_sec(3,1)+8;% add the z axis distorion
        elseif ismember(mm,[47 49 55 ])&&ismember(nn,1:15)
            Q_sec(3,1)=Q_sec(3,1)+8;% add the z axis distorion
        elseif ismember(mm,[58:60])&&ismember(nn,1:15)
            Q_sec(3,1)=Q_sec(3,1)-2;% add the z axis distorion
        elseif ismember(mm,[61:63])&&ismember(nn,1:length(disp_y))
            Q_sec(3,1)=Q_sec(3,1)+2;% add the z axis distorion
        elseif ismember(mm,[64:69])&&ismember(nn,1:20)
            Q_sec(3,1)=Q_sec(3,1)+2;% add the z axis distorion
        elseif ismember(mm,[70:71])&&ismember(nn,[1:20 53:length(disp_y)])
            Q_sec(3,1)=Q_sec(3,1)+5;% add the z axis distorion
        elseif ismember(mm,[72:73])&&ismember(nn,[1:20 56:length(disp_y) ])
            Q_sec(3,1)=Q_sec(3,1)+8;% add the z axis distorion
        elseif ismember(mm,[74])&&ismember(nn,[1:20, 50:length(disp_y)  ])
            Q_sec(3,1)=Q_sec(3,1)+2;% add the z axis distorion
        elseif ismember(mm,[75:77])            
            if ismember(nn,[1:20])
                Q_sec(3,1)=Q_sec(3,1)+2;
            elseif ismember(nn,[51:length(disp_y)])
                Q_sec(3,1)=Q_sec(3,1)-2;
            end% add the z axis distorion    
        elseif ismember(mm,[78])            
            if ismember(nn,[1:40])
                Q_sec(3,1)=Q_sec(3,1)-5;
            elseif ismember(nn,[51:length(disp_y)])
                Q_sec(3,1)=Q_sec(3,1)-0;
            end% add the z axis distorion
        elseif ismember(mm,[79])
            if ismember(nn,[1:55])
                Q_sec(3,1)=Q_sec(3,1)+2;
            elseif ismember(nn,[56:length(disp_y)])
                Q_sec(3,1)=Q_sec(3,1)-2;
            end% add the z axis distorion
        elseif ismember(mm,[80:83])
            if ismember(nn,[1:20])
                Q_sec(3,1)=Q_sec(3,1)-5;
            elseif ismember(nn,[50:length(disp_y)])
                Q_sec(3,1)=Q_sec(3,1)-5;
            end% add the z axis distorion
        elseif ismember(mm,[86])
            if ismember(nn,[1:20])
                Q_sec(3,1)=Q_sec(3,1)-5;
            elseif ismember(nn,[50:length(disp_y)])
                Q_sec(3,1)=Q_sec(3,1)-5;
            end
        elseif ismember(mm,[84 85 87:88 ])
            if ismember(nn,[1:20])
                Q_sec(3,1)=Q_sec(3,1)-8;
            elseif ismember(nn,[50:length(disp_y)])
                Q_sec(3,1)=Q_sec(3,1)-5;
            end% add the z axis distorion
        elseif ismember(mm,[89:90])
            if ismember(nn,[1:20])
                Q_sec(3,1)=Q_sec(3,1)-5;
            elseif ismember(nn,[51:length(disp_y)])
                Q_sec(3,1)=Q_sec(3,1)-0;
            end
        elseif ismember(mm,[90])
            if ismember(nn,[1:20])
                Q_sec(3,1)=Q_sec(3,1)-5;
            elseif ismember(nn,[51:length(disp_y)])
                Q_sec(3,1)=Q_sec(3,1)-5;
            end
        elseif ismember(mm,[91])
            if ismember(nn,[1:20])
                Q_sec(3,1)=Q_sec(3,1)-5;
            elseif ismember(nn,[45:length(disp_y)])
                Q_sec(3,1)=Q_sec(3,1)-5;
            end
        elseif ismember(mm,[92])
            if ismember(nn,[1:20])
                Q_sec(3,1)=Q_sec(3,1)+8;
            elseif ismember(nn,[45:length(disp_y)])
                Q_sec(3,1)=Q_sec(3,1)-2;
            end
        elseif ismember(mm,[93])
            if ismember(nn,[1:20])
                Q_sec(3,1)=Q_sec(3,1)+8;
            elseif ismember(nn,[45:length(disp_y)])
                Q_sec(3,1)=Q_sec(3,1)-2;
            end
        elseif ismember(mm,[94])
            if ismember(nn,[1:20])
                Q_sec(3,1)=Q_sec(3,1)-5;
            elseif ismember(nn,[50:125])
                Q_sec(3,1)=Q_sec(3,1)-1;
            elseif ismember(nn,[126:length(disp_y)])
                Q_sec(3,1)=Q_sec(3,1)-0;
            end
        elseif ismember(mm,[95])
            if ismember(nn,[1:20])
                Q_sec(3,1)=Q_sec(3,1)-2;
            elseif ismember(nn,[50:120])
                Q_sec(3,1)=Q_sec(3,1)-0;
            elseif ismember(nn,[121:length(disp_y)])
                Q_sec(3,1)=Q_sec(3,1)+1;
            end
        elseif ismember(mm,[96])
            if ismember(nn,[1:20])
                Q_sec(3,1)=Q_sec(3,1)-2;
            elseif ismember(nn,[50:120])
                Q_sec(3,1)=Q_sec(3,1)-0;
            elseif ismember(nn,[121:length(disp_y)])
                Q_sec(3,1)=Q_sec(3,1)+1;
            end
        elseif ismember(mm,[97])
            if ismember(nn,[1:20])
                Q_sec(3,1)=Q_sec(3,1)-2;
            elseif ismember(nn,[36:120])
                Q_sec(3,1)=Q_sec(3,1)-1;
            elseif ismember(nn,[121:length(disp_y)])
                Q_sec(3,1)=Q_sec(3,1)+1;
            end
        elseif ismember(mm,[98])
            if ismember(nn,[1:20])
                Q_sec(3,1)=Q_sec(3,1)-2;
            elseif ismember(nn,[36:120])
                Q_sec(3,1)=Q_sec(3,1)-0;
            elseif ismember(nn,[133:length(disp_y)])
                Q_sec(3,1)=Q_sec(3,1)-0;
            end
        elseif ismember(mm,[99])
            if ismember(nn,[1:20])
                Q_sec(3,1)=Q_sec(3,1)-2;
            elseif ismember(nn,[36:120])
                Q_sec(3,1)=Q_sec(3,1)-0;
            elseif ismember(nn,[133:length(disp_y)])
                Q_sec(3,1)=Q_sec(3,1)-0.5;
            end
        elseif ismember(mm,[100])
            if ismember(nn,[1:20])
                Q_sec(3,1)=Q_sec(3,1)-5;
            elseif ismember(nn,[58:120])
                Q_sec(3,1)=Q_sec(3,1)-1;
            elseif ismember(nn,[130:132])
                Q_sec(3,1)=Q_sec(3,1)+0.1;
            elseif ismember(nn,[133:length(disp_y)])
                Damp_2=0.5;
                Q_sec(3,1)=Q_sec(3,1)-0;
            end
        elseif ismember(mm,[101])
            if ismember(nn,[1:4])
                Q_sec(3,1)=Q_sec(3,1)-2;
            elseif ismember(nn,[21:30])
                Q_sec(3,1)=Q_sec(3,1)-1;
            elseif ismember(nn,[31:41])                  
                Q_sec(3,1)=Q_sec(3,1)-0;
            elseif ismember(nn,[42:49])                  
                Q_sec(3,1)=Q_sec(3,1)-0.1;
            elseif ismember(nn,[133:length(disp_y)])
                Damp_2=0.5;
                Q_sec(3,1)=Q_sec(3,1)-0.1;
            end
            %         elseif ismember(mm,41:80)&&ismember(nn,1:10)
            %         elseif ismember(mm,41:80)&&ismember(nn,1:10)
%             Q_sec(3,1)=Q_sec(3,1)+6;% add the z axis distorion
        end              

        %define the end points and frames
        kk=0;
        flag1=1;
        flag2=1;
        Damp=Damp_1;
        while (kk<=100||flag1>0.000001*R0_ring)&&flag2&&(kk<5/Damp||kk<2000)
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
                fprintf('mm=%u,nn=%u,kk=%u, flag1=%f \n',mm,nn,kk,flag1)
            end
%             if mod(kk,100)==2
%                 figure(3);                
%                 clf                
%                 hold on
%                 Rod_ploting1(q_node2{mm}(:,:,:,nn),Duplication_M,L_co,C_coord,U_E_ele2{mm}(nn,:),Dens_E_node2{mm}(:,:,nn))
%                 plot3(Q_sec(1,1)*[1,-1],Q_sec(2,1)*[1,1],Q_sec(3,1)*[1,1],'r-')
%                 plot3(Q_sec(1,end)*[1,1],Q_sec(2,end)*[1,-1],Q_sec(3,end)*[1,1],'g-')                
%                 c=colorbar;
%                 c.Label.String = 'Energy Density (mJ/mm)';
%                 hold off
%                 box on
%                 grid on
%                 xlabel('X')
%                 ylabel('Y')
%                 zlabel('Z')                
%                 axis([-600 600 -700 700 -700 700])
%                 axis equal
%                 view([1 1 1])
%                 title(['X_1=',num2str(Q_sec(1,1),'%.0f'), ...
%                     ', X_2=', num2str(Q_sec(2,end),'%.0f'), ...
%                     ', Z_1=', num2str(Q_sec(3,1),'%.0f')])
%                 figure(4);
%                 subplot(2,2,1)
%                 plot(U_M2_flag(1:kk))
%                 subplot(2,2,2)
%                 if kk>10/Damp_2
%                     plot(U_M2_flag(kk-10/Damp_2:kk)-U_M2_flag(kk-10/Damp_2))
%                 end
%                 subplot(2,2,3)
%                 X_1=reshape(Q_sec_M2{mm}(1,1,1:nn),1,[]);
%                 Y_2=reshape(Q_sec_M2{mm}(2,end,1:nn),1,[]);
%                 Z_1=reshape(Q_sec_M2{mm}(3,1,1:nn),1,[]);                
%                 hline1=plot(Y_2(1)-Y_2(1:nn),-F_tendon2{mm}(1,1:nn));
%                 ylim([80 inf])
%                 subplot(2,2,4)
%                 hline2=plot(Y_2(1)-Y_2(1:nn),F_tendon2{mm}(2,1:nn));
%                 ylim([0 inf])
%                 pause(0.01)
%             end
        end
        save('Rod_Outplane.mat',"Q_sec_M2",'q_node2',"U_Cal_sec2","F_tendon2",'U_E_ele2','Dens_E_node2')      
        if nn<8||mod(nn,1)==0||nn==length(disp_y)
            figure(3);
            hold on
            Rod_ploting1(q_node2{mm}(:,:,:,nn),Duplication_M,L_co*0.3,C_coord,U_E_ele2{mm}(nn,:),Dens_E_node2{mm}(:,:,nn))
            plot3(Q_sec(1,1)*[1,-1],Q_sec(2,1)*[1,1],Q_sec(3,1)*[1,1],'r-')
            plot3(Q_sec(1,end)*[1,1],Q_sec(2,end)*[1,-1],Q_sec(3,end)*[1,1],'g-')
            c=colorbar;
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
                ', Y_2=', num2str(Q_sec(2,end),'%.0f'), ...
                ', Z_1=', num2str(Q_sec(3,1),'%.0f')])           
            set(gcf, 'Position', [305 50 800 750]);
            drawnow
        end
        if nn<8||mod(nn,2)==1||nn==length(disp_y)
             X_1=reshape(Q_sec_M2{mm}(1,1,1:nn),1,[]);
            Y_2=reshape(Q_sec_M2{mm}(2,end,1:nn),1,[]);
            Z_1=reshape(Q_sec_M2{mm}(3,1,1:nn),1,[]);
            figure(4);
            subplot(4,1,1)
            plot(U_Cal_sec2{mm}(1:nn))            
            ylabel('U')
            subplot(4,1,2)
            plot(Z_1(1:nn))            
            box on
            grid on
            xlabel('Disp Y')            
            ylabel('Z')
            subplot(4,1,3)
            hline1=plot(disp_y(1:nn),-F_tendon2{mm}(1,1:nn));
            xlabel('Disp Y')
            set (gca,'XDir','reverse');
            ylabel('F_{T1}')
            subplot(4,1,4)
            hline2=plot(disp_y(1:nn),F_tendon2{mm}(2,1:nn));
            ylim([0 inf])
            ylabel('F_{T2}')
            xlabel('Disp Y')
            set (gca,'XDir','reverse');
            set(gcf, 'Position', [0 50 300 750]);            
            drawnow
        end
    end
    figure(6)
    clf
    mm
    for ii=1:mm
        L_tendon=[];
        Disp_Z=[];
        L_tendon(1,:)=Q_sec_M2{ii}(1,1,:);
        L_tendon(2,:)=Q_sec_M2{ii}(2,end,:);
        F_tendon=F_tendon2{ii}(1:2,:);
        Disp_Z(1,:)=Q_sec_M2{ii}(3,1,:);
        figure(6)
        subplot(2,1,1)
        hold on
        plot3 (L_tendon(1,:),L_tendon(2,:),F_tendon(1,:),'LineWidth',2,'color',[.0 .45 .74])
        plot3 (L_tendon(1,:),L_tendon(2,:),F_tendon(2,:),'LineWidth',2,'color',[.85 .33 .10])
        legend('F_{T1}','F_{T2}')
        hold off
        grid on
        xlabel('L_{t1}')
        ylabel('L_{t2}')
        zlabel('F_{T} (N)');
        view([1 -1 1])
        set (gca,'XDir','reverse')
        subplot(2,1,2)
        hold on
        plot3(ii+L_tendon(1,:)*0,L_tendon(2,:),Disp_Z,'LineWidth',2,'color',[.0 .45 .74])
        hold off
        %     axis equal
        grid on
        xlabel('L_{t1}')
        ylabel('L_{t2}')
        zlabel('Z1');
        view([1 -1 1])
        set (gca,'XDir','reverse')
        set(gcf, 'Position', [1110 50 430 750]);
    end

end

%%
%stability analysis for twisting from inplane config

load('Rod_Inplane.mat')
load('Rod_Outplane.mat')
load('Rod_Anti_twist_Outplane.mat')


mm=91
nn=1
Q_sec=Q_sec_M2{mm}(:,:,nn);
Frame_end_0=reshape(Q_sec_M2{mm}(4:12,end,nn),3,[]);
% Build the boundary conditions for the boundary
% start point x_disp fixed and end point y_disp in control
%Using the de_demension matrix to define the constrain of the structure
Proj_M=eye(12*(N_e_1+1),12*(N_e_1+1));%dedemesionalize the variables

%rod 1
D_fix=0*Q_sec;

% start point
% constrain the start point
%  constrain  X Y Z disp
D_fix(1,1)=1;%fix  the X disp
D_fix(2,1)=1;%fix the Y disp
D_fix(3,1)=1;%free the  Z disp
%  constrain/fix the local frame
D_fix(4:6,1)=0;%release X axis rotate at Y axis direction
% D_fix(5,1)=1;%
D_fix(8:9,1)=1;%Y axis 
D_fix(10:12,1)=0;%release Z axis rotate at Y axis direction
% D_fix(11,1)=1;%

%end point
% constrain the end-middle point of rd1 and rd2
%  constrain  X Y Z disp 
D_fix(1,end)=1;%fix the X disp
D_fix(2,end)=0;%free the Y disp
D_fix(3,end)=1;%fix the Z disp
%  constrain/Fix local frame
D_fix(4:6,end)=1;%control the X axis 
D_fix(7:9,end)=1;%free the Y axis
D_fix(10:12,end)=1;%free the Z axis 



B_index=D_fix==0;% find the fixed demensions
Proj_M=Proj_M(:,B_index);% extract the projection matrix of free paramters

T_angle=[0:pi/180:90*pi/180]


Damp_1=0.001
Damp_2=0.5
A=0.75e3;
Duplication_M=[1,1,1]%;-1,1,-1;-1,-1,-1;1,-1,1;]
figure(3)
clf
for jj=1:length(T_angle)
    jj
    if jj~=1
        Q_sec=Q_sec_M3{mm}(:,:,jj-1);
    end
    R_M_z=[cos(T_angle(jj))  0  -sin(T_angle(jj)) ;
        0                 1               0;
        sin(T_angle(jj))  0  cos(T_angle(jj)) ;];
    Q_sec(4:6,end)=R_M_z*Frame_end_0(:,1);

    kk=0;
    flag1=1;
    flag2=1;
    L_step=0.05;
    while (kk<=100||flag1>0.000001*R0_ring)&&flag2&&(kk<5/Damp||kk<2000)
        kk=kk+1;
        Damp=0.995^(kk)*Damp_1+(1-0.995^(kk))*Damp_2;
        [dFq_dq_sec3,Fq_sec_M3(:,jj),q_node3{mm}(:,:,:,jj),U_Cal_sec3{mm}(jj),U_v_sec3{mm}(:,jj),U_E_ele3{mm}(jj,:),Dens_E_node3{mm}(:,:,jj)] =Jocob_rod_sec(Q_sec,N_e_1,L_e_1,N_node,Par_E,A);
        U_M2_flag(kk)=U_Cal_sec3{mm}(jj);
        F_tendon3{mm}(:,jj)=Fq_sec_M3(1:2,jj);
        dFq_dq_sec_2=Proj_M'*dFq_dq_sec3*Proj_M;
        Fq_sec_M1_2=Proj_M'*Fq_sec_M3(:,jj);
        dQ_sec_1=-pinv(dFq_dq_sec_2)*Proj_M'*Fq_sec_M3(:,jj);
        dQ_sec=Proj_M*dQ_sec_1;
        dQ_sec=reshape(dQ_sec,12,[]);
        %the distance between two mid points of rod1 and rod2
        Q_sec=Q_sec+Damp*dQ_sec;
        P_wheel(:,1)=Q_sec(1:3,N_e_0+1);
        Q_sec_M3{mm}(:,:,jj)=Q_sec;
        flag1=sqrt(sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e_0+1));
        flag2=1;
        if kk<=10/Damp_2
            flag2=1;
        elseif (max(U_M2_flag((kk-10/Damp_2):kk))-min(U_M2_flag((kk-10/Damp_2):kk)))<1e1;%energy variation
            flag2=0;
        end
        if mod(kk,20)==1
            fprintf('mm=%u,jj=%u,kk=%u, flag1=%f \n',mm,jj,kk,flag1)
        end
    end
    save('Rod_Anti_twist_Outplane.mat',"Q_sec_M3",'q_node3',"U_Cal_sec3","F_tendon3",'U_E_ele3','Dens_E_node3')
    if jj<8||mod(jj,5)==0||jj==length(T_angle)
        figure(3);
        hold on
        Rod_ploting1(q_node3{mm}(:,:,:,jj),Duplication_M,L_co*0.3,C_coord,U_E_ele3{mm}(jj,:),Dens_E_node3{mm}(:,:,jj))
        plot3(Q_sec(1,1)*[1,-1],Q_sec(2,1)*[1,1],Q_sec(3,1)*[1,1],'r-')
        plot3(Q_sec(1,end)*[1,1],Q_sec(2,end)*[1,-1],Q_sec(3,end)*[1,1],'g-')
        c=colorbar;
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
            ', Y_2=', num2str(Q_sec(2,end),'%.0f'), ...
            ', Z_1=', num2str(Q_sec(3,1),'%.0f')])
        set(gcf, 'Position', [305 50 800 750]);
        drawnow
    end
    if jj<8||mod(jj,1)==0||jj==length(T_angle)
        X_1=reshape(Q_sec_M3{mm}(1,1,1:jj),1,[]);
        Y_2=reshape(Q_sec_M3{mm}(2,end,1:jj),1,[]);
        Z_1=reshape(Q_sec_M3{mm}(3,1,1:jj),1,[]);
        figure(4);
        subplot(4,1,1)
        plot(U_v_sec3{mm}(2,1:jj))
        ylabel('U')
        subplot(4,1,2)
        plot(Y_2(1:jj))
        box on
        grid on
        xlabel('T_{angle}')
        ylabel('Y_2')
        subplot(4,1,3)
        hline1=plot(T_angle(1:jj)/pi*180,-F_tendon3{mm}(1,1:jj));
        xlabel('T_{angle}')
%         set (gca,'XDir','reverse');
        ylabel('F_{T1}')
        subplot(4,1,4)
        hline2=plot(T_angle(1:jj)/pi*180,F_tendon3{mm}(2,1:jj));
        ylim([0 inf])
        ylabel('F_{T2}')
        xlabel('T_{angle}')
%         set (gca,'XDir','reverse');
        set(gcf, 'Position', [0 50 300 750]);
        drawnow
    end


end



% %%
% % build the half rod structure
% Q_sec_1=Q_sec_M2{mm}(:,:,nn);
% N_e_2=2*(size(Q_sec_1,2)-1);% get the element number from the data
% Q_sec_2(1:3,:)=diag([1,-1,1])*Q_sec_1(1:3,:);
% Q_sec_2(4:6,:)=diag([-1,-1,-1])*diag([1,-1,1])*Q_sec_1(4:6,:);
% Q_sec_2(7:9,:)=diag([1,1,1])*diag([1,-1,1])*Q_sec_1(7:9,:);
% Q_sec_2(10:12,:)=diag([1,1,1])*diag([1,-1,1])*Q_sec_1(10:12,:);
% Q_sec=[flip(Q_sec_2(:,2:end),2),Q_sec_1];
% 
% L_sec_2=L_sec_1*2;
% L_e_2=L_sec_2/N_e_2;
% L_e_1,L_e_0
% 
% %plot the start inplane config
% figure(4);
% C_coord=[.0 .45 .74; .74 .45 .0];
% L_co=20;
% Duplication_M=[1 1 1;-1 -1 1]
% hold on
% Rod_ploting2(Q_sec,Duplication_M,L_co,C_coord)
% hold off
% axis equal
% 
% % Set the boundary conditions for controlled twisting deformation
% 
% %Using the de_demension matrix to define the constrain of the structure
% Proj_M_0=eye(12*(N_e_1+1),12*(N_e_1+1));%dedemesionalize the variables
% 
% 
% % based on anti sym condition
% % the start and end point 
% 
% % start point
% D_fix=0*Q_sec;%fixed dimensions/varibales
% D_fix(1,1)=1;%fix the X disp
% D_fix(2,1)=1;%fix the Y disp
% D_fix(3,1)=1;%fix the  Z disp
% 
% %control the local rotation by constrain the local z axis
% D_fix(10:12,1)=1;%local z axis in control 
% r1_frame_0=reshape(Q_sec(4:12,1),3,3);
% 
% % the local x y z axis in control
% % because the equalism relationship
% %D_fix(1:3,1)=[-1 -1 1].*D_fix(4:6,end);%x axis in control 
% %D_fix(4:6,1)=[-1 -1 1].*D_fix(4:6,end);%x axis in control 
% %D_fix(7:9,1)=[-1 -1 1].*D_fix(4:6,end);%x axis in control 
% %D_fix(10:12,1)=[-1 -1 1].*D_fix(4:6,end);%x axis in control 
% % sus
% Con_M_1=diag(repmat([-1 -1 1],1,4));
% % Proj_M_0(1:12,end-11:end)=0
% Proj_M_0(end-11:end,1:12)=Con_M_1;
% 
% % end point in constrain
% D_fix(1,end)=1;%fix the X disp
% D_fix(2,end)=1;%control the Y disp
% D_fix(3,end)=1;%control the Z disp
% 
% % local x y z aixs in control
% D_fix([4 5 6],end)=1;%fix the rotation
% D_fix([7 8 9],end)=1;%fix the rotation
% D_fix([10 11 12],end)=1;%fix the rotation
% 
% n_m=N_e_1/2+1;%the seq numb of middle point
% 
% % the middle point position
% % based on anti sym condition
% % the distance of middle point to the z axis is constrained: x^2+y^2=C
% % introduce the R and theta to replace the x y disp of middle point
% r_m_0=Q_sec(1,n_m);
% theta_m_0=0;
% % the present First-order expansion relationship between d(theta,R) and d(x,y)
% Con_M_2_0=[-sin(theta_m_0)*r_m_0 cos(theta_m_0); cos(theta_m_0)*r_m_0 sin(theta_m_0)];
% Con_M_2_2=[-sin(theta_m_0)/r_m_0 cos(theta_m_0)/r_m_0; cos(theta_m_0) sin(theta_m_0)];
% Proj_M_0((n_m-1)*12+[1:2],(n_m-1)*12+[1:2])=Con_M_2_0;
% %fix R: the distance to the z axis
% D_fix(2,n_m)=1;%fix the R
% D_fix(3,n_m)=1;%fix the Z
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 











