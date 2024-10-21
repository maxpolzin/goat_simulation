%%

clc
clear
close all
%%
load('Band_inplane3.mat')
load('Band_outplane3.mat')


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
Duplication_M=[1 1 1];
C_coord=[.0 .45 .74; .74 .45 .0];%the color of the coordinates
R_wheel=120;% The radius of the wheel
N_spoke=12;% The spoke num of the wheel


mm=51
Q_sec=Q_sec_M1(:,:,mm)
Q_sec_rd1=Q_sec(:,1:(N_e_1+1));
Q_sec_rd2=Q_sec(:,(N_e_1+2):(2*(N_e_1+1)));
X0=Q_sec_rd1(1,1);
disp_x=[0:5:X0,X0];% sampling displacement x
figure(3);
% clf
hold on
Rod_ploting3(q_node1_rd1(:,:,:,mm),Duplication_M,3*L_co,2);
Rod_ploting3(q_node1_rd2(:,:,:,mm),Duplication_M,3*L_co,2);



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
% Duplication_M=[1 1 1;-1 1 1;1 -1 1;-1 -1 1]
Wheel_Connect_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M,1)
% plot3(P_wheel(1,6)*[1,-1],P_wheel(2,6)*[1,1],P_wheel(3,6)*[1,1],'r-','linewidth',2)
% plot3(P_wheel(1,9)*[1,1],P_wheel(2,9)*[1,-1],P_wheel(3,9)*[1,1],'g-','linewidth',2)
% l = light;
% l.Color = [1 1 1];
% l.Position = [1 0 1];
hold off
box off
grid off
xlabel('X')
ylabel('Y')
zlabel('Z')
% set(gca,'Xdir','reverse')
% set(gca,'Ydir','reverse')
axis off
axis equal
axis([-600 600 -700 700 -150 150])
view([1 -1 1])
set(gca,'FontSize',16,'fontname','Time','linewidth',1)
% set(gcf, 'Position', [100 100 750 560]*0.7);
set(gca, 'Position', [0  0  1 1])   %[x y width height]

%%
Duplication_M=[1 1 1];
mm=52
nn=26%81
Q_sec=Q_sec_M2{mm}(:,:,nn);
Q_sec_rd1=Q_sec(:,1:(N_e_1+1));
Q_sec_rd2=Q_sec(:,(N_e_1+2):(2*(N_e_1+1)));
figure(6);
% clf
hold on
Rod_ploting3(q_node2_rd1{mm}(:,:,:,nn),Duplication_M,3*L_co,2)
Rod_ploting3(q_node2_rd2{mm}(:,:,:,nn),Duplication_M,3*L_co,2)
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
% Duplication_M=[1 1 1;-1 1 1;1 -1 1;-1 -1 1]
Wheel_Connect_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M,1)

% plot3(P_wheel(1,6)*[1,-1],P_wheel(2,6)*[1,1],P_wheel(3,6)*[1,1],'r-','color','r','linewidth',2)
% plot3(P_wheel(1,9)*[1,1],P_wheel(2,9)*[1,-1],P_wheel(3,9)*[1,1],'-','color','g','linewidth',2)

%%plotting wheels
P_wheel=Q_sec(1:3,N_e_0+1);
RM_wheel=reshape(Q_sec(4:12,N_e_0+1),3,3);
% Wheel_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M)
hold off
box off
grid off
% l = light;
% l.Color = [1 1 1];
% l.Position = [1 0 1];

axis equal
axis([-310.7814  340.6787 -644.7385  674.7385 -322.3305  105.0000])

axis off
view([1 -1 1])
% title(['\itL_{T1}\rm=',num2str(2*(X0-disp_x(mm)+5),'%.0f'),', \itL_{T2}\rm=',num2str(2*Q_sec(2,N_e_1+1)+9,'%.0f')])
xlabel('X')
ylabel('Y')
zlabel('Z')
set(gca,'FontSize',16,'fontname','Time','linewidth',1)
% set(gcf, 'Position', [100 100 380 280]*0.7);
% set(gca, 'Position', [0.0  0.0  1 0.8])   %[x y width height]
% set(gca, 'Position', [0.0  0.1  0.65 0.9] )  %[x y width height]
% set(gca, 'Position', [0.0  0.0  0.7 0.9])   %[x y width height]








%%

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
Duplication_M=[1 1 1];
C_coord=[.0 .45 .74; .74 .45 .0];%the color of the coordinates
R_wheel=120;% The radius of the wheel
N_spoke=12;% The spoke num of the wheel

mm=51
Q_sec=Q_sec_M1(:,:,mm)
Q_sec_rd1=Q_sec(:,1:(N_e_1+1));
Q_sec_rd2=Q_sec(:,(N_e_1+2):(2*(N_e_1+1)));
X0=Q_sec_rd1(1,1);
disp_x=[0:5:X0,X0];% sampling displacement x
figure(3);
% clf
hold on
Rod_ploting3(q_node1_rd1(:,:,:,mm),Duplication_M,3*L_co,2);
Rod_ploting3(q_node1_rd2(:,:,:,mm),Duplication_M,3*L_co,2);
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
% Duplication_M=[1 1 1;-1 1 1;-1 -1 1;1 -1 1;];
Wheel_Connect_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M,1)
% plot3(P_wheel(1,6)*[1,-1],P_wheel(2,6)*[1,1],P_wheel(3,6)*[1,1],'r-','linewidth',2)
% plot3(P_wheel(1,9)*[1,1],P_wheel(2,9)*[1,-1],P_wheel(3,9)*[1,1],'g-','linewidth',2)
% l = light;
% l.Color = [1 1 1];
% l.Position = [1 0 1];
hold off
box off
grid off
xlabel('X')
ylabel('Y')
zlabel('Z')
% axis([-600 600 -700 700 -150 150])
axis off
axis equal
view([1 -1 1])

%%
%%
Duplication_M=[1 1 1];
mm=50
nn=41%81
Q_sec=Q_sec_M2{mm}(:,:,nn);
Q_sec_rd1=Q_sec(:,1:(N_e_1+1));
Q_sec_rd2=Q_sec(:,(N_e_1+2):(2*(N_e_1+1)));
figure(6);
% clf
hold on
Rod_ploting3(q_node2_rd1{mm}(:,:,:,nn),Duplication_M,3*L_co,2)
Rod_ploting3(q_node2_rd2{mm}(:,:,:,nn),Duplication_M,3*L_co,2)
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
% Duplication_M=[1 1 1;-1 1 1;1 -1 1;-1 -1 1]
Wheel_Connect_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M,1)

% plot3(P_wheel(1,6)*[1,-1],P_wheel(2,6)*[1,1],P_wheel(3,6)*[1,1],'r-','color','r','linewidth',2)
% plot3(P_wheel(1,9)*[1,1],P_wheel(2,9)*[1,-1],P_wheel(3,9)*[1,1],'-','color','g','linewidth',2)

%%plotting wheels
P_wheel=Q_sec(1:3,N_e_0+1);
RM_wheel=reshape(Q_sec(4:12,N_e_0+1),3,3);
% Wheel_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M)
hold off
box off
grid off
% l = light;
% l.Color = [1 1 1];
% l.Position = [1 0 1];

axis equal
axis([ -302.3881  332.3778 -631.0813  661.0761 -416.6656  104.3675])

axis off
view([1 -1 1])
% title(['\itL_{T1}\rm=',num2str(2*(X0-disp_x(mm)+5),'%.0f'),', \itL_{T2}\rm=',num2str(2*Q_sec(2,N_e_1+1)+9,'%.0f')])
xlabel('X')
ylabel('Y')
zlabel('Z')




%%
mm=50
nn=79
Q_sec=Q_sec_M2{mm}(:,:,nn);
Q_sec_rd1=Q_sec(:,1:(N_e_1+1));
Q_sec_rd2=Q_sec(:,(N_e_1+2):(2*(N_e_1+1)));
figure(5);
clf
hold on
Rod_ploting3(q_node2_rd1{mm}(:,:,:,nn),Duplication_M,L_co,1)
Rod_ploting3(q_node2_rd2{mm}(:,:,:,nn),Duplication_M,L_co,1)
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

axis equal
axis off
view([1 1 1])
title(['\itL_{T1}\rm=',num2str(2*(X0-disp_x(mm)+4.5),'%.0f'),', \itL_{T2}\rm=',num2str(2*Q_sec(2,N_e_1+1)+9,'%.0f')])
xlabel('X')
ylabel('Y')
zlabel('Z')


%%
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
% axis([-250 250 -250 250 -500 10])
axis equal
axis off
view([1 1 1])
% title(['\itL_{T1}\rm=',num2str(2*(X0-disp_x(mm)+4.5),'%.0f'),', \itL_{T2}\rm=',num2str(2*Q_sec(2,N_e_1+1),'%.0f')])
xlabel('X')
ylabel('Y')
zlabel('Z')



