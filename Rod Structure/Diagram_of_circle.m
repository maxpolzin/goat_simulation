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






Duplication_M=[1,1,1]
mm=51
Q_sec=Q_sec_M1(:,:,mm)
X0=Q_sec(1,1);
disp_x=[0:5:X0,X0];% sampling displacement x
figure(1);
clf
hold on
Rod_ploting3(q_node1(:,:,:,mm),Duplication_M,3*L_co,2);
q_node=q_node1(:,:,:,mm);
for ii=[1 size(q_node,3)]
    t = linspace(0, 2*pi,36);
    r = 1;
    x = 3*L_co*cos(t);
    y = 3*L_co*sin(t);
    if ii==1
        P(1:3,1)=q_node(1:3,1,ii);
        RM=reshape(q_node(4:12,1,ii),3,3);
    else
        P(1:3,1)=q_node(1:3,end,ii);
        RM=reshape(q_node(4:12,end,ii),3,3);
    end
    TM=[RM,P;[0 0 0 1]];
    xyz_patch=TM*[0*x;x;y;1+0*x];
    patch(xyz_patch(1,:), xyz_patch(2,:),xyz_patch(3,:), [1 1 1]);
end
axis equal
l = light;
l.Color = [1 1 1];
l.Position = [1 0 1];
%%plotting wheels
P_wheel=Q_sec(1:3,N_e_0+1);
RM_wheel=reshape(Q_sec(4:12,N_e_0+1),3,3);
% Wheel_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M)
% plot3(Q_sec(1,1)*[1,-1],Q_sec(2,1)*[1,1],Q_sec(3,1)*[1,1],'r-',LineWidth=2)
% plot3(Q_sec(1,end)*[1,1],Q_sec(2,end)*[1,-1],Q_sec(3,end)*[1,1],'g-',LineWidth=2)
% c=colorbar
% caxis([0 350])
% c.Label.String = 'Energy Density (mJ/mm)';
% colorbar('off')
hold off
box off
grid off
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
axis([-600 600 -700 700 -30 30])
axis off

view([0 0 1])
% title(['\itL_{T1}\rm=',num2str(2*(X0-disp_x(mm)),'%.0f'),', \itL_{T2}\rm=',num2str(2*Q_sec_M1(2,N_e_1+1,mm),'%.0f')])
set(gca,'FontSize',16,'fontname','Time','linewidth',1)
set(gcf, 'Position', [100 100 750 560]*0.7);
set(gca, 'Position', [0  0  1 1])   %[x y width height]
%%
mm=51
Q_sec=Q_sec_M1(:,:,mm)
X0=Q_sec(1,1);
disp_x=[0:5:X0,X0];% sampling displacement x
figure(1);
clf
hold on


Rod_ploting3(q_node1(:,:,:,mm),Duplication_M,3*L_co,C_coord);
q_node=q_node1(:,:,:,mm);
for ii=[1 size(q_node,3)]
    t = linspace(0, 2*pi,36);
    r = 1;
    x = 3*L_co*cos(t);
    y = 3*L_co*sin(t);
    if ii==1
        P(1:3,1)=q_node(1:3,1,ii);
        RM=reshape(q_node(4:12,1,ii),3,3);
    else
        P(1:3,1)=q_node(1:3,end,ii);
        RM=reshape(q_node(4:12,end,ii),3,3);
    end
    TM=[RM,P;[0 0 0 1]];
    xyz_patch=TM*[0*x;x;y;1+0*x];
    patch(xyz_patch(1,:), xyz_patch(2,:),xyz_patch(3,:), [1 1 1]);
end
axis equal
l = light;
l.Color = [1 1 1];
l.Position = [1 0 1];
%%plotting wheels
P_wheel=Q_sec(1:3,N_e_0+1);
RM_wheel=reshape(Q_sec(4:12,N_e_0+1),3,3);
% Wheel_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M)
% plot3(Q_sec(1,1)*[1,-1],Q_sec(2,1)*[1,1],Q_sec(3,1)*[1,1],'r-',LineWidth=2)
% plot3(Q_sec(1,end)*[1,1],Q_sec(2,end)*[1,-1],Q_sec(3,end)*[1,1],'g-',LineWidth=2)
% c=colorbar
% caxis([0 350])
% c.Label.String = 'Energy Density (mJ/mm)';
% colorbar('off')
hold off
box off
grid off
xlabel('X')
ylabel('Y')
zlabel('Z')
% axis([-600 600 -700 700 -30 30])
axis off
axis equal
view([1 1 1])
% title(['\itL_{T1}\rm=',num2str(2*(X0-disp_x(mm)),'%.0f'),', \itL_{T2}\rm=',num2str(2*Q_sec_M1(2,N_e_1+1,mm),'%.0f')])
set(gca,'FontSize',16,'fontname','Time','linewidth',1)
set(gcf, 'Position', [100 100 750 560]*0.7);
set(gca, 'Position', [0.1  0.1  1 0.8])   %[x y width height]

%%
mm=51
nn=1

Duplication_M=[1 1 1]%;-1 1 1;-1 -1 1;1 -1 1;];
Q_sec=Q_sec_M2{mm}(:,:,nn);
figure(3);
clf
hold on

Rod_ploting3(q_node2{mm}(:,:,:,nn),Duplication_M,3*L_co,1)
q_node=q_node2{mm}(:,:,:,nn);

for ii=[1 size(q_node,3)]
    t = linspace(0, 2*pi,36);
    r = 1;
    x = 3*L_co*cos(t);
    y = 3*L_co*sin(t);
    if ii==1
        P(1:3,1)=q_node(1:3,1,ii);
        RM=reshape(q_node(4:12,1,ii),3,3);
    else
        P(1:3,1)=q_node(1:3,end,ii);
        RM=reshape(q_node(4:12,end,ii),3,3);
    end
    TM=[RM,P;[0 0 0 1]];
    xyz_patch=TM*[0*x;x;y;1+0*x];
    patch(xyz_patch(1,:), xyz_patch(2,:),xyz_patch(3,:),'w');
end


P_wheel=Q_sec(1:3,N_e_0+1);
RM_wheel=reshape(Q_sec(4:12,N_e_0+1),3,3);
% Wheel_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M)
hold off
box off
grid off
h=gca;
%                 set(h, 'ZDir', 'reverse');
l = light;
l.Color = [1 1 1];
l.Position = [1 0 1];

axis equal
axis([-500 500 -800 800 -500 50])
axis off
view([1 -1 1])
% title(['\itL_{T1}\rm=',num2str(2*(X0-disp_x(mm)),'%.0f'),', \itL_{T2}\rm=',num2str(2*Q_sec(2,N_e_1+1),'%.0f')])
xlabel('X')
ylabel('Y')
zlabel('Z')
set(gca,'FontSize',16,'fontname','Time','linewidth',1)
set(gcf, 'Position', [100 100 750 560]*0.7);
set(gca, 'Position', [0.0  0.0  1 1])   %[x y width height]
% set(gcf, 'Color', 'none');

%%

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
