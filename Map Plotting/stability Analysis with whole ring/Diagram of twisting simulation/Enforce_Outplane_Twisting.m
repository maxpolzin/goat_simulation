%% Evaluate the stability of every bisym config by applying twisting
% trying to analysis the multistability via applying a enforced displace
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
L_co=2;%the length of the coordinates
C_coord=[.0 .45 .74; .74 .45 .0];%the color of the coordinates
Duplication_M=[1 1 1;-1 -1 1];



% Evaluate the stability of every anti-sym config by applying twisting

% double the section to generate the half ring configuration and make it as anti-sym
% the origin config
load("Rod_BiSym_Folding.mat","Q_sec_M2")% load the twist or flat configuration as the original configuration
load("Rod_BiSym_Folding.mat","Q_sec_M3")% load the twist or flat configuration as the original configuration
Q_sec=[];
Q_sec_1=[];
Q_sec_2=[];
Fq_sec=[];
U=[];

mm=41
nn=5

%generate the half ring configuration
Q_sec_1=Q_sec_M2{mm}(:,:,nn);
N_e_1=2*(size(Q_sec_1,2)-1);% get the element number from the data
Q_sec_2(1:3,:)=diag([1,-1,1])*Q_sec_1(1:3,:);
Q_sec_2(4:6,:)=diag([-1,-1,-1])*diag([1,-1,1])*Q_sec_1(4:6,:);
Q_sec_2(7:9,:)=diag([1,1,1])*diag([1,-1,1])*Q_sec_1(7:9,:);
Q_sec_2(10:12,:)=diag([1,1,1])*diag([1,-1,1])*Q_sec_1(10:12,:);
Q_sec=[flip(Q_sec_2(:,2:end),2),Q_sec_1];

L_sec_1=L_sec_0*2;
L_e_1=L_sec_1/N_e_1;


 
%%
load("Rod_Enforced_Twisting.mat")% load the twist or flat configuration as the original configuration

fig=figure(1)
% set(gca,'Position',[0.5 0.2 0.8 0.8])
% ax = gca;
% ax.YAxis.Exponent = 0;
hold on
mm=38
nn=6
for jj=1:size(E_elastic_tw{mm}(1,:,nn),2)
    TAngle2(jj)=-atan2(Q_sec_M3{mm}(6,1,jj,nn),norm(Q_sec_M3{mm}(4:5,1,jj,nn)))*180/pi;
end
[a,b]=min(E_elastic_tw{mm}(1,:,nn));
TAngle2_min=TAngle2(b);
subplot(2,2,1)
plot(TAngle2(1:31),E_elastic_tw{mm}(1,1:31,nn)*1e-4,LineWidth=1.5)
[a b]=min(E_elastic_tw{mm}(1,1:31,nn))
hold on
plot(TAngle2(b),E_elastic_tw{mm}(1,b,nn)*1e-4,'ok')
hold off
% ytickformat('%.3e')
set(gca,'FontSize',10,'fontname','Times')

mm=38
nn=28
for jj=1:size(E_elastic_tw{mm}(1,:,nn),2)
    TAngle2(jj)=-atan2(Q_sec_M3{mm}(6,1,jj,nn),norm(Q_sec_M3{mm}(4:5,1,jj,nn)))*180/pi;
end
[a,b]=min(E_elastic_tw{mm}(1,:,nn));
TAngle2_min=TAngle2(b);
subplot(2,2,2)
plot(TAngle2(1:30),E_elastic_tw{mm}(1,1:30,nn)*1e-4,LineWidth=1.5)
[a b]=min(E_elastic_tw{mm}(1,1:31,nn))
hold on
plot(TAngle2(b),E_elastic_tw{mm}(1,b,nn)*1e-4,'ok')
hold off
% ytickformat('%.3e')
set(gca,'FontSize',10,'fontname','Times')

mm=38
nn=30
for jj=1:size(E_elastic_tw{mm}(1,:,nn),2)
    TAngle2(jj)=-atan2(Q_sec_M3{mm}(6,1,jj,nn),norm(Q_sec_M3{mm}(4:5,1,jj,nn)))*180/pi;
end
[a,b]=min(E_elastic_tw{mm}(1,:,nn));
TAngle2_min=TAngle2(b);
subplot(2,2,3)
plot(TAngle2(1:31),E_elastic_tw{mm}(1,1:31,nn)*1e-4,LineWidth=1.5)
[a b]=min(E_elastic_tw{mm}(1,1:31,nn))
hold on
plot(TAngle2(b),E_elastic_tw{mm}(1,b,nn)*1e-4,'ok')
hold off
% ytickformat('%.3e')
set(gca,'FontSize',10,'fontname','Times')

mm=38
nn=32
for jj=1:size(E_elastic_tw{mm}(1,:,nn),2)
    TAngle2(jj)=-atan2(Q_sec_M3{mm}(6,1,jj,nn),norm(Q_sec_M3{mm}(4:5,1,jj,nn)))*180/pi;
end

T1=Q_sec_M3{mm}(1,9,1,nn)*2
T2=Q_sec_M3{mm}(2,end,1,nn)*2
[a,b]=min(E_elastic_tw{mm}(1,:,nn));
TAngle2_min=TAngle2(b);
subplot(2,2,4)
plot(TAngle2(1:31),E_elastic_tw{mm}(1,1:31,nn)*1e-4,LineWidth=1.5)
[a b]=min(E_elastic_tw{mm}(1,1:31,nn))
hold on
plot(TAngle2(b),E_elastic_tw{mm}(1,b,nn)*1e-4,'ok')
hold off
% ytickformat('%.3e')
hold off
set(gca,'FontSize',10,'fontname','Times')
han = axes(fig, 'visible', 'off');
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';

ylabel(han, 'Elastic energy (mJ)','FontSize',12,'fontname','Times');
xlabel(han, 'Twisting angle (°)','FontSize',12,'fontname','Times');
set(gca,'Position',[0.12 0.12 0.8 0.8])

% title(han, '共用标题');

%%


figure(3)

mm=38
nn=32

v = VideoWriter('Outplane_twisting_Conf_m38_n32.mp4', 'MPEG-4');
v.FrameRate = 30;  % Set the frame rate (frames per second)
v.Quality = 90;    % Set the video quality (0 to 100)
open(v);

for jj=1:size(E_elastic_tw{mm}(1,:,nn),2)
    TAngle2(jj)=-atan2(Q_sec_M3{mm}(6,1,jj,nn),norm(Q_sec_M3{mm}(4:5,1,jj,nn)))*180/pi;
end

T1=Q_sec_M3{mm}(1,9,1,nn)*2;
T2=Q_sec_M3{mm}(2,end,1,nn)*2;
Z1=Q_sec_M3{mm}(3,9,1,nn);
[a,b]=min(E_elastic_tw{mm}(1,:,nn));
TAngle2_min=TAngle2(b);

subplot(2,1,1)

plot(TAngle2(1:30),E_elastic_tw{mm}(1,1:30,nn)*1e-4,LineWidth=1.5)
[a b]=min(E_elastic_tw{mm}(1,1:31,nn))
hold on
plot(TAngle2(b),E_elastic_tw{mm}(1,b,nn)*1e-4,'ok')
hold off
ylabel('Elastic energy (mJ)','FontSize',12,'fontname','Times');
xlabel( 'Twisting angle (°)','FontSize',12,'fontname','Times');
title(['D_1=',num2str(T1,'%.0f'),', D_2=',num2str(T2),', Z_1=',num2str(Z1,'%.1f')])




for ii=[1:b b*ones(1,40-b)]
    figure(3);
    subplot(2,1,1)
    hold on
    L1=plot(TAngle2(ii),E_elastic_tw{mm}(1,ii,nn)*1e-4,'or',LineWidth=1.5)
    hold off
    subplot(2,1,2)
    
    n_last=size(Q_sec_M3{mm}(:,:,:,nn),3);
    n_select=ii;%b
    % n_select=b
    Q_sec=Q_sec_M3{mm}(:,:,n_select,nn);
    plot3(Q_sec(1,[1,N_e_1+1]),Q_sec(2,[1,N_e_1+1]),Q_sec(3,[1,N_e_1+1]),'-','color','g','linewidth',1)
    hold on
    plot3(Q_sec(1,N_e_1/2+1)*[1,-1],Q_sec(2,N_e_1/2+1)*[1,-1],Q_sec(3,N_e_1/2+1)*[1,1],'-','color','r','linewidth',1)
    hold off
    hold on
    for kk=1:(size(Q_sec,2)-1)
        q_node=Curve_interp(Q_sec(:,kk),Q_sec(:,kk+1),L_e_0,4);
        Rod_ploting3(q_node,Duplication_M,L_co,2)
    end
    hold off 
    l = light;
    l.Color = [1 1 1];
    l.Position = [1 0 -1];
    box on
    grid on

    xlabel('X (cm)')
    ylabel('Y (cm)')
    zlabel('Z (cm)')
    set(gca,'ZDir','reverse')
    axis equal
    ax=[-50 50 -70 70 -20 60];
%     ax=[-40 40 -50 50 -20 60];
    axis(ax)
%     axis off
    % for ii=1:3
    % %     xyzt(1:2,ii)=fix(ax([2*(ii-1)+1 2*(ii-1)+2])*4/pi/5)*5/(4/pi)
    %     xyzt(1:2,ii)=ax([2*(ii-1)+1 2*(ii-1)+2])
    %     xyzt(1:3,ii)=[xyzt(1,ii) 0 xyzt(2,ii)]
    %     if ii==1
    %         xticks(xyzt(:,ii))
    %     elseif ii==2
    %         yticks(xyzt(:,ii))
    %     elseif ii==3
    %         zticks(xyzt(:,ii))
    %     end
    %     for jj=1:length(xyzt)
    % %         XYZlb{jj,ii}=num2str(xyzt(jj,ii)*4/pi)
    %         XYZlb{jj,ii}=num2str(xyzt(jj,ii))
    %     end
    %     if ii==1
    %         xticklabels(XYZlb(:,ii))
    %     elseif ii==2
    %         yticklabels(XYZlb(:,ii))
    %     elseif ii==3
    %         zticklabels(XYZlb(:,ii))
    %     end
    % end
    view([1 0.5 0.5])
    drawnow;% 立即刷新当前绘图窗口，这是matlab绘图中动态展示的关键
    set(gcf, 'Position', [305 50 400 360]);
    set(gca,'FontSize',10,'fontname','Times')
    set(gcf, 'Position', [100 100 378 378].*[1 1 0.7 1.2]);
    pause(0.1)
    F = getframe(gcf);  % 获取当前绘图窗口的图片
    writeVideo(v, F);    
    delete(L1)
end
close(v)




