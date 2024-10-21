% twisting ploting
clc
clear
close all
load('Rod_twist_2.mat')
for ii=1:size(U_M,2)
    for jj=1:size(U_M,1)
        TAngle2_M(jj,ii)=-atan2(Q_sec_MM{ii}(6,1,jj),norm(Q_sec_MM{ii}(4:5,1,jj)))*180/pi;        
    end
    [a,b]=min(U_M(:,ii));
    TAngle2_min(ii)=TAngle2_M(b,ii);
end







figure(4)
% for ii=1:51
%     U_M_2(:,nn)=U_M{ii};
%     TAngle_M_2(:,nn)=TAngle_M{ii};
%     disp_x_M_2(:,nn)=disp_x_M{ii};
% end
hold on
surf(TAngle2_M,100-disp_x_M*2,U_M*1e-4,'FaceAlpha',0.7)
% conversion U_e=E*R^4/L  (20e6*5^4/50)/(2e4*5^4/500)

shading interp
% contour3(TAngle_M,disp_x_M,U_M,15,'LineWidth',1,'LineStyle','--','LineColor',[0.5 0.5 0.5])
plot3(TAngle2_M(:,1:5:51),100-disp_x_M(:,1:5:51)*2,U_M(:,1:5:51)*1e-4,'LineWidth',1,'LineStyle','--','Color',[0.5 0.5 0.5])
colormap turbo
grid on
box on
figure (4)
hold on
plot3(TAngle2_M(1,:),100-disp_x_M(1,:)*2,U_M(1,:)*1e-4,'-','LineWidth',2,'Color',[0.5 0.5 0.5],'MarkerFaceColor','w')
plot3(TAngle2_min,100-disp_x_min*2,U_min*1e-4,'-','LineWidth',2,'Color',[0 0 0],'MarkerFaceColor','w')
grid on
set(gca,'YDir','reverse')
xlim([0 45])
zlim([900 inf])
xlabel('Twisting Angle(°)','Rotation', -5)
ylabel('L_{T1} (cm)','Rotation', 60)

zlabel('Elastic Energy(mJ)')
set(gca,'FontSize',16,'fontname','Times')
view([0.3,-1,0.5])
set(gcf, 'Position', [100 100 378 378].*[1 1 1.4 1.2]);
v = VideoWriter('Inplane_twisting_E_Surf.mp4', 'MPEG-4');
v.FrameRate = 30;  % Set the frame rate (frames per second)
v.Quality = 90;    % Set the video quality (0 to 100)
open(v);
for ii=1:length(TAngle2_min)
    hold on
    L1=plot3(TAngle2_min(ii),100-disp_x_min(ii)*2,U_min(ii)*1e-4,'o','LineWidth',2,'Color',[0 0 0],'MarkerFaceColor','w','MarkerSize',10)
    hold off
    F = getframe(gcf);  % 获取当前绘图窗口的图片
    writeVideo(v, F);
% close(v)
    pause(0.1)
    delete(L1)
end
close(v)
% set(gca, 'Position', [0.2  0.15  0.6  0.6])   %[x y width height]
%%
figure(2)
Color_L=linspecer(3)
Style_L={'--','-.','-'}

hold on
kk=0
for nn=[31 41 51]
    kk=kk+1
Line(kk)=plot3(TAngle2_M(1:5:end,nn),(50-disp_x_M(1:5:end,nn))*2,U_M(1:5:end,nn)*1e-4,'-','LineWidth',1,'Color',Color_L(kk,:))%
Line(kk).LineStyle=Style_L{kk};
plot3(TAngle2_min(nn),disp_x_min(nn),U_min(nn)*1e-4,'o-','LineWidth',2,'Color',Color_L(kk,:),'MarkerFaceColor','w')
end
hold off
view([0 -1 0])
% legend(Line,{'L_{T1}=40','L_{T1}=20','L_{T1}=0'},'Location','best')
xlabel('Twisting Angle(°)')
ylabel('L_{T1}  (cm)','Rotation', 60)
zlabel('Elastic Energy(mJ)')
set(gca,'FontSize',16,'fontname','Times')
set(gcf, 'Position', [100 100 378 378].*[1 1 0.8 0.8]);


%%

v = VideoWriter('Inplane_twisting_Conf.mp4', 'MPEG-4');
v.FrameRate = 30;  % Set the frame rate (frames per second)
v.Quality = 90;    % Set the video quality (0 to 100)
open(v);
for mm=1:51
    figure(3);
    clf
    hold on
    [a b]=min(U_M(:,mm));
    Q_sec=Q_sec_MM{mm}(:,:,b);

    L_e=50*pi/2/(size(Q_sec,2)-1);
    Duplication_M=[1,1,1; -1,1,-1; -1,-1,1; 1,-1,-1];
    L_co=3;%the length of the coordinates

    for ii=1:size(Q_sec,2)-1
        q_node=Curve_interp(Q_sec(:,ii),Q_sec(:,ii+1),L_e,5);
        Rod_ploting3(q_node,Duplication_M,L_co,2)
    end
    plot3(Q_sec(1,1)*[1,-1],Q_sec(2,1)*[1,-1],Q_sec(3,1)*[1,-1],'-','color','r','linewidth',2)
    plot3(Q_sec(1,end)*[1,-1],Q_sec(2,end)*[1,-1],Q_sec(3,end)*[1,-1],'-','color','g','linewidth',2)
    % scatter3(Q_sec(1,1),Q_sec(2,1),Q_sec(3,1),'r')
    % scatter3(Q_sec(1,n_2q),Q_sec(2,n_2q),Q_sec(3,n_2q),'g')
    box on
    grid on
    xlabel('X (cm)')
    ylabel('Y (cm)')
    zlabel('Z (cm)')
    l = light;
    l.Color = [1 1 1];
    l.Position = [1 0 1];
    axis equal
    axis([-55 55 -80 80 -25 25])
    % axis off
    % set(gca,'Zdir','reverse')
    view([1 1 0.5])
    title(['D_1=',num2str(2*Q_sec(1,1),'%.0f'),', D_2=',num2str(2*Q_sec(2,end),'%.0f')])
    hold off
    set(gcf, 'Position', [305 50 400 360]);
    set(gca,'FontSize',16,'fontname','Times')
    set(gcf, 'Position', [100 100 378 378].*[1 1 1.44 1.2]);
      F = getframe(gcf);  % 获取当前绘图窗口的图片
    writeVideo(v, F);
    pause(0.1)
end
close(v)



































