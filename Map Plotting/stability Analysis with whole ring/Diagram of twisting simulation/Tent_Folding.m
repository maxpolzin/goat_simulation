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
L_co=1.5;%the length of the coordinates
C_coord=[.0 .45 .74; .74 .45 .0];%the color of the coordinates
Duplication_M=[1 1 1;-1 -1 1];



%% build the whole rod ring to caculate the twisting folding

load("Rod_BiSym_Folding.mat","Q_sec_M3")% load the twist or flat configuration as the original configuration
Q_sec=[];
Q_sec_1=[];
Q_sec_2=[];
Fq_sec=[];
U=[];

mm=41
nn=5

v = VideoWriter('Tent Folding_3.mp4', 'MPEG-4');
v.FrameRate = 30;  % Set the frame rate (frames per second)
v.Quality = 90;    % Set the video quality (0 to 100)
open(v);


for ll=[1:2:24 24]
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

    r_node_1=Q_sec(1:3,:);
    dr_dx_1=Q_sec(4:6,:);
    dr_dy_1=Q_sec(7:9,:);
    dr_dz_1=Q_sec(10:12,:);

    Duplication_M=[1 1 1;-1,-1,1]

    %plot the origin config
    figure(7);
    plot3(r_node_1(1,[1,n_2q]),r_node_1(2,[1,n_2q]),r_node_1(3,[1,n_2q]),'-','color','g','linewidth',1)
    hold on
    plot3(r_node_1(1,[n_1q,n_3q]),r_node_1(2,[n_1q,n_3q]),r_node_1(3,[n_1q,n_3q]),'-','color','r','linewidth',1)
    hold off
    hold on
    for ii=1:size(Q_sec,2)-1
        q_node = Curve_interp(Q_sec(:,ii),Q_sec(:,ii+1),L_e_2,5);
        Rod_ploting3(q_node,Duplication_M,2*L_co,2)
    end    
    scatter3(r_node_1(1,1),r_node_1(2,1),r_node_1(3,1),'r')
    scatter3(r_node_1(1,n_2q),r_node_1(2,n_2q),r_node_1(3,n_2q),'g')
    l = light;
    l.Color = [1 1 1];
    l.Position = [1 0 -1];
    box on
    grid on
    xlabel('X (cm)')
    ylabel('Y (cm)')
    zlabel('Z (cm)')
    V_13q=Q_sec(1:3,n_3q)-Q_sec(1:3,n_1q);
     title(['D_1=',num2str(norm(V_13q),'%.0f'), ...
         ', D_2=', num2str(-2*Q_sec(2,1),'%.0f')],FontSize=18)
     set(gca,'FontSize',16,'fontname','Times')
    set(gca,'Zdir','reverse')
%     set(gca,'Ydir','reverse')
    axis equal
    % axis off
    axis([-50 50 -70 70 -30 70])
    view([1 1 1])
    hold off
    pause(0.1)
    F = getframe(gcf);  % 获取当前绘图窗口的图片
    writeVideo(v, F);  
end
%%

load('Rod_Free_Folding.mat')
L_co=1.5
Duplication_M=[1,1,1]
for ii=[1:1:107 107]
%     q_node=q_node4{mm}(:,:,:,ii,ll,nn)
    Q_sec=Q_sec_M4{mm}(:,:,ii,ll,nn)
    

    figure(11);
    %  clf
    plot3(Q_sec(1,[1,n_2q]),Q_sec(2,[1,n_2q]),Q_sec(3,[1,n_2q]),'-','color','g','linewidth',2)    
    hold on
    plot3(Q_sec(1,[n_1q,n_3q]),Q_sec(2,[n_1q,n_3q]),Q_sec(3,[n_1q,n_3q]),'-','color','r','linewidth',2)
    %     Rod_ploting3(q_node,Duplication_M,2*L_co,2)
    for ii=1:size(Q_sec,2)-1
        q_node = Curve_interp(Q_sec(:,ii),Q_sec(:,ii+1),L_e_2,5);
        Rod_ploting3(q_node,Duplication_M,2*L_co,2)
    end    
    scatter3(Q_sec(1,1),Q_sec(2,1),Q_sec(3,1),'r')
    scatter3(Q_sec(1,n_2q),Q_sec(2,n_2q),Q_sec(3,n_2q),'g')
    box on
    grid on
    xlabel('X (cm)')
    ylabel('Y (cm)')
    zlabel('Z (cm)')
    l = light;
    l.Color = [1 1 1];
    l.Position = [1 0 -1];
    %     title(['X_1=',num2str(Q_sec(1,n_1q),'%.0f'), ...
    %             ', Y_2=', num2str(-Q_sec(2,1),'%.0f'), ...
    %             ', Z_1=', num2str(Q_sec(3,n_1q),'%.0f')],FontSize=18)
    V_13q=Q_sec(1:3,n_3q)-Q_sec(1:3,n_1q);
     title(['D_1=',num2str(norm(V_13q),'%.0f'), ...
         ', D_2=', num2str(-2*Q_sec(2,1),'%.0f')],FontSize=18)
    axis equal
    axis([-50 50 -70 70 -30 70])
%     axis off
    set(gca,'Zdir','reverse')
%     set(gca,'Ydir','reverse')
    set(gca,'FontSize',16,'fontname','Times')
    view([1 1 1])
    hold off
    pause(0.1)
    F = getframe(gcf);  % 获取当前绘图窗口的图片
    writeVideo(v, F);  
end
close(v)












