
%% caculate bisymmetric outplane configurations in the range inside the boundary based on inplane deformation
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
L_co=2;%the length of the coordinates
C_coord=[.0 .45 .74; .74 .45 .0];%the color of the coordinates

Duplication_M=[1 1 1;-1 1 1;-1 -1 1;1 -1 1;];
%% caculate bisymmetric configurations in the range inside the boundary
load('Rod_BiSym_Folding.mat',"Q_sec_M1")
load("Rod_BiSym_Folding.mat","Q_sec_M2")
A=0.02e8;
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
D_fix=0*Q_sec_M1(:,:,1);%fixed dimensions/varibales
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
for mm=8%N_x+1  
    figure(3);
    clf
    Q_sec=Q_sec_M1(:,:,mm);
    % caculate the range of disp y based on the boundary 
    disp_y_0=Q_sec(2,end);
    N_y=fix((Q_sec(2,end)-Q_sec(1,1))/1)
    disp_y=flip([Q_sec(1,1):1:disp_y_0, disp_y_0]);% sampling displacement y
    for n_last=1:length(disp_y)
        mm
        n_last
        if n_last~=1
            Q_sec=Q_sec_M2{mm}(:,:,n_last-1);
        end
        Q_sec(2,end)=disp_y(n_last);
        if ismember(n_last,[1:20])
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
            Q_sec_M2{mm}(:,:,n_last)=Q_sec;
            flag1=sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e_0+1);
            flag2=mean(abs(Proj_M'*Fq_sec(:,mm)));
            if kk==5000
                flag3=1;
            end
        end
        if mod(n_last,5)==1
            figure(3);
            hold on
            Rod_ploting(q_node,Duplication_M,L_co,C_coord)
            hold off
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
pic_num=1
for mm=1:51
    figure(4);
    clf
    n_last=size(Q_sec_M2{mm}(:,:,:),3);
    Q_sec=Q_sec_M2{mm}(:,:,n_last);
    hold on
    for nn=1:(size(Q_sec,2)-1)
        q_node=Curve_interp(Q_sec(:,nn),Q_sec(:,nn+1),L_e_0,4);
        Rod_ploting(q_node,Duplication_M,L_co,C_coord)
    end
    hold off
    box on
    grid on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    axis equal
    axis([-60 60 -60 60 -60 60])
    view([1 1 1])
    drawnow;% 立即刷新当前绘图窗口，这是matlab绘图中动态展示的关键
    F = getframe(gcf);  % 获取当前绘图窗口的图片
    Im = frame2im(F);   % 返回与电影帧相关的图片数据
    [Amap, map] = rgb2ind(Im, 256); % 将RGB图片转化为索引图片
    if pic_num == 1
        imwrite(Amap, map, 'Outplane_1.gif', 'gif', 'Loopcount', Inf, 'DelayTime', 0.1);
        % 将第一张图片写入‘sin.gif’文件中，并且将GIF播放次数设置成无穷，即保存的GIF图会一直动下去
    else
        imwrite(Amap, map,'Outplane_1.gif','gif','WriteMode','append','DelayTime',0.1);
        % 依次将其他的图片写入到GIF文件当中
    end
    pic_num = pic_num + 1;
    pause(0.1)
end













