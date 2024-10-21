%% Evaluate the stability of every bisym config by applying twisting
% trying to analysis the multistability via applying a enforced displace
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
Duplication_M=[1 1 1;-1 -1 1];



%% Evaluate the stability of every anti-sym config by applying twisting

% double the section to generate the half ring configuration and make it as anti-sym
% the origin config
load("Rod_BiSym_Folding.mat","Q_sec_M2")% load the twist or flat configuration as the original configuration
% load("Rod_BiSym_Folding.mat","Q_sec_M3")% load the twist or flat configuration as the original configuration
load("Rod_Enforced_Twisting.mat")
Q_sec=[];
Q_sec_1=[];
Q_sec_2=[];
Fq_sec=[];
U=[];

mm=38
nn=30

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

%plot the doubled origin config
figure(4);
C_coord=[.0 .45 .74; .74 .45 .0];
L_co=2;
hold on
Rod_ploting(Q_sec,Duplication_M,L_co,C_coord)
hold off
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

% the local x y z axis of start and end point in control
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
%
n_m=N_e_1/2+1;%the seq numb of middle point

% the middle point position
% based on anti sym condition
% the distance of middle point to the z axis is constrained: x^2+y^2=C
% introduce the R and theta to replace the x y disp of middle point
r_m_0=Q_sec(1,n_m);
theta_m_0=0;
% the present First-order expansion relationship between d(theta,R) and d(x,y)
Con_M_2_2=[-sin(theta_m_0)/r_m_0 cos(theta_m_0)/r_m_0; cos(theta_m_0) sin(theta_m_0)];
Con_M_2_1=inv(Con_M_2_2);
Proj_M_0((n_m-1)*12+[1:2],(n_m-1)*12+[1:2])=Con_M_2_1;
%fix R: the distance to the z axis
% D_fix(1,n_m)=1;% theta
D_fix(2,n_m)=1;%fix the R
D_fix(3,n_m)=1;%fix the Z
%the angle between local x axis and the plane X-Y in control
% D_fix(4,n_m)=1;
D_fix(6,n_m)=1;
%local y axis in the plane X-Y
D_fix(9,n_m)=1;

%
T_angle=linspace(0,15*pi/180,31);

% Q_sec(:,1)
% Q_sec(:,end);
figure(4);
clf
pic_num =1
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
%     L_step=0.1;
   L_step_1=0.05
   T_angle_M{mm}(1:length(T_angle),nn)=T_angle;
    while (kk==0||flag1>0.0001*L_step||kk<2/L_step)&&kk<800
        kk=kk+1;
        L_step=0.01*0.99^kk+L_step_1*(1-0.99^kk);
        [dFq_dq_sec,Fq_sec(:,ii),q_node,U(ii),E_elastic_tw{mm}(:,ii,nn)] = Jocob_rod_sec(Q_sec,N_e_1,L_e_1,N_node,Par_E,A);
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
            Con_M_2_1=inv(Con_M_2_2);
            Proj_M_0((n_m-1)*12+[1:2],(n_m-1)*12+[1:2])=Con_M_2_1;
        end
        B_index=find(D_fix==0);% find the fixed demensions
        Proj_M=Proj_M_0(:,B_index);% extract the projection matrix of free paramters
        dFq_dq_sec_2=Proj_M'*dFq_dq_sec*Proj_M;
        dQ_sec_2=-pinv(dFq_dq_sec_2)*Proj_M'*Fq_sec(:,ii);
        dQ_sec=L_step*Proj_M*dQ_sec_2;
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
            subplot(2,1,1)
            hold on
            Rod_ploting(q_node,Duplication_M,L_co,C_coord)
            hold off
            axis equal
            box on
            grid on
            view([1 1 1])
            subplot(2,1,2)
            hold on
            plot(E_elastic_tw{mm}(1,1:ii,nn),'-*')
            hold off            
            set(gcf, 'Position', [0 50 300 750]);
        end

    end


    figure(4);
    clf
    subplot(2,1,1)
    plot(T_angle(1:ii)*180/pi,E_elastic_tw{mm}(1,1:ii,nn),'-*')
    hold on
    plot(T_angle(ii)*180/pi,E_elastic_tw{mm}(1,ii,nn),'o')
    hold off
    ylabel('U')
    xlabel("Enforced tilt angle at middle point (°)")
    subplot(2,1,2)
    hold on
    Rod_ploting(q_node,Duplication_M,L_co,C_coord)    
    plot3(Q_sec(1,[1,N_e_1+1]),Q_sec(2,[1,N_e_1+1]),Q_sec(3,[1,N_e_1+1]),'-','color','b','linewidth',1)
    plot3([1,-1]*Q_sec(1,n_m),[1,-1]*Q_sec(2,n_m),[1,1]*Q_sec(3,n_m),'-','color','r','linewidth',1)
    title(['X_1=',num2str(Q_sec(1,n_m),'%.0f'), ...
                    ', Y_2=', num2str(-Q_sec(2,1),'%.0f'), ...
                    ', Z_1=', num2str(Q_sec(3,n_m),'%.0f')],FontSize=18)
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
    set(gcf, 'Position', [305 50 400 750]);
        drawnow;% 立即刷新当前绘图窗口，这是matlab绘图中动态展示的关键
    F = getframe(gcf);  % 获取当前绘图窗口的图片
    Im = frame2im(F);   % 返回与电影帧相关的图片数据
    [Amap, map] = rgb2ind(Im, 256); % 将RGB图片转化为索引图片
    if pic_num == 1
        imwrite(Amap, map,['Enforced_Twisting_mm_',num2str(mm),'_nn',num2str(nn),'.gif'],'gif', 'Loopcount', Inf, 'DelayTime', 0.1);
        % 将第一张图片写入‘sin.gif’文件中，并且将GIF播放次数设置成无穷，即保存的GIF图会一直动下去
    else
        imwrite(Amap, map,['Enforced_Twisting_mm_',num2str(mm),'_nn',num2str(nn),'.gif'],'gif','WriteMode','append','DelayTime',0.1);
        % 依次将其他的图片写入到GIF文件当中
    end
    pic_num = pic_num + 1;
    pause(0.1)

end

ii=1
jj=length(T_angle)
figure(6);
hold on
plot3(nn*ones(jj-ii+1),ii:jj,E_elastic_tw{mm}(1,ii:jj,nn),'-*')
% ylim([1,31])
hold off
pause(0.1)
view([1 1 1])
box on
grid on

save('Rod_Enforced_Twisting.mat','E_elastic_tw','Q_sec_M3','-append')
 
%%
mm=41
nn=5
pic_num=1
ii=24
figure(4);
clf
n_last=size(Q_sec_M3{mm}(:,:,:,nn),3);
Q_sec=Q_sec_M3{mm}(:,:,n_last,nn);

for kk=1:(size(Q_sec,2)-1)
    hold on
    q_node=Curve_interp(Q_sec(:,kk),Q_sec(:,kk+1),L_e_0,4);
    Rod_ploting(q_node,Duplication_M,L_co,C_coord)
    hold off
end
% hold on
% l = light;
% l.Color = [1 1 1];
% l.Position = [1 0 -1];
box on
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
axis([-60 60 -70 70 -60 60])
view([1 1 1])
drawnow;% 立即刷新当前绘图窗口，这是matlab绘图中动态展示的关键
    






