%% Evaluate the stability of every bisym config by applying twisting
% try to free the twisting freedom to see if the structure will stablize at
% a twisted state
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
load("Rod_BiSym_Folding.mat","Q_sec_M3")% load the twist or flat configuration as the original configuration
load('Rod_Free_Twisting.mat')
Q_sec=[];
Q_sec_1=[];
Q_sec_2=[];
Fq_sec=[];
U_flag=[];


for mm=41
    
    for nn=2%:size(Q_sec_M2{mm}(:,:,:),3)
        pic_num_Xdisp=nn;
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
        D_fix(2,1)=0;%fix the Y disp
        D_fix(3,1)=1;%fix the  Z disp

        % the local x y z axis of start and end point in control
        % because the equalism relationship
        %D_fix(1:3,1)=[-1 -1 1].*D_fix(4:6,end);%x axis in control
        %D_fix(4:6,1)=[-1 -1 1].*D_fix(4:6,end);%x axis in control
        %D_fix(7:9,1)=[-1 -1 1].*D_fix(4:6,end);%x axis in control
        %D_fix(10:12,1)=[-1 -1 1].*D_fix(4:6,end);%x axis in control
        % sus
        Con_M_1=diag(repmat([-1 -1 1],1,4));
        % Proj_M_0(1:12,end-11:end)=0
        Proj_M_0(end-11:end,1:12)=Con_M_1;
        D_fix([4 5 6],1)=0;%free x axis
        D_fix(5,1)=0;% contrain to XZ plane
        D_fix([7 8 9],1)=0;%fix y axis
        D_fix([10 11 12],1)=0;%free z axis
        D_fix(11,1)=0;% contrain to XZ plane

        % end point in constrain
        D_fix(1,end)=1;%fix the X disp
        D_fix(2,end)=1;%control the Y disp
        D_fix(3,end)=1;%control the Z disp

        % local x y z aixs in control
        D_fix([4 5 6],end)=1;%free x axis
        D_fix(5,end)=1;% contrain to XZ plane
        D_fix([7 8 9],end)=1;%fix y axis
        D_fix([10 11 12],end)=1;%free z axis
        D_fix(11,end)=1;% contrain to XZ plane

        n_m=N_e_1/2+1;%the seq numb of middle point

        % the middle point position
        % based on anti sym condition
        % the distance of middle point to the z axis is constrained: x^2+y^2=C
        % introduce the R and theta to replace the x y disp of middle point
        r_m_0=Q_sec(1,n_m);
        theta_m_0=0;
        % the present First-order expansion relationship between d(theta,R) and d(x,y)
        Con_M_2_2=[-sin(theta_m_0)/r_m_0 cos(theta_m_0)/r_m_0; cos(theta_m_0) sin(theta_m_0)];
        Proj_M_0((n_m-1)*12+(1:2),(n_m-1)*12+(1:2))=Con_M_2_2;
        %fix R: the distance to the z axis
        D_fix(2,n_m)=1;%fix the R
        D_fix(3,n_m)=1;%fix the Z
        % %the angle between local x axis and the plane X-Y in control
        % D_fix(6,n_m)=1;

        %local y axis in the plane X-Y
        D_fix(9,n_m)=1;


        kk=0;
        flag1=1;
        flag2=1;
        L_step=0.2;
        Fq_sec=[];
        U_flag=[];

        %add distortion at the middle point
        Q_sec(6,n_m)=sin(-15*pi/180);
        if ismember(mm,30:36)
            Q_sec(6,n_m)=sin(-10*pi/180);
        elseif ismember(mm,37:39)
            Q_sec(6,n_m)=sin(-5*pi/180);
        elseif ismember(mm,40:44)
            Q_sec(6,n_m)=sin(-10*pi/180);
        elseif ismember(mm,45:51)
            Q_sec(6,n_m)=sin(-15*pi/180);
        end
        pic_num=1;
        while (kk<=50||flag1>0.000001*R0_ring)||(flag2&&(kk<5/L_step||kk<2000))
            kk=kk+1;
            [dFq_dq_sec,Fq_sec(:,kk,nn),q_node{mm}(:,:,:,kk,nn),U_flag(kk),E_elastic_tw{mm}(:,kk,nn)] = Jocob_rod_sec(Q_sec,N_e_1,L_e_1,N_node,Par_E,A);
            % update the boundary information
            if r_m_0<10e-4
                r_m=0;
                theta_m=theta_m_0;
                D_fix(1:2,n_m)=1;
            else
                D_fix(1,n_m)=0;
                D_fix(2,n_m)=1;
                r_m=sqrt(sum(Q_sec(1:2,n_m).^2));
                theta_m=atan2(Q_sec(2,n_m),Q_sec(1,n_m));
                % the present First-order expansion relationship between d(theta,R) and d(x,y)
                Con_M_2_2=[-sin(theta_m)/r_m cos(theta_m)/r_m; cos(theta_m) sin(theta_m)];
                Proj_M_0((n_m-1)*12+(1:2),(n_m-1)*12+(1:2))=Con_M_2_2;
            end
            B_index=find(D_fix==0);% find the fixed demensions
            Proj_M=Proj_M_0(:,B_index);% extract the projection matrix of free paramters
            F_tendon2{mm}(:,kk)=Fq_sec(:,kk);
            dFq_dq_sec_2=Proj_M'*dFq_dq_sec*Proj_M;
            dQ_sec_2=-pinv(dFq_dq_sec_2)*Proj_M'*Fq_sec(:,kk);
            dQ_sec=L_step*Proj_M*dQ_sec_2;
            dQ_sec=reshape(dQ_sec,12,[]);
            Q_sec=Q_sec+dQ_sec;
            % canonicalize the diplacement based on the varition of angle
            if r_m_0<10e-4
                Q_sec(1:2,n_m)=[0,0]';
            else
                d_theta=Con_M_2_2(1,:)*dQ_sec(1:2,n_m);
                theta_m=theta_m+d_theta;
                Q_sec(1:2,n_m)=r_m_0*[cos(theta_m);sin(theta_m)];
            end
            Q_sec_M3{mm}(:,:,kk,nn)=Q_sec;
            flag1=sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e_1+1);
            flag2=mean(abs(Proj_M'*Fq_sec(:,kk)));
            if kk<=10/L_step
                flag2=1;
            elseif (max(U_flag((kk-10/L_step):kk))-min(U_flag((kk-10/L_step):kk)))<0.1e7%energy variation
                flag2=0;
            end

            if kk==1||mod(kk,10)==0
                fprintf('mm=%u,nn=%u,kk=%u, flag1=%f \n',mm,nn,kk,flag1)
                figure(4);
                clf
                hold on
                Rod_ploting(q_node{mm}(:,:,:,kk,nn),Duplication_M,L_co,C_coord)
                plot3(Q_sec(1,[1,N_e_1+1]),Q_sec(2,[1,N_e_1+1]),Q_sec(3,[1,N_e_1+1]),'-','color','b','linewidth',1)
                plot3([1,-1]*Q_sec(1,n_m),[1,-1]*Q_sec(2,n_m),[1,1]*Q_sec(3,n_m),'-','color','r','linewidth',1)
                box on
                grid on
                xlabel('X')
                ylabel('Y')
                zlabel('Z')
                axis equal
                axis([-70 70 -70 70 -50 70])
                set(gca,'Zdir','reverse')
                view([1 1 1])
                hold off
                title(['X_1=',num2str(Q_sec(1,n_m),'%.0f'), ...
                    ', Y_2=', num2str(-Q_sec(2,1),'%.0f'), ...
                    ', Z_1=', num2str(Q_sec(3,n_m),'%.0f')],FontSize=18)
                drawnow;% 立即刷新当前绘图窗口，这是matlab绘图中动态展示的关键
%                 F = getframe(gcf);  % 获取当前绘图窗口的图片
%                 Im = frame2im(F);   % 返回与电影帧相关的图片数据
%                 [Amap, map] = rgb2ind(Im, 256); % 将RGB图片转化为索引图片
%                 if pic_num== 1
%                     imwrite(Amap, map,['Free_Twisting_mm_',num2str(mm),'_nn',num2str(nn),'.gif'],'gif', 'Loopcount', Inf, 'DelayTime', 0.1);
%                     % 将第一张图片写入‘sin.gif’文件中，并且将GIF播放次数设置成无穷，即保存的GIF图会一直动下去
%                 else
%                     imwrite(Amap, map,['Free_Twisting_mm_',num2str(mm),'_nn',num2str(nn),'.gif'],'gif','WriteMode','append','DelayTime',0.1);
%                     % 依次将其他的图片写入到GIF文件当中
%                 end
%                 pic_num = pic_num + 1;
%                 set(gcf, 'Position', [305 50 400 360]);
                X_1=reshape(Q_sec_M3{mm}(1,1,1:kk,nn),1,[]);
                Y_2=reshape(Q_sec_M3{mm}(2,end,1:kk,nn),1,[]);
                Z_1=reshape(Q_sec_M3{mm}(3,1,1:kk,nn),1,[]);
                figure(5);
                subplot(2,1,1)
                plotyy(1:kk,U_flag(1:kk),1:kk,E_elastic_tw{mm}(1,1:kk,nn))
                ylabel('U')
                subplot(2,1,2)
                plot(Y_2(1:kk))
                box on
                grid on
                ylabel('Y_2')

                set(gcf, 'Position', [0 50 300 750]);
                drawnow
            end

        end
        figure(4);
        clf
        hold on
        Rod_ploting(q_node{mm}(:,:,:,kk,nn),Duplication_M,L_co,C_coord)
        plot3(Q_sec(1,[1,N_e_1+1]),Q_sec(2,[1,N_e_1+1]),Q_sec(3,[1,N_e_1+1]),'-','color','b','linewidth',1)
        plot3([1,-1]*Q_sec(1,n_m),[1,-1]*Q_sec(2,n_m),[1,1]*Q_sec(3,n_m),'-','color','r','linewidth',1)
        box on
        grid on
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        axis equal
        axis([-70 70 -70 70 -50 70])
        set(gca,'Zdir','reverse')
        view([1 1 1])
        hold off
        title(['X_1=',num2str(Q_sec(1,n_m),'%.0f'), ...
            ', Y_2=', num2str(-Q_sec(2,1),'%.0f'), ...
            ', Z_1=', num2str(Q_sec(3,n_m),'%.0f')],FontSize=18)
        drawnow;% 立即刷新当前绘图窗口，这是matlab绘图中动态展示的关键
%         F = getframe(gcf);  % 获取当前绘图窗口的图片
%         Im = frame2im(F);   % 返回与电影帧相关的图片数据
%         [Amap, map] = rgb2ind(Im, 256); % 将RGB图片转化为索引图片
%         if pic_num_Xdisp == 1
%             imwrite(Amap, map,['Free_Twisting_Ydis','_mm',num2str(mm) '.gif'],'gif', 'Loopcount', Inf, 'DelayTime', 0.1);
%             % 将第一张图片写入‘sin.gif’文件中，并且将GIF播放次数设置成无穷，即保存的GIF图会一直动下去
%         else
%             imwrite(Amap, map,['Free_Twisting_Ydis','_mm',num2str(mm) '.gif'],'gif','WriteMode','append','DelayTime',0.1);
%             % 依次将其他的图片写入到GIF文件当中
%         end
        pic_num_Xdisp = pic_num_Xdisp + 1;
        set(gcf, 'Position', [305 50 400 360]);
    end
   



% figure(6);
% hold on
% plot3(nn*ones(jj-ii+1),ii:jj,E_elastic_tw{mm}(1,ii:jj,nn),'-*')
% hold off
% pause(0.1)
% view([1 1 1])
% box on
% grid on
save('Rod_Free_Twisting.mat','E_elastic_tw','Q_sec_M3')
end


figure(7)
clf
for mm=1:51
    nn=1
    V1=Q_sec_M3{mm}([4:6],1,:,nn)
    V1=reshape(V1,3,[])
    V2=Q_sec_M3{mm}([4:6],n_m,:,nn)
    V2=reshape(V2,3,[])

    Twist_angle1{mm}=atan2(V1(3,:),V1(1,:))*180/pi
    Twist_angle2{mm}=atan2(V1(3,:),sqrt(sum(V2(1:2,:).^2,1)))*180/pi

    % figure(7)
    % subplot(2,1,1)
    % plot(Twist_angle1{mm})
    % hold on
    % plot(Twist_angle2{mm})
    % subplot(2,1,2)    
    hold on
    plot3(mm+0*abs(Twist_angle2{mm}),abs(Twist_angle2{mm}),abs(Twist_angle1{mm}))
    labels = cellstr(num2str((50:-1:0)'));
    set(gcf, 'Position', [0 50 400 350]);
    hold off
    grid on
end
figure(7)
legend(labels)

%%
mm=41;
nn=5;
pic_num=1;
ii=24;
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

box on
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
axis([-60 60 -70 70 -60 60])
view([1 1 1])
drawnow;% 立即刷新当前绘图窗口，这是matlab绘图中动态展示的关键
    






