%% Evaluate the stability of every bisym config by applying twisting
% trying to analysis the multistability via applying a enforced displace
%% start the simulation
clear
clc
close all
%%
% input mechinical information
R0_ring=50;%radius of the ring cm
L_sec_0=R0_ring*pi/2;% set as the 1/4circle 
R_rod=5;%radius of the rod cm
E=20e6;%elastic modulus 10^4 Pa
G=E/2*(1+0.4);%shearing modulus
Par_E=[E*pi*R_rod^2, E*pi*R_rod^4/64 0.2*G*pi*R_rod^4/32]';% the axis 
% bending and torsion stiffness
% unit of evergy N*cm=10^-2 J
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
load("Rod_BiSym_Folding.mat","Q_sec_M2")% load the Bisys_folding configuration as the original configuration
load('Rod_Enforced_Twisting_1.mat')
Q_sec=[];
Q_sec_1=[];
Q_sec_2=[];
Fq_sec=[];
U_flag=[];

for mm=41
    for nn=20%:length(Q_sec_M2{mm}(1,1,:));

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
        Con_M_1=diag(repmat([-1 -1 1],1,4));
        % Proj_M_0(1:12,end-11:end)=0
        Proj_M_0(end-11:end,1:12)=Con_M_1;

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
        D_fix(1,n_m)=1;%fix the R
        D_fix(2,n_m)=1;%fix the R
        D_fix(3,n_m)=1;%fix the Z
        %the angle between local x axis and the plane X-Y in control
        D_fix(4,n_m)=1;
        D_fix(6,n_m)=1;
        %local y axis in the plane X-Y
        D_fix(9,n_m)=1;

        %
        Num_smp=9
        T_angle{mm}=linspace(0,2*pi/180,Num_smp);

        % Q_sec(:,1)
        % Q_sec(:,end);
        figure(4);
        clf
        pic_num =1
        for ii=1:length(T_angle{mm})
            ii
            figure(3);
            clf
            if ii~=1
                if nn>=40&mm>=42
                    Q_sec=Q_sec_M3_tw{mm}(:,:,1,nn);
                else
                    Q_sec=Q_sec_M3_tw{mm}(:,:,ii-1,nn);
                end
            end
            
            %redefine the middle points and frames
            Q_sec(6,n_m)=sin(T_angle{mm}(ii));

            kk=0;
            flag1=1;
            
            L_step_1=0.05;
            if nn>=54&&nn<56&mm==45
                L_step_1=0.1;
            end

            if nn>=56&mm==45
                L_step_1=0.1;
                if ii>=1
                    L_step_1=0.5/2;
                end
            end

            if nn>=42&mm==46
                L_step_1=0.125;            
            end

            if nn>=36&mm==47
                L_step_1=0.125;            
            end
            if mm==48
                L_step_1=0.25;
                if  nn>=36&&nn<=45
                    L_step_1=0.05;
                elseif nn>=45
                    L_step_1=0.25;
                end
            end
            
 
            if mm==49
                L_step_1=0.25;
                if  nn>=36&&nn<=43
                    L_step_1=0.1;
                elseif nn<=80
                    L_step_1=0.25;
                end
            end

            if mm==50
                L_step_1=0.25;
                if  nn>=33&&nn<=41
                    L_step_1=0.05;
                elseif nn>41
                    L_step_1=0.25;
                end
            end

            if mm==51
                L_step_1=0.25;
                if  nn>=33&&nn<=40
                    L_step_1=0.05;
                elseif nn>40
                    L_step_1=0.25;
                end
            end

            while (kk==0||kk<2/L_step_1)||flag2&&kk<2000
                flag2=1;
                kk=kk+1;
                L_step=0.0001*0.99^kk+L_step_1*(1-0.99^kk);
                [dFq_dq_sec,Fq_sec(:,ii),q_node{mm}(:,:,:,ii,nn),U_flag(kk),E_elastic_tw{mm}(:,ii,nn)] = Jocob_rod_sec(Q_sec,N_e_1,L_e_1,N_node,Par_E,A);
                % update the boundary information
                if r_m_0<10e-4
                    r_m=0;
                    theta_m=theta_m_0;
                    D_fix(1:2,n_m)=1;
                else
%                     D_fix(1,n_m)=0;
                    D_fix(1:3,n_m)=1;
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
                V_13q=2*Q_sec(1:3,n_m);
                if r_m_0<10e-4
                    Q_sec(1:2,n_m)=[0,0]';
                else
                    d_theta=Con_M_2_2(1,:)*dQ_sec([1:2],n_m);
                    theta_m=theta_m+d_theta;
                    Q_sec(1:2,n_m)=r_m_0*[cos(theta_m);sin(theta_m)];
                end
                Q_sec_M3_tw{mm}(:,:,ii,nn)=Q_sec;
                flag1=sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e_1+1);
                if mm<=11
%                     E_lim=1e1;
                    E_lim=1e2;
                elseif mm<=21
                    E_lim=4e1;
                elseif mm<=31
                    E_lim=1e2;
                elseif mm<=41
                    E_lim=1e2;
                elseif mm<=44
                    E_lim=2e2;
                elseif mm<=49
                    E_lim=2e2;
                else
                    E_lim=4e2;
                end
                
                if kk<=0.5/L_step_1
                    flag2=1;
                elseif (max(U_flag((kk-0.5/L_step_1):kk))-min(U_flag((kk-0.5/L_step_1):kk)))<E_lim%energy variation
                    flag2=0;
                end
                
                if kk==1||mod(kk,100)==0
                    fprintf('mm=%u,nn=%u,ii=%u, kk=%u, flag2=%f, D13=%f \n',mm,nn,ii,kk,flag2,norm(V_13q))
                    figure(3);
                    subplot(3,1,3)
                    hold on
                    Rod_ploting(q_node{mm}(:,:,:,ii,nn),Duplication_M,L_co,C_coord)
                    plot3(Q_sec(1,[1,end]),Q_sec(2,[1,end]),Q_sec(3,[1,end]),'-','color','b','linewidth',1)
                    plot3(Q_sec(1,n_m)*[1,-1],Q_sec(2,n_m)*[1,-1],Q_sec(3,n_m)*[1,1],'-','color','r','linewidth',1)
                    scatter3(Q_sec(1,1),Q_sec(2,1),Q_sec(3,1),'r')
                    scatter3(Q_sec(1,end),Q_sec(2,end),Q_sec(3,end),'g')
                    hold off
                    title(['X_1=',num2str(Q_sec(1,n_m),'%.0f'), ...
                        ', Y_2=', num2str(-Q_sec(2,1),'%.0f'), ...
                        ', Z_1=', num2str(Q_sec(3,n_m),'%.0f')],FontSize=18)
                    set(gca,'Zdir','reverse')
                    axis equal
                    box on
                    grid on
                    view([1 1 1])
                end

                if mod(kk,2)==1
                    figure(3);
                    if kk==1
                        clf
                    end
                    subplot(3,1,1)
                    plot(1:kk,U_flag(1:kk))
                    hold on
                    plot(kk,E_elastic_tw{mm}(1,ii,nn),'-*')
                    hold off
                    ylabel('U_{flag}&U_E')

                    subplot(3,1,2)
                    if kk>0.5/L_step_1
                        plot((kk-0.5/L_step_1):kk,U_flag((kk-0.5/L_step_1):kk)-U_flag((kk-0.5/L_step_1)))
                    end
                    box on
                    grid on
                    ylabel('U_{flag}')
                    set(gcf, 'Position', [0 50 300 750]);
                    drawnow                    
                end

            end


            figure(4);
            clf
            subplot(2,1,1)
            plot(0:ii-1,E_elastic_tw{mm}(1,1:ii,nn),'-*')
            hold on
            plot(ii-1,E_elastic_tw{mm}(1,ii,nn),'o')
            hold off
            ylabel('U')
            xlabel("Enforced tilt angle at middle point (°)")
            subplot(2,1,2)
            hold on
            Rod_ploting(q_node{mm}(:,:,:,ii,nn),Duplication_M,L_co,C_coord)
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
        end


        [E_TW_min{mm}(nn), Num_min{mm}(nn)]=min(E_elastic_tw{mm}(1,:,nn));
        save('Rod_Enforced_Twisting_1_2.mat','E_elastic_tw','Q_sec_M3_tw','E_TW_min','Num_min','T_angle')
        

        
        N_excl=2;
        f3 = @(p,x) p(1)*x + p(2);   
        p3=[];
        for jj=1:Num_smp
            x=0:Num_smp-1;
            x(jj)=[];
            y=E_elastic_tw{mm}(1,:,nn);
            y(jj)=[];
           p3(jj,:) =  polyfit(x,y,1); % coefficients
        end
        Stiff=p3(:,1)
        [a b]=sort(abs(Stiff-mean(Stiff)));        
        x=0:Num_smp-1;
        x(b(end-N_excl+1:end))=[];
        y=E_elastic_tw{mm}(1,:,nn);
        y(b(end-N_excl+1:end))=[];
        p3 = polyfit(x,y,1);
        Stiff_3{mm}(nn)=p3(1);

        f2 = @(p,x) p(1)*x.^2 + p(2); %p(1)*x.^3 + p(2)*x.^2 + p(3);
        p0 = [1,1]; % initial guess
        lb = [-inf,-inf]; % lower bounds
        ub = [inf,inf]; % upper bounds
        options = optimoptions('lsqcurvefit','Display','off'); % options
        for jj=1:Num_smp
            x=0:Num_smp-1;
            x(jj)=[];
            y=E_elastic_tw{mm}(1,:,nn);
            y(jj)=[];
            p2(jj,:) = lsqcurvefit(f2,p0,x,y,lb,ub,options); % coefficients
        end
        Stiff=p2(:,1)
        [a b]=sort(abs(Stiff-mean(Stiff)));        
        x=0:Num_smp-1;
        x(b(end-N_excl+1:end))=[];
        y=E_elastic_tw{mm}(1,:,nn);
        y(b(end-N_excl+1:end))=[];            
        p2 = lsqcurvefit(f2,p0,x,y,lb,ub,options); % coefficients
        Stiff_2{mm}(nn)=p2(1);

        f1 = @(p,x) p(1)*x.^3 + p(2)*x.^2 + p(3);
        p0 = [1,1,1]; % initial guess
        lb = [-inf,-inf,-inf]; % lower bounds
        ub = [inf,inf,inf]; % upper bounds
        options = optimoptions('lsqcurvefit','Display','off'); % options
        p1= lsqcurvefit(f1,p0,x,y,lb,ub,options); % coefficients
        Stiff_1{mm}(nn)=p1(2); 

        save('Rod_Enforced_Twisting_1_2.mat','Stiff_1','Stiff_2','Stiff_3','-append')

        figure(5)        
        plot(0:Num_smp-1,E_elastic_tw{mm}(1,:,nn),'sk')
        hold on
        plot(x,y,'+k',Markersize=15)
        plot(0:Num_smp-1,f1(p1,0:Num_smp-1),'-or');
        plot(0:Num_smp-1,f2(p2,0:Num_smp-1),'-og');
        plot(0:Num_smp-1,f3(p3,0:Num_smp-1),'-ob');
        hold off
        set(gcf, 'Position', [705 50 400 350]);

        figure(6)
        subplot(3,1,1)
        plot(Stiff_1{mm}(1:nn))
        subplot(3,1,2)
        plot(Stiff_2{mm}(1:nn))        
        subplot(3,1,3)
        plot(Stiff_3{mm}(1:nn))
        set(gcf, 'Position', [1105 50 400 750]);
               
    end
    figure(7)
    %     subplot(3,1,1)
    %     hold on
    %     scatter3(ii+0*Stiff_1{ii},0:length(Stiff_1{ii})-1,Stiff_1{ii},[],Stiff_2{ii})
    %     hold off
    %     view([1 1 1])

    subplot(2,1,1)
    hold on
    plot3(51-(mm+0*Stiff_2{mm}),reshape(Q_sec_M2{mm}(2,end,:),1,[]),Stiff_2{mm})
    hold off
    view([1 0 0])
    axis square
    xlabel('Tendon 2')
    ylabel('Tendon 1')
    xlim([0 70])
    ylim([0 70])
    grid on    
    hold off
    subplot(2,1,2)
    hold on
    plot3(51-(mm+0*Stiff_2{mm}),[0:(length(Stiff_2{mm})-1)],Stiff_2{mm})
    hold off
    view([1 0 0])
    axis square
    xlabel('Tendon 2')
    ylabel('Tendon 1')
    xlim([0 70])
    ylim([0 70])
    grid on    
    hold off


    set(gcf, 'Position', [1105 50 400 750]);
end

%%
load('Rod_Enforced_Twisting_1_2.mat')
load("Rod_BiSym_Folding.mat","Q_sec_M2")
mm=51

for ii=1:mm
    figure(7)
    %  subplot(2,1,1)
    hold on
    plot3(51-(ii+0*Stiff_2{ii}),reshape(Q_sec_M2{ii}(2,end,:),1,[]),Stiff_2{ii})
    hold off
    view([0.1 0 1])
    axis square
    set(gcf, 'Position', [1105 500 400 300]);
%         subplot(2,1,2)
    figure(8)
        hold on
        scatter3(51-(ii+0*Stiff_2{ii}),reshape(Q_sec_M2{ii}(2,end,:),1,[]),Stiff_2{ii},[],1e-4*180/pi*Stiff_2{ii},'filled')
        scatter3(reshape(Q_sec_M2{ii}(2,end,:),1,[]),51-(ii+0*Stiff_2{ii}),Stiff_2{ii},[],1e-4*180/pi*Stiff_2{ii},'filled')
        hold off
        view([0 0 1])
   
% %         subplot(3,1,3)
%         hold on
%         scatter3(ii+0*Stiff_3{ii},0:length(Stiff_3{ii})-1,Stiff_3{ii},[],Stiff_2{ii})
%         view([1 1 1])
    hold off
    
end
% % B_stiff
% clear a b AT_stiff
% for ii=1:51
%      [a{ii} b{ii}]=find(Stiff_2{ii}<-70);
%      AT_stiff{ii,1}=51-(ii+0*Stiff_2{ii}(b{ii}));
%      AT_stiff{ii,2}=reshape(Q_sec_M2{ii}(2,end,b{ii}),1,[]);  
%      AT_stiff{ii,3}=Stiff_2{ii}(b{ii});
% end
% Color(1,:)=[0.0 0.45 0.74];
% Color(2,:)=[0.80,0.00,0.20];
% Color(3,:)=[0.93,0.69,0.13];
% figure(8)
% hold on
% for ii=1:51
%     scatter3(AT_stiff{ii,1},AT_stiff{ii,2},AT_stiff{ii,3},'o',MarkeredgeColor=Color(2,:),MarkerfaceColor=Color(2,:))%,LineWidth=2)
%     scatter3(AT_stiff{ii,2},AT_stiff{ii,1},AT_stiff{ii,3},'o',MarkeredgeColor=Color(2,:),MarkerfaceColor=Color(2,:))%,LineWidth=2)
% end
% hold off
% 
% clear a b AIn_stiff
% for ii=1:51
%      [a{ii} b{ii}]=find(reshape(Q_sec_M2{ii}(3,1,:)-Q_sec_M2{ii}(3,5,:),1,[])<(21*pi/4)&Stiff_2{ii}>-70);
%      AIn_stiff{ii,1}=51-(ii+0*Stiff_2{ii}(b{ii}));
%      AIn_stiff{ii,2}=reshape(Q_sec_M2{ii}(2,end,b{ii}),1,[]);  
%      AIn_stiff{ii,3}=Stiff_2{ii}(b{ii});
% end
% 
% figure(8)
% hold on
% for ii=1:51
%     scatter3(AIn_stiff{ii,1},AIn_stiff{ii,2},AIn_stiff{ii,3},'s',MarkeredgeColor=Color(1,:),MarkerfaceColor=Color(1,:))%,LineWidth=2)
%     scatter3(AIn_stiff{ii,2},AIn_stiff{ii,1},AIn_stiff{ii,3},'s',MarkeredgeColor=Color(1,:),MarkerfaceColor=Color(1,:))%,LineWidth=2)
% end
% hold off
% 
% clear a b AOut_stiff
% for ii=1:51
%      [a{ii} b{ii}]=find(reshape(Q_sec_M2{ii}(3,1,:)-Q_sec_M2{ii}(3,5,:),1,[])>(21*pi/4)&Stiff_2{ii}>-70);
%      AOut_stiff{ii,1}=51-(ii+0*Stiff_2{ii}(b{ii}));
%      AOut_stiff{ii,2}=reshape(Q_sec_M2{ii}(2,end,b{ii}),1,[]);  
%      AOut_stiff{ii,3}=Stiff_2{ii}(b{ii});
% end
% 
figure(8)
% hold on
% for ii=1:51
%     scatter3(AOut_stiff{ii,1},AOut_stiff{ii,2},AOut_stiff{ii,3},'s',MarkeredgeColor=Color(3,:),MarkerfaceColor=Color(3,:))%,LineWidth=2)
%     scatter3(AOut_stiff{ii,2},AOut_stiff{ii,1},AOut_stiff{ii,3},'s',MarkeredgeColor=Color(3,:),MarkerfaceColor=Color(3,:))%,LineWidth=2)
% end
% hold off
 axis equal
 axis('off')
    xlabel('Tendon 2')
    ylabel('Tendon 1')
    xlim([0 70])
    ylim([0 70])
    grid on
    c=colorbar 
    c.Label.String = 'Torsion stiffness (N·mm/rad)';
set(gcf, 'Position', [1105 50 400 400]);

%%
nn_L=68

f1 = @(p,x) p(1)*x.^3 + p(2)*x.^2 + p(3);
p0 = [1,1,1]; % initial guess
lb = [-inf,-inf,-inf]; % lower bounds
ub = [inf,inf,inf]; % upper bounds
options = optimoptions('lsqcurvefit','Display','off'); % options
for nn=1:nn_L
    p1 = lsqcurvefit(f1,p0,0:6,E_elastic_tw{mm}(1,:,nn),lb,ub,options); % coefficients
    Stiff_1{mm}(nn)=p1(2);
end

f2 = @(p,x) p(1)*x.^2 + p(2); %p(1)*x.^3 + p(2)*x.^2 + p(3);
p0 = [1,1]; % initial guess
lb = [-inf,-inf]; % lower bounds
ub = [inf,inf]; % upper bounds
options = optimoptions('lsqcurvefit','Display','off'); % options
for nn=1:nn_L
    p2 = lsqcurvefit(f2,p0,0:6,E_elastic_tw{mm}(1,:,nn),lb,ub,options); % coefficients
    Stiff_2{mm}(nn)=p2(1);
end

f3 = @(p,x) p(1)*x + p(2); %p(1)*x.^3 + p(2)*x.^2 + p(3);
for nn=1:nn_L
    p3 = polyfit(0:6,E_elastic_tw{mm}(1,:,nn),1);
    Stiff_3{mm}(nn)=p3(1);
end


figure(6)
subplot(3,1,1)
plot(Stiff_1{mm}(1:nn))
subplot(3,1,2)
plot(Stiff_1{mm}(1:nn))
subplot(3,1,3)
plot(Stiff_3{mm}(1:nn))
set(gcf, 'Position', [1105 50 400 750]);



%% %%
% pic_num = 1
% for ii=1:Num_min{mm}(nn)
%     figure(4);
%     clf
%     subplot(3,1,1)
%     plot(0:length(E_elastic_tw{mm}(1,:,nn))-1,1e-2*E_elastic_tw{mm}(1,:,nn),'-',LineWidth=1.5)
%     hold on
%     plot(ii-1,1e-2*E_elastic_tw{mm}(1,ii,nn),'o',MarkerFaceColor=[0.85 0.33 0.1])
%     hold off
%     ylabel('U_E [J]')
%     xlabel("Twisting angle to tendon2 (°)")
%     subplot(3,1,2)
%     hold on
%     Rod_ploting(q_node{mm}(:,:,:,ii,nn),Duplication_M,L_co,C_coord)
%     Line(1)=plot3(Q_sec(1,[1,N_e_1+1]),Q_sec(2,[1,N_e_1+1]),Q_sec(3,[1,N_e_1+1]),'-','color','b','linewidth',1);
%     Line(2)=plot3([1,-1]*Q_sec(1,n_m),[1,-1]*Q_sec(2,n_m),[1,1]*Q_sec(3,n_m),'-','color','r','linewidth',1);
%     title(['X_1=',num2str(Q_sec(1,n_m),'%.0f'), ...
%         ', Y_2=', num2str(-Q_sec(2,1),'%.0f'), ...
%         ', Z_1=', num2str(Q_sec(3,n_m),'%.0f')],FontSize=18)
%     legend(Line,{'Tendon 1','Tendon 2'},FontSize=12)
%     box on
%     grid on
%     xlabel('X')
%     ylabel('Y')
%     zlabel('Z')
%     axis equal
%     axis([-70 70 -70 70 -50 50])
%     set(gca,'Zdir','reverse')
%     view([1 1 1])
%     hold off
%     set(gcf, 'Position', [305 50 400 750]);
%     drawnow;% 立即刷新当前绘图窗口，这是matlab绘图中动态展示的关键
% %     F = getframe(gcf);  % 获取当前绘图窗口的图片
% %     Im = frame2im(F);   % 返回与电影帧相关的图片数据
% %     [Amap, map] = rgb2ind(Im, 256); % 将RGB图片转化为索引图片
% %     if pic_num == 1
% %         imwrite(Amap, map,['Enforced_Twisting_mm_',num2str(mm),'_nn',num2str(nn),'.gif'],'gif', 'Loopcount', Inf, 'DelayTime', 0.1);
% %         % 将第一张图片写入‘sin.gif’文件中，并且将GIF播放次数设置成无穷，即保存的GIF图会一直动下去
% %     else
% %         imwrite(Amap, map,['Enforced_Twisting_mm_',num2str(mm),'_nn',num2str(nn),'.gif'],'gif','WriteMode','append','DelayTime',0.1);
% %         % 依次将其他的图片写入到GIF文件当中
% %     end
% 
%     if nn==1
%         subplot(3,1,3)
%         Twist_angle=atan2(-Q_sec_M3_tw{mm}([6],1,:,nn),Q_sec_M3_tw{mm}([4],1,:,nn))*180/pi
%         Twist_angle=reshape(Twist_angle,1,[])
%         plot(Twist_angle,1e-2*E_elastic_tw{mm}(1,:,nn),'-',LineWidth=1.5)
%         hold on
%         plot(Twist_angle(ii),1e-2*E_elastic_tw{mm}(1,ii,nn),'o',MarkerFaceColor=[0.85 0.33 0.1])
%         hold off
%         ylabel('U_E [J]')
%         xlabel("Twisting angle to tendon1 (°)")
%         drawnow;% 立即刷新当前绘图窗口，这是matlab绘图中动态展示的关键
% %         F = getframe(gcf);  % 获取当前绘图窗口的图片
% %         Im = frame2im(F);   % 返回与电影帧相关的图片数据
% %         [Amap, map] = rgb2ind(Im, 256); % 将RGB图片转化为索引图片
% %         if pic_num == 1
% %             imwrite(Amap, map,['Enforced_Twisting_mm_',num2str(mm),'_inplane','.gif'],'gif', 'Loopcount', Inf, 'DelayTime', 0.1);
% %             % 将第一张图片写入‘sin.gif’文件中，并且将GIF播放次数设置成无穷，即保存的GIF图会一直动下去
% %         else
% %             imwrite(Amap, map,['Enforced_Twisting_mm_',num2str(mm),'_inplane','.gif'],'gif','WriteMode','append','DelayTime',0.1);
% %             % 依次将其他的图片写入到GIF文件当中
% %         end
%     end
%     
%     pic_num = pic_num + 1;
% end



 
%%
% mm=41
% nn=5
% pic_num=1
% ii=24
% figure(4);
% clf
% n_last=size(Q_sec_M3{mm}(:,:,:,nn),3);
% Q_sec=Q_sec_M3{mm}(:,:,n_last,nn);
% 
% for kk=1:(size(Q_sec,2)-1)
%     hold on
%     q_node=Curve_interp(Q_sec(:,kk),Q_sec(:,kk+1),L_e_0,4);
%     Rod_ploting(q_node,Duplication_M,L_co,C_coord)
%     hold off
% end
% 
% box on
% grid on
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% axis equal
% axis([-60 60 -70 70 -60 60])
% view([1 1 1])
% drawnow;% 立即刷新当前绘图窗口，这是matlab绘图中动态展示的关键
    






