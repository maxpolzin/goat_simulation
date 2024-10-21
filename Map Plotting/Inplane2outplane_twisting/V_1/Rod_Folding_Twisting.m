% the script simulate the rod ring folding and twisting deformation 
% author: Qinghua Guan
% test remote pc
%% input mechinical information
% define the caculation parameters of one section
%the element number for caculation
% the length of the element
% the initial solution of the section
clear
clc

% load set its original state
load("Rod_twist_2.mat","Q_sec_MM","U_M")% load the twist or flat configuration as the original configuration
N_disp=48%the disp number of the state 1:51 related to 0:50
[a b]=min(U_M(:,N_disp))
Q_sec_1=Q_sec_MM{N_disp}(:,:,b)
N_e=2*(size(Q_sec_1,2)-1);% get the element number from the data
Q_sec_2(1:3,:)=diag([-1,1,-1])*Q_sec_1(1:3,:)
Q_sec_2(4:6,:)=diag([-1,-1,-1])*diag([-1,1,-1])*Q_sec_1(4:6,:)
Q_sec_2(7:9,:)=diag([1,1,1])*diag([-1,1,-1])*Q_sec_1(7:9,:)
Q_sec_2(10:12,:)=diag([-1,-1,-1])*diag([-1,1,-1])*Q_sec_1(10:12,:)
Q_sec=[Q_sec_1,flip(Q_sec_2(:,1:end-1),2)]

L_sec=100*pi/2;
L_sec_1=L_sec;
L_sec_2=L_sec;
L_e=L_sec_1/N_e;

R=5;
E=20e6
G=E/2*(1+0.4)
[Par_E]=[E*pi*R^2, E*pi*R^4/64 0.2*G*pi*R^4/32]';
A=0.04e8;%can be to high or too low, depending on overall elatic energy level. the constrain parameter for 

% Q_sec=interp1(0:N_e,Q_sec',0:0.5:N_e)'
% 
% N_e=2*N_e
% L_e=L_e/2

% Q_sec=Q_sec_M(:,:,66); 
figure(1)
C_coord=[.0 .45 .74; .74 .45 .0];
L_co=5;

r_node_1=Q_sec(1:3,:);
dr_dx_1=Q_sec(4:6,:);
dr_dy_1=Q_sec(7:9,:);
dr_dz_1=Q_sec(10:12,:);


hold on
plot3(r_node_1(1,:),r_node_1(2,:),r_node_1(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dy_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dy_1(2,:)], ...
    [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dy_1(3,:)],'-','color',C_coord(1,:),'linewidth',1)
plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dz_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dz_1(2,:)], ...
    [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dz_1(3,:)],'-','color',C_coord(2,:),'linewidth',1)


r_node_2=[1 -1 -1]'.*Q_sec(1:3,:);
dr_dx_2=[-1 1 1]'.*Q_sec(4:6,:);
dr_dy_2=[1 -1 -1]'.*Q_sec(7:9,:);
dr_dz_2=[-1 1 1]'.*Q_sec(10:12,:);

plot3(r_node_2(1,:),r_node_2(2,:),r_node_2(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dy_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dy_2(2,:)], ...
    [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dy_2(3,:)],'-','color',C_coord(1,:),'linewidth',1)
plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dz_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dz_2(2,:)], ...
    [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dz_2(3,:)],'-','color',C_coord(2,:),'linewidth',1)

box on
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
axis([-70 70 -70 70 -50 50])
view([1 1 1])
hold off





%plot its original configuration
N_node=3;

%%
%set the boundary conditional

%Using the de_demension matrix to define the constrain of the structure
Proj_M=eye(12*(N_e+1),12*(N_e+1));%dedemesionalize the variables

D_fix=0*Q_sec;%fixed dimensions/varibales
%BC1 Antisymmetric twisted folding

% r1 start point

% position:
% control: X
% Fixed: Y Z
D_fix(1,1)=1;%control the X disp
D_fix(2,1)=1;%fix the Y disp
D_fix(3,1)=1;%fix the  Z disp

% Rotation
% fix the local y axis
D_fix(7:9,1)=1;%fix the rotation

% r2
% position:
% Free: X
% Fixed: Y Z
D_fix(1,end)=1;%control the Y disp
D_fix(2,end)=1;%fix the X disp
D_fix(3,end)=1;%fix the  Z disp

% Rotation
% Fix: local y axis 
D_fix(7:9,end)=1;%fix the rotation

%r3
% position: 
N_m=N_e/2+1
% control Y disp: 
% fixed Z disp:
D_fix(2,N_m)=1;% Y disp in control
D_fix(3,N_m)=1;% Z disp fixed
B_index=find(D_fix==0);
%
Proj_M=Proj_M(:,B_index);

% displacement x

% Angle_y=linspace(0,90,19)*pi/180;
%define the start points and frames
P_m_0=Q_sec(1:3,N_m);
disp_x=linspace(0,50,101);
U=0*disp_x;
for mm=1:101
    mm
    if mm~=1
       Q_sec=Q_sec_M(:,:,mm-1);    
    end
    Q_sec(2,N_m)=P_m_0(2)-disp_x(mm);
    if ismember(mm,[1:10])
        Q_sec(1,N_m)=Q_sec(1,N_m)+2;
%         Q_sec(3,N_m)=Q_sec(3,N_m)+1;
    else
        Q_sec(1,N_m)=Q_sec(1,N_m);
    end
    %define the end points and frames
    kk=0;
    flag1=1;
    flag3=0;
    while (kk==0||flag1>0.00001||flag2>0.10)&&kk<4000
        kk=kk+1;
        [dFq_dq_sec,Fq_sec(:,mm),q_node,U(mm)] = Jocob_rod_sec(Q_sec,N_e,L_e,N_node,Par_E,A);
        dFq_dq_sec_2=Proj_M'*dFq_dq_sec*Proj_M;
        dQ_sec_2=-pinv(dFq_dq_sec_2)*Proj_M'*Fq_sec(:,mm);
        dQ_sec=Proj_M*dQ_sec_2;
        dQ_sec=reshape(dQ_sec,12,[]);
        Q_sec=Q_sec+0.5*dQ_sec;
        flag1=sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e+1);
        flag2=mean(abs(Proj_M'*Fq_sec(:,mm)));
        if kk==4000
            flag3=1
        elseif mod(kk,1000)==1
            kk
        end
    end
    
    Angle_T_1(mm)=atan2d(Q_sec(9,1),Q_sec(7,1));
    Angle_T_2(mm)=atan2d(Q_sec(9,end),-Q_sec(8,end));
    Q_sec_M(:,:,mm)=Q_sec;
    
    % cla
    if 1==1 %mod(mm,1)==0
        figure(1)
        C_coord=[.0 .45 .74; .74 .45 .0];
        for ii=1:N_e
            r_node_1=q_node(1:3,:,ii);
            dr_dx_1=q_node(4:6,:,ii);
            dr_dy_1=q_node(7:9,:,ii);
            dr_dz_1=q_node(10:12,:,ii);
            L_co=2;
            figure(1)
            hold on
            plot3(r_node_1(1,:),r_node_1(2,:),r_node_1(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
            plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dy_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dy_1(2,:)], ...
                [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dy_1(3,:)],'-','color',C_coord(1,:),'linewidth',1)
            plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dz_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dz_1(2,:)], ...
                [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dz_1(3,:)],'-','color',C_coord(2,:),'linewidth',1)
            r_node_2=[1 -1 -1]'.*r_node_1;
            dr_dx_2=[-1 1 1]'.*dr_dx_1;
            dr_dy_2=[1 -1 -1]'.*dr_dy_1;
            dr_dz_2=[-1 1 1]'.*dr_dz_1;
            plot3(r_node_2(1,:),r_node_2(2,:),r_node_2(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
            plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dy_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dy_2(2,:)], ...
                [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dy_2(3,:)],'-','color',C_coord(1,:),'linewidth',1)
            plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dz_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dz_2(2,:)], ...
                [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dz_2(3,:)],'-','color',C_coord(2,:),'linewidth',1)

        end
        box on
        grid on
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        axis equal
        view([1 1 1])
        axis([-70 70 -100 100 -50 50])
        title(num2str(mm))
        pause(0.1)
    end

end


save("Rod_Folding_Twisting_2.mat","Q_sec_M",'U_M')
% reshape(Q_sec_M([10 12],1,:),2,[])
% for ii=1:51
%     Angle_T_1(ii)=atan2d(Q_sec_M(10,1,ii),Q_sec_M(12,1,ii));
%     Angle_T_2(ii)=atan2d(-Q_sec_M(10,end,ii),Q_sec_M(12,end,ii));
% end
% Angle_T=Angle_T_1-Angle_T_2
% figure (2)
% hold on
% % plot3(disp_x(1:7)*0+Angle_T(1:7),disp_x(1:7),U(1:7),'-o')
% plot3(Angle_T,disp_x,U,'-o')
% grid on
% box on
% xlabel('TAngle')
% ylabel('disp_x')
% zlabel('U')
% 
% figure (3)
% hold on
% % plot3(disp_x(1:7)*0+Angle_T(1:7),disp_x(1:7),-Fq_sec(1,1:7),'-o')
% plot3(Angle_T,disp_x,-Fq_sec(1,:),'-o')
% grid on
% box on
% xlabel('TAngle')
% ylabel('disp_x')
% zlabel('F_x')



%%








